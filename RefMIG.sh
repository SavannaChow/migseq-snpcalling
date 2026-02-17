#!/bin/bash
# ==============================================================================
# MIG Analysis Full Pipeline - Smart Hybrid Mode (v1.7) -- RefMIG
# ==============================================================================

# --- 自動偵測「目前執行 shell」設定檔路徑 ---
# 注意：腳本由 bash 執行時不可 source zshrc，否則 zsh 指令會報錯。
if [ -n "${BASH_VERSION:-}" ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        CONF_FILE="$HOME/.bash_profile"
    else
        CONF_FILE="$HOME/.bashrc"
    fi
elif [ -n "${ZSH_VERSION:-}" ]; then
    CONF_FILE="$HOME/.zshrc"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    CONF_FILE="$HOME/.bash_profile"
else
    CONF_FILE="$HOME/.bashrc"
fi
[ ! -f "$CONF_FILE" ] && touch "$CONF_FILE"
source "$CONF_FILE"

# ------------------------------------------------------------------------------
# 環境依賴檢查 (Dependency Check)
# ------------------------------------------------------------------------------
ENV_CHECK_FILE=".pipeline_env_ready"

check_dependencies() {
    if [ -f "$ENV_CHECK_FILE" ]; then
        echo "偵測env檢查標記,跳過驗證。"
        return 0
    fi

    echo "chech env dependencies"

    local tools=("fastp" "bwa" "samtools" "parallel" "angsd" "ngsLD" "prune_graph" "bcftools" "Rscript")

    local missing_tools=()
    local manual_install=()

    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [ ${#missing_tools[@]} -eq 0 ]; then
        echo "所有套件已就緒。"
        touch "$ENV_CHECK_FILE"
    else
        echo "缺少以下套件: ${missing_tools[*]}"
        for m_tool in "${missing_tools[@]}"; do
            local action=""
            case "$m_tool" in
                fastp) action="conda install -c bioconda fastp" ;;
                bwa) action="sudo apt-get install bwa" ;;
                samtools) action="sudo apt-get install samtools" ;;
                parallel) action="sudo apt-get install parallel" ;;
                angsd) action="http://www.popgen.dk/angsd/index.php/Installation" ;;
                ngsLD) action="https://github.com/fgvieira/ngsLD" ;;
                prune_graph) action="https://github.com/fgvieira/ngsLD (Included in ngsLD)" ;;
                bcftools) action="sudo apt-get install bcftools" ;;
                Rscript) action="sudo apt-get install r-base" ;;
                *) action="";;
            esac
            if [[ "$action" == http* ]]; then
                echo "  - $m_tool: 請手動編譯安裝，參考網址: $action"
                manual_install+=("$m_tool")
            else
                read -p "  - 偵測到 $m_tool 沒有安裝，是否嘗試自動安裝看看? (y/n): " do_install
                if [[ "$do_install" == "y" ]]; then
                    eval "$action"
                else
                    manual_install+=("$m_tool")
                fi
            fi
        done
    fi

    if [ ${#manual_install[@]} -ne 0 ]; then
        echo "錯誤：請先解決上述手動安裝套件後再運行。"
        exit 1
    fi
}

# ------------------------------------------------------------------------------
# 常數與共用設定
# ------------------------------------------------------------------------------
HEADER_TEXT="===========================================================================
   MIG-seq Genomic Analysis Pipeline (Refgenome mapping)
   開發者：Savanna Chow (savanna201@gmail.com) 使用 Gemini 協助開發
   分析邏輯基於 https://github.com/jamesfifer/JapanRE
   Credit: AllenChen's lab, Biodiversity Research Center, Academia Sinica
==========================================================================="

STAGE1="01_Trimming"
STAGE2="02_Alignment"
STAGE3="03_PCA_Analysis"
STAGE4="04_Clone_Detection"
STAGE5="05_SNP_Calling"

THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu)
JOBS=$(( THREADS / 4 )); [ "$JOBS" -lt 1 ] && JOBS=1

RUN_S1=n; RUN_S2=n; RUN_S3=n; RUN_S4=n; RUN_S5=n; RUN_S6=n
RUN_MODE="2"
AUTO_PCA_CHOICE="2"
AUTO_CLONE_CHOICE="2"
CHAIN_STAGES=false
RUN_SCOPE=""

PROJECT_NAME=""
BASE_PROJECT_NAME=""
RAW_PATH=""
REF_GENOME=""
BAM_LIST=""
SUMMARY_INPUT_DIR=""
MAPPED_BAM_DIR=""
TRIM_INPUT_DIR=""
LD_SITES_INPUT=""

# ------------------------------------------------------------------------------
# 內嵌模組: Genome Fetcher (原邏輯保留)
# ------------------------------------------------------------------------------
run_fetch_genome_module() {
    echo "======================================================="
    echo "啟動 NCBI Genome 下載模組"
    echo "======================================================="

    BASE_DIR=~/RefGenome
    mkdir -p "$BASE_DIR"

    ENV_CHECK_FILE_FG="$BASE_DIR/.env_verified"
    REQUIRED_CMDS=("esearch" "bwa" "samtools" "curl" "gunzip")

    if [[ ! -f "$ENV_CHECK_FILE_FG" ]]; then
        echo "Initializing environment check for Genome Fetcher..."

        EXTRA_PATHS=("/usr/local/bin" "/opt/bin" "$HOME/bin" "$HOME/edirect")
        for p in "${EXTRA_PATHS[@]}"; do
            [[ -d "$p" ]] && export PATH="$p:$PATH"
        done

        for cmd in "${REQUIRED_CMDS[@]}"; do
            if ! command -v "$cmd" &> /dev/null; then
                if [[ "$cmd" == "esearch" ]]; then
                    echo "EDirect not found. Installing..."
                    sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
                    export PATH="$HOME/edirect:$PATH"
                else
                    echo "Error: Required command '$cmd' is missing."
                    return 1
                fi
            fi
        done

        touch "$ENV_CHECK_FILE_FG"
        echo "Environment verified and cached."
    fi

    while true; do
        read -p "請輸入搜尋關鍵字 (物種名或 BioProject, 輸入 q 離開): " QUERY

        if [[ "$QUERY" == "q" || "$QUERY" == "Q" ]]; then
            echo "取消下載。"
            return 0
        fi

        TEMP_DATA=$(mktemp)

        TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
        HISTORY_FILE="$BASE_DIR/SearchRecord_${QUERY// /_}_${TIMESTAMP}.txt"

        echo "正在檢索 NCBI 資料庫並分析數據結構..."
        esearch -db assembly -query "$QUERY" \
        | esummary \
        | xtract -pattern DocumentSummary -def "NA" \
          -element AssemblyAccession AssemblyName AssemblyStatus AssemblyType \
                   ScaffoldN50 Coverage Isolate RefSeq_category \
                    SubmissionDate LastUpdateDate SubmitterOrganization FtpPath_GenBank \
          -block Stat -if @category -equals total_length -element Stat > "$TEMP_DATA"

        if [ ! -s "$TEMP_DATA" ]; then
            echo "找不到符合結果。"
            rm "$TEMP_DATA"
            continue
        fi

        awk -F'\t' '{
            total_mb = $13 / 1000000
            printf "-------------------- [ Index: %-4d ] --------------------\n", NR
            printf "Submitter: %s(%s)\n", $11, $3
            printf "Assembly Accession | Assembly Name : %s | %s\n", $1, $2
            printf "Assembly Status    | Assembly Type : %s | %s\n", $3, $4
            printf "Scaffold N50       | Coverage      : %s | %s X\n", $5, $6
            printf "Isolate            | RefSeq Class  : %s | %s\n", $7, $8
            printf "Submission Date    | Update Date   : %s | %s\n", $9, $10
            printf "Total Genome Size  : %.0f MB\n", total_mb
            printf "NCBI FTP 路徑: %s\n\n", $12
        }' "$TEMP_DATA" | tee "$HISTORY_FILE"

        echo "搜尋結果已存檔至: $HISTORY_FILE"

        read -p "選擇要下載的Genome編號 (輸入 r 重新搜尋, q 離開): " INDEX

        if [[ "$INDEX" == "q" ]]; then
            rm "$TEMP_DATA"
            return 0
        fi

        if [[ "$INDEX" == "r" ]]; then
            rm "$TEMP_DATA"
            echo "重新開始搜尋..."
            continue
        fi

        SELECTED_LINE=$(sed -n "${INDEX}p" "$TEMP_DATA")

        if [ -z "$SELECTED_LINE" ]; then
            echo "Invalid selection."
            rm "$TEMP_DATA"
            continue
        fi

        ACCESSION=$(echo "$SELECTED_LINE" | cut -f1)
        rm "$TEMP_DATA"
        break
    done

    echo "擷取 $ACCESSION FTP位置..."
    FTP_BASE=$(esearch -db assembly -query "$ACCESSION" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | tr ' ' '\n' | grep "$ACCESSION" | head -n 1 | sed 's|^ftp://|https://|')

    if [ -z "$FTP_BASE" ]; then
        echo "FTP path not found."
        return 1
    fi

    TARGET_DIR="$BASE_DIR/$ACCESSION"
    mkdir -p "$TARGET_DIR"

    FILE_NAME=$(basename "$FTP_BASE")
    DOWNLOAD_FILE="${FILE_NAME}_genomic.fna.gz"
    FULL_URL="${FTP_BASE}/$DOWNLOAD_FILE"

    echo "下載Genome $ACCESSION 至 $TARGET_DIR..."
    wget -c --tries=0 -P "$TARGET_DIR" "$FULL_URL"

    echo "--------------------------------------------------"
    echo "基因組 $ACCESSION 下載完成至 $TARGET_DIR。"
    read -p "是否繼續設定環境變數與indexing？(y/n): " CONTINUE_PROC
    if [[ "$CONTINUE_PROC" != "y" ]]; then
        echo "已終止。檔案保留在 $TARGET_DIR。"
        return 0
    fi

    clear

    if [ -n "$SELECTED_LINE" ]; then
        echo "=================================================="
        echo "           [ 已下載基因組詳細資訊 ]"
        echo "--------------------------------------------------"
        RAW_ACC=$(echo "$SELECTED_LINE" | awk -F'\t' '{print $1}')
        RAW_NAME=$(echo "$SELECTED_LINE" | awk -F'\t' '{print $2}')

        CLEAN_ACC=$(echo "$RAW_ACC" | sed 's/[^a-zA-Z0-9]/_/g' | sed 's/_\+/_/g')
        CLEAN_NAME=$(echo "$RAW_NAME" | sed 's/[^a-zA-Z0-9]/_/g' | sed 's/_\+/_/g')
        SUGGESTED_INPUT="${CLEAN_NAME}_${CLEAN_ACC}"

        echo "$SELECTED_LINE" | awk -F'\t' '{
            total_mb = $13 / 1000000
            printf "Submitter          : %s(%s)\n", $11, $3
            printf "Assembly Accession : %s\n", $1
            printf "Assembly Name      : %s\n", $2
            printf "Assembly Status    : %s | Type: %s\n", $3, $4
            printf "Scaffold N50       : %s\n", $5
            printf "Coverage           : %s X\n", $6
            printf "Isolate            : %s\n", $7
            printf "RefSeq Class       : %s\n", $8
            printf "Submission Date    : %s\n", $9
            printf "Total Genome Size  : %.0f MB\n", total_mb
        }'
        echo "=================================================="
    fi

    source "$CONF_FILE"
    MAPFILE=()
    MAPVAL=()
    while IFS='=' read -r name value; do
        if [[ "$value" =~ \.(fa|fasta|fna)$ ]]; then
            MAPFILE+=("$name")
            MAPVAL+=("$value")
        fi
    done < <(env)

    echo ""
    echo "--- 可用的參考基因組 ---"
    for i in "${!MAPFILE[@]}"; do
        echo "$((i+1))) \$${MAPFILE[$i]} (${MAPVAL[$i]})"
    done
    echo ""
    echo "--------------------------------------------------"
    echo "              請輸入參考基因組名稱                   "
    echo "--------------------------------------------------"
    echo "不要跟現有的重複，也不可使用空白或特殊字元"
    echo "嚴禁包含小數點（.）、破折號（-）或其他標點符號"
    echo "只能包含字母（a-z, A-Z）、數字（0-9）以及下底線（_）。"
    echo "--------------------------------------------------"
    echo "請命名剛剛下載好的參考基因組名稱"
    echo "或是按Enter直接無腦使用建議值-> ${SUGGESTED_INPUT}"
    read -p "請輸入：" USER_INPUT

    if [ -z "$USER_INPUT" ]; then
        FINAL_INPUT="$SUGGESTED_INPUT"
    else
        FINAL_INPUT=$(echo "$USER_INPUT" | sed 's/[^a-zA-Z0-9]/_/g' | sed 's/_\+/_/g')
    fi

    ENV_VAR="Ref_${FINAL_INPUT}"
    FNA_FILE="$TARGET_DIR/${FILE_NAME}_genomic.fna"
    ABS_PATH=$(realpath "$FNA_FILE")

    echo "最終變數名稱將設定為: $ENV_VAR"

    echo "解壓縮..."
    gunzip -f "$TARGET_DIR/$DOWNLOAD_FILE"

    echo "建立索引 Building BWA index (this may take a while)..."
    echo "bwa index $ABS_PATH"
    bwa index "$ABS_PATH"
    echo "寫入環境變數$CONF_FILE"
    echo "export $ENV_VAR=\"$ABS_PATH\"" >> "$CONF_FILE"
    source "$CONF_FILE"
    echo "環境變數 '$ENV_VAR' 已加入 $CONF_FILE。"
    echo "--------------------------------------------------"
    echo "Process complete."
    echo "Genome path: $ABS_PATH"
    echo "處理成功。環境變數已寫入 $CONF_FILE"
    echo "請關閉視窗，重新載入分析才能讀到新的Ref Genome位置"
    exec $SHELL
}

# ------------------------------------------------------------------------------
# 輔助函式
# ------------------------------------------------------------------------------
ask_to_run() {
    local step_name=$1
    local check_target=$2
    local skip_var_name=$3
    local exists=false
    if [ -f "$check_target" ]; then
        exists=true
    elif [ -d "$check_target" ] && [ "$(ls -A "$check_target" 2>/dev/null)" ]; then
        exists=true
    fi

    if [ "$exists" = true ]; then
        echo "-------------------------------------------------------"
        echo "[偵測到已存在的分析結果]: $step_name"
        read -p "是否跳過此步驟，並使用現有結果？ (y:跳過 / n:重新分析): " choice < /dev/tty
        [[ "$choice" == "y" || "$choice" == "Y" ]] && eval "$skip_var_name=true" || eval "$skip_var_name=false"
    else
        eval "$skip_var_name=false"
    fi
}

setup_output_dirs() {
    mkdir -p "$STAGE1/trim" "$STAGE1/fastp_report" "$STAGE2/bam" "$STAGE2/mapped_bam" "$STAGE2/mapping_results" "$STAGE3" "$STAGE4" "$STAGE5"
}

start_logging() {
    local current_time
    current_time=$(date +"%Y%m%d_%H%M%S")
    LOG_DIR="00_Logs"
    mkdir -p "$LOG_DIR"
    LOG_FILE="$LOG_DIR/${PROJECT_NAME}_${current_time}.log"

    echo "$HEADER_TEXT" > "$LOG_FILE"
    echo "專案名稱: $PROJECT_NAME" >> "$LOG_FILE"
    echo "啟動時間: $(date)" >> "$LOG_FILE"
    echo "-------------------------------------------------------" >> "$LOG_FILE"

    exec > >(tee -i -a "$LOG_FILE") 2>&1
}

print_runtime_config() {
    echo "======================================================="
    echo "  分析配置"
    echo "-------------------------------------------------------"
    echo ""
    echo "  啟動時間   : $(date)"
    echo "  專案名稱   : $PROJECT_NAME"
    echo "  執行執行緒 : $THREADS (Jobs: $JOBS)"
    [ -n "$RAW_PATH" ] && echo "  原始路徑   : $RAW_PATH"
    [ -n "$REF_GENOME" ] && echo "  參考基因組 : $REF_GENOME"
    echo "  日誌檔案   : $LOG_FILE"
    echo "  流程串接   : $([[ "$CHAIN_STAGES" == true ]] && echo "啟用(完整流程)" || echo "停用(單獨分析)")"
    echo "======================================================="
    echo ""
}

select_analysis_scope() {
    clear
    echo "$HEADER_TEXT"
    echo ""
    echo "======================================================="
    echo "  請先選擇要做哪種分析"
    echo "======================================================="
    echo "1) 執行完整分析 (Stage 1-6, 會串接上一階段輸出)"
    echo "2) 只跑 Stage 1: Fastp Trimming"
    echo "3) 只跑 Stage 2: BWA Alignment"
    echo "4) 只跑 Stage 3: PCA Outlier 分析"
    echo "5) 只跑 Stage 4: Clone Detection"
    echo "6) 只跑 Stage 5: LD Pruning"
    echo "7) 只跑 Stage 6: Final SNP Calling"
    echo "8) 自定義多階段 (不自動串接)"
    echo "q) 離開"
    echo ""
    read -p "選擇分析範圍: " RUN_SCOPE

    if [[ "$RUN_SCOPE" == "q" || "$RUN_SCOPE" == "Q" ]]; then
        echo "使用者取消操作，程式結束。"
        exit 0
    fi

    case "$RUN_SCOPE" in
        1)
            RUN_S1=y; RUN_S2=y; RUN_S3=y; RUN_S4=y; RUN_S5=y; RUN_S6=y
            CHAIN_STAGES=true
            ;;
        2) RUN_S1=y ;;
        3) RUN_S2=y ;;
        4) RUN_S3=y ;;
        5) RUN_S4=y ;;
        6) RUN_S5=y ;;
        7) RUN_S6=y ;;
        8)
            read -p "執行 Stage 1 Trimming? (y/n): " RUN_S1
            read -p "執行 Stage 2 Alignment? (y/n): " RUN_S2
            read -p "執行 Stage 3 PCA Outlier Filtering? (y/n): " RUN_S3
            read -p "執行 Stage 4 Clone Filtering? (y/n): " RUN_S4
            read -p "執行 Stage 5 LD Pruning? (y/n): " RUN_S5
            read -p "執行 Stage 6 Final SNP Calling? (y/n): " RUN_S6
            CHAIN_STAGES=false
            ;;
        *)
            echo "錯誤：無效選項。"
            exit 1
            ;;
    esac
}

configure_project_name() {
    echo "請輸入專案名稱,不要有特殊或空白字元"
    echo "*分析產生的檔案都將以專案名稱為開頭"
    read -p "請輸入: " BASE_PROJECT_NAME

    if [[ "$CHAIN_STAGES" == true ]]; then
        PROJECT_NAME="$BASE_PROJECT_NAME"
    else
        PROJECT_NAME="${BASE_PROJECT_NAME}_standalone_$(date +%Y%m%d_%H%M%S)"
        echo "[單獨分析模式] 為避免覆蓋既有輸出，本次專案名稱自動調整為: $PROJECT_NAME"
    fi
}

select_ref_genome() {
    while true; do
        source "$CONF_FILE"
        MAPFILE=()
        MAPVAL=()

        while IFS='=' read -r name value; do
            if [[ "$value" =~ \.(fa|fasta|fna)$ ]]; then
                MAPFILE+=("$name")
                MAPVAL+=("$value")
            fi
        done < <(env)

        echo ""
        echo "--- 可用的參考基因組 ---"
        for i in "${!MAPFILE[@]}"; do
            echo "$((i+1))) \$${MAPFILE[$i]} (${MAPVAL[$i]})"
        done

        MANUAL_OPTION=$(( ${#MAPFILE[@]} + 1 ))
        echo "------------------------"
        echo "d) 下載/新增參考基因組 (呼叫 NCBI Fetcher)"
        echo "$MANUAL_OPTION) 手動輸入絕對路徑"
        echo "q) 離開程式"
        echo "------------------------"

        read -p "請選擇參考基因組 (1-$MANUAL_OPTION, d, 或 q): " REF_CHOICE

        if [[ "$REF_CHOICE" == "d" || "$REF_CHOICE" == "D" ]]; then
            run_fetch_genome_module
            echo ""
            echo ">>> 返回主選單，正在刷新參考基因組清單..."
            continue
        elif [[ "$REF_CHOICE" == "q" || "$REF_CHOICE" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        fi

        if [[ "$REF_CHOICE" =~ ^[0-9]+$ ]] && [ "$REF_CHOICE" -ge 1 ] && [ "$REF_CHOICE" -le "${#MAPFILE[@]}" ]; then
            REF_GENOME="${MAPVAL[$((REF_CHOICE-1))]}"
            break
        elif [[ "$REF_CHOICE" == "$MANUAL_OPTION" ]]; then
            read -e -p "請輸入絕對路徑 (或輸入 'b' 返回選單): " REF_GENOME
            if [[ "$REF_GENOME" == "b" || "$REF_GENOME" == "B" ]]; then
                echo "返回選單..."
                continue
            fi
            break
        else
            echo "錯誤：無效的選擇。"
        fi
    done

    REF_GENOME=$(echo "$REF_GENOME" | sed "s/['\"]//g")

    if [ -z "$REF_GENOME" ] || [ ! -f "$REF_GENOME" ]; then
        echo "錯誤：檔案不存在於指定路徑：$REF_GENOME"
        exit 1
    fi

    if [ ! -f "${REF_GENOME}.bwt" ]; then
        echo "建立 BWA index..."
        bwa index "$REF_GENOME" || { echo "BWA index 失敗"; exit 1; }
    fi

    if [ ! -f "${REF_GENOME}.fai" ]; then
        echo "建立 Samtools index..."
        samtools faidx "$REF_GENOME" || { echo "Samtools faidx 失敗"; exit 1; }
    fi

    echo "使用參考基因組: $REF_GENOME"
}

collect_inputs() {
    if [[ "$RUN_S3" == "y" || "$RUN_S4" == "y" ]]; then
        echo ""
        echo "-------------------------------------------------------"
        echo "請選擇分析運行模式："
        echo "1) 自動模式 (預設決策)"
        echo "2) 互動模式 (每次詢問)"
        echo ""
        read -p "請輸入選項 (1, 2, 或 q 離開): " RUN_MODE
        [[ "$RUN_MODE" == "q" || "$RUN_MODE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }

        if [ "$RUN_MODE" == "1" ]; then
            read -p "  > 偵測到 PCA Outlier 時處理方式 (1:移除, 2:保留, q:離開): " AUTO_PCA_CHOICE
            [[ "$AUTO_PCA_CHOICE" == "q" || "$AUTO_PCA_CHOICE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }
            read -p "  > 偵測到 Clone 樣本時處理方式 (1:移除, 2:保留, q:離開): " AUTO_CLONE_CHOICE
            [[ "$AUTO_CLONE_CHOICE" == "q" || "$AUTO_CLONE_CHOICE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }
        fi
    fi

    if [[ "$RUN_S1" == "y" ]]; then
        read -e -p "請輸入 Stage1 原始序列 (raw data) 資料夾路徑: " RAW_PATH
        [ ! -d "$RAW_PATH" ] && { echo "錯誤：找不到路徑 $RAW_PATH"; exit 1; }
        RAW_PATH=$(realpath "$RAW_PATH")
    fi

    if [[ "$RUN_S2" == "y" ]]; then
        if [[ "$CHAIN_STAGES" == true ]]; then
            TRIM_INPUT_DIR="$STAGE1/trim"
        else
            read -e -p "請輸入 Stage2 要用的 trimmed fastq 資料夾路徑: " TRIM_INPUT_DIR
            [ ! -d "$TRIM_INPUT_DIR" ] && { echo "錯誤：找不到路徑 $TRIM_INPUT_DIR"; exit 1; }
            TRIM_INPUT_DIR=$(realpath "$TRIM_INPUT_DIR")
        fi
    fi

    if [[ "$RUN_S2" == "y" || "$RUN_S5" == "y" || "$RUN_S6" == "y" ]]; then
        select_ref_genome
    fi

    if [[ "$RUN_S3" == "y" ]]; then
        if [[ "$CHAIN_STAGES" == true ]]; then
            SUMMARY_INPUT_DIR="$STAGE2/mapping_results"
            MAPPED_BAM_DIR="$STAGE2/mapped_bam"
        else
            read -e -p "請輸入 Stage3 mapping_results 目錄路徑: " SUMMARY_INPUT_DIR
            read -e -p "請輸入 Stage3 mapped_bam 目錄路徑: " MAPPED_BAM_DIR
            [ ! -d "$SUMMARY_INPUT_DIR" ] && { echo "錯誤：找不到路徑 $SUMMARY_INPUT_DIR"; exit 1; }
            [ ! -d "$MAPPED_BAM_DIR" ] && { echo "錯誤：找不到路徑 $MAPPED_BAM_DIR"; exit 1; }
            SUMMARY_INPUT_DIR=$(realpath "$SUMMARY_INPUT_DIR")
            MAPPED_BAM_DIR=$(realpath "$MAPPED_BAM_DIR")
        fi
    fi

    if [[ "$RUN_S4" == "y" ]]; then
        if [[ "$CHAIN_STAGES" == true ]]; then
            :
        else
            read -e -p "請輸入 Stage4 clone detection 要使用的 BAM list (.bamfile) 路徑: " BAM_LIST_CLONE_INPUT
            [ ! -s "$BAM_LIST_CLONE_INPUT" ] && { echo "錯誤：找不到或空檔 $BAM_LIST_CLONE_INPUT"; exit 1; }
            BAM_LIST_CLONE_INPUT=$(realpath "$BAM_LIST_CLONE_INPUT")
        fi
    fi

    if [[ "$RUN_S5" == "y" && "$CHAIN_STAGES" != true ]]; then
        read -e -p "請輸入 Stage5 LD pruning 要使用的 BAM list (.bamfile) 路徑: " BAM_LIST_LD_INPUT
        [ ! -s "$BAM_LIST_LD_INPUT" ] && { echo "錯誤：找不到或空檔 $BAM_LIST_LD_INPUT"; exit 1; }
        BAM_LIST_LD_INPUT=$(realpath "$BAM_LIST_LD_INPUT")
    fi

    if [[ "$RUN_S6" == "y" && "$CHAIN_STAGES" != true ]]; then
        read -e -p "請輸入 Stage6 Final SNP 要使用的 BAM list (.bamfile) 路徑: " BAM_LIST_FINAL_INPUT
        [ ! -s "$BAM_LIST_FINAL_INPUT" ] && { echo "錯誤：找不到或空檔 $BAM_LIST_FINAL_INPUT"; exit 1; }
        BAM_LIST_FINAL_INPUT=$(realpath "$BAM_LIST_FINAL_INPUT")

        read -e -p "請輸入 Stage6 要使用的 LD pruned sites 檔案路徑: " LD_SITES_INPUT
        [ ! -f "$LD_SITES_INPUT" ] && { echo "錯誤：找不到檔案 $LD_SITES_INPUT"; exit 1; }
        LD_SITES_INPUT=$(realpath "$LD_SITES_INPUT")
    fi
}

confirm_run() {
    clear
    [[ "$RUN_MODE" == "1" ]] && MODE_STR="自動模式 (Auto)" || MODE_STR="互動模式 (Manual)"

    echo "                      執行前最終確認"
    echo ""
    echo "======================================================="
    echo "  分析參數確認"
    echo "-------------------------------------------------------"
    echo ""
    printf "  %-15s : %s\n" "專案名稱" "$PROJECT_NAME"
    [ -n "$RAW_PATH" ] && printf "  %-15s : %s\n" "原始資料路徑" "$RAW_PATH"
    [ -n "$REF_GENOME" ] && printf "  %-15s : %s\n" "參考基因組" "$REF_GENOME"
    [[ "$RUN_S3" == "y" || "$RUN_S4" == "y" ]] && printf "  %-15s : %s\n" "執行模式" "$MODE_STR"
    printf "  %-15s : %s\n" "流程串接" "$([[ "$CHAIN_STAGES" == true ]] && echo "啟用" || echo "停用")"
    echo "  運行範圍       :"
    [[ "$RUN_S1" == "y" ]] && echo "    - Stage 1 Fastp Trimming"
    [[ "$RUN_S2" == "y" ]] && echo "    - Stage 2 BWA Alignment"
    [[ "$RUN_S3" == "y" ]] && echo "    - Stage 3 PCA Outlier Filtering"
    [[ "$RUN_S4" == "y" ]] && echo "    - Stage 4 Clone Identification"
    [[ "$RUN_S5" == "y" ]] && echo "    - Stage 5 LD Pruning"
    [[ "$RUN_S6" == "y" ]] && echo "    - Stage 6 Final SNP Calling"

    if [[ "$RUN_MODE" == "1" && ( "$RUN_S3" == "y" || "$RUN_S4" == "y" ) ]]; then
        printf "  %-15s : %s\n" "PCA Outlier 決策" "$([[ "$AUTO_PCA_CHOICE" == "1" ]] && echo "移除" || echo "保留")"
        printf "  %-15s : %s\n" "Clone 樣本決策" "$([[ "$AUTO_CLONE_CHOICE" == "1" ]] && echo "移除" || echo "保留")"
    fi

    printf "  %-15s : %s (Jobs: %s)\n" "執行緒" "$THREADS" "$JOBS"
    printf "  %-15s : %s\n" "日誌檔案" "$LOG_FILE"
    echo ""
    echo "-------------------------------------------------------"
    read -p "  以上是否正確？ (y: 開始執行 / n: 結束): " CONFIRM

    if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
        echo "[取消] 使用者中止程式。"
        exit 0
    fi
}

# ------------------------------------------------------------------------------
# Stage function 定義
# ------------------------------------------------------------------------------
run_stage1_fastp() {
    ask_to_run "Fastp Quality Control" "$STAGE1/trim" SKIP_FASTP
    if [[ "$SKIP_FASTP" == true ]]; then
        return
    fi

    echo "執行 Fastp..."
    ls "${RAW_PATH}"/*_R1_001.fast* > raw_list.txt
    parallel -j "$JOBS" --bar "
      r1=\"{}\"
      r2=\$(echo \"\$r1\" | sed 's/_R1_/_R2_/')
      ext=\$(basename \"\$r1\" | sed 's/.*_R1_001//')
      base=\$(basename \"\$r1\" _R1_001\$ext)
      fastp -i \"\$r1\" -I \"\$r2\" -o \"$STAGE1/trim/\${base}_R1_001.fastq.gz\" -O \"$STAGE1/trim/\${base}_R2_001.fastq.gz\" \
        --thread 2 --qualified_quality_phred 30 --length_required 80 \
        --html \"$STAGE1/fastp_report/\${base}.html\" --json \"$STAGE1/fastp_report/\${base}.json\"
    " < raw_list.txt

    N_TRIM=$(wc -l < raw_list.txt)
    echo "-------------------------------------------------------"
    echo "[Stage 1 完成回報]"
    echo "處理樣本總數: $N_TRIM"
    echo "清理後的檔案目錄: $STAGE1/trim/"
    echo "Fastp報告目錄: $STAGE1/fastp_report/"
    echo "-------------------------------------------------------"
}

run_stage2_alignment() {
    ask_to_run "BWA Alignment" "$STAGE2/mapped_bam" SKIP_BWA
    if [[ "$SKIP_BWA" == true ]]; then
        return
    fi

    echo "執行 BWA Mapping..."
    while read -r r1; do
        base=$(basename "$r1" _R1_001.fastq.gz)
        r2="$(dirname "$r1")/${base}_R2_001.fastq.gz"
        bwa mem -t "$THREADS" "$REF_GENOME" "$r1" "$r2" > "$STAGE2/bam/${base}.sam"
        samtools view -Sb "$STAGE2/bam/${base}.sam" > "$STAGE2/bam/${base}.bam"
        samtools view -bF4 -@ "$THREADS" "$STAGE2/bam/${base}.bam" > "$STAGE2/mapped_bam/${base}.bam"
        samtools sort -@ "$THREADS" -o "$STAGE2/mapped_bam/${base}_sorted.bam" "$STAGE2/mapped_bam/${base}.bam"
        mv "$STAGE2/mapped_bam/${base}_sorted.bam" "$STAGE2/mapped_bam/${base}.bam"
        samtools index "$STAGE2/mapped_bam/${base}.bam"
        samtools flagstat "$STAGE2/bam/${base}.bam" > "$STAGE2/mapping_results/${base}.txt"
        rm "$STAGE2/bam/${base}.sam"
    done < <(ls "$TRIM_INPUT_DIR"/*_R1_001.fastq.gz)

    N_MAPPED=$(ls "$STAGE2/mapped_bam/"*.bam | wc -l)
    echo "-------------------------------------------------------"
    echo "[Stage 2 完成回報]"
    echo "使用的參考基因組: $REF_GENOME"
    echo "比對指令參數: bwa mem -t $THREADS"
    echo "完成步驟: SAM轉換、BAM過濾(F4)、排序與建立索引"
    echo "產出的比對檔案(BAM): $STAGE2/mapped_bam/"
    echo "比對率統計結果: $STAGE2/mapping_results/"
    echo "共計完成比對樣本數: $N_MAPPED"
    echo "-------------------------------------------------------"
}

run_stage3_pca() {
    local summary_csv bam_list_local
    summary_csv="$STAGE3/${PROJECT_NAME}_mapping_summary.csv"
    bam_list_local="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"

    echo "[Stage 3] 生成比對報表與 PCA 品質檢測..."

    ls -d "$(realpath "$MAPPED_BAM_DIR")"/*.bam > "$bam_list_local"
    ask_to_run "PCA Analysis & Outlier Detection" "$summary_csv" SKIP_PCA

    if [ "$SKIP_PCA" = true ]; then
        echo ">>> 跳過 PCA 分析，嘗試載入先前的過濾結果..."
        if [ -f "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
            echo "[!] 載入既有的 PCA 過濾後清單: $(wc -l < "$STAGE3/${PROJECT_NAME}_after_pca.bamfile") 個樣本"
            if [[ "$CHAIN_STAGES" == true ]]; then
                BAM_LIST="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
            fi
        else
            echo "[!] 未發現過濾後清單，假設上次選擇保留所有樣本 (或無 Outlier)。"
            if [[ "$CHAIN_STAGES" == true ]]; then
                BAM_LIST="$bam_list_local"
            fi
        fi
        return
    fi

    echo "輸出mapping表頭與數據"
    echo "Sample,Total,Mapped,Properly_Paired,With_Mate_Mapped,Singletons,Mate_Diff_Chr,Mate_Diff_Chr_MapQ5,Secondary,Supplementary,Duplicates,Paired_in_Seq,Read1,Read2" > "$summary_csv"

    for f in "${SUMMARY_INPUT_DIR}"/*.txt; do
        sample=$(basename "$f" .txt)
        total=$(grep "in total" "$f" | head -1 | awk '{print $1}')
        mapped=$(grep " mapped (" "$f" | awk '{print $1}')
        properly_paired=$(grep "properly paired" "$f" | awk '{print $1}')
        with_mate=$(grep "with itself and mate mapped" "$f" | awk '{print $1}')
        singletons=$(grep "singletons" "$f" | awk '{print $1}')
        mate_diff_chr=$(grep "with mate mapped to a different chr$" "$f" | awk '{print $1}')
        mate_diff_chr_q5=$(grep "with mate mapped to a different chr (mapQ>=5)" "$f" | awk '{print $1}')
        secondary=$(grep " secondary" "$f" | awk '{print $1}')
        supplementary=$(grep " supplementary" "$f" | awk '{print $1}')
        duplicates=$(grep " duplicates" "$f" | awk '{print $1}')
        paired_in_seq=$(grep "paired in sequencing" "$f" | awk '{print $1}')
        read1=$(grep " read1" "$f" | awk '{print $1}')
        read2=$(grep " read2" "$f" | awk '{print $1}')
        echo "$sample,$total,$mapped,$properly_paired,$with_mate,$singletons,$mate_diff_chr,$mate_diff_chr_q5,$secondary,$supplementary,$duplicates,$paired_in_seq,$read1,$read2" >> "$summary_csv"
    done

    echo "-------------------------------------------------------"
    echo "[Stage 3 進度回報]"
    echo "Mapping Summary 資料表已生成: $summary_csv"
    echo "對應統計來源: $SUMMARY_INPUT_DIR"
    echo "-------------------------------------------------------"

    ABS_STAGE3=$(realpath "$STAGE3")
    ABS_SUMMARY_CSV=$(realpath "$summary_csv")
    ABS_BAM_LIST=$(realpath "$bam_list_local")

    echo "輸出PCA.r，使用mapping_rate, pairing_efficiency, singleton_rate三大數據做分析"
    echo "如果有疑問可以自已使用R-studio調整目錄內的PCA.r後再做plot"
    PCA_SCRIPT="$STAGE3/${PROJECT_NAME}_PCA.r"
    cat << R_CODE > "$PCA_SCRIPT"
# --- [自動環境檢查] ---
required_packages <- c("readr", "dplyr", "ggplot2", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='https://cran.csie.ntu.edu.tw/')
lapply(required_packages, library, character.only = TRUE)

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
# 1. 讀取數據 (使用絕對路徑)
df <- read_csv("${ABS_SUMMARY_CSV}", show_col_types = FALSE)
original_bams <- read.table("${ABS_BAM_LIST}", header = FALSE, stringsAsFactors = FALSE)
colnames(original_bams) <- c("FilePath")
original_bams\$Sample <- gsub(".bam$", "", basename(original_bams\$FilePath))

df_pca <- df %>%
  transmute(
    Sample = Sample,
    Total = as.numeric(Total),
    Mapped = as.numeric(Mapped),
    Properly_Paired = as.numeric(Properly_Paired),
    Singletons = as.numeric(Singletons),
    mapping_rate = ifelse(Total > 0, Mapped / Total, NA_real_),
    properly_paired_rate = ifelse(Total > 0, Properly_Paired / Total, NA_real_),
    pairing_efficiency = ifelse(Mapped > 0, Properly_Paired / Mapped, NA_real_),
    singleton_rate = ifelse(Total > 0, Singletons / Total, NA_real_),
    logTotal = log10(Total + 1)
  ) %>%
  filter(if_all(c(mapping_rate, pairing_efficiency, singleton_rate, properly_paired_rate, logTotal), ~ !is.na(.)))

# 執行 PCA
X <- df_pca %>% select(mapping_rate, pairing_efficiency, singleton_rate) %>% as.data.frame()
pca_cor <- prcomp(X, center = TRUE, scale. = TRUE)

# 提取特徵統計量
eig_val <- pca_cor\$sdev^2
prop_var <- eig_val / sum(eig_val)
cum_var <- cumsum(prop_var)

pca_stat <- data.frame(
  PC = paste0("PC", 1:length(eig_val)),
  Eigenvalue = eig_val,
  Variance_Percent = prop_var * 100,
  Cumulative_Percent = cum_var * 100
)

# 輸出詳細統計表 (絕對路徑)
write.table(pca_stat, "${ABS_STAGE3}/${PROJECT_NAME}_PCA_statistics.txt", sep="\t", quote=FALSE, row.names=FALSE)

# 2. 提取特徵向量
loadings <- as.data.frame(pca_cor\$rotation)
loadings\$Variable <- rownames(loadings)

# 3. 離群值判定
scores <- as.data.frame(pca_cor\$x[, 1:2])
scores\$Sample <- df_pca\$Sample
n <- nrow(scores); p_vars <- 2
cutoff <- qf(0.95, p_vars, n - 1) * (p_vars * (n - 1) / (n - p_vars))
scores\$dist_sq <- mahalanobis(scores[, 1:2], center = colMeans(scores[, 1:2]), cov = cov(scores[, 1:2]))
scores\$is_outlier <- scores\$dist_sq > cutoff

# 4. 繪製 PDF (包含數據點、橢圓與特徵向量)
# 計算向量縮放比例以適應數據分布
scaling_factor <- max(abs(scores[, 1:2])) / max(abs(loadings[, 1:2])) * 0.8
# 繪製 PDF
p <- ggplot() +
# 繪製樣本點與橢圓
  stat_ellipse(data = scores, aes(x = PC1, y = PC2), type = "norm", level = 0.95, linetype = "dashed", color = "darkgreen") +
  geom_point(data = scores, aes(x = PC1, y = PC2, color = is_outlier), size = 2.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_text_repel(data = subset(scores, is_outlier), aes(x = PC1, y = PC2, label = Sample), color = "red") +
# 繪製特徵向量 (特徵向量方向)
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * scaling_factor, yend = PC2 * scaling_factor),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue", alpha = 0.7) +
  geom_text(data = loadings, aes(x = PC1 * scaling_factor * 1.1, y = PC2 * scaling_factor * 1.1, label = Variable),
            color = "blue", fontface = "bold") +
# 標籤加上貢獻率
  labs(title = paste0(df_pca\$Sample[1], "... PCA Analysis"),
       x = paste0("PC1 (", round(prop_var[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(prop_var[2]*100, 1), "%)")) +
  theme_classic()

ggsave("${ABS_STAGE3}/${PROJECT_NAME}_PCA_Mapping_Quality.pdf", plot = p, width = 8, height = 7)

# 產出過濾後的結果 (絕對路徑)
valid_samples <- scores\$Sample[!scores\$is_outlier]
write.table(scores\$Sample[scores\$is_outlier], "${ABS_STAGE3}/${PROJECT_NAME}_outliers.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
to_keep_bams <- original_bams %>% filter(Sample %in% valid_samples) %>% select(FilePath)
write.table(to_keep_bams, "${ABS_STAGE3}/${PROJECT_NAME}_after_pca.bamfile", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 於 Terminal 輸出簡單統計摘要
cat("\n[PCA 統計摘要]\n")
print(pca_stat)
R_CODE

    Rscript "$PCA_SCRIPT"

    N_ORIG=$(wc -l < "$bam_list_local")
    N_AFTER=$(wc -l < "$STAGE3/${PROJECT_NAME}_after_pca.bamfile")
    N_DIFF=$((N_ORIG - N_AFTER))

    echo "-------------------------------------------------------"
    echo "[Mapping Outlier偵測結果]"
    echo "原始輸入樣本總數: $N_ORIG"
    echo "PCA 檢測建議保留樣本數: $N_AFTER"
    echo "被剔除的樣本數: $N_DIFF"
    echo "-------------------------------------------------------"

    OUTLIER_FILE="$STAGE3/${PROJECT_NAME}_outliers.txt"
    if [ -s "$OUTLIER_FILE" ]; then
        echo "偵測到 $N_DIFF 個 PCA Outliers:"
        cat "$OUTLIER_FILE"

        if [ "$RUN_MODE" == "1" ]; then
            PCA_DECISION=$AUTO_PCA_CHOICE
            echo "自動模式：套用預設選項 ($PCA_DECISION)"
        else
            read -p "是否移除Outlier樣本？(1:移除, 2:保留, q:離開): " PCA_DECISION < /dev/tty
        fi

        if [[ "$PCA_DECISION" == "q" || "$PCA_DECISION" == "Q" ]]; then
            echo "使用者選擇離開，分析中止。"
            exit 0
        fi

        if [[ "$CHAIN_STAGES" == true ]]; then
            if [ "$PCA_DECISION" == "1" ]; then
                BAM_LIST="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
                echo "[!] 已套用過濾後的 BAM 清單，當前分析為: $(wc -l < "$BAM_LIST") 個樣本。"
            else
                BAM_LIST="$bam_list_local"
                echo "[+] 已選擇保留outlier樣本，維持原始分析樣本數。"
            fi
        fi
    else
        [[ "$CHAIN_STAGES" == true ]] && BAM_LIST="$bam_list_local"
    fi
}

run_stage4_clone() {
    local clone_bam_list
    if [[ "$CHAIN_STAGES" == true ]]; then
        [ -z "$BAM_LIST" ] && BAM_LIST="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"
        clone_bam_list="$BAM_LIST"
    else
        clone_bam_list="$BAM_LIST_CLONE_INPUT"
    fi

    echo "[Stage 4] 執行 Clone 偵測分析 (ANGSD IBS)..."
    ask_to_run "Clone Detection (IBS Matrix)" "$STAGE4/${PROJECT_NAME}_clone_identification.ibsMat" SKIP_CLONE_CALC

    if [ "$SKIP_CLONE_CALC" = true ]; then
        echo ">>> 跳過 Clone 計算，嘗試載入先前去重結果..."
        if [ -f "$STAGE4/${PROJECT_NAME}_after_clones.bamfile" ]; then
            echo "[!] 載入既有的 Clone 去重後清單: $(wc -l < "$STAGE4/${PROJECT_NAME}_after_clones.bamfile") 個樣本"
            [[ "$CHAIN_STAGES" == true ]] && BAM_LIST="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
        else
            echo "[!] 未發現過濾後清單，假設上次選擇保留 Clone 或無 Clone。"
        fi
        return
    fi

    N_IND_TMP=$(wc -l < "$clone_bam_list")
    MIN_IND_TMP=$(echo "$N_IND_TMP * 0.7 / 1" | bc)

    angsd -bam "$clone_bam_list" -GL 1 -P 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -minInd "$MIN_IND_TMP" \
          -snp_pval 1e-5 -minMaf 0.05 -doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 \
          -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -out "$STAGE4/${PROJECT_NAME}_clone_identification"

    gzip -kfd "$STAGE4"/*.gz 2>/dev/null || true
    ABS_STAGE4=$(realpath "$STAGE4")
    ABS_BAM_LIST=$(realpath "$clone_bam_list")

    cat << R_CODE > "$STAGE4/${PROJECT_NAME}_identify_clones.r"
# --- [自動環境檢查] 確保 RStudio 環境具備必要套件 ---
required_packages <- c("readr", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='https://cran.csie.ntu.edu.tw/')
lapply(required_packages, library, character.only = TRUE)

library(readr)
library(dplyr)
# 1. 讀取遺傳距離數據與樣本路徑 (已轉換為絕對路徑)
# ma: 從 ANGSD 產出的 .ibsMat 讀取樣本間的遺傳距離矩陣
ma <- as.matrix(read.table("${ABS_STAGE4}/${PROJECT_NAME}_clone_identification.ibsMat"))
bams <- read.table("${ABS_BAM_LIST}", header = FALSE, stringsAsFactors = FALSE)
colnames(bams) <- c("FilePath")
bams\$Sample <- gsub(".bam$", "", basename(bams\$FilePath))
dimnames(ma) <- list(bams\$Sample, bams\$Sample)

# 2. 執行層次聚類 (Hierarchical Clustering)
# 使用平均連鎖法 (Average Linkage) 計算樣本間的親緣關係
hc <- hclust(as.dist(ma), "ave")

# 根據樣本數量動態調整 PDF 寬度
dynamic_width <- max(12, nrow(bams) * 0.2)

# 輸出原始聚類圖，不含任何門檻標記
pdf("${ABS_STAGE4}/${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf", width = dynamic_width, height = 10)
plot(hc, cex=0.7, main="Raw Clustering of Samples (IBS Distance)")
dev.off()

# 統計跳躍點偵測邏輯：自動判定潛在的 Clone 門檻
h <- hc\$height
gaps <- diff(h) # 計算分支高度間的間隙

# 定義背景雜訊區域 (取前 5 個分支高度作為底噪參考)
ref_idx <- 1:min(5, length(gaps))
noise_mean <- mean(gaps[ref_idx])
noise_sd <- sd(gaps[ref_idx])

# 尋找第一個顯著超過背景雜訊 (3倍標準差) 的跳躍點
jump_idx <- which(gaps > (noise_mean + 3 * noise_sd))[1]

# 判斷分支邏輯：若偵測到早期顯著跳躍點則作為門檻，否則使用物理約束區間
if (!is.na(jump_idx) && jump_idx <= 10) {
  threshold_h <- (h[jump_idx] + h[jump_idx + 1]) / 2
} else {
  threshold_h <- max(0.05, min(h[1] * 1.2, 0.2))
}

# 3. 執行樣本分群並建立結果表格
clusters <- cutree(hc, h = threshold_h)
df_clusters <- data.frame(Sample = bams\$Sample, FilePath = bams\$FilePath, ClusterID = clusters)

# 識別包含複本樣本的群組 (一個 ClusterID 中若有多個樣本則視為互為 Clone)
clone_cluster_ids <- which(table(clusters) > 1)

# 寫出建議移除的名單 (每個 Clone 群組僅保留一個代表樣本)
write.table(df_clusters %>% 
              filter(ClusterID %in% names(clone_cluster_ids)) %>% 
              group_by(ClusterID) %>% 
              slice(-1) %>% 
              ungroup() %>% 
              select(Sample), 
            "${ABS_STAGE4}/${PROJECT_NAME}_clones_to_review.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 寫出保留下來的樣本路徑清單 (作為後續分析用的新 BAM List)
write.table(df_clusters[!duplicated(df_clusters\$ClusterID), "FilePath"], 
            "${ABS_STAGE4}/${PROJECT_NAME}_after_clones.bamfile", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 4. 繪製帶有判定門檻（紅線）的診斷圖
pdf("${ABS_STAGE4}/${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf", width = dynamic_width, height = 10)
plot(hc, cex=0.7, main=paste("Expert-Logic Gap Detection (h =", round(threshold_h, 4), ")"))
abline(h = threshold_h, col = "red", lty = 2, lwd = 2) 
dev.off()

# 輸出統計診斷數值至終端機
cat(paste("[R] 統計診斷：底噪高度", round(h[1], 4), "| 自動門檻", round(threshold_h, 4), "\n"))
R_CODE

    Rscript "$STAGE4/${PROJECT_NAME}_identify_clones.r"

    CLONE_REVIEW_FILE="$STAGE4/${PROJECT_NAME}_clones_to_review.txt"
    N_CLONE_BEFORE=$(wc -l < "$clone_bam_list")
    N_CLONE_AFTER=$(wc -l < "$STAGE4/${PROJECT_NAME}_after_clones.bamfile")
    N_CLONE_COUNT=$(wc -l < "$CLONE_REVIEW_FILE")

    echo "-------------------------------------------------------"
    echo "[Clone樣本偵測報告]"
    echo "1. 原始 PDF (無門檻): $STAGE4/${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf"
    echo "2. 過濾 PDF (含紅線): $STAGE4/${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf"
    echo "進入偵測之樣本總數: $N_CLONE_BEFORE"
    echo "剔除Clone後樣本數: $N_CLONE_AFTER"
    echo "偵測到 $N_CLONE_COUNT 個潛在Clone樣本"
    echo "-------------------------------------------------------"

    if [ -s "$CLONE_REVIEW_FILE" ]; then
        echo "偵測到以下潛在Clone樣本 (建議移除名單):"
        cat "$CLONE_REVIEW_FILE"

        if [ "$RUN_MODE" == "1" ]; then
            CLONE_DECISION=$AUTO_CLONE_CHOICE
            echo "自動模式：套用預設選項 ($CLONE_DECISION)"
        else
            read -p "是否移除上述Clone樣本？(1:移除, 2:保留, q:離開): " CLONE_DECISION < /dev/tty
        fi

        if [[ "$CLONE_DECISION" == "q" || "$CLONE_DECISION" == "Q" ]]; then
            echo "使用者選擇離開，分析中止。"
            exit 0
        fi

        if [[ "$CHAIN_STAGES" == true ]]; then
            if [ "$CLONE_DECISION" == "1" ]; then
                BAM_LIST="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
                echo "[!] 已套用去除clones後的最終 BAM 清單，當前分析樣本數: $(wc -l < "$BAM_LIST") 個樣本。"
            else
                echo "[+] 已選擇保留Clone樣本，維持樣本數: $N_CLONE_BEFORE"
            fi
        fi
    fi
}

resolve_bam_list_for_stage5_or_6() {
    if [[ "$CHAIN_STAGES" == true ]]; then
        LIST_FINAL="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
        LIST_PCA_ONLY="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        LIST_FULL="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"

        if [ -f "$LIST_FINAL" ]; then
            BAM_LIST="$LIST_FINAL"
        elif [ -f "$LIST_PCA_ONLY" ]; then
            BAM_LIST="$LIST_PCA_ONLY"
        else
            BAM_LIST="$LIST_FULL"
        fi
    else
        if [[ "$RUN_S5" == "y" ]]; then
            BAM_LIST="$BAM_LIST_LD_INPUT"
        elif [[ "$RUN_S6" == "y" ]]; then
            BAM_LIST="$BAM_LIST_FINAL_INPUT"
        fi
    fi

    if [ ! -s "$BAM_LIST" ]; then
        echo "錯誤：指定的清單檔案「$BAM_LIST」不存在或為空。"
        exit 1
    fi

    N_IND=$(wc -l < "$BAM_LIST")
    MIN_IND=$(echo "$N_IND * 0.7 / 1" | bc)
    echo "當前分析清單：$BAM_LIST"
    echo "樣本總數：$N_IND，SNP Calling 門檻 (70%)：$MIN_IND"
}

run_stage5_ld_pruning() {
    ask_to_run "LD Pruning (Site Map)" "$STAGE5/LDpruned_snp.sites" SKIP_S5
    if [[ "$SKIP_S5" == true ]]; then
        return
    fi

    echo "[Stage 5] 執行 LD Pruning 產生連鎖不平衡位點表..."
    angsd -b "$BAM_LIST" -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -baq 1 -setMinDepth 5 -SNP_pval 1e-6 -skipTriallelic 1 -doHWE 1 -Hetbias_pval 0.00001 -minInd "$MIN_IND" -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 32 -doCounts 1 -ref "$REF_GENOME" -P 1 -out "$STAGE5/allsnps"

    gzip -kfd "$STAGE5"/*.gz
    gunzip -c "$STAGE5/allsnps.mafs.gz" | tail -n +2 | cut -f 1,2 > "$STAGE5/mc1.sites"
    N_SITES=$(wc -l < "$STAGE5/mc1.sites")

    ngsLD --geno "$STAGE5/allsnps.geno" --verbose 1 --probs 1 --n_ind "$N_IND" --n_sites "$N_SITES" --max_kb_dist 50 --pos "$STAGE5/mc1.sites" --n_threads "$THREADS" --extend_out 1 --out "$STAGE5/allsnpsites.LD"
    prune_graph --header -v -n "$THREADS" --in "$STAGE5/allsnpsites.LD" --weight-field "r2" --weight-filter "dist <=10000 && r2 >= 0.5" --out "$STAGE5/allsnpsites.pos"

    sed 's/:/\t/g' "$STAGE5/allsnpsites.pos" | awk '$2!=""' | sort -k1 > "$STAGE5/LDpruned_snp.sites"
    angsd sites index "$STAGE5/LDpruned_snp.sites"

    N_SITES_AFTER=$(wc -l < "$STAGE5/LDpruned_snp.sites")
    echo "-------------------------------------------------------"
    echo "[Stage 5 完成回報]"
    echo "LD Pruning 前的初始位點數: $N_SITES"
    echo "LD Pruning 後保留的位點數: $N_SITES_AFTER"
    echo "去LD後位點索引檔路徑: $STAGE5/LDpruned_snp.sites"
    echo "-------------------------------------------------------"
}

run_stage6_final_snp() {
    local ld_sites
    if [[ "$CHAIN_STAGES" == true ]]; then
        ld_sites="$STAGE5/LDpruned_snp.sites"
    else
        ld_sites="$LD_SITES_INPUT"
    fi

    if [ ! -f "$ld_sites" ]; then
        echo "錯誤：未發現 LD 位點表 ($ld_sites)。"
        exit 1
    fi

    ask_to_run "Final VCF" "$STAGE5/${PROJECT_NAME}_snps_final.vcf" SKIP_S6
    if [[ "$SKIP_S6" == true ]]; then
        return
    fi

    echo "[Stage 6] 執行最終 SNP Calling..."
    angsd -sites "$ld_sites" -b "$BAM_LIST" -GL 1 -P 1 -minInd "$MIN_IND" -minMapQ 20 -minQ 25 -sb_pval 1e-5 -Hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGlf 2 -dosnpstat 1 -doPost 1 -doGeno 8 -doBcf 1 --ignore-RG 0 -doHWE 1 -ref "$REF_GENOME" -out "$STAGE5/${PROJECT_NAME}_snps_final"

    bcftools view -O v -o "$STAGE5/${PROJECT_NAME}_snps_final.vcf" "$STAGE5/${PROJECT_NAME}_snps_final.bcf"
    FINAL_SNPS=$(bcftools view -H "$STAGE5/${PROJECT_NAME}_snps_final.vcf" | wc -l)

    echo "-------------------------------------------------------"
    echo "[Stage 6 完成]"
    echo "最終產出的 SNP 總數量: $FINAL_SNPS"
    echo "最終 VCF 檔案路徑: $STAGE5/${PROJECT_NAME}_snps_final.vcf"
    echo "-------------------------------------------------------"
}

# ------------------------------------------------------------------------------
# 主程序
# ------------------------------------------------------------------------------
main() {
    check_dependencies
    select_analysis_scope
    configure_project_name
    setup_output_dirs
    start_logging
    collect_inputs
    print_runtime_config
    confirm_run

    [[ "$RUN_S1" == "y" ]] && run_stage1_fastp
    [[ "$RUN_S2" == "y" ]] && run_stage2_alignment

    if [[ "$RUN_S3" == "y" ]]; then
        run_stage3_pca
    fi

    if [[ "$RUN_S4" == "y" ]]; then
        run_stage4_clone
    fi

    if [[ "$RUN_S5" == "y" || "$RUN_S6" == "y" ]]; then
        resolve_bam_list_for_stage5_or_6
    fi

    [[ "$RUN_S5" == "y" ]] && run_stage5_ld_pruning
    [[ "$RUN_S6" == "y" ]] && run_stage6_final_snp

    echo "======================================================="
    echo "分析結束: $(date)"
    echo "產出 VCF: $STAGE5/${PROJECT_NAME}_snps_final.vcf"
    echo "日誌位置: $LOG_FILE"
    echo "VCF可用PGDSpider做轉換"
    echo "https://software.bioinformatics.unibe.ch/pgdspider/"
    echo "此分析參考 James Fifer https://github.com/jamesfifer/JapanRE"
    echo "======================================================="
}

main "$@"
