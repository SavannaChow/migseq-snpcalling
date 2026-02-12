#!/bin/bash
# ==============================================================================
# MIG Analysis Full Pipeline - Smart Hybrid Mode (v8.12)
# ==============================================================================

# --- [新增] 自動偵測 Shell 設定檔路徑 ---
if [[ "$SHELL" == */zsh ]]; then
    CONF_FILE="$HOME/.zshrc"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    CONF_FILE="$HOME/.bash_profile"
else
    CONF_FILE="$HOME/.bashrc"
fi
[ ! -f "$CONF_FILE" ] && touch "$CONF_FILE"
# ---------------------------------------

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
    
    local -A tools=(
        ["fastp"]="conda install -c bioconda fastp"
        ["bwa"]="sudo apt-get install bwa"
        ["samtools"]="sudo apt-get install samtools"
        ["parallel"]="sudo apt-get install parallel"
        ["angsd"]="http://www.popgen.dk/angsd/index.php/Installation"
        ["ngsLD"]="https://github.com/fgvieira/ngsLD"
        ["prune_graph"]="https://github.com/fgvieira/ngsLD (Included in ngsLD)"
        ["bcftools"]="sudo apt-get install bcftools"
        ["Rscript"]="sudo apt-get install r-base"
    )

    local missing_tools=()
    local manual_install=()

    for tool in "${!tools[@]}"; do
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
            local action="${tools[$m_tool]}"
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

check_dependencies

# ------------------------------------------------------------------------------
# 內嵌模組: Genome Fetcher (原 fetch_genome.sh 邏輯)
# ------------------------------------------------------------------------------
run_fetch_genome_module() {
    echo "======================================================="
    echo "啟動 NCBI Genome 下載模組"
    echo "======================================================="

    # 設定基礎下載目錄
    BASE_DIR=~/RefGenome
    mkdir -p "$BASE_DIR"

    # 設定隱藏標記檔路徑
    ENV_CHECK_FILE_FG="$BASE_DIR/.env_verified"
    REQUIRED_CMDS=("esearch" "bwa" "samtools" "curl" "gunzip")

    # 環境偵測邏輯
    if [[ ! -f "$ENV_CHECK_FILE_FG" ]]; then
        echo "Initializing environment check for Genome Fetcher..."

        # 擴展可能的軟體路徑
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
                    return 1 # 改為 return 以免結束主程式
                fi
            fi
        done

        # 建立隱藏標記檔
        touch "$ENV_CHECK_FILE_FG"
        echo "Environment verified and cached."
    fi

    # 進入搜尋與選擇迴圈
    while true; do
        # 1. 使用者輸入關鍵字
        # 不使用 clear 以保留主程式上下文
        read -p "請輸入搜尋關鍵字 (物種名或 BioProject, 輸入 q 離開): " QUERY
        
        if [[ "$QUERY" == "q" || "$QUERY" == "Q" ]]; then
            echo "取消下載。"
            return 0
        fi

        # 2. 執行檢索並暫存結果
        TEMP_DATA=$(mktemp)
        
        TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
        HISTORY_FILE="$BASE_DIR/SearchRecord_${QUERY// /_}_${TIMESTAMP}.txt"

        echo "正在檢索 NCBI 資料庫並解析數據結構..."
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

        # 3. 視覺化格式輸出
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

        # 4. 互動式選擇
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

    # 5.獲取 FTP 路徑並下載
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
    echo "基因組檔案下載完成。"
    read -p "是否繼續設定環境變數並執行解壓縮與索引建置？(y/n): " CONTINUE_PROC
    if [[ "$CONTINUE_PROC" != "y" ]]; then
        echo "程序已終止。檔案保留在 $TARGET_DIR。"
        return 0
    fi

    # 6. 設定環境變數
    clear
    
    # --- 新增功能：顯示剛才下載的基因組資訊 ---
    if [ -n "$SELECTED_LINE" ]; then
        echo "=================================================="
        echo "           [ 已下載基因組詳細資訊 ]"
        echo "--------------------------------------------------"
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
        echo ""
    fi
    # --------------------------------------------------------

    echo "目前系統中已定義的變數名稱與路徑 (Ref_XXX):"
    if grep -q "^export Ref_" "$CONF_FILE"; then
        grep "^export Ref_" "$CONF_FILE" | sed 's/export //g' | sed 's/=/  -->  /g'
    else
        echo "(目前尚無設定任何 Ref_ 變數)"
    fi
    echo "--------------------------------------------------"
    echo "              請輸入參考基因組名稱                   "
    echo "--------------------------------------------------"
    echo "不要跟現有的重複，也不可使用空白或特殊字元"
    echo "嚴禁包含小數點（.）、破折號（-）或其他標點符號"
    echo "只能包含字母（a-z, A-Z）、數字（0-9）以及下底線（_）。"
    echo "--------------------------------------------------"
    read -p "REF_" USER_INPUT

    if [ -z "$USER_INPUT" ]; then
        echo "錯誤：未輸入名稱，停止執行後續解壓縮與索引程序。"
        return 1
    fi
    
    ENV_VAR="Ref_${USER_INPUT}"
    FNA_FILE="$TARGET_DIR/${FILE_NAME}_genomic.fna"
    ABS_PATH=$(realpath "$FNA_FILE")
    echo "自動刷新環境變數"
    exec $SHELL
    # 寫入 .bashrc
    echo "export $ENV_VAR=\"$ABS_PATH\"" >> "$CONF_FILE"
    # 注意：這裡 source 只對當前 shell 有效，主程式需在外部再次 source
    source "$CONF_FILE"
    echo "環境變數 '$ENV_VAR' 已加入 "$CONF_FILE"。"


    # 7. 解壓縮與建立索引
    echo "解壓縮..."
    gunzip -f "$TARGET_DIR/$DOWNLOAD_FILE"

    echo "建立索引 Building BWA index (this may take a while)..."
    echo "bwa index $ABS_PATH"
    bwa index "$ABS_PATH"

    echo "--------------------------------------------------"
    echo "Process complete."
    echo "Genome path: $ABS_PATH"
    echo "處理成功。環境變數已寫入 $CONF_FILE"
    echo "為使新設定生效，將重新載入 Shell。"
    echo "自動刷新環境變數"
    exec $SHELL
}


# ------------------------------------------------------------------------------
# 0. 輔助函式定義
# ------------------------------------------------------------------------------
ask_to_run() {
    local step_name=$1
    local check_target=$2
    local skip_var_name=$3
    local exists=false
    if [ -f "$check_target" ]; then exists=true; elif [ -d "$check_target" ] && [ "$(ls -A "$check_target" 2>/dev/null)" ]; then exists=true; fi
    if [ "$exists" = true ]; then
        echo "-------------------------------------------------------"
        echo "[偵測到已存在的分析結果]: $step_name"
        read -p "是否跳過此步驟，並使用現有結果？ (y:跳過 / n:重新分析): " choice < /dev/tty
        [[ "$choice" == "y" || "$choice" == "Y" ]] && eval "$skip_var_name=true" || eval "$skip_var_name=false"
    else
        eval "$skip_var_name=false"
    fi
}

# ------------------------------------------------------------------------------
# 1. 優先參數輸入與流程選擇
# ------------------------------------------------------------------------------
# 1.1 基本資訊
# ==============================================================================
#   MIG-seq Genomic Analysis Pipeline 2026
#   Biodiversity Research Center, Academia Sinica
# ==============================================================================
clear
echo "==========================================================================="
echo "   MIG-seq Genomic Analysis Pipeline (Refgenome mapping)"
echo "   開發者：Savanna Chow (savanna201@gmail.com) 使用 Gemini 協助開發"
echo "   分析邏輯基於 https://github.com/jamesfifer/JapanRE"
echo "   Credit: AllenChen's lab, Biodiversity Research Center, Academia Sinica"
echo "==========================================================================="
echo "此Script會執行以下六個分析："
echo "  序列修剪 (Fastp Trimming)"
echo "  基因組比對 (BWA Alignment)"
echo "  樣本品質過濾 (PCA Outlier Filtering)"
echo "  複本樣本鑑定 (Clone Identification)"
echo "  連鎖不平衡過濾 (LD Pruning & Site Map Generation)"
echo "  最終變異位點標定 (Final SNP Calling & VCF Output)"
echo "**僅在Ubuntu測試過。路徑與檔案名稱嚴禁空格或特殊字元**"
echo "============================================================================"
echo "                                                       "
echo "                                                       "
read -p "請輸入專案名稱(分析產生的檔案都將以專案名稱為開頭,不要有特殊或空白字元) " PROJECT_NAME
read -e -p "請輸入原始序列 (raw data) 資料夾路徑: " RAW_PATH
[ ! -d "$RAW_PATH" ] && { echo "錯誤：找不到路徑 $RAW_PATH"; exit 1; }
RAW_PATH=$(realpath "$RAW_PATH")

# 1.2 參考基因組配置 (整合下載選單)
# ==============================================================================
while true; do
    # 每次迴圈重新載入 .bashrc 以獲取最新的 Ref_ 變數
    source "$CONF_FILE"
    MAPFILE=()
    MAPVAL=()

    # 過濾環境變數中以 .fa, .fasta, .fna 結尾的路徑
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

    # 分支邏輯
    if [[ "$REF_CHOICE" == "d" || "$REF_CHOICE" == "D" ]]; then
        # 呼叫內嵌的下載模組
        run_fetch_genome_module
        
        echo ""
        echo ">>> 返回主選單，正在刷新參考基因組清單..."
        # 不 break，continue 回到 while 開頭
        continue
    elif [[ "$REF_CHOICE" == "q" || "$REF_CHOICE" == "Q" ]]; then
        echo "使用者取消操作，程式結束。"
        exit 0
    fi

    # 驗證輸入為純數字且在選項範圍內
    if [[ "$REF_CHOICE" =~ ^[0-9]+$ ]] && [ "$REF_CHOICE" -ge 1 ] && [ "$REF_CHOICE" -le "${#MAPFILE[@]}" ]; then
        REF_GENOME="${MAPVAL[$((REF_CHOICE-1))]}"
        break # 選擇成功，跳出迴圈
    elif [ "$REF_CHOICE" -eq "$MANUAL_OPTION" ]; then
        read -e -p "請輸入絕對路徑 (或輸入 'b' 返回選單): " REF_GENOME
        if [[ "$REF_GENOME" == "b" || "$REF_GENOME" == "B" ]]; then
            echo "返回選單..."
            continue
        fi
        break # 選擇成功，跳出迴圈
    else
        echo "錯誤：無效的選擇。"
        # 輸入錯誤，不 break，自動重跑迴圈
    fi
done

# 移除潛在的引號或空白，確保路徑純淨
REF_GENOME=$(echo "$REF_GENOME" | sed "s/['\"]//g")

# 核心檔案存在檢查
if [ -z "$REF_GENOME" ] || [ ! -f "$REF_GENOME" ]; then
    echo "錯誤：檔案不存在於指定路徑：$REF_GENOME"
    exit 1
fi

# 索引狀態檢查與建立
if [ ! -f "${REF_GENOME}.bwt" ]; then
    echo "建立 BWA index..."
    bwa index "$REF_GENOME" || { echo "BWA index 失敗"; exit 1; }
fi

if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "建立 Samtools index..."
    samtools faidx "$REF_GENOME" || { echo "Samtools faidx 失敗"; exit 1; }
fi

echo "使用參考基因組: $REF_GENOME"
# ==============================================================================



# 1.3 模式選擇
echo ""
echo "-------------------------------------------------------"
echo "請選擇分析運行模式："
echo "1) 自動模式 (預先設定篩選策略，遇到決策點不中斷)"
echo "2) 互動模式 (各階段結束後，手動決定是否移除樣本)"
echo ""
read -p "請輸入選項 (1, 2, 或 q 離開): " RUN_MODE

if [[ "$RUN_MODE" == "q" || "$RUN_MODE" == "Q" ]]; then
    echo "使用者取消操作，程式結束。"
    exit 0
fi

if [ "$RUN_MODE" == "1" ]; then
    echo "[模式選擇：自動模式]"
    read -p "  > 偵測到 PCA Outlier 時處理方式 (1:移除, 2:保留, q:離開): " AUTO_PCA_CHOICE
    [[ "$AUTO_PCA_CHOICE" == "q" || "$AUTO_PCA_CHOICE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }

    read -p "  > 偵測到 Clone 樣本時處理方式 (1:移除, 2:保留, q:離開): " AUTO_CLONE_CHOICE
    [[ "$AUTO_CLONE_CHOICE" == "q" || "$AUTO_CLONE_CHOICE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }

    echo "預設 PCA 決策: $AUTO_PCA_CHOICE"
    echo "預設 Clone 決策: $AUTO_CLONE_CHOICE"
else
    echo "[模式選擇：互動模式]"
fi
echo ""

# 1.4 階段選擇 (模組化重構版本)
echo "======================================================="
echo "  1.4 分析流程選擇 (Stage Selector)"
echo "-------------------------------------------------------"
echo ""
echo "1) 執行完整分析 (Stage 1-6)"
echo "2) 僅執行序列清理與比對參考基因組 (Trim & Align)"
echo "3) 僅執行潛在問題樣本過濾分析 (Outlier & Clone Filtering)"
echo "4) 僅執行 LD Pruning (分析連鎖不平衡位點)"
echo "5) 執行最終 SNP Calling (基於現有無LD位點)"
echo "6) 自定義流程"
echo ""
read -p "選擇運行範圍 (或輸入 q 離開): " RUN_CHOICE

if [[ "$RUN_CHOICE" == "q" || "$RUN_CHOICE" == "Q" ]]; then
    echo "使用者取消操作，程式結束。"
    exit 0
fi

echo ""
case "$RUN_CHOICE" in
    1) RUN_S1=y; RUN_S2=y; RUN_S3=y; RUN_S4=y; RUN_S5=y; RUN_S6=y ;;
    2) RUN_S1=y; RUN_S2=y; RUN_S3=n; RUN_S4=n; RUN_S5=n; RUN_S6=n ;;
    3) RUN_S1=n; RUN_S2=n; RUN_S3=y; RUN_S4=y; RUN_S5=n; RUN_S6=n ;;
    4) RUN_S1=n; RUN_S2=n; RUN_S3=n; RUN_S4=n; RUN_S5=y; RUN_S6=n ;;
    5) RUN_S1=n; RUN_S2=n; RUN_S3=n; RUN_S4=n; RUN_S5=n; RUN_S6=y ;;
    *) 
       read -p "執行 Stage 1 Trimming? (y/n, 或 q 離開): " RUN_S1
       [[ "$RUN_S1" == "q" ]] && { echo "使用者取消操作。"; exit 0; }
       read -p "執行 Stage 2 Alignment? (y/n, 或 q 離開): " RUN_S2
       [[ "$RUN_S2" == "q" ]] && { echo "使用者取消操作。"; exit 0; }
       read -p "執行 Stage 3 PCA Outlier Filtering? (y/n, 或 q 離開): " RUN_S3
       [[ "$RUN_S3" == "q" ]] && { echo "使用者取消操作。"; exit 0; }
       read -p "執行 Stage 4 Clone Filtering? (y/n, 或 q 離開): " RUN_S4
       [[ "$RUN_S4" == "q" ]] && { echo "使用者取消操作。"; exit 0; }
       read -p "執行 Stage 5 LD Pruning (Generate Site Map)? (y/n, 或 q 離開): " RUN_S5
       [[ "$RUN_S5" == "q" ]] && { echo "使用者取消操作。"; exit 0; }
       read -p "執行 Stage 6 Final SNP Calling (Apply Map)? (y/n, 或 q 離開): " RUN_S6
       [[ "$RUN_S6" == "q" ]] && { echo "使用者取消操作。"; exit 0; }
       ;;
esac
# ------------------------------------------------------------------------------
# 2. 目錄與日誌初始化
# ------------------------------------------------------------------------------
CURRENT_TIME=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="00_Logs"; STAGE1="01_Trimming"; STAGE2="02_Alignment"; STAGE3="03_PCA_Analysis"
STAGE4="04_Clone_Detection"; STAGE5="05_SNP_Calling"
mkdir -p "$LOG_DIR" "$STAGE1/trim" "$STAGE1/fastp_report" "$STAGE2/bam" "$STAGE2/mapped_bam" "$STAGE2/mapping_results" "$STAGE3" "$STAGE4" "$STAGE5"

LOG_FILE="$LOG_DIR/${PROJECT_NAME}_${CURRENT_TIME}.log"
exec > >(tee -i "$LOG_FILE") 2>&1

THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu)
JOBS=$(( THREADS / 4 )); [ "$JOBS" -lt 1 ] && JOBS=1

# 2. 目錄與日誌初始化後
echo "======================================================="
echo "  分析配置總結"
echo "-------------------------------------------------------"
echo ""
echo "  啟動時間   : $(date)"
echo "  專案名稱   : $PROJECT_NAME"
echo "  執行執行緒 : $THREADS (Jobs: $JOBS)"
echo "  原始路徑   : $RAW_PATH"
echo "  參考基因組 : $REF_GENOME"
echo ""
echo "  日誌檔案   : $LOG_FILE"
echo "======================================================="
echo ""


# ------------------------------------------------------------------------------
# 2.5 執行前最終確認 (Final Confirmation)
# ------------------------------------------------------------------------------
# 轉換模式與流程的顯示字串
clear
[[ "$RUN_MODE" == "1" ]] && MODE_STR="自動模式 (Auto)" || MODE_STR="互動模式 (Manual)"

case "$RUN_CHOICE" in
    1) RANGE_STR="完整分析 (Stage 1-6)" ;;
    2) RANGE_STR="序列清理與比對 (Stage 1-2)" ;;
    3) RANGE_STR="問題樣本過濾 (Stage 3-4)" ;;
    4) RANGE_STR="LD Pruning (Stage 5)" ;;
    5) RANGE_STR="最終 SNP Calling (Stage 6)" ;;
    *) RANGE_STR="自定義流程" ;;
esac
echo "                      執行前最終確認                      "
echo ""
echo "======================================================="
echo "  分析配置總結 (Analysis Configuration Summary)"
echo "-------------------------------------------------------"
echo ""
printf "  %-15s : %s\n" "專案名稱" "$PROJECT_NAME"
printf "  %-15s : %s\n" "原始資料路徑" "$RAW_PATH"
printf "  %-15s : %s\n" "參考基因組" "$REF_GENOME"
printf "  %-15s : %s\n" "執行模式" "$MODE_STR"
printf "  %-15s : %s\n" "運行範圍" "$RANGE_STR"

if [ "$RUN_MODE" == "1" ]; then
    printf "  %-15s : %s\n" "PCA Outlier 決策" "$([[ "$AUTO_PCA_CHOICE" == "1" ]] && echo "移除" || echo "保留")"
    printf "  %-15s : %s\n" "Clone 樣本決策" "$([[ "$AUTO_CLONE_CHOICE" == "1" ]] && echo "移除" || echo "保留")"
fi

printf "  %-15s : %s (Jobs: %s)\n" "並行執行緒" "$THREADS" "$JOBS"
printf "  %-15s : %s\n" "日誌檔案" "$LOG_FILE"
echo ""
echo "-------------------------------------------------------"
read -p "  以上配置是否正確？ (y: 開始執行 / n: 結束並重新設定): " CONFIRM
echo ""

if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
    echo "  [取消] 使用者中止程式，請重新啟動腳本進行設定。"
    # 移除剛建立但尚未使用的日誌檔案與目錄（選擇性）
    exit 0
fi

echo "  [確認] 配置正確，啟動分析流程..."
echo ""



# ------------------------------------------------------------------------------
# 3. 執行分析流程
# ------------------------------------------------------------------------------

# --- [Stage 1: Fastp Trimming] ---
if [[ "$RUN_S1" == "y" ]]; then
    ask_to_run "Fastp Quality Control" "$STAGE1/trim" SKIP_FASTP
    if [[ "$SKIP_FASTP" != true ]]; then
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

    fi
fi

# --- [Stage 2: BWA Alignment] ---
if [[ "$RUN_S2" == "y" ]]; then
    ask_to_run "BWA Alignment" "$STAGE2/mapped_bam" SKIP_BWA
    if [[ "$SKIP_BWA" != true ]]; then
        echo "執行 BWA Mapping..."
        while read -r r1; do
            base=$(basename "$r1" _R1_001.fastq.gz)
            r2="$STAGE1/trim/${base}_R2_001.fastq.gz"
            bwa mem -t "$THREADS" "$REF_GENOME" "$r1" "$r2" > "$STAGE2/bam/${base}.sam"
            samtools view -Sb "$STAGE2/bam/${base}.sam" > "$STAGE2/bam/${base}.bam"
            samtools view -bF4 -@ "$THREADS" "$STAGE2/bam/${base}.bam" > "$STAGE2/mapped_bam/${base}.bam"
            samtools sort -@ "$THREADS" -o "$STAGE2/mapped_bam/${base}_sorted.bam" "$STAGE2/mapped_bam/${base}.bam"
            mv "$STAGE2/mapped_bam/${base}_sorted.bam" "$STAGE2/mapped_bam/${base}.bam"
            samtools index "$STAGE2/mapped_bam/${base}.bam"
            samtools flagstat "$STAGE2/bam/${base}.bam" > "$STAGE2/mapping_results/${base}.txt"
            rm "$STAGE2/bam/${base}.sam"
        done < <(ls "$STAGE1/trim"/*_R1_001.fastq.gz)
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
    fi
fi

# ------------------------------------------------------------------------------
# 6. Stage 3: 生成 Mapping Summary 與 PCA 篩選
# ------------------------------------------------------------------------------
if [[ "$RUN_S3" == "y" ]]; then
    echo "[5/8] 生成比對報表與 PCA 品質檢測..."
    SUMMARY_CSV="$STAGE3/${PROJECT_NAME}_mapping_summary.csv"
    BAM_LIST="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"

    # 確保基礎 BAM LIST 存在 (對應分類資料夾路徑)
    ls -d "$PWD/$STAGE2/mapped_bam/"*.bam > "$BAM_LIST"

    ask_to_run "PCA Analysis & Outlier Detection" "$SUMMARY_CSV" SKIP_PCA

    if [ "$SKIP_PCA" = true ]; then
        echo ">>> 跳過 PCA 分析，嘗試載入先前的過濾結果..."
        if [ -f "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
            BAM_LIST="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
            echo "[!] 載入既有的 PCA 過濾後清單: $(wc -l < "$BAM_LIST") 個樣本"
        else
            echo "[!] 未發現過濾後清單，假設上次選擇保留所有樣本 (或無 Outlier)。"
        fi
    else
        # --- PCA 流程 ---
        echo "Sample,Total,Mapped,Properly_Paired,With_Mate_Mapped,Singletons,Mate_Diff_Chr,Mate_Diff_Chr_MapQ5,Secondary,Supplementary,Duplicates,Paired_in_Seq,Read1,Read2" > "$SUMMARY_CSV"

        for f in "$STAGE2/mapping_results"/*.txt; do
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
            echo "$sample,$total,$mapped,$properly_paired,$with_mate,$singletons,$mate_diff_chr,$mate_diff_chr_q5,$secondary,$supplementary,$duplicates,$paired_in_seq,$read1,$read2" >> "$SUMMARY_CSV"
        done

        echo "-------------------------------------------------------"
        echo "[Stage 3 進度回報]"
        echo "Mapping Summary 資料表已生成: $SUMMARY_CSV"
        echo "資料表頭包含: Sample, Total, Mapped, Properly_Paired, Singletons 等 14 項指標"
        echo "對應統計來源: $STAGE2/mapping_results/"
        echo "-------------------------------------------------------"

# --- 在產生 R script 之前，先將路徑轉為絕對路徑 ---
        ABS_STAGE3=$(realpath "$STAGE3")
        ABS_SUMMARY_CSV=$(realpath "$SUMMARY_CSV")
        ABS_BAM_LIST=$(realpath "$BAM_LIST")

# --- PCA R 腳本 (增強統計資訊與向量顯示版) ---
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

        # 統計過濾結果
        N_ORIG=$(wc -l < "$BAM_LIST")
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
                read -p "是否移除離群樣本？(1:移除, 2:保留, q:離開): " PCA_DECISION < /dev/tty
            fi

            if [[ "$PCA_DECISION" == "q" || "$PCA_DECISION" == "Q" ]]; then
                echo "使用者選擇離開，分析中止。"
                exit 0
            fi

            if [ "$PCA_DECISION" == "1" ]; then
                BAM_LIST="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
                echo "[!] 已套用過濾後的 BAM 清單，當前分析為: $(wc -l < "$BAM_LIST") 個樣本。"
            else
                echo "[+] 已選擇保留outlier樣本，維持原始分析樣本數。"
            fi
        fi
    fi
fi

# ------------------------------------------------------------------------------
# 7. 執行 Clone 偵測 
# ------------------------------------------------------------------------------
if [[ "$RUN_S4" == "y" ]]; then
    # 若上一步沒跑 S3，則需要確保有基礎的 BAM LIST，路徑對接至 Stage 3 的產出
    [ -z "$BAM_LIST" ] && BAM_LIST="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"

    echo "[6/8] 執行 Clone 偵測分析 (ANGSD IBS)..."
    # 呼叫斷點續傳函式，檢查是否存在 IBS 矩陣檔案
    ask_to_run "Clone Detection (IBS Matrix)" "$STAGE4/${PROJECT_NAME}_clone_identification.ibsMat" SKIP_CLONE_CALC

    if [ "$SKIP_CLONE_CALC" = true ]; then
        echo ">>> 跳過 Clone 計算，嘗試載入先前的去重結果..."
        # 檢查是否存在去clone後的 BAM 清單，若有則直接讀取
        if [ -f "$STAGE4/${PROJECT_NAME}_after_clones.bamfile" ]; then
            BAM_LIST="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
            echo "[!] 載入既有的 Clone 去重後清單: $(wc -l < "$BAM_LIST") 個樣本"
        else
            echo "[!] 未發現過濾後清單，假設上次選擇保留 Clone 或無 Clone。"
        fi
    else
        # --- 原本的 Clone 流程 ---
        N_IND_TMP=$(wc -l < "$BAM_LIST")
        # 設定最小樣本覆蓋門檻，此處維持 70% 樣本數 
        MIN_IND_TMP=$(echo "$N_IND_TMP * 0.7 / 1" | bc)

        # 使用 ANGSD 計算 IBS (Identity by State) 矩陣，用於評估樣本遺傳一致性
        echo "使用 ANGSD 計算 IBS (Identity by State) 矩陣"
        angsd -bam "$BAM_LIST" -GL 1 -P 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -minInd "$MIN_IND_TMP" \
              -snp_pval 1e-5 -minMaf 0.05 -doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 \
              -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -out "$STAGE4/${PROJECT_NAME}_clone_identification"

        # 解壓 ANGSD 產出的中間檔以供 R 讀取
        gzip -kfd "$STAGE4"/*.gz 2>/dev/null || true
        # --- 在進入 R 前，將路徑轉換為絕對路徑以確保 RStudio 兼容性 ---
        ABS_STAGE4=$(realpath "$STAGE4")
        ABS_BAM_LIST=$(realpath "$BAM_LIST")
        # --- Clone R 腳本：執行層次聚類與門檻判定 ---
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

        # Clone 統計與決策傳遞
        CLONE_REVIEW_FILE="$STAGE4/${PROJECT_NAME}_clones_to_review.txt"
        N_CLONE_BEFORE=$(wc -l < "$BAM_LIST")
        N_CLONE_AFTER=$(wc -l < "$STAGE4/${PROJECT_NAME}_after_clones.bamfile")
        N_CLONE_COUNT=$(wc -l < "$CLONE_REVIEW_FILE")

        echo "-------------------------------------------------------"
        echo "[Clone樣本偵測報告]"
        echo "1. 原始 PDF (無門檻): $STAGE4/${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf"
        echo "2. 決策 PDF (含紅線): $STAGE4/${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf"
        echo "進入偵測之樣本總數: $N_CLONE_BEFORE"
        echo "剔除Clone樣本後建議樣本數: $N_CLONE_AFTER"
        echo "偵測到 $N_CLONE_COUNT 個潛在Clone樣本"
        echo "-------------------------------------------------------"

        if [ -s "$CLONE_REVIEW_FILE" ]; then
            echo "偵測到以下潛在Clone樣本 (建議移除名單):"
            cat "$CLONE_REVIEW_FILE"
            echo "-------------------------------------------------------"

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

            if [ "$CLONE_DECISION" == "1" ]; then
                BAM_LIST="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
                echo "[!] 已套用去除clones後的最終 BAM 清單，當前分析樣本數: $(wc -l < "$BAM_LIST") 個樣本。"
            else
                echo "[+] 已選擇保留Clone樣本，維持樣本數: $N_CLONE_BEFORE"
            fi
        fi
    fi
fi

# --- [共用邏輯：BAM 清單定位與樣本數計算] ---
# 若執行 S5 或 S6，皆需確認輸入樣本清單與動態門檻
# BAM 清單定位 (區分「完整自動」與「手動自選」)
# 動態掃描與樣本清單選取
# ==============================================================================
if [[ "$RUN_S5" == "y" || "$RUN_S6" == "y" ]]; then
    # 預設清單定義 (用於自動模式或快速參考)
    LIST_FINAL="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
    LIST_PCA_ONLY="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
    LIST_FULL="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"

    if [ "$RUN_CHOICE" == "1" ]; then
        # 完整流程自動模式：優先級判定
        if [ -f "$LIST_FINAL" ]; then BAM_LIST="$LIST_FINAL"
        elif [ -f "$LIST_PCA_ONLY" ]; then BAM_LIST="$LIST_PCA_ONLY"
        else BAM_LIST="$LIST_FULL"
        fi
        echo "[自動流程] 選取當前最完善清單: $(basename "$BAM_LIST")"
    else
        # 手動/自選模式：自動掃描目錄下所有 .bamfile
        echo "-------------------------------------------------------"
        echo "[手動模式] 正在掃描目錄下可用的 .bamfile 清單..."
        
        # 搜尋當前目錄與二層子目錄內所有的 .bamfile 並存入陣列
        mapfile -t FOUND_LISTS < <(find . -maxdepth 3 -name "*.bamfile" | sort)
        
        if [ ${#FOUND_LISTS[@]} -eq 0 ]; then
            echo "警告:未在專案目錄中偵測到任何 .bamfile。"
            read -e -p "請手動輸入自訂清單的完整路徑: " BAM_LIST < /dev/tty
        else
            echo "偵測到以下樣本清單，請選擇欲使用的檔案："
            for i in "${!FOUND_LISTS[@]}"; do
                # 標註哪些是腳本預設產出的檔案，方便識別
                note=""
                [[ "${FOUND_LISTS[$i]}" == *"$LIST_FINAL" ]] && note="(清理過PCA Outlier與Clone的清單)"
                [[ "${FOUND_LISTS[$i]}" == *"$LIST_PCA_ONLY" ]] && note="(僅經過 PCA outlier 過濾清單)"
                [[ "${FOUND_LISTS[$i]}" == *"$LIST_FULL" ]] && note="(預設什麼都沒清理的清單)"
                
                printf "%2d) %s %s\n" "$((i+1))" "${FOUND_LISTS[$i]}" "$note"
            done
            echo " m) 手動輸入其他路徑"
            echo " q) 離開程式"
            
            read -p "請輸入選項 (1-${#FOUND_LISTS[@]}, m, 或 q): " FILE_CHOICE < /dev/tty
            
            if [[ "$FILE_CHOICE" == "q" || "$FILE_CHOICE" == "Q" ]]; then
                echo "使用者取消操作，程式結束。"
                exit 0
            elif [[ "$FILE_CHOICE" =~ ^[0-9]+$ ]] && [ "$FILE_CHOICE" -ge 1 ] && [ "$FILE_CHOICE" -le "${#FOUND_LISTS[@]}" ]; then
                BAM_LIST="${FOUND_LISTS[$((FILE_CHOICE-1))]}"
            elif [[ "$FILE_CHOICE" == "m" || "$FILE_CHOICE" == "M" ]]; then
                read -e -p "請輸入自訂清單路徑: " BAM_LIST < /dev/tty
            else
                echo "錯誤：無效的選擇。程式中止。"
                exit 1
            fi
        fi
    fi

    # 驗證最終選擇的檔案
    if [ ! -s "$BAM_LIST" ]; then
        echo "錯誤：指定的清單檔案「$BAM_LIST」不存在或為空。"
        exit 1
    fi

    # 根據選定的清單內容，動態更新樣本總數與 70% 門檻 (MIN_IND)
    N_IND=$(wc -l < "$BAM_LIST")
    MIN_IND=$(echo "$N_IND * 0.7 / 1" | bc)
    echo "當前分析清單：$BAM_LIST"
    echo "樣本總數：$N_IND，SNP Calling 門檻 (70%)：$MIN_IND"
fi


# --- [Stage 5: LD Pruning & Site Map Generation] ---
if [[ "$RUN_S5" == "y" ]]; then
    echo "[7/8] 執行 LD Pruning 產生連鎖不平衡位點表..."
    ask_to_run "LD Pruning (Site Map)" "$STAGE5/LDpruned_snp.sites" SKIP_S5
    
    if [[ "$SKIP_S5" != true ]]; then
        # 執行初步 SNP Calling 以獲取全位點資訊
        echo "執行初步 SNP Calling 以獲取全位點資訊"
        angsd -b "$BAM_LIST" -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -baq 1 -setMinDepth 5 -SNP_pval 1e-6 -skipTriallelic 1 -doHWE 1 -Hetbias_pval 0.00001 -minInd "$MIN_IND" -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 32 -doCounts 1 -ref "$REF_GENOME" -P 1 -out "$STAGE5/allsnps"
        
        gzip -kfd "$STAGE5"/*.gz
        gunzip -c "$STAGE5/allsnps.mafs.gz" | tail -n +2 | cut -f 1,2 > "$STAGE5/mc1.sites"
        N_SITES=$(wc -l < "$STAGE5/mc1.sites")
        
        # 計算 LD 矩陣與執行位點修剪 (Pruning)
        echo "# 計算 LD 矩陣"
        ngsLD --geno "$STAGE5/allsnps.geno" --verbose 1 --probs 1 --n_ind "$N_IND" --n_sites "$N_SITES" --max_kb_dist 50 --pos "$STAGE5/mc1.sites" --n_threads "$THREADS" --extend_out 1 --out "$STAGE5/allsnpsites.LD"
        echo "# 執行位點修剪 (Pruning)"
        prune_graph --header -v -n "$THREADS" --in "$STAGE5/allsnpsites.LD" --weight-field "r2" --weight-filter "dist <=10000 && r2 >= 0.5" --out "$STAGE5/allsnpsites.pos"
        
        # 轉換格式並建立索引
        echo "# 轉換LD pruned檔案格式給angsd並建立索引"
        sed 's/:/\t/g' "$STAGE5/allsnpsites.pos" | awk '$2!=""' | sort -k1 > "$STAGE5/LDpruned_snp.sites"
        echo "# angsd建立無LD的loci索引"
        angsd sites index "$STAGE5/LDpruned_snp.sites"
        
        echo "[完成] LD Pruned Site Map 已產出: $STAGE5/LDpruned_snp.sites"
        N_SITES_AFTER=$(wc -l < "$STAGE5/LDpruned_snp.sites")
        echo "-------------------------------------------------------"
        echo "[Stage 5 完成回報]"
        echo "LD Pruning 前的初始位點數: $N_SITES"
        echo "LD Pruning 後保留的位點數: $N_SITES_AFTER"
        echo "位點索引檔路徑: $STAGE5/LDpruned_snp.sites"
        echo "-------------------------------------------------------"
    fi
fi

# --- [Stage 6: Final SNP Calling] ---
if [[ "$RUN_S6" == "y" ]]; then
    echo "[8/8] 執行最終 SNP Calling (基於 LD Pruned Sites)..."
    
    # 強制檢查上游位點表是否存在
    if [ ! -f "$STAGE5/LDpruned_snp.sites" ]; then
        echo "錯誤：未發現 LD Pruning 產出的位點表 ($STAGE5/LDpruned_snp.sites)。"
        echo "請先執行 Stage 5 或確保該路徑下檔案完整。"
        exit 1
    fi

    ask_to_run "Final VCF" "$STAGE5/${PROJECT_NAME}_snps_final.vcf" SKIP_S6
    
    if [[ "$SKIP_S6" != true ]]; then
        # 套用 LD-pruned sites 進行高精準度 SNP Calling
        echo "# 套用 LD-pruned sites 進行高精準度 SNP Calling"
        angsd -sites "$STAGE5/LDpruned_snp.sites" -b "$BAM_LIST" -GL 1 -P 1 -minInd "$MIN_IND" -minMapQ 20 -minQ 25 -sb_pval 1e-5 -Hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGlf 2 -dosnpstat 1 -doPost 1 -doGeno 8 -doBcf 1 --ignore-RG 0 -doHWE 1 -ref "$REF_GENOME" -out "$STAGE5/${PROJECT_NAME}_snps_final"
        
        echo "# 轉換 BCF 為 VCF"
        # 轉換 BCF 為 VCF
        bcftools view -O v -o "$STAGE5/${PROJECT_NAME}_snps_final.vcf" "$STAGE5/${PROJECT_NAME}_snps_final.bcf"
        
        FINAL_SNPS=$(bcftools view -H "$STAGE5/${PROJECT_NAME}_snps_final.vcf" | wc -l)
        echo "-------------------------------------------------------"
        echo "[Stage 6 完成]"
        echo "最終產出的 SNP 總數量: $FINAL_SNPS"
        echo "最終 VCF 檔案路徑: $STAGE5/${PROJECT_NAME}_snps_final.vcf"
        echo "-------------------------------------------------------"
        echo "[完成] 最終 SNP Call Set 已產出: $STAGE5/${PROJECT_NAME}_snps_final.vcf"
    fi
fi
echo "======================================================="
echo "分析結束: $(date)"
echo "產出 VCF: $STAGE5/${PROJECT_NAME}_snps_final.vcf"
echo "日誌位置: $LOG_FILE"
echo "請檢查Log檔案是否有潛在分析問題"
echo "======================================================="