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
PROJECT_NAME_FILE=".project_name"

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

STAGE1_BASE="01_Trimming"
STAGE2_BASE="02_Alignment"
STAGE3_BASE="03_PCA_Analysis"
STAGE4_BASE="04_Clone_Detection"
STAGE5_BASE="05_AllSNP_SiteMap"
STAGE6_BASE="06_LD_Pruning"
STAGE7_BASE="07_Final_SNP_Calling"
STAGE8_BASE="08_STRUCTURE_Auto"
STAGE9_BASE="09_Genetic_Divergence"

STAGE1="$STAGE1_BASE"
STAGE2="$STAGE2_BASE"
STAGE3="$STAGE3_BASE"
STAGE4="$STAGE4_BASE"
STAGE5="$STAGE5_BASE"
STAGE6="$STAGE6_BASE"
STAGE7="$STAGE7_BASE"
STAGE8="$STAGE8_BASE"
STAGE9="$STAGE9_BASE"

THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu)
JOBS=$(( THREADS / 4 )); [ "$JOBS" -lt 1 ] && JOBS=1

RUN_S1=n; RUN_S2=n; RUN_S3=n; RUN_S4=n; RUN_S5=n; RUN_S6=n; RUN_S7=n; RUN_S8=n; RUN_S9=n
RUN_S7_WITH_LD="y"
RUN_S7_SKIP_LD="n"
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
ALL_SITES_INPUT=""
BAM_LIST_PCA_INPUT=""
STR_INPUT=""
S7_STR_FILE=""
BAM_LIST_DIV_ALL_INPUT=""
STAGE9_SKIP_LD_VCF_INPUT=""
STAGE9_LAST_RUN_DIR=""
S9_RUN_STATS2="n"

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

build_mapping_summary_csv() {
    local source_dir="$1"
    local out_csv="$2"
    local txt_files=()
    local f sample total mapped properly_paired with_mate singletons
    local mate_diff_chr mate_diff_chr_q5 secondary supplementary duplicates paired_in_seq read1 read2

    txt_files=("$source_dir"/*.txt)
    if [ "${#txt_files[@]}" -eq 0 ] || [ ! -e "${txt_files[0]}" ]; then
        echo "錯誤：找不到 mapping 統計檔案（$source_dir/*.txt）。"
        return 1
    fi

    echo "Sample,Total,Mapped,Properly_Paired,With_Mate_Mapped,Singletons,Mate_Diff_Chr,Mate_Diff_Chr_MapQ5,Secondary,Supplementary,Duplicates,Paired_in_Seq,Read1,Read2" > "$out_csv"

    for f in "${txt_files[@]}"; do
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
        echo "$sample,$total,$mapped,$properly_paired,$with_mate,$singletons,$mate_diff_chr,$mate_diff_chr_q5,$secondary,$supplementary,$duplicates,$paired_in_seq,$read1,$read2" >> "$out_csv"
    done
}

setup_output_dirs() {
    if [[ "$RUN_S1" == "y" ]]; then
        mkdir -p "$STAGE1/trim" "$STAGE1/fastp_report"
    fi

    if [[ "$RUN_S2" == "y" ]]; then
        mkdir -p "$STAGE2/bam" "$STAGE2/mapped_bam" "$STAGE2/mapping_results"
    fi

    if [[ "$RUN_S3" == "y" ]]; then
        mkdir -p "$STAGE3"
    fi

    if [[ "$RUN_S4" == "y" ]]; then
        mkdir -p "$STAGE4"
    fi

    if [[ "$RUN_S5" == "y" ]]; then
        mkdir -p "$STAGE5"
    fi

    if [[ "$RUN_S6" == "y" ]]; then
        mkdir -p "$STAGE6"
    fi

    if [[ "$RUN_S7" == "y" ]]; then
        mkdir -p "$STAGE7"
    fi

    if [[ "$RUN_S8" == "y" ]]; then
        mkdir -p "$STAGE8"
    fi

    if [[ "$RUN_S9" == "y" ]]; then
        mkdir -p "$STAGE9"
    fi
}

configure_stage_paths() {
    if [[ "$CHAIN_STAGES" == true ]]; then
        STAGE1="$STAGE1_BASE"
        STAGE2="$STAGE2_BASE"
        STAGE3="$STAGE3_BASE"
        STAGE4="$STAGE4_BASE"
        STAGE5="$STAGE5_BASE"
        STAGE6="$STAGE6_BASE"
        STAGE7="$STAGE7_BASE"
        STAGE8="$STAGE8_BASE"
        STAGE9="$STAGE9_BASE"
    else
        STAGE1="$STAGE1_BASE"
        STAGE2="$STAGE2_BASE"
        STAGE3="$STAGE3_BASE"
        STAGE4="$STAGE4_BASE"
        STAGE5="$STAGE5_BASE"
        STAGE6="$STAGE6_BASE"
        STAGE7="$STAGE7_BASE"
        STAGE8="$STAGE8_BASE"
        STAGE9="$STAGE9_BASE"
    fi
}

normalize_bamfile_to_absolute() {
    local input_bamfile="$1"
    local output_bamfile="$2"
    local input_dir raw_line abs_line

    input_dir=$(dirname "$input_bamfile")
    : > "$output_bamfile"
    while IFS= read -r raw_line; do
        [ -z "$raw_line" ] && continue
        if [[ "$raw_line" = /* ]]; then
            abs_line="$raw_line"
        else
            abs_line="$input_dir/$raw_line"
        fi
        abs_line=$(realpath "$abs_line")
        echo "$abs_line" >> "$output_bamfile"
    done < "$input_bamfile"
}

check_structure_harvester_ready() {
    local sh_path sh_dir
    sh_path=$(command -v structureHarvester.py 2>/dev/null || true)
    [ -z "$sh_path" ] && return 1
    sh_dir=$(dirname "$sh_path")
    [ -f "$sh_dir/harvesterCore.py" ] || return 1
    return 0
}

ensure_structure_harvester() {
    local install_dir tmp_dir repo_dir
    install_dir="/usr/local/bin"

    if check_structure_harvester_ready; then
        return 0
    fi

    echo "[Stage 7] 偵測不到可用的 structureHarvester.py / harvesterCore.py。"
    echo "可自動執行以下步驟："
    echo "1) git clone https://github.com/dentearl/structureHarvester.git"
    echo "2) 複製 structureHarvester.py 與 harvesterCore.py 到 $install_dir"
    read -p "是否嘗試自動安裝？(y/n): " auto_install_harvester
    if [[ "$auto_install_harvester" != "y" && "$auto_install_harvester" != "Y" ]]; then
        return 1
    fi

    if ! command -v git >/dev/null 2>&1; then
        echo "錯誤：找不到 git，請手動安裝 structureHarvester。"
        return 1
    fi

    tmp_dir=$(mktemp -d)
    repo_dir="$tmp_dir/structureHarvester"
    if ! git clone --depth 1 https://github.com/dentearl/structureHarvester.git "$repo_dir"; then
        echo "錯誤：無法下載 structureHarvester repository。"
        rm -rf "$tmp_dir"
        return 1
    fi

    if [ ! -f "$repo_dir/structureHarvester.py" ] || [ ! -f "$repo_dir/harvesterCore.py" ]; then
        echo "錯誤：下載內容不完整，缺少 structureHarvester.py 或 harvesterCore.py。"
        rm -rf "$tmp_dir"
        return 1
    fi

    if [ -w "$install_dir" ]; then
        cp "$repo_dir/structureHarvester.py" "$install_dir/" || { rm -rf "$tmp_dir"; return 1; }
        cp "$repo_dir/harvesterCore.py" "$install_dir/" || { rm -rf "$tmp_dir"; return 1; }
        chmod 755 "$install_dir/structureHarvester.py" 2>/dev/null || true
    else
        echo "需要管理員權限寫入 $install_dir，將嘗試使用 sudo。"
        sudo cp "$repo_dir/structureHarvester.py" "$install_dir/" || { rm -rf "$tmp_dir"; return 1; }
        sudo cp "$repo_dir/harvesterCore.py" "$install_dir/" || { rm -rf "$tmp_dir"; return 1; }
        sudo chmod 755 "$install_dir/structureHarvester.py" || true
    fi

    rm -rf "$tmp_dir"
    hash -r

    if check_structure_harvester_ready; then
        echo "[完成] structureHarvester 已安裝到 $install_dir。"
        return 0
    fi

    echo "警告：已複製檔案但尚未在 PATH 中偵測到 structureHarvester.py。"
    echo "請確認 $install_dir 在 PATH，或手動執行："
    echo "$install_dir/structureHarvester.py --help"
    return 1
}

select_bamfile_input() {
    local prompt="$1"
    local out_var="$2"
    local found_lists=()
    local choice selected_path

    while true; do
        found_lists=()
        while IFS= read -r f; do
            found_lists+=("$f")
        done < <(find . -maxdepth 3 -type f -name "*.bamfile" | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_lists[@]}" -gt 0 ]; then
            echo "偵測到以下 bamfile："
            for i in "${!found_lists[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_lists[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到任何 bamfile。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_lists[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_lists[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入 bamfile 完整路徑: " selected_path
            [ ! -s "$selected_path" ] && { echo "錯誤：找不到或空檔 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_lists[@]}" ]; then
            selected_path=$(realpath "${found_lists[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

select_stage34_bamfile_input() {
    local prompt="$1"
    local out_var="$2"
    local found_lists=()
    local choice selected_path f

    while true; do
        found_lists=()
        while IFS= read -r f; do
            found_lists+=("$f")
        done < <(find . -maxdepth 5 -type f \( -name "*after_pca.bamfile" -o -name "*after_clones.bamfile" \) | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_lists[@]}" -gt 0 ]; then
            echo "偵測到以下 Stage3/Stage4 bamfile："
            for i in "${!found_lists[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_lists[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到 Stage3/Stage4 輸出的 bamfile。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_lists[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_lists[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入 Stage3/Stage4 bamfile 完整路徑: " selected_path
            [ ! -s "$selected_path" ] && { echo "錯誤：找不到或空檔 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_lists[@]}" ]; then
            selected_path=$(realpath "${found_lists[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

select_trim_dir_input() {
    local prompt="$1"
    local out_var="$2"
    local found_dirs=()
    local choice selected_path

    while true; do
        found_dirs=()
        while IFS= read -r d; do
            found_dirs+=("$d")
        done < <(find . -maxdepth 3 -type d -name "trim" | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_dirs[@]}" -gt 0 ]; then
            echo "偵測到以下 trim 資料夾："
            for i in "${!found_dirs[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_dirs[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到 trim 資料夾。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_dirs[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_dirs[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入 trimmed fastq 資料夾完整路徑: " selected_path
            [ ! -d "$selected_path" ] && { echo "錯誤：找不到路徑 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_dirs[@]}" ]; then
            selected_path=$(realpath "${found_dirs[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

select_directory_input() {
    local prompt="$1"
    local out_var="$2"
    local found_dirs=()
    local choice selected_path

    while true; do
        found_dirs=()
        while IFS= read -r d; do
            found_dirs+=("$d")
        done < <(find . -maxdepth 2 -mindepth 1 -type d | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_dirs[@]}" -gt 0 ]; then
            echo "目前可選資料夾："
            for i in "${!found_dirs[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_dirs[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到可用資料夾。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_dirs[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_dirs[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入資料夾完整路徑: " selected_path
            [ ! -d "$selected_path" ] && { echo "錯誤：找不到路徑 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_dirs[@]}" ]; then
            selected_path=$(realpath "${found_dirs[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

select_ld_sites_input() {
    local prompt="$1"
    local out_var="$2"
    local found_sites=()
    local choice selected_path f

    while true; do
        found_sites=()
        while IFS= read -r f; do
            # 只列出已完成 angsd sites index 的 .sites 檔案
            if [ -f "${f}.idx" ]; then
                found_sites+=("$f")
            fi
        done < <(find . -maxdepth 3 -type f -name "*.sites" | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_sites[@]}" -gt 0 ]; then
            echo "偵測到以下已建立 index 的 sites 檔案："
            for i in "${!found_sites[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_sites[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到已建立 index 的 .sites 檔案。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_sites[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_sites[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入 LD pruned sites 檔案完整路徑: " selected_path
            [ ! -f "$selected_path" ] && { echo "錯誤：找不到檔案 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_sites[@]}" ]; then
            selected_path=$(realpath "${found_sites[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

select_str_input() {
    local prompt="$1"
    local out_var="$2"
    local found_str=()
    local choice selected_path

    while true; do
        found_str=()
        while IFS= read -r f; do
            found_str+=("$f")
        done < <(find . -maxdepth 3 -type f -name "*.str" | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_str[@]}" -gt 0 ]; then
            echo "偵測到以下 .str 檔案："
            for i in "${!found_str[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_str[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到 .str 檔案。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_str[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_str[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入 .str 檔案完整路徑: " selected_path
            [ ! -f "$selected_path" ] && { echo "錯誤：找不到檔案 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_str[@]}" ]; then
            selected_path=$(realpath "${found_str[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

select_vcf_input() {
    local prompt="$1"
    local out_var="$2"
    local found_vcf=()
    local choice selected_path

    while true; do
        found_vcf=()
        while IFS= read -r f; do
            found_vcf+=("$f")
        done < <(find . -maxdepth 4 -type f -name "*.vcf" | sort)

        echo "-------------------------------------------------------"
        echo "$prompt"
        if [ "${#found_vcf[@]}" -gt 0 ]; then
            echo "偵測到以下 .vcf 檔案："
            for i in "${!found_vcf[@]}"; do
                printf "%2d) %s\n" "$((i+1))" "${found_vcf[$i]}"
            done
        else
            echo "目前目錄下尚未偵測到 .vcf 檔案。"
        fi
        echo " r) 重新掃描"
        echo " m) 手動輸入完整路徑"
        echo " q) 離開程式"

        if [ "${#found_vcf[@]}" -gt 0 ]; then
            read -p "請選擇 (1-${#found_vcf[@]}, r, m, q): " choice
        else
            read -p "請選擇 (r, m, q): " choice
        fi

        if [[ "$choice" == "q" || "$choice" == "Q" ]]; then
            echo "使用者取消操作，程式結束。"
            exit 0
        elif [[ "$choice" == "r" || "$choice" == "R" ]]; then
            continue
        elif [[ "$choice" == "m" || "$choice" == "M" ]]; then
            read -e -p "請輸入 .vcf 檔案完整路徑: " selected_path
            [ ! -f "$selected_path" ] && { echo "錯誤：找不到檔案 $selected_path"; continue; }
            selected_path=$(realpath "$selected_path")
            eval "$out_var=\"$selected_path\""
            return 0
        elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "${#found_vcf[@]}" ]; then
            selected_path=$(realpath "${found_vcf[$((choice-1))]}")
            eval "$out_var=\"$selected_path\""
            return 0
        else
            echo "錯誤：無效選擇。"
        fi
    done
}

set_stage7_default_values() {
    local dataset_default="$1"
    local numind_default="$2"
    local numloci_default="$3"
    local detected_cores="${THREADS:-1}"
    [ -z "$numind_default" ] && numind_default=129
    [ -z "$numloci_default" ] && numloci_default=1289
    [[ "$detected_cores" =~ ^[0-9]+$ ]] || detected_cores=1
    [ "$detected_cores" -lt 1 ] && detected_cores=1
    S7_MAXPOPS=10
    S7_BURNIN=50000
    S7_MCMC=500000
    S7_DATASET="$dataset_default"
    S7_KRUNS=10
    S7_NUMIND="$numind_default"
    S7_NUMLOCI="$numloci_default"
    S7_PLOIDY=2
    S7_MISSING="-9"
    S7_ONEROWPERIND=0
    S7_LABEL=True
    S7_POPDATA=False
    S7_POPFLAG=False
    S7_LOCDATA=True
    S7_PHENO=False
    S7_EXTRACOLS=False
    S7_MARKERS=True
    S7_DOMINANT=False
    S7_MAPDIST=False
    S7_PHASE=False
    S7_PHASEINFO=False
    S7_MARKOV=False
    S7_PARALLEL=True
    S7_CORE_DEF=number
    S7_CORES="$detected_cores"
    S7_HARVEST=True
}

infer_stage7_str_metrics() {
    local str_file="$1"
    local numind_var="$2"
    local numloci_var="$3"
    local total_lines data_lines header numind numloci

    [ ! -f "$str_file" ] && return 1
    normalize_structure_str_header "$str_file"

    total_lines=$(wc -l < "$str_file")
    [ -z "$total_lines" ] && total_lines=0
    if [ "$total_lines" -le 0 ]; then
        return 1
    fi

    header=$(head -n1 "$str_file")
    numloci=$(echo "$header" | awk '{print NF}')
    [ -z "$numloci" ] && return 1

    data_lines=$((total_lines - 1))
    if [ "$data_lines" -lt 0 ]; then
        return 1
    fi
    numind=$((data_lines / 2))

    eval "$numind_var=\"$numind\""
    eval "$numloci_var=\"$numloci\""
    return 0
}

normalize_structure_str_header() {
    local str_file="$1"
    local first_line tmp_file

    [ ! -f "$str_file" ] && return 1
    first_line=$(head -n1 "$str_file")
    if [ -n "${first_line//[[:space:]]/}" ]; then
        return 0
    fi

    tmp_file=$(mktemp)
    awk 'NR==1 && $0 ~ /^[[:space:]]*$/ {next} {print}' "$str_file" > "$tmp_file"
    mv "$tmp_file" "$str_file"
    echo "[修正] 偵測到 STR 首行空白，已移除：$str_file"
    return 0
}

sanitize_vcf_sample_ids_inplace() {
    local vcf_file="$1"
    local map_file="$2"
    local tmp_file

    [ ! -f "$vcf_file" ] && return 1
    : > "$map_file"
    tmp_file=$(mktemp)

    awk -v OFS="\t" -v MAP="$map_file" '
        function clean_name(raw, base) {
            base = raw
            gsub(/^.*\//, "", base)
            sub(/\.(bam|cram|sam)$/, "", base)
            gsub(/[^A-Za-z0-9_-]/, "_", base)
            gsub(/^_+/, "", base)
            gsub(/_+$/, "", base)
            if (base == "") base = "sample"
            return base
        }

        /^#CHROM\t/ {
            for (i = 10; i <= NF; i++) {
                old = $i
                new = clean_name(old)
                print old "\t" new >> MAP
                $i = new
            }
            print
            next
        }
        { print }
    ' "$vcf_file" > "$tmp_file"

    mv "$tmp_file" "$vcf_file"
    return 0
}

prompt_stage7_cores() {
    local value detected_cores
    detected_cores="${THREADS:-1}"
    [[ "$detected_cores" =~ ^[0-9]+$ ]] || detected_cores=1
    [ "$detected_cores" -lt 1 ] && detected_cores=1

    while true; do
        read -p "25. cores [$S7_CORES] (偵測本機核心: $detected_cores): " value
        [ -z "$value" ] && value="$S7_CORES"
        if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -ge 1 ]; then
            if [ "$value" -gt "$detected_cores" ]; then
                echo "警告：輸入值大於本機核心數 $detected_cores。"
                read -p "仍要使用 $value 嗎？(y/n): " confirm_high_core
                if [[ "$confirm_high_core" != "y" && "$confirm_high_core" != "Y" ]]; then
                    continue
                fi
            fi
            S7_CORES="$value"
            return 0
        fi
        echo "錯誤：請輸入 >=1 的整數。"
    done
}

prompt_stage7_dataset_from_str() {
    local choice tmp_str tmp_numind tmp_numloci

    while true; do
        echo "4. dataset 與 .str 輸入來源 (dataset 固定等於 .str 檔名)"
        echo "  目前 .str: $S7_STR_FILE"
        echo "  目前 dataset: $S7_DATASET"
        echo "1) 使用目前 .str"
        echo "2) 掃描並選擇其他 .str"
        echo "3) 手動輸入 .str 完整路徑"
        echo "q) 返回"
        read -p "請選擇: " choice

        case "$choice" in
            1)
                if infer_stage7_str_metrics "$S7_STR_FILE" tmp_numind tmp_numloci; then
                    S7_DATASET=$(basename "$S7_STR_FILE" .str)
                    S7_NUMIND="$tmp_numind"
                    S7_NUMLOCI="$tmp_numloci"
                    echo "已帶入: dataset=$S7_DATASET, numind=$S7_NUMIND, numloci=$S7_NUMLOCI"
                else
                    echo "警告：無法由目前 .str 自動推算 numind/numloci。"
                fi
                return 0
                ;;
            2)
                select_str_input "請選擇 Stage7 要使用的 STRUCTURE .str 檔案" tmp_str
                S7_STR_FILE="$tmp_str"
                if infer_stage7_str_metrics "$S7_STR_FILE" tmp_numind tmp_numloci; then
                    S7_DATASET=$(basename "$S7_STR_FILE" .str)
                    S7_NUMIND="$tmp_numind"
                    S7_NUMLOCI="$tmp_numloci"
                    echo "已帶入: dataset=$S7_DATASET, numind=$S7_NUMIND, numloci=$S7_NUMLOCI"
                else
                    S7_DATASET=$(basename "$S7_STR_FILE" .str)
                    echo "警告：無法由新 .str 自動推算 numind/numloci，僅更新 dataset=$S7_DATASET"
                fi
                return 0
                ;;
            3)
                read -e -p "請輸入 .str 檔案完整路徑: " tmp_str
                if [ ! -f "$tmp_str" ]; then
                    echo "錯誤：找不到檔案 $tmp_str"
                    continue
                fi
                S7_STR_FILE=$(realpath "$tmp_str")
                if infer_stage7_str_metrics "$S7_STR_FILE" tmp_numind tmp_numloci; then
                    S7_DATASET=$(basename "$S7_STR_FILE" .str)
                    S7_NUMIND="$tmp_numind"
                    S7_NUMLOCI="$tmp_numloci"
                    echo "已帶入: dataset=$S7_DATASET, numind=$S7_NUMIND, numloci=$S7_NUMLOCI"
                else
                    S7_DATASET=$(basename "$S7_STR_FILE" .str)
                    echo "警告：無法由新 .str 自動推算 numind/numloci，僅更新 dataset=$S7_DATASET"
                fi
                return 0
                ;;
            q|Q)
                return 0
                ;;
            *)
                echo "錯誤：無效選擇。"
                ;;
        esac
    done
}

save_stage7_params() {
    local out_file="$1"
    {
        printf 'S7_MAXPOPS=%q\n' "$S7_MAXPOPS"
        printf 'S7_BURNIN=%q\n' "$S7_BURNIN"
        printf 'S7_MCMC=%q\n' "$S7_MCMC"
        printf 'S7_STR_FILE=%q\n' "$S7_STR_FILE"
        printf 'S7_DATASET=%q\n' "$S7_DATASET"
        printf 'S7_KRUNS=%q\n' "$S7_KRUNS"
        printf 'S7_NUMIND=%q\n' "$S7_NUMIND"
        printf 'S7_NUMLOCI=%q\n' "$S7_NUMLOCI"
        printf 'S7_PLOIDY=%q\n' "$S7_PLOIDY"
        printf 'S7_MISSING=%q\n' "$S7_MISSING"
        printf 'S7_ONEROWPERIND=%q\n' "$S7_ONEROWPERIND"
        printf 'S7_LABEL=%q\n' "$S7_LABEL"
        printf 'S7_POPDATA=%q\n' "$S7_POPDATA"
        printf 'S7_POPFLAG=%q\n' "$S7_POPFLAG"
        printf 'S7_LOCDATA=%q\n' "$S7_LOCDATA"
        printf 'S7_PHENO=%q\n' "$S7_PHENO"
        printf 'S7_EXTRACOLS=%q\n' "$S7_EXTRACOLS"
        printf 'S7_MARKERS=%q\n' "$S7_MARKERS"
        printf 'S7_DOMINANT=%q\n' "$S7_DOMINANT"
        printf 'S7_MAPDIST=%q\n' "$S7_MAPDIST"
        printf 'S7_PHASE=%q\n' "$S7_PHASE"
        printf 'S7_PHASEINFO=%q\n' "$S7_PHASEINFO"
        printf 'S7_MARKOV=%q\n' "$S7_MARKOV"
        printf 'S7_PARALLEL=%q\n' "$S7_PARALLEL"
        printf 'S7_CORE_DEF=%q\n' "$S7_CORE_DEF"
        printf 'S7_CORES=%q\n' "$S7_CORES"
        printf 'S7_HARVEST=%q\n' "$S7_HARVEST"
    } > "$out_file"
}

prompt_stage7_int() {
    local var="$1" prompt="$2" default="$3" value
    while true; do
        read -p "$prompt [$default]: " value
        [ -z "$value" ] && value="$default"
        if [[ "$value" =~ ^[0-9]+$ ]]; then
            eval "$var=\"$value\""
            return 0
        fi
        echo "錯誤：請輸入整數。"
    done
}

prompt_stage7_str() {
    local var="$1" prompt="$2" default="$3" value
    read -p "$prompt [$default]: " value
    [ -z "$value" ] && value="$default"
    eval "$var=\"$value\""
}

prompt_stage7_bool() {
    local var="$1" prompt="$2" default="$3" value norm
    while true; do
        read -p "$prompt [$default]: " value
        [ -z "$value" ] && value="$default"
        norm=$(echo "$value" | tr '[:upper:]' '[:lower:]')
        case "$norm" in
            true|t|1|y|yes) eval "$var=True"; return 0 ;;
            false|f|0|n|no) eval "$var=False"; return 0 ;;
            *) echo "錯誤：請輸入 True/False (或 y/n)。" ;;
        esac
    done
}

stage7_bool_is_true() {
    local value norm
    value="$1"
    norm=$(echo "$value" | tr '[:upper:]' '[:lower:]')
    case "$norm" in
        true|t|1|y|yes) return 0 ;;
        *) return 1 ;;
    esac
}

stage7_bool_to_int() {
    if stage7_bool_is_true "$1"; then
        echo 1
    else
        echo 0
    fi
}

prompt_stage7_core_def() {
    local value
    while true; do
        read -p "24. core_def (number/percent) [$S7_CORE_DEF]: " value
        [ -z "$value" ] && value="$S7_CORE_DEF"
        case "$value" in
            number|percent)
                S7_CORE_DEF="$value"
                return 0
                ;;
            *)
                echo "錯誤：請輸入 number 或 percent。"
                ;;
        esac
    done
}

show_stage7_params() {
    echo "---------------- Stage7 參數 ----------------"
    echo " 1) maxpops=$S7_MAXPOPS (最大K值)"
    echo " 2) burnin=$S7_BURNIN"
    echo " 3) mcmc=$S7_MCMC"
    echo " 4) dataset=$S7_DATASET (.str: $S7_STR_FILE)"
    echo " 5) kruns=$S7_KRUNS (每個K重複次數)"
    echo " 6) numind=$S7_NUMIND"
    echo " 7) numloci=$S7_NUMLOCI"
    echo " 8) ploidy=$S7_PLOIDY"
    echo " 9) missing=$S7_MISSING"
    echo "10) onerowperind=$S7_ONEROWPERIND (每個個體一行或兩行)"
    echo "11) label=$S7_LABEL (樣本ID欄位)"
    echo "12) popdata=$S7_POPDATA (是否提供既有族群編號欄位)"
    echo "13) popflag=$S7_POPFLAG (USEPOPINFO時, 是否提供強制指定族群欄位)"
    echo "14) locdata=$S7_LOCDATA (是否提供location ID欄位; 對應LOCPRIOR)"
    echo "15) pheno=$S7_PHENO"
    echo "16) extracols=$S7_EXTRACOLS"
    echo "17) markers=$S7_MARKERS"
    echo "18) dominant=$S7_DOMINANT"
    echo "19) mapdist=$S7_MAPDIST"
    echo "20) phase=$S7_PHASE"
    echo "21) phaseinfo=$S7_PHASEINFO"
    echo "22) markov=$S7_MARKOV (是否使用linkage model)"
    echo "23) parallel=$S7_PARALLEL"
    echo "24) core_def=$S7_CORE_DEF"
    echo "25) cores=$S7_CORES"
    echo "26) harvest=$S7_HARVEST"
    echo "---------------------------------------------"
}

print_stage7_parameter_notes() {
    echo "---------------- Stage7 參數註解 ----------------"
    echo "[LABEL]"
    echo "  - LABEL 代表資料檔是否包含樣本ID欄位。"
    echo ""
    echo "[POPDATA]"
    echo "  - 在資料檔中提供「事先定義的族群編號」。放在每個個體基因型前面的一欄，通常在 sample_ID 之後。"
    echo "  - popdata=1: 想檢驗既有族群是否合理"
    echo "  - popdata=0: 純探索式 STRUCTURE」（無先驗族群）"
    echo ""
    echo "[POPFLAG]"
    echo "  - popflag 用於 USEPOPINFO 模式，控制個體是否強制指定族群。"
    echo "  - 啟用時需在 popdata 後再加一欄。"
    echo ""
    echo "[LOCPRIOR / LOCDATA]"
    echo "  - LOCPRIOR 用於族群訊號弱時，加入 sampling location 的輕量先驗。"
    echo "  - 它不是把地理位置當成強制分群，而只是提供機率偏好"
    echo "  - location ID 放在 sample_ID 後；若有 popdata，則放在 popdata 後。"
    echo "  - location ID 必須是整數(1,2,3...)，同地點個體用同數字。"
    echo "  - diploid 個體的兩行資料，location ID 必須完全一致。"
    echo ""
    echo "[LOCDATA 補充]"
    echo "  - locdata 也可提供 loci 額外資訊(如染色體或 map position)。"
    echo "  - 此類資訊只在 linkage model 下需要。"
    echo ""
    echo "[ONEROWPERIND]"
    echo "  - onerowperind=0: 每個 diploid 個體以兩行呈現，同一個體的兩條同源染色體分別佔一行，每一行對應一套等位基因"
    echo "  - onerowperind=1: 每個 diploid 個體以一行呈現，每個 marker 需在同一行中連續給出兩個等位基因值作為一對 genotype"
    echo ""
    echo "[MARKOV]"
    echo "  - markov 表示是否使用 linkage model。"
    echo "  - 非連鎖模型通常關閉。"
    echo ""
    echo "[MAXPOPS / KRUNS]"
    echo "  - maxpops: 最大 K 值(必填)，決定要測試到幾個族群假說。"
    echo "  - kruns: 每個 K 重複跑幾次，建議 >=10 以評估穩定性。"
    echo "-------------------------------------------------"
}

prompt_stage7_param_by_index() {
    local idx="$1"
    case "$idx" in
        1) prompt_stage7_int S7_MAXPOPS "1. maxpops" "$S7_MAXPOPS" ;;
        2) prompt_stage7_int S7_BURNIN "2. burnin" "$S7_BURNIN" ;;
        3) prompt_stage7_int S7_MCMC "3. mcmc" "$S7_MCMC" ;;
        4) prompt_stage7_dataset_from_str ;;
        5) prompt_stage7_int S7_KRUNS "5. kruns" "$S7_KRUNS" ;;
        6) prompt_stage7_int S7_NUMIND "6. numind" "$S7_NUMIND" ;;
        7) prompt_stage7_int S7_NUMLOCI "7. numloci" "$S7_NUMLOCI" ;;
        8) prompt_stage7_int S7_PLOIDY "8. ploidy" "$S7_PLOIDY" ;;
        9) prompt_stage7_str S7_MISSING "9. missing code" "$S7_MISSING" ;;
        10) prompt_stage7_int S7_ONEROWPERIND "10. onerowperind (0/1)" "$S7_ONEROWPERIND" ;;
        11) prompt_stage7_bool S7_LABEL "11. label (True/False)" "$S7_LABEL" ;;
        12) prompt_stage7_bool S7_POPDATA "12. popdata (True/False)" "$S7_POPDATA" ;;
        13) prompt_stage7_bool S7_POPFLAG "13. popflag (True/False)" "$S7_POPFLAG" ;;
        14) prompt_stage7_bool S7_LOCDATA "14. locdata (True/False)" "$S7_LOCDATA" ;;
        15) prompt_stage7_bool S7_PHENO "15. pheno (True/False)" "$S7_PHENO" ;;
        16) prompt_stage7_bool S7_EXTRACOLS "16. extracols (True/False)" "$S7_EXTRACOLS" ;;
        17) prompt_stage7_bool S7_MARKERS "17. markers (True/False)" "$S7_MARKERS" ;;
        18) prompt_stage7_bool S7_DOMINANT "18. dominant (True/False)" "$S7_DOMINANT" ;;
        19) prompt_stage7_bool S7_MAPDIST "19. mapdist (True/False)" "$S7_MAPDIST" ;;
        20) prompt_stage7_bool S7_PHASE "20. phase (True/False)" "$S7_PHASE" ;;
        21) prompt_stage7_bool S7_PHASEINFO "21. phaseinfo (True/False)" "$S7_PHASEINFO" ;;
        22) prompt_stage7_bool S7_MARKOV "22. markov (True/False)" "$S7_MARKOV" ;;
        23) prompt_stage7_bool S7_PARALLEL "23. parallel (True/False)" "$S7_PARALLEL" ;;
        24) prompt_stage7_core_def ;;
        25) prompt_stage7_cores ;;
        26) prompt_stage7_bool S7_HARVEST "26. harvest (True/False)" "$S7_HARVEST" ;;
        *) echo "無效參數編號: $idx" ;;
    esac
}

prompt_stage7_all() {
    local i
    for i in $(seq 1 26); do
        prompt_stage7_param_by_index "$i"
    done
}

prompt_stage7_selected() {
    local choice idx
    show_stage7_params
    read -p "請輸入要修改的編號（可空白/逗號分隔，輸入 all 代表全部）: " choice
    if [[ "$choice" == "all" || "$choice" == "ALL" || -z "$choice" ]]; then
        prompt_stage7_all
        return
    fi

    choice=${choice//,/ }
    for idx in $choice; do
        if [[ "$idx" =~ ^[0-9]+$ ]] && [ "$idx" -ge 1 ] && [ "$idx" -le 26 ]; then
            prompt_stage7_param_by_index "$idx"
        else
            echo "跳過無效編號: $idx"
        fi
    done
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
    echo "  流程串接   : $([[ "$CHAIN_STAGES" == true ]] && echo "啟用" || echo "停用")"
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
    echo "1) 只跑 Stage 1: Fastp Trimming"
    echo "2) 只跑 Stage 2: BWA Alignment"
    echo "3) 只跑 Stage 3: PCA Outlier 分析"
    echo "4) 只跑 Stage 4: Clone Detection"
    echo "5) 只跑 Stage 5: All SNP Site Map"
    echo "6) 只跑 Stage 6: LD Pruned SNP Site Map"
    echo "7) 只跑 Stage 7: Final SNP Calling (with LD pruning / skip LD pruning)"
    echo "8) 只跑 Stage 8: Structure Auto Generator"
    echo "9) 只跑 Stage 9: Analysis of Genetic Divergence"
    echo "10) 自定義多階段 (不自動串接)"
    echo "q) 離開"
    echo ""
    read -p "選擇分析範圍: " RUN_SCOPE

    if [[ "$RUN_SCOPE" == "q" || "$RUN_SCOPE" == "Q" ]]; then
        echo "使用者取消操作，程式結束。"
        exit 0
    fi

    case "$RUN_SCOPE" in
        1) RUN_S1=y ;;
        2) RUN_S2=y ;;
        3) RUN_S3=y ;;
        4) RUN_S4=y ;;
        5) RUN_S5=y ;;
        6) RUN_S6=y ;;
        7) RUN_S7=y ;;
        8) RUN_S8=y ;;
        9) RUN_S9=y ;;
        10)
            read -p "執行 Stage 1 Trimming? (y/n): " RUN_S1
            read -p "執行 Stage 2 Alignment? (y/n): " RUN_S2
            read -p "執行 Stage 3 PCA Outlier Filtering? (y/n): " RUN_S3
            read -p "執行 Stage 4 Clone Filtering? (y/n): " RUN_S4
            read -p "執行 Stage 5 All SNP Site Map? (y/n): " RUN_S5
            read -p "執行 Stage 6 LD Pruned SNP Site Map? (y/n): " RUN_S6
            read -p "執行 Stage 7 Final SNP Calling? (y/n): " RUN_S7
            read -p "執行 Stage 8 Structure Auto Generator? (y/n): " RUN_S8
            read -p "執行 Stage 9 Genetic Divergence? (y/n): " RUN_S9
            CHAIN_STAGES=false
            ;;
        *)
            echo "錯誤：無效選項。"
            exit 1
            ;;
    esac
}

configure_project_name() {
    local stored_name
    if [ -s "$PROJECT_NAME_FILE" ]; then
        stored_name=$(head -n1 "$PROJECT_NAME_FILE" | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
        if [[ "$stored_name" =~ ^[A-Za-z0-9._-]+$ ]]; then
            BASE_PROJECT_NAME="$stored_name"
            PROJECT_NAME="$stored_name"
            echo "專案名稱已由 $PROJECT_NAME_FILE 載入：$PROJECT_NAME"
            return
        fi
        echo "警告：$PROJECT_NAME_FILE 內容格式不正確，將重新設定專案名稱。"
    fi

    while true; do
        echo "首次執行：尚未偵測到專案名稱記錄檔 ($PROJECT_NAME_FILE)"
        echo "請輸入專案名稱（英數字、點、底線、連字號）"
        read -p "請輸入: " BASE_PROJECT_NAME
        if [[ -z "$BASE_PROJECT_NAME" ]]; then
            echo "錯誤：專案名稱不可為空。"
            continue
        fi
        if [[ ! "$BASE_PROJECT_NAME" =~ ^[A-Za-z0-9._-]+$ ]]; then
            echo "錯誤：專案名稱僅允許英數字、點(.)、底線(_)與連字號(-)。"
            continue
        fi
        PROJECT_NAME="$BASE_PROJECT_NAME"
        printf "%s\n" "$PROJECT_NAME" > "$PROJECT_NAME_FILE"
        echo "已建立 $PROJECT_NAME_FILE，專案名稱：$PROJECT_NAME"
        break
    done
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
        echo "$MANUAL_OPTION) 手動輸入絕對路徑"
        echo "q) 離開程式"
        echo "------------------------"

        read -p "請選擇參考基因組 (1-$MANUAL_OPTION, 或 q): " REF_CHOICE

        if [[ "$REF_CHOICE" == "q" || "$REF_CHOICE" == "Q" ]]; then
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
            if [[ "$RUN_S3" == "y" ]]; then
                read -p "  > 偵測到 PCA Outlier 時處理方式 (1:移除, 2:保留, q:離開): " AUTO_PCA_CHOICE
                [[ "$AUTO_PCA_CHOICE" == "q" || "$AUTO_PCA_CHOICE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }
            fi
            if [[ "$RUN_S4" == "y" ]]; then
                read -p "  > 偵測到 Clone 樣本時處理方式 (1:移除, 2:保留, q:離開): " AUTO_CLONE_CHOICE
                [[ "$AUTO_CLONE_CHOICE" == "q" || "$AUTO_CLONE_CHOICE" == "Q" ]] && { echo "使用者取消操作。"; exit 0; }
            fi
        fi
    fi

    if [[ "$RUN_S1" == "y" ]]; then
        select_directory_input "請選擇 Stage1 原始序列 (raw data) 資料夾路徑" RAW_PATH
    fi

    if [[ "$RUN_S2" == "y" ]]; then
        if [[ "$CHAIN_STAGES" == true || "$RUN_S1" == "y" ]]; then
            TRIM_INPUT_DIR="$STAGE1/trim"
        else
            select_trim_dir_input "請選擇 Stage2 要用的 trimmed fastq 資料夾路徑" TRIM_INPUT_DIR
        fi
    fi

    if [[ "$RUN_S2" == "y" || "$RUN_S5" == "y" || "$RUN_S6" == "y" || "$RUN_S7" == "y" || "$RUN_S9" == "y" ]]; then
        select_ref_genome
    fi

    if [[ "$RUN_S3" == "y" ]]; then
        if [[ "$CHAIN_STAGES" == true || "$RUN_S2" == "y" ]]; then
            :
        else
            select_bamfile_input "請選擇 Stage3 要使用的 BAM list (.bamfile)" BAM_LIST_PCA_INPUT
        fi
    fi

    if [[ "$RUN_S4" == "y" ]]; then
        if [[ "$CHAIN_STAGES" == true || "$RUN_S3" == "y" || "$RUN_S2" == "y" ]]; then
            :
        else
            select_bamfile_input "請選擇 Stage4 clone detection 要使用的 BAM list (.bamfile)" BAM_LIST_CLONE_INPUT
        fi
    fi

    if [[ "$RUN_S5" == "y" && "$CHAIN_STAGES" != true && "$RUN_S4" != "y" && "$RUN_S3" != "y" && "$RUN_S2" != "y" ]]; then
        select_bamfile_input "請選擇 Stage5 All SNP sites 要使用的 BAM list (.bamfile)" BAM_LIST_LD_INPUT
    fi

    if [[ "$RUN_S6" == "y" && "$CHAIN_STAGES" != true ]]; then
        if [[ "$RUN_S5" != "y" && "$RUN_S4" != "y" && "$RUN_S3" != "y" && "$RUN_S2" != "y" ]]; then
            select_bamfile_input "請選擇 Stage6 LD pruning 要使用的 BAM list (.bamfile)" BAM_LIST_FINAL_INPUT
        fi
        if [[ "$RUN_S5" != "y" ]]; then
            select_ld_sites_input "請選擇 Stage6 要使用的 All SNP sites 檔案" ALL_SITES_INPUT
        fi
    fi

    if [[ "$RUN_S7" == "y" && "$CHAIN_STAGES" != true && "$RUN_S6" != "y" && "$RUN_S5" != "y" && "$RUN_S4" != "y" && "$RUN_S3" != "y" && "$RUN_S2" != "y" ]]; then
        select_bamfile_input "請選擇 Stage7 Final SNP 要使用的 BAM list (.bamfile)" BAM_LIST_FINAL_INPUT
    fi

    if [[ "$RUN_S7" == "y" ]]; then
        read -p "Stage7 是否要做 LD pruning（產生無 linkage SNP）? (y/n) [y]: " RUN_S7_WITH_LD
        read -p "Stage7 是否要跳過 LD pruning（保留潛在 linkage SNP）? (y/n) [n]: " RUN_S7_SKIP_LD
        [ -z "$RUN_S7_WITH_LD" ] && RUN_S7_WITH_LD="y"
        [ -z "$RUN_S7_SKIP_LD" ] && RUN_S7_SKIP_LD="n"
        if [[ "$RUN_S7_WITH_LD" != "y" && "$RUN_S7_SKIP_LD" != "y" ]]; then
            echo "錯誤：Stage7 至少需選擇一種輸出模式。"
            exit 1
        fi
        if [[ "$RUN_S7_WITH_LD" == "y" && "$RUN_S6" != "y" ]]; then
            select_ld_sites_input "請選擇 Stage7(with LD pruning) 要使用的 LD pruned sites 檔案" LD_SITES_INPUT
        fi
        if [[ "$RUN_S7_SKIP_LD" == "y" && "$RUN_S5" != "y" ]]; then
            select_ld_sites_input "請選擇 Stage7(skip LD pruning) 要使用的 All SNP sites 檔案" ALL_SITES_INPUT
        fi
    fi

    if [[ "$RUN_S8" == "y" ]]; then
        if [[ "$RUN_S7" != "y" || "$RUN_S7_WITH_LD" != "y" ]]; then
            select_str_input "請選擇 Stage8 要使用的 STRUCTURE .str 檔案（建議 Stage7 with LD pruning 輸出）" STR_INPUT
        fi
    fi

    if [[ "$RUN_S9" == "y" ]]; then
        echo "-------------------------------------------------------"
        echo "[Stage 9 提醒]"
        echo "1) 你需要一個『所有族群』bamfile 來建立 AllSites。"
        echo "2) 你需要多個 population.bamfile 做 divergence 比較。"
        echo "3) 你需要參考基因組。"
        echo ""
        echo "Stage9 Fst 統計選項說明："
        echo "- stats : 輸出整體 FST.Unweight 與 FST.Weight（建議必跑）"
        echo "- stats2: 輸出視窗/區段層級統計（檔案較大、耗時較長）"
        if [[ "$RUN_S7" != "y" ]]; then
            read -p "是否指定 Skip LD-pruning VCF 路徑供 Stage9 紀錄？(y/n) [n]: " use_stage9_vcf
            [ -z "$use_stage9_vcf" ] && use_stage9_vcf="n"
            if [[ "$use_stage9_vcf" == "y" || "$use_stage9_vcf" == "Y" ]]; then
                select_vcf_input "請選擇 Stage9 對應的 Skip LD-pruning VCF（僅作紀錄）" STAGE9_SKIP_LD_VCF_INPUT
            fi
        fi
        read -p "是否執行 stats2（較耗時）? (y/n) [n]: " S9_RUN_STATS2
        [ -z "$S9_RUN_STATS2" ] && S9_RUN_STATS2="n"
        if [[ "$RUN_S4" != "y" && "$RUN_S3" != "y" ]]; then
            select_stage34_bamfile_input "請選擇 Stage9 要使用的『所有族群』BAM list (.bamfile)" BAM_LIST_DIV_ALL_INPUT
        fi
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
    [[ "$RUN_S5" == "y" ]] && echo "    - Stage 5 All SNP Site Map"
    [[ "$RUN_S6" == "y" ]] && echo "    - Stage 6 LD Pruned SNP Site Map"
    [[ "$RUN_S7" == "y" ]] && echo "    - Stage 7 Final SNP Calling (with LD pruning / skip LD pruning)"
    [[ "$RUN_S8" == "y" ]] && echo "    - Stage 8 Structure Auto Generator"
    [[ "$RUN_S9" == "y" ]] && echo "    - Stage 9 Analysis of Genetic Divergence"
    [[ "$RUN_S7" == "y" ]] && printf "  %-15s : %s\n" "Stage7 做LD pruning" "$RUN_S7_WITH_LD"
    [[ "$RUN_S7" == "y" ]] && printf "  %-15s : %s\n" "Stage7 跳過LD pruning" "$RUN_S7_SKIP_LD"
    [[ "$RUN_S9" == "y" ]] && printf "  %-15s : %s\n" "Stage9 全族群BAM" "$BAM_LIST_DIV_ALL_INPUT"
    [[ "$RUN_S9" == "y" ]] && printf "  %-15s : %s\n" "Stage9 參考No-LD" "$STAGE9_SKIP_LD_VCF_INPUT"
    [[ "$RUN_S9" == "y" ]] && printf "  %-15s : %s\n" "Stage9 跑stats2" "$S9_RUN_STATS2"

    if [[ "$RUN_MODE" == "1" && "$RUN_S3" == "y" ]]; then
        printf "  %-15s : %s\n" "PCA Outlier 決策" "$([[ "$AUTO_PCA_CHOICE" == "1" ]] && echo "移除" || echo "保留")"
    fi
    if [[ "$RUN_MODE" == "1" && "$RUN_S4" == "y" ]]; then
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

    echo "執行 Fastp (增量模式：只處理新增樣本)..."
    local raw_list_all raw_list_todo n_all n_done n_todo
    raw_list_all="$STAGE1/raw_list_all.txt"
    raw_list_todo="$STAGE1/raw_list_todo.txt"
    : > "$raw_list_all"
    : > "$raw_list_todo"

    ls "${RAW_PATH}"/*_R1_001.fast* > "$raw_list_all"

    while read -r r1; do
        base=$(basename "$r1" | sed 's/_R1_001\.fastq\.gz$//; s/_R1_001\.fastq$//; s/_R1_001\.fq\.gz$//; s/_R1_001\.fq$//')
        out_r1="$STAGE1/trim/${base}_R1_001.fastq.gz"
        out_r2="$STAGE1/trim/${base}_R2_001.fastq.gz"

        # 已有完整 trim 輸出就跳過；缺任何一個則視為待處理
        if [[ ! -s "$out_r1" || ! -s "$out_r2" ]]; then
            echo "$r1" >> "$raw_list_todo"
        fi
    done < "$raw_list_all"

    n_all=$(wc -l < "$raw_list_all")
    n_todo=$(wc -l < "$raw_list_todo")
    n_done=$((n_all - n_todo))

    echo "原始樣本總數: $n_all"
    echo "已完成樣本數: $n_done"
    echo "本次需新增 trimming 樣本數: $n_todo"

    if [[ "$n_todo" -eq 0 ]]; then
        echo "[Stage 1] 無新增樣本，略過 trimming。"
        return
    fi

    parallel -j "$JOBS" --bar "
      r1=\"{}\"
      r2=\$(echo \"\$r1\" | sed 's/_R1_/_R2_/')
      ext=\$(basename \"\$r1\" | sed 's/.*_R1_001//')
      base=\$(basename \"\$r1\" _R1_001\$ext)
      fastp -i \"\$r1\" -I \"\$r2\" -o \"$STAGE1/trim/\${base}_R1_001.fastq.gz\" -O \"$STAGE1/trim/\${base}_R2_001.fastq.gz\" \
        --thread 2 --qualified_quality_phred 30 --length_required 80 \
        --html \"$STAGE1/fastp_report/\${base}.html\" --json \"$STAGE1/fastp_report/\${base}.json\"
    " < "$raw_list_todo"

    N_TRIM=$(wc -l < "$raw_list_todo")
    echo "-------------------------------------------------------"
    echo "[Stage 1 完成回報]"
    echo "本次新增處理樣本數: $N_TRIM"
    echo "累積原始樣本總數: $n_all"
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

    # 固定輸出下游 Stage 3/4/5/6 會使用的核心 bam list 介面
    ls -d "$(realpath "$STAGE2/mapped_bam")"/*.bam > "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile"
    BAM_LIST="$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile"

    # Stage2 固定輸出 mapping summary CSV，供 Stage3 直接複製使用
    if ! build_mapping_summary_csv "$STAGE2/mapping_results" "$STAGE2/${PROJECT_NAME}_mapping_summary.csv"; then
        echo "錯誤：Stage2 無法產生 mapping summary CSV。"
        exit 1
    fi

    N_MAPPED=$(ls "$STAGE2/mapped_bam/"*.bam | wc -l)
    echo "-------------------------------------------------------"
    echo "[Stage 2 完成回報]"
    echo "使用的參考基因組: $REF_GENOME"
    echo "比對指令參數: bwa mem -t $THREADS"
    echo "完成步驟: SAM轉換、BAM過濾(F4)、排序與建立索引"
    echo "產出的比對檔案(BAM): $STAGE2/mapped_bam/"
    echo "比對率統計結果: $STAGE2/mapping_results/"
    echo "Mapping Summary CSV: $STAGE2/${PROJECT_NAME}_mapping_summary.csv"
    echo "共計完成比對樣本數: $N_MAPPED"
    echo "-------------------------------------------------------"
}

run_stage3_pca() {
    local summary_csv summary_csv_stage2 bam_list_local summary_source_dir
    summary_csv_stage2="$STAGE2/${PROJECT_NAME}_mapping_summary.csv"
    summary_csv="$STAGE3/${PROJECT_NAME}_mapping_summary.csv"
    bam_list_local="$STAGE3/${PROJECT_NAME}_bwa_mapped.bamfile"

    echo "[Stage 3] 生成比對報表與 PCA 品質檢測..."

    if [[ "$CHAIN_STAGES" == true ]]; then
        if [ ! -s "$bam_list_local" ]; then
            if [ -s "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile" ]; then
                cp "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile" "$bam_list_local"
            else
                ls -d "$(realpath "$STAGE2/mapped_bam")"/*.bam > "$bam_list_local"
            fi
        fi
        summary_source_dir="$STAGE2/mapping_results"
    else
        if [ -n "$BAM_LIST_PCA_INPUT" ]; then
            normalize_bamfile_to_absolute "$BAM_LIST_PCA_INPUT" "$bam_list_local"
        elif [ -s "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile" ]; then
            cp "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile" "$bam_list_local"
        elif [ -n "$BAM_LIST" ] && [ -s "$BAM_LIST" ]; then
            normalize_bamfile_to_absolute "$BAM_LIST" "$bam_list_local"
        else
            echo "錯誤：Stage3 需要 BAM list，但尚未提供且前序 Stage2 也未產生。"
            exit 1
        fi
        summary_source_dir="$STAGE3/mapping_results_from_bam"
        mkdir -p "$summary_source_dir"

        while IFS= read -r bam_path; do
            [ ! -f "$bam_path" ] && { echo "錯誤：BAM 不存在: $bam_path"; exit 1; }
            sample=$(basename "$bam_path" .bam)
            samtools flagstat "$bam_path" > "$summary_source_dir/${sample}.txt"
        done < "$bam_list_local"
    fi

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

    # mapping summary CSV 的主檔固定在 Stage2，Stage3 複製一份給 R 使用。
    mkdir -p "$STAGE2"
    if ! build_mapping_summary_csv "$summary_source_dir" "$summary_csv_stage2"; then
        echo "錯誤：無法建立 Stage2 mapping summary CSV：$summary_csv_stage2"
        exit 1
    fi
    cp "$summary_csv_stage2" "$summary_csv"
    if [ ! -s "$summary_csv" ]; then
        echo "錯誤：Stage3 無法複製 mapping summary CSV：$summary_csv"
        exit 1
    fi

    echo "-------------------------------------------------------"
    echo "[Stage 3 進度回報]"
    echo "Mapping Summary 資料表已生成: $summary_csv"
    echo "Stage2 Mapping Summary 主檔: $summary_csv_stage2"
    echo "對應統計來源: $summary_source_dir"
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
    PCA_R_STATUS=$?
    OUTLIER_FILE="$STAGE3/${PROJECT_NAME}_outliers.txt"

    # Stage3 輸入/輸出資料契約保護：
    # input:  X_bwa_mapped.bamfile
    # output: X_after_pca.bamfile
    if [ ! -s "$bam_list_local" ]; then
        echo "錯誤：Stage3 輸入 bamfile 不存在或為空：$bam_list_local"
        exit 1
    fi
    mkdir -p "$STAGE3"
    if [ "$PCA_R_STATUS" -ne 0 ]; then
        echo "[Stage 3] PCA 分析未完整完成，改用保底策略：保留全部樣本。"
        : > "$OUTLIER_FILE"
        cp "$bam_list_local" "$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
    else
        # 不改 R code，直接以 R 產出的 outliers.txt 重新建立 after_pca.bamfile，
        # 確保一定從 X_bwa_mapped.bamfile 正確移除 outlier 對應樣本。
        if [ -s "$OUTLIER_FILE" ]; then
            awk '
            NR==FNR { out[$1]=1; next }
            {
              bam=$0
              sample=bam
              sub(/^.*\//, "", sample)
              sub(/\.bam$/, "", sample)
              if (!(sample in out)) print bam
            }' "$OUTLIER_FILE" "$bam_list_local" > "$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        else
            cp "$bam_list_local" "$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        fi
    fi

    if [ ! -s "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
        echo "[Stage 3] after_pca 為空或不存在，啟用保底：保留全部樣本。"
        cp "$bam_list_local" "$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
    fi
    if [ ! -s "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
        echo "錯誤：無法建立 $STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        exit 1
    fi

    N_ORIG=$(wc -l < "$bam_list_local")
    N_AFTER=$(wc -l < "$STAGE3/${PROJECT_NAME}_after_pca.bamfile")
    N_DIFF=$((N_ORIG - N_AFTER))

    echo "-------------------------------------------------------"
    echo "[Mapping Outlier偵測結果]"
    echo "原始輸入樣本總數: $N_ORIG"
    echo "PCA 檢測建議保留樣本數: $N_AFTER"
    echo "被剔除的樣本數: $N_DIFF"
    echo "-------------------------------------------------------"

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

        if [ "$PCA_DECISION" == "1" ]; then
            BAM_LIST="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
            echo "[!] 已套用過濾後的 BAM 清單，當前分析為: $(wc -l < "$BAM_LIST") 個樣本。"
        else
            BAM_LIST="$bam_list_local"
            echo "[+] 已選擇保留outlier樣本，維持原始分析樣本數。"
        fi
    else
        BAM_LIST="$bam_list_local"
    fi
}

run_stage4_clone() {
    local clone_bam_list
    if [[ "$CHAIN_STAGES" == true ]]; then
        [ -z "$BAM_LIST" ] && BAM_LIST="$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile"
        clone_bam_list="$BAM_LIST"
    else
        if [ -n "$BAM_LIST_CLONE_INPUT" ]; then
            clone_bam_list="$BAM_LIST_CLONE_INPUT"
        elif [ -n "$BAM_LIST" ] && [ -s "$BAM_LIST" ]; then
            clone_bam_list="$BAM_LIST"
        elif [ -s "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
            clone_bam_list="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        elif [ -s "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile" ]; then
            clone_bam_list="$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile"
        else
            echo "錯誤：Stage4 需要 BAM list，但尚未提供且前序 Stage2/3 未產生。"
            exit 1
        fi
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

        if [ "$CLONE_DECISION" == "1" ]; then
            BAM_LIST="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
            echo "[!] 已套用去除clones後的最終 BAM 清單，當前分析樣本數: $(wc -l < "$BAM_LIST") 個樣本。"
        else
            BAM_LIST="$clone_bam_list"
            echo "[+] 已選擇保留Clone樣本，維持樣本數: $N_CLONE_BEFORE"
        fi
    fi

    if [ ! -s "$CLONE_REVIEW_FILE" ]; then
        BAM_LIST="$clone_bam_list"
    fi
}

resolve_bam_list_for_stage5_to_7() {
    if [[ "$CHAIN_STAGES" == true ]]; then
        LIST_FINAL="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
        LIST_PCA_ONLY="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        LIST_FULL="$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile"

        if [ -f "$LIST_FINAL" ]; then
            BAM_LIST="$LIST_FINAL"
        elif [ -f "$LIST_PCA_ONLY" ]; then
            BAM_LIST="$LIST_PCA_ONLY"
        else
            BAM_LIST="$LIST_FULL"
        fi
    else
        if [ -n "$BAM_LIST" ] && [ -s "$BAM_LIST" ]; then
            :
        elif [ -s "$STAGE4/${PROJECT_NAME}_after_clones.bamfile" ]; then
            BAM_LIST="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
        elif [ -s "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
            BAM_LIST="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        elif [ -s "$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile" ]; then
            BAM_LIST="$STAGE2/${PROJECT_NAME}_bwa_mapped.bamfile"
        elif [[ "$RUN_S5" == "y" ]]; then
            BAM_LIST="$BAM_LIST_LD_INPUT"
        elif [[ "$RUN_S6" == "y" || "$RUN_S7" == "y" ]]; then
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

run_stage5_all_snp_sites() {
    ask_to_run "All SNP Site Map" "$STAGE5/all_snp.sites" SKIP_S5
    if [[ "$SKIP_S5" == true ]]; then
        return
    fi

    echo "[Stage 5] 建立 All SNP sites map..."
    angsd -b "$BAM_LIST" -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 30 -baq 1 -setMinDepth 5 -SNP_pval 1e-6 -skipTriallelic 1 -doHWE 1 -Hetbias_pval 0.00001 -minInd "$MIN_IND" -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 32 -doCounts 1 -ref "$REF_GENOME" -P 1 -out "$STAGE5/allsnps"

    gzip -kfd "$STAGE5"/*.gz
    gunzip -c "$STAGE5/allsnps.mafs.gz" | tail -n +2 | cut -f 1,2 > "$STAGE5/all_snp.sites"
    N_SITES=$(wc -l < "$STAGE5/all_snp.sites")
    angsd sites index "$STAGE5/all_snp.sites"

    echo "-------------------------------------------------------"
    echo "[Stage 5 完成回報]"
    echo "All SNP 位點數: $N_SITES"
    echo "All SNP sites map: $STAGE5/all_snp.sites"
    echo "-------------------------------------------------------"
}

run_stage6_ld_pruning_sites() {
    local all_sites_file n_sites_after
    ask_to_run "LD Pruned SNP Site Map" "$STAGE6/LD_pruned_snp.sites" SKIP_S6
    if [[ "$SKIP_S6" == true ]]; then
        return
    fi

    if [[ "$RUN_S5" == "y" ]]; then
        all_sites_file="$STAGE5/all_snp.sites"
    else
        all_sites_file="$ALL_SITES_INPUT"
    fi

    if [ ! -s "$all_sites_file" ]; then
        echo "錯誤：Stage6 需要 All SNP sites map。"
        echo "請先執行 Stage5，或提供有效的 All SNP .sites。"
        exit 1
    fi
    if [ ! -s "$STAGE5/allsnps.geno" ]; then
        echo "錯誤：缺少 $STAGE5/allsnps.geno。請先執行 Stage5 建立 All SNP 基礎檔案。"
        exit 1
    fi

    N_SITES=$(wc -l < "$all_sites_file")
    ngsLD --geno "$STAGE5/allsnps.geno" --verbose 1 --probs 1 --n_ind "$N_IND" --n_sites "$N_SITES" --max_kb_dist 50 --pos "$all_sites_file" --n_threads "$THREADS" --extend_out 1 --out "$STAGE6/allsnpsites.LD"
    prune_graph --header -v -n "$THREADS" --in "$STAGE6/allsnpsites.LD" --weight-field "r2" --weight-filter "dist <=10000 && r2 >= 0.5" --out "$STAGE6/allsnpsites.pos"
    sed 's/:/\t/g' "$STAGE6/allsnpsites.pos" | awk '$2!=""' | sort -k1 > "$STAGE6/LD_pruned_snp.sites"
    angsd sites index "$STAGE6/LD_pruned_snp.sites"
    n_sites_after=$(wc -l < "$STAGE6/LD_pruned_snp.sites")
    echo "-------------------------------------------------------"
    echo "[Stage 6 完成回報]"
    echo "LD Pruning 前位點數: $N_SITES"
    echo "LD Pruning 後位點數: $n_sites_after"
    echo "LD pruned sites map: $STAGE6/LD_pruned_snp.sites"
    echo "-------------------------------------------------------"
}

run_stage7_final_snp_with_mode() {
    local target_mode="$1"
    local sites_file="$2"
    local output_prefix="$3"
    local mode_label="$4"
    local vcf_sample_map_file

    if [ ! -f "$sites_file" ]; then
        echo "錯誤：$mode_label 需要位點表，但找不到：$sites_file"
        exit 1
    fi
    if [ ! -f "${sites_file}.idx" ]; then
        echo "錯誤：未發現 sites index 檔案 (${sites_file}.idx)。"
        exit 1
    fi

    echo "[Stage 7] 執行 $mode_label ..."
    angsd -sites "$sites_file" -b "$BAM_LIST" -GL 1 -P 1 -minInd "$MIN_IND" -minMapQ 20 -minQ 25 -sb_pval 1e-5 -Hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGlf 2 -dosnpstat 1 -doPost 1 -doGeno 8 -doBcf 1 --ignore-RG 0 -doHWE 1 -ref "$REF_GENOME" -out "$output_prefix"

    bcftools view -O v -o "${output_prefix}.vcf" "${output_prefix}.bcf"
    vcf_sample_map_file="${output_prefix}_vcf_sample_rename_map.tsv"
    sanitize_vcf_sample_ids_inplace "${output_prefix}.vcf" "$vcf_sample_map_file"
    FINAL_SNPS=$(bcftools view -H "${output_prefix}.vcf" | wc -l)

    # --------------------------------------------------------------------------
    # VCF -> STRUCTURE(.str) 轉檔 (PGDSpider3-cli)
    # --------------------------------------------------------------------------
    local final_vcf final_str spid_file pgdspider_jar pgd_cmd
    local final_vcf_abs final_str_abs spid_file_abs jar_cmd_path
    final_vcf="${output_prefix}.vcf"
    final_str="${output_prefix}.str"
    spid_file="$STAGE7/VCF2STR.spid"

    cat << 'SPID_CODE' > "$spid_file"
# spid-file generated: Wed Feb 18 21:28:41 CST 2026

# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Include a file with population definitions:
VCF_PARSER_POP_FILE_QUESTION=
# What is the ploidy of the data? (DIPLOID/HAPLOID)
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include INDELS as STANDARD genetic markers? (TRUE/FALSE)
VCF_PARSER_INDEL_QUESTION=false
# Do you want to include non-polymorphic SNPs? (TRUE/FALSE)
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field? (TRUE/FALSE)
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data? (TRUE/FALSE)
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false

# STRUCTURE / fastSTRUCTURE Writer questions
WRITER_FORMAT=STRUCTURE

# Do you want to include inter-marker distances? (TRUE/FALSE)
STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
# Specify which data type should be included in the STRUCTURE / fastSTRUCTURE file  (STRUCTURE / fastSTRUCTURE can only analyze one data type per file): (MICROSAT/SNP/STANDARD/DNA)
STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
# If SNP data are encoded as nucleotides,  enter the integers that code for nucleotide A, T, C, G (comma separated, e.g 1,2,3,4):
STRUCTURE_WRITER_SNP_CODE_QUESTION=
# Save more specific fastSTRUCTURE format? (TRUE/FALSE)
STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
SPID_CODE

    final_vcf_abs=$(realpath "$final_vcf")
    final_str_abs="$(realpath "$(dirname "$final_str")")/$(basename "$final_str")"
    spid_file_abs=$(realpath "$spid_file")

    pgdspider_jar=""
    if [ -n "${PGDSPIDER_JAR:-}" ] && [ -f "${PGDSPIDER_JAR}" ]; then
        pgdspider_jar=$(realpath "${PGDSPIDER_JAR}")
    elif jar_cmd_path=$(command -v PGDSpider3-cli.jar 2>/dev/null); then
        pgdspider_jar=$(realpath "$jar_cmd_path")
    elif [ -f "$PWD/PGDSpider3-cli.jar" ]; then
        pgdspider_jar=$(realpath "$PWD/PGDSpider3-cli.jar")
    elif [ -f "PGDSpider3-cli.jar" ]; then
        pgdspider_jar=$(realpath "PGDSpider3-cli.jar")
    fi

    local jar_show vcf_show str_show spid_show
    if [ -n "$pgdspider_jar" ]; then
        printf -v jar_show '%q' "$pgdspider_jar"
    else
        jar_show="/ABSOLUTE/PATH/TO/PGDSpider3-cli.jar"
    fi
    printf -v vcf_show '%q' "$final_vcf_abs"
    printf -v str_show '%q' "$final_str_abs"
    printf -v spid_show '%q' "$spid_file_abs"
    pgd_cmd="java -Xmx1024m -Xms512m -jar $jar_show -inFile $vcf_show -inFormat VCF -outFile $str_show -outFormat STRUCTURE -spid $spid_show"

    if ! command -v java >/dev/null 2>&1; then
        echo "[警告] 未安裝 java，無法自動轉換 VCF -> STR。"
        echo "請安裝 Java 後手動執行："
        echo "$pgd_cmd"
    elif [ -z "$pgdspider_jar" ]; then
        echo "[警告] 找不到 PGDSpider3-cli.jar，無法自動轉換 VCF -> STR。"
        echo "已輸出 spid 檔案：$spid_file"
        echo "請準備好 PGDSpider3-cli.jar 後手動執行："
        echo "$pgd_cmd"
    else
        echo "[Stage 7] 轉換 VCF -> STRUCTURE(.str): $target_mode"
        java -Xmx1024m -Xms512m -jar "$pgdspider_jar" \
          -inFile "$final_vcf_abs" \
          -inFormat VCF \
          -outFile "$final_str_abs" \
          -outFormat STRUCTURE \
          -spid "$spid_file_abs"

        if [ $? -eq 0 ]; then
            normalize_structure_str_header "$final_str_abs"
            echo "[完成] STRUCTURE 檔案已產出: $final_str"
            echo "spid 設定檔: $spid_file"
        else
            echo "[警告] PGDSpider 轉檔失敗。"
            echo "已輸出 spid 檔案：$spid_file"
            echo "請檢查 PGDSpider 安裝與參數後手動執行："
            echo "$pgd_cmd"
        fi
    fi

    echo "[$mode_label] SNP 數量: $FINAL_SNPS"
    echo "[$mode_label] VCF: $final_vcf"
    echo "[$mode_label] STR: $final_str"
}

run_stage7_final_snp() {
    local ld_sites all_sites
    ask_to_run "Stage7 Final SNP Calling" "$STAGE7/${PROJECT_NAME}_snps_final_with_LD_Pruning.vcf" SKIP_S7
    if [[ "$SKIP_S7" == true ]]; then
        return
    fi

    if [[ "$RUN_S7_WITH_LD" == "y" ]]; then
        if [[ "$RUN_S6" == "y" ]]; then
            ld_sites="$STAGE6/LD_pruned_snp.sites"
        else
            ld_sites="$LD_SITES_INPUT"
        fi
        run_stage7_final_snp_with_mode "with_LD_Pruning" "$ld_sites" "$STAGE7/${PROJECT_NAME}_snps_final_with_LD_Pruning" "Final SNP with LD pruning"
    fi

    if [[ "$RUN_S7_SKIP_LD" == "y" ]]; then
        if [[ "$RUN_S5" == "y" ]]; then
            all_sites="$STAGE5/all_snp.sites"
        else
            all_sites="$ALL_SITES_INPUT"
        fi
        run_stage7_final_snp_with_mode "Skip_LD_Pruning" "$all_sites" "$STAGE7/${PROJECT_NAME}_snps_final_Skip_LD_Pruning" "Skip LD-pruning SNP Calling"
    fi
}

run_stage9_genetic_divergence() {
    local stage9_dir stage9_pop_dir stage9_all_bamfile stage9_all_bam_abs
    local stage9_fst_dir stage9_stats2_dir stage9_matrix_dir
    local stage9_all_n stage9_all_minind
    local pop_files=() pop_norm_files=() pop_names=()
    local choice input_path target_name target_path
    local f n_pop minind_pop pop_name pop_norm
    local i j p1 p2 pair_tag pair_sfs fst_prefix fst_idx
    local fst_stats_file fst_stats2_file fst_summary_file
    local stats_line raw_line uw_val wt_val val

    if ! command -v realSFS >/dev/null 2>&1; then
        echo "錯誤：找不到 realSFS，無法執行 Stage9。"
        return 1
    fi

    if [ -z "$BAM_LIST_DIV_ALL_INPUT" ] || [ ! -s "$BAM_LIST_DIV_ALL_INPUT" ]; then
        if [ -n "$BAM_LIST" ] && [ -s "$BAM_LIST" ]; then
            BAM_LIST_DIV_ALL_INPUT="$BAM_LIST"
        elif [ -s "$STAGE4/${PROJECT_NAME}_after_clones.bamfile" ]; then
            BAM_LIST_DIV_ALL_INPUT="$STAGE4/${PROJECT_NAME}_after_clones.bamfile"
        elif [ -s "$STAGE3/${PROJECT_NAME}_after_pca.bamfile" ]; then
            BAM_LIST_DIV_ALL_INPUT="$STAGE3/${PROJECT_NAME}_after_pca.bamfile"
        else
            echo "錯誤：Stage9 全族群 bamfile 不存在或為空：$BAM_LIST_DIV_ALL_INPUT"
            return 1
        fi
    fi

    if [ -z "$STAGE9_SKIP_LD_VCF_INPUT" ] && [ -f "$STAGE7/${PROJECT_NAME}_snps_final_Skip_LD_Pruning.vcf" ]; then
        STAGE9_SKIP_LD_VCF_INPUT="$STAGE7/${PROJECT_NAME}_snps_final_Skip_LD_Pruning.vcf"
    fi

    stage9_dir="$STAGE9/divergence"
    stage9_pop_dir="$stage9_dir/population_bamfiles"
    stage9_fst_dir="$stage9_dir/fst_results"
    stage9_matrix_dir="$stage9_dir/fst_matrices"
    mkdir -p "$stage9_pop_dir" "$stage9_fst_dir" "$stage9_matrix_dir"
    if [[ "$S9_RUN_STATS2" == "y" || "$S9_RUN_STATS2" == "Y" ]]; then
        stage9_stats2_dir="$stage9_dir/fst_stats2"
        mkdir -p "$stage9_stats2_dir"
    fi
    STAGE9_LAST_RUN_DIR="$stage9_dir"

    stage9_all_bam_abs="$stage9_dir/all_populations.bamfile"
    normalize_bamfile_to_absolute "$BAM_LIST_DIV_ALL_INPUT" "$stage9_all_bam_abs"
    stage9_all_bamfile="$stage9_all_bam_abs"

    echo "-------------------------------------------------------"
    echo "[Stage 9] 分析資料夾已建立：$stage9_dir"
    [ -n "$STAGE9_SKIP_LD_VCF_INPUT" ] && echo "[Stage 9] 參考 Skip LD VCF: $STAGE9_SKIP_LD_VCF_INPUT"
    echo "請準備每個族群的 population.bamfile 放入：$stage9_pop_dir"
    echo "檔名建議：pop1.bamfile, pop2.bamfile ..."
    echo "-------------------------------------------------------"

    while true; do
        pop_files=()
        while IFS= read -r f; do
            pop_files+=("$f")
        done < <(find "$stage9_pop_dir" -maxdepth 1 -type f -name "*.bamfile" ! -name "all_populations.bamfile" | sort)

        echo "目前偵測到 ${#pop_files[@]} 個 population.bamfile："
        for i in "${!pop_files[@]}"; do
            printf "%2d) %s\n" "$((i+1))" "$(basename "${pop_files[$i]}")"
        done
        echo "r) 重新掃描"
        echo "m) 從其他位置匯入 bamfile"
        echo "c) 準備完成，繼續執行 Stage9"
        echo "q) 取消 Stage9"
        read -p "請選擇: " choice

        case "$choice" in
            r|R) ;;
            m|M)
                read -e -p "請輸入要匯入的 bamfile 路徑: " input_path
                if [ ! -s "$input_path" ]; then
                    echo "錯誤：找不到或空檔 $input_path"
                    continue
                fi
                input_path=$(realpath "$input_path")
                target_name=$(basename "$input_path")
                target_path="$stage9_pop_dir/$target_name"
                cp "$input_path" "$target_path"
                echo "已匯入：$target_path"
                ;;
            c|C)
                if [ "${#pop_files[@]}" -lt 2 ]; then
                    echo "錯誤：至少需要 2 個 population.bamfile 才能計算族群間 divergence。"
                    continue
                fi
                break
                ;;
            q|Q)
                echo "使用者取消 Stage9。"
                return 0
                ;;
            *)
                echo "錯誤：無效選擇。"
                ;;
        esac
    done

    stage9_all_n=$(wc -l < "$stage9_all_bamfile")
    stage9_all_minind=$(( stage9_all_n * 8 / 10 ))
    [ "$stage9_all_minind" -lt 1 ] && stage9_all_minind=1

    echo "[Stage 9 - Step 1] 建立 AllSites..."
    angsd -b "$stage9_all_bamfile" -GL 1 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -Hetbias_pval 1e-5 -minInd "$stage9_all_minind" -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8 -P 1 -out "$stage9_dir/AllSites"
    gunzip -c "$stage9_dir/AllSites.mafs.gz" | tail -n +2 | cut -f1,2 > "$stage9_dir/AllSites.sites"
    angsd sites index "$stage9_dir/AllSites.sites"

    echo "[Stage 9 - Step 2] 逐族群計算 SAF..."
    pop_norm_files=()
    pop_names=()
    for f in "${pop_files[@]}"; do
        pop_name=$(basename "$f" .bamfile)
        pop_norm="$stage9_dir/${pop_name}.normalized.bamfile"
        normalize_bamfile_to_absolute "$f" "$pop_norm"
        pop_norm_files+=("$pop_norm")
        pop_names+=("$pop_name")

        n_pop=$(wc -l < "$pop_norm")
        minind_pop=$(( n_pop * 8 / 10 ))
        [ "$minind_pop" -lt 1 ] && minind_pop=1

        angsd -sites "$stage9_dir/AllSites.sites" -b "$pop_norm" -GL 1 -P 1 -minInd "$minind_pop" -doSaf 1 -anc "$REF_GENOME" -ref "$REF_GENOME" -out "$stage9_dir/${pop_name}"
    done

    echo "[Stage 9 - Step 3] 產生每族群 SFS..."
    for pop_name in "${pop_names[@]}"; do
        realSFS -cores "$THREADS" "$stage9_dir/${pop_name}.saf.idx" > "$stage9_dir/${pop_name}.sfs"
    done

    echo "[Stage 9 - Step 4] 計算 pairwise Fst..."
    fst_summary_file="$stage9_matrix_dir/fst_pairwise_summary.tsv"
    echo -e "Population1\tPopulation2\tFST.Unweight\tFST.Weight" > "$fst_summary_file"

    for ((i=0; i<${#pop_names[@]}; i++)); do
        for ((j=i+1; j<${#pop_names[@]}; j++)); do
            p1="${pop_names[$i]}"
            p2="${pop_names[$j]}"
            pair_tag="${p1}_${p2}"
            pair_sfs="$stage9_fst_dir/${p1}.${p2}.sfs"
            fst_prefix="$stage9_fst_dir/$pair_tag"
            realSFS -cores "$THREADS" "$stage9_dir/${p1}.saf.idx" "$stage9_dir/${p2}.saf.idx" > "$pair_sfs"
            realSFS fst index "$stage9_dir/${p1}.saf.idx" "$stage9_dir/${p2}.saf.idx" -sfs "$pair_sfs" -fstout "$fst_prefix" -cores "$THREADS"

            fst_idx="${fst_prefix}.fst.idx"
            fst_stats_file="$stage9_fst_dir/${pair_tag}.fst.txt"
            {
                echo "# Pair: $p1 vs $p2"
                echo "# Columns: FST.Unweight<TAB>FST.Weight"
                realSFS fst stats "$fst_idx" -cores "$THREADS" 2>&1
            } > "$fst_stats_file"

            if [[ "$S9_RUN_STATS2" == "y" || "$S9_RUN_STATS2" == "Y" ]]; then
                fst_stats2_file="$stage9_stats2_dir/${pair_tag}.fst.stats2.txt"
                {
                    echo "# Pair: $p1 vs $p2"
                    echo "# Output columns from realSFS fst stats2"
                    echo "# Typically includes: region  chr  midPos  Nsites (and related window stats)"
                    realSFS fst stats2 "$fst_idx" -cores "$THREADS" 2>&1
                } > "$fst_stats2_file"
            fi

            stats_line=$(grep -E 'FST\.Unweight.*Fst\.Weight' "$fst_stats_file" | tail -n1)
            if [ -n "$stats_line" ]; then
                uw_val=$(echo "$stats_line" | sed -E 's/.*FST\.Unweight[^:]*:([0-9eE+.-]+).*/\1/')
                wt_val=$(echo "$stats_line" | sed -E 's/.*Fst\.Weight:([0-9eE+.-]+).*/\1/')
            else
                raw_line=$(grep -E '^[[:space:]]*[0-9eE+.-]+[[:space:]]+[0-9eE+.-]+[[:space:]]*$' "$fst_stats_file" | tail -n1)
                uw_val=$(echo "$raw_line" | awk '{print $1}')
                wt_val=$(echo "$raw_line" | awk '{print $2}')
            fi
            [ -z "$uw_val" ] && uw_val="NA"
            [ -z "$wt_val" ] && wt_val="NA"
            echo -e "$p1\t$p2\t$uw_val\t$wt_val" >> "$fst_summary_file"
        done
    done

    {
        printf "Population"
        for p1 in "${pop_names[@]}"; do printf "\t%s" "$p1"; done
        printf "\n"
        for p1 in "${pop_names[@]}"; do
            printf "%s" "$p1"
            for p2 in "${pop_names[@]}"; do
                if [ "$p1" = "$p2" ]; then
                    val="0"
                else
                    val=$(awk -F'\t' -v a="$p1" -v b="$p2" 'NR>1 && (($1==a && $2==b) || ($1==b && $2==a)) {print $3; exit}' "$fst_summary_file")
                    [ -z "$val" ] && val="NA"
                fi
                printf "\t%s" "$val"
            done
            printf "\n"
        done
    } > "$stage9_matrix_dir/fst_unweight_matrix.tsv"

    {
        printf "Population"
        for p1 in "${pop_names[@]}"; do printf "\t%s" "$p1"; done
        printf "\n"
        for p1 in "${pop_names[@]}"; do
            printf "%s" "$p1"
            for p2 in "${pop_names[@]}"; do
                if [ "$p1" = "$p2" ]; then
                    val="0"
                else
                    val=$(awk -F'\t' -v a="$p1" -v b="$p2" 'NR>1 && (($1==a && $2==b) || ($1==b && $2==a)) {print $4; exit}' "$fst_summary_file")
                    [ -z "$val" ] && val="NA"
                fi
                printf "\t%s" "$val"
            done
            printf "\n"
        done
    } > "$stage9_matrix_dir/fst_weight_matrix.tsv"

    echo "-------------------------------------------------------"
    echo "[Stage 9 完成回報]"
    echo "AllSites: $stage9_dir/AllSites.sites"
    echo "Population 檔案資料夾: $stage9_pop_dir"
    echo "FST 結果資料夾: $stage9_fst_dir"
    if [[ "$S9_RUN_STATS2" == "y" || "$S9_RUN_STATS2" == "Y" ]]; then
        echo "FST stats2 資料夾: $stage9_stats2_dir"
    else
        echo "FST stats2: 已略過（可下次選 y 啟用）"
    fi
    echo "FST matrix 資料夾: $stage9_matrix_dir"
    echo "輸出資料夾: $stage9_dir"
    echo "-------------------------------------------------------"
}

run_stage8_structure_auto() {
    local stage7_str_input stage7_str_abs stage7_str_base
    local stage7_numind_default stage7_numloci_default
    local latest_params reuse_prev edit_mode
    local run_label run_dir
    local stage7_str_filename
    local run_dir_show
    local strprefix harprefix postfix
    local label_i popdata_i popflag_i locdata_i pheno_i extracols_i
    local markers_i dominant_i mapdist_i phase_i markov_i
    local parallel_enabled harvest_enabled
    local myK run seed

    if [ -n "$STR_INPUT" ]; then
        stage7_str_input="$STR_INPUT"
    elif [ -f "$STAGE7/${PROJECT_NAME}_snps_final_with_LD_Pruning.str" ]; then
        stage7_str_input="$STAGE7/${PROJECT_NAME}_snps_final_with_LD_Pruning.str"
    else
        stage7_str_input="$STR_INPUT"
    fi
    if [ ! -f "$stage7_str_input" ]; then
        echo "錯誤：找不到 Stage 8 輸入 .str：$stage7_str_input"
        echo "請先執行 Stage7 with LD pruning SNP Calling，或手動指定 .str。"
        return 1
    fi

    stage7_str_abs=$(realpath "$stage7_str_input")
    normalize_structure_str_header "$stage7_str_abs"
    stage7_str_base=$(basename "$stage7_str_abs" .str)
    S7_STR_FILE="$stage7_str_abs"

    if infer_stage7_str_metrics "$S7_STR_FILE" stage7_numind_default stage7_numloci_default; then
        :
    else
        stage7_numind_default=129
        stage7_numloci_default=1289
        echo "[警告] 無法由 .str 自動推算 numind/numloci，將使用預設值。"
    fi

    set_stage7_default_values "$stage7_str_base" "$stage7_numind_default" "$stage7_numloci_default"
    print_stage7_parameter_notes

    latest_params=$(find "$STAGE8" -maxdepth 2 -type f -name "stage7_params.env" | sort | tail -n1)
    if [ -n "$latest_params" ] && [ -f "$latest_params" ]; then
        echo "偵測到上次 Stage7 參數檔：$latest_params"
        read -p "是否沿用上次參數？(y/n): " reuse_prev
        if [[ "$reuse_prev" == "y" || "$reuse_prev" == "Y" ]]; then
            # shellcheck disable=SC1090
            source "$latest_params"
        else
            # shellcheck disable=SC1090
            source "$latest_params"
            read -p "參數要如何調整？(m:修改部分 / a:全部重填) [m]: " edit_mode
            [ -z "$edit_mode" ] && edit_mode="m"
            if [[ "$edit_mode" == "a" || "$edit_mode" == "A" ]]; then
                prompt_stage7_all
            else
                prompt_stage7_selected
            fi
        fi
    else
        echo "未發現舊參數檔，請輸入 Stage7 參數。"
        prompt_stage7_all
    fi

    if [ ! -f "$S7_STR_FILE" ]; then
        echo "[警告] 參數中的 .str 檔案不存在，恢復使用當前 Stage7 輸入檔案。"
        S7_STR_FILE="$stage7_str_abs"
    fi
    S7_DATASET=$(basename "$S7_STR_FILE" .str)
    if infer_stage7_str_metrics "$S7_STR_FILE" stage7_numind_default stage7_numloci_default; then
        S7_NUMIND="$stage7_numind_default"
        S7_NUMLOCI="$stage7_numloci_default"
    fi

    show_stage7_params
    run_label="STRUCTURE"
    run_dir="$STAGE8/$run_label"
    mkdir -p "$run_dir"

    if [ ! -f "$S7_STR_FILE" ]; then
        echo "錯誤：Stage7 .str 檔案不存在：$S7_STR_FILE"
        return 1
    fi

    normalize_structure_str_header "$S7_STR_FILE"
    stage7_str_filename=$(basename "$S7_STR_FILE")
    cp "$S7_STR_FILE" "$run_dir/$stage7_str_filename"

    label_i=$(stage7_bool_to_int "$S7_LABEL")
    popdata_i=$(stage7_bool_to_int "$S7_POPDATA")
    popflag_i=$(stage7_bool_to_int "$S7_POPFLAG")
    locdata_i=$(stage7_bool_to_int "$S7_LOCDATA")
    pheno_i=$(stage7_bool_to_int "$S7_PHENO")
    extracols_i=$(stage7_bool_to_int "$S7_EXTRACOLS")
    markers_i=$(stage7_bool_to_int "$S7_MARKERS")
    dominant_i=$(stage7_bool_to_int "$S7_DOMINANT")
    mapdist_i=$(stage7_bool_to_int "$S7_MAPDIST")
    phase_i=$(stage7_bool_to_int "$S7_PHASE")
    markov_i=$(stage7_bool_to_int "$S7_MARKOV")

    if stage7_bool_is_true "$S7_PARALLEL"; then
        parallel_enabled=true
    else
        parallel_enabled=false
    fi

    if stage7_bool_is_true "$S7_HARVEST"; then
        harvest_enabled=true
    else
        harvest_enabled=false
    fi

    if $harvest_enabled; then
        if ! ensure_structure_harvester; then
            echo "[警告] structureHarvester 未就緒，runstructure 將保留指令但可能無法成功執行 harvest。"
        fi
    fi

    save_stage7_params "$run_dir/stage7_params.env"

    cat << EOF > "$run_dir/mainparams"
Basic program parameters
#define MAXPOPS 	 $S7_MAXPOPS
#define BURNIN  	 $S7_BURNIN
#define NUMREPS 	 $S7_MCMC

Input file
#define INFILE 	 ${S7_DATASET}.str

Data file format
#define NUMINDS 	 $S7_NUMIND
#define NUMLOCI 	 $S7_NUMLOCI
#define PLOIDY  	 $S7_PLOIDY
#define MISSING 	 $S7_MISSING
#define ONEROWPERIND 	 $S7_ONEROWPERIND

#define LABEL   	 $label_i
#define POPDATA 	 $popdata_i
#define POPFLAG 	 $popflag_i
#define LOCDATA 	 $locdata_i
#define PHENOTYPE 	 $pheno_i
#define EXTRACOLS 	 $extracols_i
#define MARKERNAMES 	 $markers_i
#define RECESSIVEALLELES 	 $dominant_i
#define MAPDISTANCES 	 $mapdist_i

Advanced data file options
#define PHASED    	 $phase_i
#define MARKOVPHASE 	 $markov_i
#define NOTAMBIGUOUS 	 -999
EOF

    cat << 'EOF' > "$run_dir/extraparams"
#define NOADMIX 0
#define LINKAGE 0
#define USEPOPINFO 0
#define LOCPRIOR 0
#define FREQSCORR 1
#define ONEFST 0
#define INFERALPHA 1
#define POPALPHAS 0
#define ALPHA 1.0
#define INFERLAMBDA 0
#define POPSPECIFICLAMBDA 0
#define LAMBDA 1.0
#define FPRIORMEAN 0.01
#define FPRIORSD 0.05
#define UNIFPRIORALPHA 1
#define ALPHAMAX 10.0
#define ALPHAPRIORA 1.0
#define ALPHAPRIORB 2.0
#define LOG10RMIN -4.0
#define LOG10RMAX 1.0
#define LOG10RPROPSD 0.1
#define LOG10RSTART -2.0
#define GENSBACK 2
#define MIGRPRIOR 0.01
#define PFROMPOPFLAGONLY 0
#define LOCISPOP 0
#define LOCPRIORINIT 1.0
#define MAXLOCPRIOR 20.0
#define PRINTNET 1
#define PRINTLAMBDA 1
#define PRINTQSUM 1
#define SITEBYSITE 0
#define PRINTQHAT 0
#define UPDATEFREQ 100
#define PRINTLIKES 0
#define INTERMEDSAVE 0
#define ECHODATA 1
#define ANCESTDIST 0
#define NUMBOXES 1000
#define ANCESTPINT 0.90
#define COMPUTEPROB 1
#define ADMBURNIN 500
#define ALPHAPROPSD 0.025
#define STARTATPOPINFO 0
#define RANDOMIZE 0
#define SEED 2245
#define METROFREQ 10
#define REPORTHITRATE 0
EOF

    if command -v structure >/dev/null 2>&1; then
        strprefix=""
    else
        strprefix="./"
        echo "[警告] structure 不在 PATH 中，runstructure 將使用 ./structure。"
    fi

    if $harvest_enabled; then
        if check_structure_harvester_ready; then
            harprefix=""
        else
            harprefix="./"
            echo "[警告] structureHarvester.py 不在 PATH 中，runstructure 將使用 ./structureHarvester.py。"
        fi
    fi

    {
        echo "#!/bin/sh"
        echo "mkdir results_f log harvester"
        for myK in $(seq 1 "$S7_MAXPOPS"); do
            echo "mkdir k${myK}"
        done
        echo ""
        echo "cd log"
        for myK in $(seq 1 "$S7_MAXPOPS"); do
            echo "mkdir k${myK}"
        done
        echo ""
        echo "cd .."
        echo ""

        if $parallel_enabled; then
            if [ "$S7_CORE_DEF" == "percent" ]; then
                postfix="%"
            elif [ "$S7_CORE_DEF" == "number" ]; then
                postfix=""
            else
                postfix=""
                echo "[警告] core_def=$S7_CORE_DEF 非法，已改用 number 模式。"
            fi
            echo "cat structureCommands | parallel -j ${S7_CORES}${postfix}"
            echo ""
        else
            for myK in $(seq 1 "$S7_MAXPOPS"); do
                for run in $(seq 1 "$S7_KRUNS"); do
                    seed=$((100000 + RANDOM % 900000))
                    echo "${strprefix}structure -D ${seed} -K ${myK} -m mainparams -o k${myK}/${S7_DATASET}_k${myK}_run${run} 2>&1 | tee log/k${myK}/${S7_DATASET}_k${myK}_run${run}.log"
                done
            done
        fi

        for myK in $(seq 1 "$S7_MAXPOPS"); do
            printf "mv k%d " "$myK"
        done
        echo "results_f/"
        echo "mkdir harvester_input"
        echo "cp results_f/k*/*_f harvester_input"
        echo "echo 'Your structure run has finished.'"

        if $harvest_enabled; then
            echo "# Run structureHarvester"
            echo "${harprefix}structureHarvester.py --dir harvester_input --out harvester --evanno --clumpp"
            echo "echo 'structureHarvester run has finished.'"
        fi

        echo "#Clean up harvester input files."
        echo "zip ${S7_DATASET}_Harvester_Upload.zip harvester_input/*"
        echo "mv ${S7_DATASET}_Harvester_Upload.zip harvester/"
        echo "rm -rf harvester_input"
    } > "$run_dir/runstructure"

    if $parallel_enabled; then
        if ! command -v parallel >/dev/null 2>&1; then
            echo "[警告] parallel 未安裝，runstructure 已產生但平行模式無法執行。"
        fi
        : > "$run_dir/structureCommands"
        for myK in $(seq 1 "$S7_MAXPOPS"); do
            for run in $(seq 1 "$S7_KRUNS"); do
                seed=$((100000 + RANDOM % 900000))
                echo "${strprefix}structure -D ${seed} -K ${myK} -m mainparams -o k${myK}/${S7_DATASET}_k${myK}_run${run} 2>&1 > log/k${myK}/${S7_DATASET}_k${myK}_run${run}.log" >> "$run_dir/structureCommands"
            done
        done
    fi

    chmod 755 "$run_dir/runstructure"
    echo "[完成] Stage 8 檔案已產生："
    echo "  - $run_dir/mainparams"
    echo "  - $run_dir/extraparams"
    [ -f "$run_dir/structureCommands" ] && echo "  - $run_dir/structureCommands"
    echo "  - $run_dir/runstructure"

    read -p "是否立刻執行 runstructure？(y/n): " RUNSTRUCTURE_NOW
    if [[ "$RUNSTRUCTURE_NOW" == "y" || "$RUNSTRUCTURE_NOW" == "Y" ]]; then
        (cd "$run_dir" && ./runstructure)
    else
        printf -v run_dir_show '%q' "$run_dir"
        echo "已略過執行。可稍後手動執行：cd $run_dir_show && ./runstructure"
    fi
}

# ------------------------------------------------------------------------------
# 主程序
# ------------------------------------------------------------------------------
main() {
    check_dependencies
    select_analysis_scope
    configure_project_name
    configure_stage_paths
    start_logging
    collect_inputs
    print_runtime_config
    confirm_run
    setup_output_dirs

    [[ "$RUN_S1" == "y" ]] && run_stage1_fastp
    [[ "$RUN_S2" == "y" ]] && run_stage2_alignment

    if [[ "$RUN_S3" == "y" ]]; then
        run_stage3_pca
    fi

    if [[ "$RUN_S4" == "y" ]]; then
        run_stage4_clone
    fi

    if [[ "$RUN_S5" == "y" || "$RUN_S6" == "y" || "$RUN_S7" == "y" ]]; then
        resolve_bam_list_for_stage5_to_7
    fi

    [[ "$RUN_S5" == "y" ]] && run_stage5_all_snp_sites
    [[ "$RUN_S6" == "y" ]] && run_stage6_ld_pruning_sites
    [[ "$RUN_S7" == "y" ]] && run_stage7_final_snp
    [[ "$RUN_S8" == "y" ]] && run_stage8_structure_auto
    [[ "$RUN_S9" == "y" ]] && run_stage9_genetic_divergence

    echo "======================================================="
    echo "分析結束: $(date)"
    [[ "$RUN_S7" == "y" ]] && [[ "$RUN_S7_WITH_LD" == "y" ]] && echo "產出 VCF(with LD pruning): $STAGE7/${PROJECT_NAME}_snps_final_with_LD_Pruning.vcf"
    [[ "$RUN_S7" == "y" ]] && [[ "$RUN_S7_SKIP_LD" == "y" ]] && echo "產出 VCF(skip LD pruning): $STAGE7/${PROJECT_NAME}_snps_final_Skip_LD_Pruning.vcf"
    [[ "$RUN_S9" == "y" ]] && echo "Stage9 輸出資料夾: $STAGE9_LAST_RUN_DIR"
    echo "日誌位置: $LOG_FILE"
    echo "VCF可用PGDSpider做轉換"
    echo "https://software.bioinformatics.unibe.ch/pgdspider/"
    echo "此分析參考 James Fifer https://github.com/jamesfifer/JapanRE"
    echo "======================================================="
}

main "$@"
