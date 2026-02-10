#!/bin/bash

BASE_REF_DIR="$HOME/RefGenome"
BASHRC="$HOME/.bashrc"
EDIRECT_DIR="$HOME/edirect"

# 1. 基礎環境偵測
export PATH="$EDIRECT_DIR:$PATH"
for cmd in esearch bwa samtools curl gunzip bc; do
    if ! command -v $cmd &> /dev/null; then
        echo "錯誤：缺少指令 '$cmd'。請確保已安裝該軟體包。"
        exit 1
    fi
done

# 2. 檢索與物種選擇
while true; do
    echo "================================================="
    read -p "請輸入搜尋關鍵字 (物種名或 BioProject): " QUERY
    [[ -z "$QUERY" ]] && continue

    echo -e "\n正在檢索 NCBI 資料庫並解析數據結構..."
    RAW_XML_FILE="assembly_raw_output.xml"
    esearch -db assembly -query "$QUERY" | esummary > "$RAW_XML_FILE"

    RAW_LIST=$(cat "$RAW_XML_FILE" | xtract -pattern DocumentSummary -def "N/A" \
        -element AssemblyAccession AssemblyName AssemblyStatus \
        -element ScaffoldN50 Coverage AssemblyType \
        -first   GB_BioProjects/Bioproj/BioprojectAccn \
        -element BioSampleAccn \
        -first   Biosource/Isolate \
        -element SubmitterOrganization \
        -element SubmissionDate \
        -first   Meta/representative-status \
        -block Stat -if category -equals total_length -element Stat)

    if [ -z "$RAW_LIST" ]; then
        echo "找不到符合結果。"
        continue
    fi

    mapfile -t DISPLAY_RESULTS <<< "$RAW_LIST"
    echo -e "檢索結果列表：\n"
    for i in "${!DISPLAY_RESULTS[@]}"; do
        IFS=$'\t' read -r ACC NAME STATUS N50 COV TYPE BPROJ SAMN ISO ORG DATE REF SIZE <<< "${DISPLAY_RESULTS[$i]}"
        
        ISO=$(echo "$ISO" | xargs); [[ -z "$ISO" || "$ISO" == " " ]] && ISO="N/A"
        BPROJ=$(echo "$BPROJ" | xargs)
        SAMN=$(echo "$SAMN" | xargs)
        ORG=$(echo "$ORG" | xargs)
        REF=$(echo "$REF" | xargs); [[ "$REF" == "na" ]] && REF=""

        [[ "$SIZE" =~ ^[0-9]+$ ]] && SIZE_DISPLAY=$(printf "%.1fMb" "$(echo "scale=2; $SIZE/1000000" | bc)") || SIZE_DISPLAY="N/A"
        [[ "$N50" =~ ^[0-9]+$ ]] && N50_DISPLAY=$(printf "%.1fMb" "$(echo "scale=2; $N50/1000000" | bc)") || N50_DISPLAY="N/A"

        if [[ "$BPROJ" != "N/A" && "$SAMN" != "N/A" ]]; then SOURCE_INFO="$BPROJ ($SAMN)";
        elif [[ "$BPROJ" != "N/A" ]]; then SOURCE_INFO="$BPROJ";
        elif [[ "$SAMN" != "N/A" ]]; then SOURCE_INFO="$SAMN";
        else SOURCE_INFO="N/A"; fi

        STATUS_HEADER="$STATUS${REF:+ ($REF)}"
        echo "[$((i+1))] $ACC | $NAME ($STATUS_HEADER)"
        echo "    Quality: N50: $N50_DISPLAY | Coverage: ${COV}x | Size: $SIZE_DISPLAY"
        echo "    Type: $TYPE | Status: $STATUS"
        echo "    Source: $SOURCE_INFO | Isolate: $ISO"
        echo "    Organization: $ORG ($(echo "$DATE" | cut -d' ' -f1))"
        echo "--------------------------------------------------------------------------------"
    done
    echo "[r] 重新搜尋 | [q] 退出"

    read -p "請選取編號 (1-${#DISPLAY_RESULTS[@]}): " CHOICE
    if [[ "$CHOICE" == "r" ]]; then continue;
    elif [[ "$CHOICE" == "q" ]]; then exit 0;
    elif [[ "$CHOICE" =~ ^[0-9]+$ ]] && [ "$CHOICE" -ge 1 ] && [ "$CHOICE" -le "${#DISPLAY_RESULTS[@]}" ]; then
        IFS=$'\t' read -r ACC NAME STATUS N50 COV TYPE BPROJ SAMN ISO ORG DATE REF SIZE <<< "${DISPLAY_RESULTS[$((CHOICE-1))]}"
        SELECTED_ACC=$ACC
        SELECTED_BPROJ=$(echo "$BPROJ" | xargs)
        break
    else echo "無效輸入。"
    fi
done

# 3. 提取 FTP 路徑
echo "正在提取 $SELECTED_ACC 的下載路徑..."
FINAL_FTP_BASE=$(esearch -db assembly -query "$SELECTED_ACC" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)

if [ -z "$FINAL_FTP_BASE" ]; then
    echo "錯誤：找不到下載路徑。"
    exit 1
fi

# 4. 路徑配置與自動重試下載
FILE_STEM=$(basename "$FINAL_FTP_BASE")
GZ_FILE="${FILE_STEM}_genomic.fna.gz"
FULL_FTP_URL="${FINAL_FTP_BASE}/${GZ_FILE}"

SUB_DIR=${SELECTED_BPROJ}
[[ "$SUB_DIR" == "N/A" || -z "$SUB_DIR" ]] && SUB_DIR=$SELECTED_ACC
TARGET_DIR="$BASE_REF_DIR/$SUB_DIR"
mkdir -p "$TARGET_DIR" && cd "$TARGET_DIR" || exit

echo "正在下載至 $TARGET_DIR (包含自動斷線重試機制)..."
RETRY_COUNT=0
MAX_RETRIES=10
SUCCESS=false

until [ $RETRY_COUNT -ge $MAX_RETRIES ]; do
    curl -L -C - --connect-timeout 30 --retry 3 --keepalive-time 60 -O "$FULL_FTP_URL"
    if [ $? -eq 0 ]; then
        SUCCESS=true
        break
    else
        RETRY_COUNT=$((RETRY_COUNT+1))
        echo "傳輸中斷 (Error 18)，正在進行第 $RETRY_COUNT 次續傳重試..."
        sleep 2
    fi
done

if [ "$SUCCESS" = false ]; then
    echo "錯誤：達到最大重試次數，下載失敗。"
    exit 1
fi

gunzip -f "$GZ_FILE"
FNA_FILE="${FILE_STEM}_genomic.fna"

# 5. 索引與環境變數 (維持原邏輯)
echo "執行索引構建 (bwa & samtools)..."
bwa index "$FNA_FILE"
samtools faidx "$FNA_FILE"

ABS_PATH=$(realpath "$FNA_FILE")
while true; do
    read -p "請輸入環境變數名稱 (必須以 Ref 開頭): " VAR_NAME
    [[ "$VAR_NAME" =~ ^Ref ]] && break
    echo "命名錯誤：必須以 Ref 開頭。"
done

echo -e "\n# Ref Genome Entry: $SELECTED_ACC" >> "$BASHRC"
echo "export $VAR_NAME=\"$ABS_PATH\"" >> "$BASHRC"

echo "------------------------------------------------"
echo "處理成功。請執行: source ~/.bashrc"