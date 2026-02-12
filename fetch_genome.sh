#!/bin/bash

# 設定基礎下載目錄
BASE_DIR=~/RefGenome

# 設定隱藏標記檔路徑
ENV_CHECK_FILE="$BASE_DIR/.env_verified"
REQUIRED_CMDS=("esearch" "bwa" "samtools" "curl" "gunzip")

# 環境偵測邏輯：若標記檔不存在則執行偵測
if [[ ! -f "$ENV_CHECK_FILE" ]]; then
    echo "Initializing environment check..."

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
                exit 1
            fi
        fi
    done

    # 建立隱藏標記檔
    touch "$ENV_CHECK_FILE"
    echo "Environment verified and cached."
fi

# 進入搜尋與選擇迴圈
while true; do
    # 1. 使用者輸入關鍵字
    read -p "請輸入搜尋關鍵字 (物種名或 BioProject): " QUERY

    # 2. 執行檢索並暫存結果
    TEMP_DATA=$(mktemp)
    
    # 產生時間戳記與完整存檔路徑
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    HISTORY_FILE="$BASE_DIR/SearchRecord_${QUERY// /_}_${TIMESTAMP}.txt"

echo "正在檢索 NCBI 資料庫並解析數據結構"
esearch -db assembly -query "\"${QUERY}\"[Assembly Accession]" \
| esummary \
| xtract -pattern DocumentSummary -def "NA" \
  -lbl "AssemblyAccession"       -element AssemblyAccession \
  -lbl "AssemblyName"            -element AssemblyName \
  -lbl "AssemblyStatus"          -element AssemblyStatus \
  -lbl "BioprojectAccn"          -block Bioproj   -sep "," -element BioprojectAccn \
  -lbl "SubmitterOrganization"   -element SubmitterOrganization \
  -lbl "Isolate"                 -block Biosource -sep "," -element Isolate \
  -lbl "Coverage"                -element Coverage \
  -lbl "ScaffoldN50"             -element ScaffoldN50 \
  -lbl "AssemblyType"            -element AssemblyType \
  -lbl "RefSeq_category"         -element RefSeq_category \
  -lbl "BioSampleAccn"           -element BioSampleAccn \
  -lbl "SubmissionDate"          -element SubmissionDate \
  -lbl "LastUpdateDate"          -element LastUpdateDate \
  -lbl "FtpPath_GenBank"         -element FtpPath_GenBank \
  -lbl "total_length"            -block Stat -if @category -equals total_length -element Stat \
> "$TEMP_DATA"

    if [ ! -s "$TEMP_DATA" ]; then
        echo "找不到符合結果。"
        rm "$TEMP_DATA"
        continue
    fi

# 1  AssemblyAccession
# 2  AssemblyName
# 3  AssemblyStatus
# 4  BioprojectAccn
# 5  SubmitterOrganization
# 6  Isolate
# 7  Coverage
# 8  ScaffoldN50
# 9  AssemblyType
# 10 RefSeq_category
# 11 BioSampleAccn
# 12 SubmissionDate
# 13 LastUpdateDate
# 14 FtpPath_GenBank
# 15 total_length


    # 3. 視覺化格式輸出並同步存檔至 BASE_DIR
awk -F'\t' '
{
    data[$1] = $2
}
END{
    total_mb = "NA"
    if (data["total_length"] != "NA" && data["total_length"] ~ /^[0-9]+$/) {
        total_mb = data["total_length"] / 1000000
    }

    printf "-------------------- [ Index: %-4d ] --------------------\n", 1
    printf "ID & ACC       | %s | %s\n", data["AssemblyName"], data["AssemblyAccession"]
    printf "SOURCE         | Project: %s | BioSample: %s | Isolate: %s\n", data["BioprojectAccn"], data["BioSampleAccn"], data["Isolate"]
    printf "SPECS          | Status: %s | Type: %s | RefSeq: %s\n", data["AssemblyStatus"], data["AssemblyType"], data["RefSeq_category"]

    if (total_mb == "NA") printf "GENOME         | NA\n"
    else printf "GENOME         | %.0f Mb (total)\n", total_mb

    printf "QUALITY        | ScaffoldN50: %s | Coverage: %s\n", data["ScaffoldN50"], data["Coverage"]
    printf "SUBMITTER      | %s\n", data["SubmitterOrganization"]
    printf "DATE           | Submitted: %s | Updated: %s\n", data["SubmissionDate"], data["LastUpdateDate"]
    printf "FTP            | %s_genomic.fna.gz\n\n", data["FtpPath_GenBank"]
}
' "$TEMP_DATA" | tee "$HISTORY_FILE"

    echo "搜尋結果已存檔至: $HISTORY_FILE"

    # 4. 互動式選擇
    read -p "選擇要下載的Genome編號 (輸入 r 重新搜尋): " INDEX

    # 判斷是否重新搜尋
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
    break # 選擇成功，跳出迴圈進入下載程序
done

# 5. 獲取 FTP 路徑並下載
echo "擷取 $ACCESSION FTP位置..."
FTP_BASE=$(esearch -db assembly -query "$ACCESSION" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)

if [ -z "$FTP_BASE" ]; then
    echo "FTP path not found."
    exit 1
fi

TARGET_DIR="$BASE_DIR/$ACCESSION"
mkdir -p "$TARGET_DIR"

FILE_NAME=$(basename "$FTP_BASE")
DOWNLOAD_FILE="${FILE_NAME}_genomic.fna.gz"
FULL_URL="${FTP_BASE}/$DOWNLOAD_FILE"

echo "下載Genome $ACCESSION 至 $TARGET_DIR..."
wget -c --tries=0 -P "$TARGET_DIR" "$FULL_URL"

# 6. 設定環境變數 (移至解壓縮之前)
echo "目前系統中已定義的變數名稱與路徑 (Ref_XXX):"
if grep -q "^export Ref_" ~/.bashrc; then
    # 使用 sed 去除 "export " 並將第一個 "=" 替換成 " -> " 方便閱讀
    grep "^export Ref_" ~/.bashrc | sed 's/export //g' | sed 's/=/  -->  /g'
else
    echo "(目前尚無設定任何 Ref_ 變數)"
fi
echo "--------------------------------------------------"
echo ""
echo "請輸入參考基因組名稱, 不要跟現有的重複，也不可使用空白或特殊字元"
read -p "REF_" USER_INPUT

# 若使用者未輸入名稱則停止
if [ -z "$USER_INPUT" ]; then
    echo "錯誤：未輸入名稱，停止執行後續解壓縮與索引程序。"
    exit 1
fi
# 自動合成變數名稱
ENV_VAR="Ref_${USER_INPUT}"

# 預先定義解壓後的絕對路徑
FNA_FILE="$TARGET_DIR/${FILE_NAME}_genomic.fna"
ABS_PATH=$(realpath "$FNA_FILE")

# 寫入 .bashrc
echo "export $ENV_VAR=\"$ABS_PATH\"" >> ~/.bashrc
source ~/.bashrc
echo "環境變數 '$ENV_VAR' 已加入 ~/.bashrc。"


# 7. 解壓縮與建立索引 (只有在環境變數設定後才跑)
echo "解壓縮..."
gunzip -f "$TARGET_DIR/$DOWNLOAD_FILE"

echo "建立索引 Building BWA index (this may take a while)..."
bwa index "$ABS_PATH"


echo "--------------------------------------------------"
echo "Process complete."
echo "Genome path: $ABS_PATH"
echo "處理成功。請執行: source ~/.bashrc 使設定生效。"