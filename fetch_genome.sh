#!/bin/bash


# 設定基礎下載目錄
BASE_DIR=~/RefGenome


# 設定隱藏標記檔路徑
ENV_CHECK_FILE="$BASE_DIR/.env_verified"
REQUIRED_CMDS=("esearch" "bwa" "samtools" "curl" "gunzip")

# 若標記檔存在則跳過偵測
if [[ -f "$ENV_CHECK_FILE" ]]; then
    exit 0
fi

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



# 1. 使用者輸入關鍵字
read -p "請輸入搜尋關鍵字 (物種名或 BioProject)" QUERY

# 2. 執行檢索並暫存結果
TEMP_DATA=$(mktemp)

echo "正在檢索 NCBI 資料庫並解析數據結構"
esearch -db assembly -query "$QUERY" -retmax 5000 \
| esummary \
| xtract -pattern DocumentSummary -def "NA" \
  -element AssemblyAccession AssemblyName AssemblyStatus \
           BioprojectAccn SubmitterOrganization Isolate \
           Coverage ScaffoldN50 AssemblyType RefSeq_category BioSampleAccn \
  -block Stat -if @category -equals total_length -element Stat > "$TEMP_DATA"

if [ ! -s "$TEMP_DATA" ]; then
    echo "找不到符合結果."
    rm "$TEMP_DATA"
    exit 1
fi

# 3. 視覺化格式輸出
awk -F'\t' '{
    total_mb     = $12 / 1000000
    ungapped_mb  = $13 / 1000000
    printf "-------------------- [ Index: %-4d ] --------------------\n", NR
    printf "ID & NAME      | %s (%s)\n", $1, $2
    printf "SOURCE         | Project: %s | BioSample: %s | Isolate: %s\n", $4, $11, $6
    printf "SPECS          | Status: %s | Type: %s | RefSeq: %s\n", $3, $9, $10
    printf "GENOME         | %.0f Mb (total) | %.1f Mb (ungapped)\n", total_mb, ungapped_mb
    printf "QUALITY        | ScaffoldN50: %s | Coverage: %s\n", $8, $7
    printf "SUBMITTER      | %s\n\n", $5
}' "$TEMP_DATA"


# 4. 互動式選擇
read -p "選擇要下載的Genome編號" INDEX
SELECTED_LINE=$(sed -n "${INDEX}p" "$TEMP_DATA")

if [ -z "$SELECTED_LINE" ]; then
    echo "Invalid selection."
    rm "$TEMP_DATA"
    exit 1
fi

ACCESSION=$(echo "$SELECTED_LINE" | cut -f1)

# 5. 獲取 FTP 路徑並下載
echo "擷取 $ACCESSION FTP位置..."
FTP_BASE=$(esearch -db assembly -query "$ACCESSION" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)

if [ -z "$FTP_BASE" ]; then
    echo "FTP path not found."
    rm "$TEMP_DATA"
    exit 1
fi

TARGET_DIR="$BASE_DIR/$ACCESSION"
mkdir -p "$TARGET_DIR"

FILE_NAME=$(basename "$FTP_BASE")
DOWNLOAD_FILE="${FILE_NAME}_genomic.fna.gz"
FULL_URL="${FTP_BASE}/$DOWNLOAD_FILE"

echo "下載Genome $ACCESSION 至 $TARGET_DIR..."
wget -c --tries=0 -P "$TARGET_DIR" "$FULL_URL"

# 6. 解壓縮與建立索引
echo "解壓縮..."
gunzip -f "$TARGET_DIR/$DOWNLOAD_FILE"
FNA_FILE="$TARGET_DIR/${FILE_NAME}_genomic.fna"
ABS_PATH=$(realpath "$FNA_FILE")

echo "建立索引 Building BWA index (this may take a while)..."
bwa index "$ABS_PATH"

# 7. 設定環境變數 (強制 Ref_ 字首)
echo "--------------------------------------------------"
read -p "請輸入參考基因組名稱, 不可使用空白或特殊字元 (e.g., Ahya_XXX): " USER_INPUT

# 自動合成變數名稱
ENV_VAR="Ref_${USER_INPUT}"

# 寫入 .bashrc (使用絕對路徑與雙引號包裹)
echo "export $ENV_VAR=\"$ABS_PATH\"" >> ~/.bashrc
source ~/.bashrc

rm "$TEMP_DATA"
echo "--------------------------------------------------"
echo "Process complete."
echo "Genome path: $ABS_PATH"
echo "Environment variable '$ENV_VAR' has been added to ~/.bashrc."
echo "處理成功。請執行: source ~/.bashrc 使設定生效。"


