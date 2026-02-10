#!/bin/bash
# ==============================================================================
# MIG Analysis Full Pipeline - Smart Hybrid Mode (v8.12)
# ==============================================================================
source ~/.bashrc

# ------------------------------------------------------------------------------
# 環境依賴檢查 (Dependency Check)
# ------------------------------------------------------------------------------
ENV_CHECK_FILE=".pipeline_env_ready"

check_dependencies() {
    if [ -f "$ENV_CHECK_FILE" ]; then
        echo "[系統] 偵測到環境檢查標記，跳過軟體驗證。"
        return 0
    fi

    echo "[系統] 正在檢查執行環境依賴項..."
    
    # 定義工具清單：指令名稱 | 安裝指令或網址
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
        echo "[系統] 所有核心套件已就緒。"
        touch "$ENV_CHECK_FILE"
    else
        echo "[警告] 缺少以下套件: ${missing_tools[*]}"
        for m_tool in "${missing_tools[@]}"; do
            local action="${tools[$m_tool]}"
            if [[ "$action" == http* ]]; then
                echo "  - $m_tool: 請手動編譯安裝，參考網址: $action"
                manual_install+=("$m_tool")
            else
                read -p "  - 偵測到 $m_tool 缺失，是否嘗試自動安裝? (y/n): " do_install
                if [[ "$do_install" == "y" ]]; then
                    eval "$action"
                else
                    manual_install+=("$m_tool")
                fi
            fi
        done
    fi

    if [ ${#manual_install[@]} -ne 0 ]; then
        echo "錯誤：請先解決上述手動安裝套件後再運行腳本。"
        exit 1
    fi
}

check_dependencies



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
        read -p "是否跳過此步驟並使用現有結果？ (y:跳過 / n:重新分析): " choice < /dev/tty
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

# 1.2 參考基因組配置
source ~/.bashrc
MAPFILE=()
MAPVAL=()

# 過濾環境變數中以 .fa, .fasta, .fna 結尾的路徑
while IFS='=' read -r name value; do
    if [[ "$value" =~ \.(fa|fasta|fna)$ ]]; then
        MAPFILE+=("$name")
        MAPVAL+=("$value")
    fi
done < <(env)

echo "--- 可用的參考基因組 ---"
for i in "${!MAPFILE[@]}"; do
    echo "$((i+1))) \$${MAPFILE[$i]} (${MAPVAL[$i]})"
done

MANUAL_OPTION=$(( ${#MAPFILE[@]} + 1 ))
echo "$MANUAL_OPTION) 手動輸入絕對路徑"

read -p "請選擇參考基因組 (1-$MANUAL_OPTION): " REF_CHOICE

# 驗證輸入為純數字且在選項範圍內
if [[ "$REF_CHOICE" =~ ^[0-9]+$ ]] && [ "$REF_CHOICE" -ge 1 ] && [ "$REF_CHOICE" -le "${#MAPFILE[@]}" ]; then
    REF_GENOME="${MAPVAL[$((REF_CHOICE-1))]}"
elif [ "$REF_CHOICE" -eq "$MANUAL_OPTION" ]; then
    read -e -p "請輸入絕對路徑: " REF_GENOME
else
    echo "錯誤：無效的選擇。"
    exit 1
fi

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

# 1.3 模式選擇
echo "-------------------------------------------------------"
echo "請選擇分析運行模式："
echo "1) 自動模式 (預先設定篩選策略，遇到決策點不中斷)"
echo "2) 互動模式 (各階段結束後，手動決定是否移除樣本)"
read -p "請輸入選項 (1 或 2): " RUN_MODE

if [ "$RUN_MODE" == "1" ]; then
    echo "[模式選擇：自動模式]"
    read -p "  > 偵測到 PCA Outlier 時處理方式 (1:移除, 2:保留): " AUTO_PCA_CHOICE
    read -p "  > 偵測到 Clone 樣本時處理方式 (1:移除, 2:保留): " AUTO_CLONE_CHOICE
    echo "預設 PCA 決策: $AUTO_PCA_CHOICE"
    echo "預設 Clone 決策: $AUTO_CLONE_CHOICE"
else
    echo "[模式選擇：互動模式]"
fi

# 1.4 階段選擇 (模組化重構版本)
echo "-------------------------------------------------------"
echo "分析流程選擇 (Stage Selector):"
echo "1) 執行完整分析 (Stage 1-6)"
echo "2) 僅執行序列清理與比對參考基因組 (1. Trimming & 2. Alignment)"
echo "3) 僅執行潛在問題樣本過濾分析 (3. PCA Outlier & 4. Clone Filtering)"
echo "4) 僅執行 LD Pruning (產生連鎖不平衡位點表)"
echo "5) 執行最終 SNP Calling (基於現有無LD位點進行)"
echo "6) 自定義流程"
read -p "選擇運行範圍: " RUN_CHOICE

case "$RUN_CHOICE" in
    1) RUN_S1=y; RUN_S2=y; RUN_S3=y; RUN_S4=y; RUN_S5=y; RUN_S6=y ;;
    2) RUN_S1=y; RUN_S2=y; RUN_S3=n; RUN_S4=n; RUN_S5=n; RUN_S6=n ;;
    3) RUN_S1=n; RUN_S2=n; RUN_S3=y; RUN_S4=y; RUN_S5=n; RUN_S6=n ;;
    4) RUN_S1=n; RUN_S2=n; RUN_S3=n; RUN_S4=n; RUN_S5=y; RUN_S6=n ;;
    5) RUN_S1=n; RUN_S2=n; RUN_S3=n; RUN_S4=n; RUN_S5=n; RUN_S6=y ;;
    *) 
       read -p "執行 Stage 1 Trimming? (y/n): " RUN_S1
       read -p "執行 Stage 2 Alignment? (y/n): " RUN_S2
       read -p "執行 Stage 3 PCA Outlier Filtering? (y/n): " RUN_S3
       read -p "執行 Stage 4 Clone Filtering? (y/n): " RUN_S4
       read -p "執行 Stage 5 LD Pruning (Generate Site Map)? (y/n): " RUN_S5
       read -p "執行 Stage 6 Final SNP Calling (Apply Map)? (y/n): " RUN_S6
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

echo "======================================================="
echo "分析啟動時間: $(date)"
echo "專案名稱: $PROJECT_NAME"
echo "======================================================="

THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu)
JOBS=$(( THREADS / 4 )); [ "$JOBS" -lt 1 ] && JOBS=1

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
                read -p "是否移除離群樣本？(1:移除, 2:保留): " PCA_DECISION < /dev/tty
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
                read -p "是否移除上述Clone樣本？(1:移除, 2:保留): " CLONE_DECISION < /dev/tty
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
            echo " q) 手動輸入其他路徑(副檔名必須為bamfile)"
            
            read -p "請輸入選項 (1-${#FOUND_LISTS[@]} 或 q): " FILE_CHOICE < /dev/tty
            
            if [[ "$FILE_CHOICE" =~ ^[0-9]+$ ]] && [ "$FILE_CHOICE" -le "${#FOUND_LISTS[@]}" ]; then
                BAM_LIST="${FOUND_LISTS[$((FILE_CHOICE-1))]}"
            else
                read -e -p "請輸入自訂清單路徑: " BAM_LIST < /dev/tty
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