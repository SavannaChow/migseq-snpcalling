# migseq-snpcalling
MIG-seq Genomic Analysis Pipeline (Refgenome mapping)
開發者：Savanna Chow (savanna201@gmail.com) 使用 Gemini 協助開發
分析邏輯基於 https://github.com/jamesfifer/JapanRE
Credit: AllenChen's lab, Biodiversity Research Center, Academia Sinica


分析流程如下：

序列修剪 (Fastp Trimming)
基因組比對 (BWA Alignment)
樣本品質過濾 (PCA Outlier Filtering)
複本樣本鑑定 (Clone Identification)
連鎖不平衡過濾 (LD Pruning & Site Map Generation)
最終變異位點標定 (Final SNP Calling & VCF Output)

==================================================================
## Stage 1 – Read Quality Control and Trimming
# 1 序列修剪 Fastp Trimming

本階段目的是在進入基因組比對之前，將原始定序資料中品質不佳或過短的 reads 移除，使後續比對與變異分析建立在較穩定且一致的資料基礎上。此步驟直接影響 mapping 效率、變異偵測的錯誤率，以及後續 PCA 與 clone 判定的可靠度。

---

## 使用工具

Fastp

---

## 核心目的

1 移除低品質鹼基
2 過濾整體品質不足的 reads
3 移除過短 reads
4 提供品質評估報告

---

## 主要參數與其功能

### 1 `--qualified_quality_phred 30`

此參數設定合格鹼基品質門檻為 Q30。

Q30 代表錯誤機率約為 0.1%。此設定可確保保留下來的 reads 其鹼基錯誤率極低，減少後續 false SNP 的產生。若門檻過低，低品質區域可能造成假性變異；若門檻過高，則可能過度刪減 reads。

---

### 2 `--length_required 80`

此參數設定 reads 修剪後的最小長度為 80 bp。

過短的 reads 容易在基因組中出現多重比對或非特異性比對，會降低比對精確度。設定最低長度可減少 ambiguous mapping，提高後續 SNP calling 的準確性。

---

### 3 `--thread 2`

設定每個樣本使用 2 執行緒處理。

配合外部平行化工具同時處理多個樣本，使整體計算資源分配更穩定，避免單一樣本過度占用 CPU。

---

### 4 `--html` 與 `--json`

Fastp 會同時輸出：

HTML 報告
JSON 數值報告

HTML 用於人工檢視品質分布與 trimming 效果，JSON 可供後續批次統計或自動分析使用。

---


原始 reads 通常包含：

1 末端品質下降區域
2 偶發性錯誤鹼基
3 部分低品質 reads

若直接進入 BWA 比對，低品質區域會導致：

1 mapping rate 不穩定
2 錯誤比對增加
3 SNP calling 產生偏誤
4 PCA 結果出現技術性分群

因此在進入 Alignment 之前先統一品質標準，是整個分析流程穩定性的前提。

---

## 與後續步驟的關聯

高品質的 trimmed reads 可：

1 提高 BWA 比對準確度
2 降低 ANGSD genotype likelihood 計算誤差
3 提升 clone identification 的一致性
4 減少 LD 計算中因低品質位點造成的雜訊

本階段為整個 RefMIG 管線的基礎品質控制步驟，其輸出品質直接決定後續所有統計分析的可信度。



==================================================================
# 2 基因組比對 BWA Alignment

本階段將修剪後的雙端序列比對至參考基因組，產生排序後的 BAM 檔，並進行基本品質與格式處理，以供後續 ANGSD 與各項群體遺傳分析使用。

---

## 使用工具

BWA MEM
SAMtools

---

## 分析目的

1 將每筆 reads 對齊至參考基因組
2 建立標準化的 BAM 檔格式
3 去除低品質比對
4 建立 index 供下游工具快速存取

---

## BWA MEM

BWA MEM 為適用於中長讀長的比對演算法，採用 seed and extend 架構，能有效處理 paired end reads 並自動判定正確配對關係。

### 指令核心概念

1 輸入為 trimming 後的 R1 與 R2 reads
2 參考基因組需事先建立 index
3 輸出為 SAM 格式比對結果

BWA MEM 會：

1 判定 reads 是否唯一比對
2 計算 mapping quality
3 標示 properly paired reads
4 標記 secondary 與 supplementary alignment

---

## SAMtools 處理流程

BWA 輸出的 SAM 檔會經由 SAMtools 進一步處理：

### 1 轉換為 BAM 格式

SAM 轉 BAM 可大幅降低檔案大小並提升處理效率。

---

### 2 排序

依照 reference 座標排序，為後續 index 及變異分析做準備。

---

### 3 建立 index

建立 `.bai` 檔，使工具可快速隨機存取特定位點。

---

### 4 去除低品質比對

通常會設定：

`-q`
最低 mapping quality 門檻

此步驟可移除低可信度比對，降低假變異。

---

### 5 移除非主要比對

可排除：

secondary alignment
supplementary alignment
unmapped reads

確保後續分析基於可靠比對。

---

## Mapping Quality 的意義

Mapping Quality 代表該 reads 比對到該位置的可信度，為對數轉換機率值。

低 mapping quality 代表該 reads 可能在基因組多處皆可比對。

過濾此類 reads 可減少後續 SNP 計算時的雜訊。

---

## Properly Paired Reads

BWA MEM 會根據：

1 insert size 分布
2 方向是否符合預期
3 兩端是否比對至合理區間

判斷 reads 是否為 properly paired。

此資訊可用於 mapping quality 統計與樣本品質控管。

---

## 本階段產出

每個樣本最終產生：

1 排序完成的 BAM 檔
2 對應 index 檔
3 基本 mapping 統計結果

此 BAM 為後續 ANGSD genotype likelihood 計算、PCA、clone identification 及 LD 計算的輸入資料。


==================================================================

# 3 樣本品質過濾 PCA Outlier Filtering

本階段用「比對統計量」做 PCA，找出 mapping 行為明顯偏離大多數樣本的 outlier，輸出建議剔除名單與過濾後樣本清單。這一階段不是用 SNP 或 genotype likelihood 做 PCA，而是用每個樣本的比對品質指標。 

---

## 使用工具

R
主要套件 readr dplyr ggplot2 ggrepel 

---

## 使用的輸入統計量與其意義

R 腳本會從 mapping summary 表中取出下列欄位並轉成數值，再衍生三個核心 PCA 變數 

1 mapping_rate
Mapped / Total
代表總 reads 中成功比對到參考基因組的比例 

2 pairing_efficiency
Properly_Paired / Mapped
代表已比對 reads 中，符合 paired end 配對條件的比例 

3 singleton_rate
Singletons / Total
代表總 reads 中，成對比對失敗而落單的比例 

腳本也會同時計算 properly_paired_rate 與 logTotal，並移除含 NA 的樣本列，確保 PCA 輸入完整。 

---

## PCA 的 R 指令與用途

### 1 PCA 計算

用三個變數建立矩陣 X，做標準化後的 PCA
center TRUE 代表置中
scale. TRUE 代表依變數標準差做縮放 

```r
X <- df_pca %>% select(mapping_rate, pairing_efficiency, singleton_rate) %>% as.data.frame()
pca_cor <- prcomp(X, center = TRUE, scale. = TRUE)
```



### 2 PCA 統計量輸出

以 sdev 計算 eigenvalue，並輸出每個 PC 的解釋變異百分比與累積百分比，供判讀 PC1 與 PC2 是否足以反映主要差異 

```r
eig_val <- pca_cor$sdev^2
prop_var <- eig_val / sum(eig_val)
cum_var <- cumsum(prop_var)
pca_stat <- data.frame(
  PC = paste0("PC", 1:length(eig_val)),
  Eigenvalue = eig_val,
  Variance_Percent = prop_var * 100,
  Cumulative_Percent = cum_var * 100
)
write.table(pca_stat, "${ABS_STAGE3}/${PROJECT_NAME}_PCA_statistics.txt", sep="\t", quote=FALSE, row.names=FALSE)
```



### 3 特徵向量 loadings

提取 rotation 作為各變數在 PC1 與 PC2 上的方向與貢獻，後續會畫成向量疊在 PCA 圖上 

```r
loadings <- as.data.frame(pca_cor$rotation)
loadings$Variable <- rownames(loadings)
```



### 4 Outlier 判定方式

取 PC1 與 PC2 的 scores，使用 Mahalanobis distance 衡量每個點到中心的距離
cutoff 用 F 分布的 95 百分位建立門檻
distance 大於 cutoff 的樣本標為 outlier 

```r
scores <- as.data.frame(pca_cor$x[, 1:2])
scores$Sample <- df_pca$Sample
n <- nrow(scores); p_vars <- 2
cutoff <- qf(0.95, p_vars, n - 1) * (p_vars * (n - 1) / (n - p_vars))
scores$dist_sq <- mahalanobis(scores[, 1:2], center = colMeans(scores[, 1:2]), cov = cov(scores[, 1:2]))
scores$is_outlier <- scores$dist_sq > cutoff
```



### 5 PCA 圖的內容

輸出 PDF 圖包含四種元素 

1 样本點
outlier 與非 outlier 用不同顏色區分 

2 95 百分位橢圓
用常態假設的 ellipse 顯示主要分布範圍 

3 outlier 標籤
只對 outlier 加上 Sample 名稱 

4 特徵向量
將 mapping_rate pairing_efficiency singleton_rate 的 loadings 畫成箭頭並標註變數名
scaling_factor 用來把向量縮放到與資料點相近的尺度  

圖上 PC1 與 PC2 的軸標題會附上解釋變異百分比 

---

## Outlier 名單與保留清單的產生邏輯

1 outliers 名單
把 is_outlier 為 TRUE 的 Sample 輸出成文字檔 

2 建議保留樣本
把 is_outlier 為 FALSE 的 Sample 轉回對應 BAM 路徑清單，輸出 after_pca bamfile 

```r
valid_samples <- scores$Sample[!scores$is_outlier]
write.table(scores$Sample[scores$is_outlier], "${ABS_STAGE3}/${PROJECT_NAME}_outliers.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
to_keep_bams <- original_bams %>% filter(Sample %in% valid_samples) %>% select(FilePath)
write.table(to_keep_bams, "${ABS_STAGE3}/${PROJECT_NAME}_after_pca.bamfile", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
---

## 終端機輸出摘要

R 腳本會在終端機印出 PCA 統計表 pca_stat，並由 shell 端統計原始樣本數與建議保留樣本數，最後列出 outlier 名單。 


==================================================================


Stage 4 做的事是用 ANGSD 先從 BAM 清單估計每對樣本之間的 IBS distance matrix，然後交給 R 用階層式聚類把距離非常接近的一群樣本抓出來，當作潛在 clone 或重複樣本，再輸出去除後的 BAM 清單供後面 Stage 5 和 Stage 6 接著跑。

---

## Stage 4 輸入與前置條件

輸入是一個文字檔 BAM_LIST，每行一個 bam 的絕對或相對路徑。程式在 Stage 4 開始時會先嘗試沿用前面 PCA 篩完的清單，若 Stage 3 沒跑或沒設定到，會退回用 Stage 3 產生的原始清單。

另外你在 Stage 4 會看到一個 minInd 的動態門檻，做法是把當前樣本數乘上 0.7 取整數，代表某個位點至少要有 70% 的樣本有資料，這個位點才會被納入 IBS 計算。

---

## ANGSD 在 Stage 4 做了什麼

Stage 4 的核心指令在這一行。

它的邏輯是先決定要用哪些 reads 和哪些位點，接著以 genotype likelihood 為基礎估計 allele frequency，最後輸出樣本兩兩之間的 IBS 矩陣。

### 1) 輸入與模型

* `-bam "$BAM_LIST"` 以 BAM 清單做多樣本分析輸入。
* `-GL 1` 指定 genotype likelihood 的計算模型，這個是 ANGSD 的選項，常見用法是選一種內建 GL model 來從 BAM 計算 likelihood。
* `-doGlf 2` 會輸出 genotype likelihood 檔案，方便後續流程需要時沿用同一套 likelihood 結果。

### 2) Read 與比對品質篩選

* `-uniqueOnly 1` 只留唯一比對 reads。
* `-remove_bads 1` 移除 ANGSD 定義的 bad reads 類別。
* `-minMapQ 20` 只留 mapping quality 至少 20 的 reads。
* `-minQ 30` 只留 base quality 至少 30 的鹼基。

這一段的效果是把 IBS 的距離估計盡量建立在高可信度的比對與鹼基上，避免低品質 reads 讓兩個不相關樣本看起來更像或更不像。

### 3) 位點篩選與等位基因頻率估計

* `-minInd "$MIN_IND_TMP"` 每個位點至少要有足夠多樣本有資料才留。
* `-snp_pval 1e-5` 以 SNP test 的 p value 做粗篩，偏向只留訊號明顯的變異位點。
* `-minMaf 0.05` 等位基因頻率門檻 0.05，會排掉非常罕見的 allele，降低單次錯誤或污染造成的假訊號。
* `-doMajorMinor 1` 推定 major allele 和 minor allele。
* `-doMaf 1` 估計 MAF。
* `-doCounts 1` 輸出 allele counts。

### 4) IBS matrix 輸出

* `-doIBS 1 -makeMatrix 1` 這兩個一起用，會把 IBS distance 做成矩陣輸出。

輸出檔案的核心是
`04_Clone_Detection/${PROJECT_NAME}_clone_identification.ibsMat` 
後面 R 就是直接讀這個矩陣。

程式也會把 ANGSD 產出的 gzip 中間檔解壓一份留在 Stage 4 目錄，目的是讓 R 和你自己人工檢查時更容易讀。

---

## R 在 Stage 4 怎麼判定 clone

R script 的輸入有兩個
第一個是 `ibsMat` 距離矩陣
第二個是 BAM_LIST 讓它把樣本名稱對回矩陣的 row name 和 column name。

### 1) 距離到樹狀圖

* `hc <- hclust(as.dist(ma), "ave")` 用 average linkage 做階層式聚類。
* 先輸出一張完全不加門檻線的原始 dendrogram，檔名是 `${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf`。

### 2) 自動門檻 threshold_h 的估計方式

R 用的是你寫的 jump detection

* `h <- hc$height` 取每次合併的高度
* `gaps <- diff(h)` 看相鄰高度差，等於找合併高度突然跳很大的一點。
* 把前 5 個 gaps 當成底噪，算 mean 和 sd。
* 找第一個超過 mean 加 3 倍 sd 的 gaps 當 jump_idx。
* 如果 jump_idx 很早出現而且在前 10 個合併內，就用 jump 前後兩個高度的中點當門檻，否則走 fallback 規則，把門檻限制在一個不會太離譜的範圍內。

這段邏輯的意義是讓門檻盡量抓到兩群距離尺度的分界點，也就是極近距離的一群樣本先聚成小團，再往上合併到其他樣本時高度突然變大。

### 3) 分群與輸出名單

* `clusters <- cutree(hc, h = threshold_h)` 直接用門檻高度切樹得到 cluster ID。
* 定義 clone cluster 是任何一個 cluster 內樣本數大於 1。
* `clones_to_review.txt` 的寫法是每個 clone cluster 保留第一個樣本，把同 cluster 的其他樣本列為建議移除。
* `after_clones.bamfile` 的寫法是每個 cluster 只保留一個代表樣本的 FilePath，這就是後續分析要用的去重清單。

### 4) 帶門檻線的診斷圖與終端輸出

* 輸出第二張有紅線門檻的 dendrogram，檔名 `${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf`。
* 終端會印一行底噪高度與自動門檻，方便你在 log 裡快速比對不同 run 的結果。

---

## Stage 4 跑完後你會得到哪些檔案

在 `04_Clone_Detection` 目錄下最重要的幾個

* `${PROJECT_NAME}_clone_identification.ibsMat` 兩兩 IBS distance matrix 
* `${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf` 未加門檻的樹 
* `${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf` 有紅線門檻的樹 
* `${PROJECT_NAME}_clones_to_review.txt` 建議移除名單 
* `${PROJECT_NAME}_after_clones.bamfile` 去重後 BAM 清單 

---

## Stage 4 的決策如何影響後續流程

如果偵測到 clones_to_review.txt 有內容，程式會在自動模式套用你一開始選的策略，或在互動模式問你要不要移除。選擇移除時，它會把 BAM_LIST 指向 after_clones.bamfile，後面 Stage 5 和 Stage 6 就會用去重後的樣本集合。

這裡需要你自己判讀的點是，Stage 4 把距離非常近的樣本當成 clone 候選，它沒有替你判定原因。真實 clone，重複上機，樣本標籤錯置，跨樣本 index hopping 的極端情形，都可能造成非常近的 IBS，最後要靠你對 sample metadata 和 dendrogram 的群內結構做判斷。


==================================================================

# 5 連鎖不平衡過濾 LD Pruning 與 Site Map 產生

本階段使用 ANGSD 產生 LD 計算所需的 genotype likelihood 與位點資訊，接著用 ngsLD 計算 SNP 間的 r²，再透過 prune_graph 過濾高度連鎖的位點，最後產生 LD pruned site list。最終 site map 由 ANGSD 輸出的位點資訊回溯生成。
---

## 分析目的

1 計算 SNP 間的 linkage disequilibrium
2 過濾高度相關之位點
3 建立獨立位點集合
4 產出 LD pruned site map 作為 Stage 6 的輸入

---

## 使用工具

ANGSD
ngsLD
prune_graph

---

## 第一步 由 ANGSD 準備 LD 計算資料

### ANGSD 核心指令與功能

* `-GL 1`
  以 genotype likelihood 模型讀取 BAM

* `-doMajorMinor 1`
  推定 major 與 minor allele

* `-doMaf 1`
  估計 allele frequency

* `-doGlf 2`
  輸出 Beagle 格式 genotype likelihood

* `-minMapQ 20`

* `-minQ 30`

* `-uniqueOnly 1`

* `-minInd`

* `-minMaf 0.05`

上述條件與前面 Stage 3 與 Stage 4 保持一致，確保使用同一品質標準。

### 輸出格式

主要輸出為：

`.beagle.gz`
包含每個樣本在每個 SNP 位點上的三種 genotype likelihood

`.mafs.gz`
包含 allele frequency 與位點資訊

`.pos.gz`
包含位點座標

這些檔案作為 ngsLD 的直接輸入。

---

## 第二步 使用 ngsLD 計算 LD

### ngsLD 計算原理

ngsLD 不需要硬性 genotype calling，而是直接基於 genotype likelihood 計算 SNP 間的 LD。

計算的是每一對 SNP 間的 r²。

在 likelihood 框架下，ngsLD 會：

1 從 Beagle likelihood 計算每個 SNP 的 genotype posterior
2 基於 posterior 推估 SNP1 與 SNP2 間的共變異程度
3 計算標準化 LD 指標 r²

### r² 的定義

r² = 0
兩位點完全不相關

r² 接近 1
兩位點高度相關，等位基因組合幾乎固定一起出現

ngsLD 輸出為每對 SNP 的：

chr1 pos1 chr2 pos2 r² 值

通常僅在距離小於指定 window 的 SNP 對中計算。

---

## 第三步 使用 prune_graph 進行 LD 過濾

### prune_graph 的邏輯

ngsLD 會輸出成對 SNP 的 r² 值。

當 r² 高於指定 threshold 時，代表兩 SNP 存在高度連鎖。

prune_graph 會：

1 將每個 SNP 視為 graph 中的一個節點
2 若兩 SNP 的 r² 超過門檻，就在兩者間建立連線
3 在圖中尋找最大獨立集合
4 從每個高度連鎖的群組中保留一個節點

最終得到彼此之間 LD 低於門檻的 SNP 集合。

---

## 為什麼需要 LD pruning

若不進行 LD pruning：

1 PCA 會被某些高 LD 區域主導
2 STRUCTURE 或 ADMIXTURE 可能放大局部 linkage 訊號
3 某些 genomic 區域權重過重

LD pruning 可使剩餘 SNP 近似統計獨立。

---

## 最終輸出

`LDpruned_snp.sites` 是 Stage 5 的最終位點清單，用來限制 Stage 6 僅在 LD pruning 後保留下來的 SNP 上進行變異計算。此檔案不包含 genotype、allele frequency 或 VCF 結構，只是一個座標索引集合，格式符合 ANGSD `-sites` 參數的讀取規格。

---

## `LDpruned_snp.sites` 的內容格式

每一行代表一個 SNP 位點，包含：

* 染色體或 scaffold 名稱
* 位置座標 pos

例如：

```
chr1    102345
chr1    203876
chr2    99812
```

ANGSD 讀取此檔時，會依照 chr 與 pos 直接定位 BAM 中對應座標，只在這些位點執行後續計算。

---

## `LDpruned_snp.sites` 的產生來源

1. ANGSD 在 Stage 5 前段產生 `.beagle.gz`、`.mafs.gz`、`.pos.gz` 等位點資料
2. ngsLD 使用 Beagle likelihood 計算 SNP 對的 r²
3. prune_graph 依據 r² 門檻挑選彼此低連鎖的一組 SNP
4. 依保留下來的 SNP ID 對照 ANGSD 的位點座標檔 `.pos.gz`
5. 輸出對應 chr 與 pos 為 `LDpruned_snp.sites`

因此座標資訊來源是 ANGSD，篩選邏輯來源是 ngsLD 與 prune_graph。

---

## Stage 6 中 ANGSD 的實際使用方式

Stage 6 會再次呼叫 ANGSD，並以 `-sites` 限定位點集合：

```bash
angsd \
  -sites "$STAGE5/LDpruned_snp.sites" \
  -bam "$BAM_LIST" \
  -GL 1 \
  -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 \
  -minInd "$MIN_IND" \
  -snp_pval 1e-5 -minMaf 0.05 \
  -doMajorMinor 1 -doMaf 1 -doCounts 1 \
  -doGeno 32 -doPost 1 -doGlf 2 \
  -out "$STAGE5/${PROJECT_NAME}_snps_final"
```

---

## `-sites` 在 Stage 6 的作用

`-sites LDpruned_snp.sites` 會使 ANGSD：

1. 僅讀取並處理 `.sites` 檔列出的 chr 與 pos
2. 僅在這些位點計算 genotype likelihood、posterior 與 SNP 結果
3. 輸出結果僅包含這份位點集合

未出現在 `.sites` 檔中的其他位點，即使 BAM 中有 reads 覆蓋，也完全不會被納入計算或輸出。

---

## 對最終 VCF 的影響

由於 Stage 6 僅在 LDpruned_snp.sites 定義的位點集合上進行變異計算，因此最終產生的 SNP 與 Stage 5 的 LD pruning 結果完全一致。

這確保輸出的變異集合經過連鎖不平衡過濾，適合用於需要統計獨立位點假設的下游分析，例如 PCA、STRUCTURE、ADMIXTURE、Fst 計算與群體遺傳推論。
