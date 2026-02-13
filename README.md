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

### Overview

This stage performs quality filtering and adapter trimming of paired end Illumina reads using **fastp**.
Processing is parallelised across samples using GNU parallel.

The goal of this stage is to:

* Remove low quality bases
* Remove short reads
* Standardise read quality prior to mapping
* Generate per sample QC reports

---

### Software Used

* fastp
* GNU parallel

---

### Input

Raw paired end reads:

```
RAW_PATH/*_R1_001.fast*
RAW_PATH/*_R2_001.fast*
```

---

### Core fastp Command

```bash
fastp \
  -i R1.fastq.gz \
  -I R2.fastq.gz \
  -o trimmed_R1.fastq.gz \
  -O trimmed_R2.fastq.gz \
  --thread 2 \
  --qualified_quality_phred 30 \
  --length_required 80 \
  --html sample.html \
  --json sample.json
```

---

### Parameter Explanation

| Parameter                      | Meaning                                                                                  |
| ------------------------------ | ---------------------------------------------------------------------------------------- |
| `--thread 2`                   | Two threads per sample. Parallelisation is handled at the sample level via GNU parallel. |
| `--qualified_quality_phred 30` | Bases below Q30 are considered low quality.                                              |
| `--length_required 80`         | Reads shorter than 80 bp after trimming are discarded.                                   |
| `--html`                       | Generates interactive QC report.                                                         |
| `--json`                       | Generates machine readable QC summary.                                                   |

---

### Parallelisation Strategy

Sample level parallelisation is used:

```
parallel -j $JOBS
```

Each fastp process runs with 2 threads.
This prevents CPU oversubscription and improves throughput when processing many samples.

---

### Output

| Directory              | Contents                            |
| ---------------------- | ----------------------------------- |
| `STAGE1/trim/`         | Trimmed paired end FASTQ files      |
| `STAGE1/fastp_report/` | Per sample HTML and JSON QC reports |

Trimmed reads:

* retain only high quality bases
* are at least 80 bp long
* are ready for downstream alignment

---

### Design Rationale

* Q30 threshold reduces sequencing error propagation into mapping and SNP calling.
* 80 bp minimum length maintains reliability of BWA MEM alignment.
* Sample level parallelisation maximises throughput while maintaining reproducible per sample processing.


==================================================================
## Stage 2 – Read Mapping to Reference Genome

### Overview

This stage performs paired-end read alignment against a reference genome using **BWA-MEM**, followed by BAM conversion, filtering, sorting, indexing, and mapping statistics summarisation.

The workflow is designed to:

1. Align paired trimmed FASTQ files
2. Remove unmapped reads
3. Produce coordinate-sorted BAM files
4. Generate mapping quality summary reports

---

### Software Used

* `bwa mem`
* `samtools`

---

### Input

* Trimmed paired-end reads

  ```
  STAGE1/trim/*_R1_001.fastq.gz
  STAGE1/trim/*_R2_001.fastq.gz
  ```

* Reference genome

  ```
  $REF_GENOME
  ```

---

### Core Alignment Command

```bash
bwa mem -t $THREADS $REF_GENOME R1.fastq.gz R2.fastq.gz
```

#### Parameters

| Parameter     | Meaning                                  |
| ------------- | ---------------------------------------- |
| `-t $THREADS` | Number of CPU threads used for alignment |
| `$REF_GENOME` | Indexed reference genome                 |
| `R1/R2`       | Paired-end trimmed reads                 |

`bwa mem` performs gapped alignment using the Burrows-Wheeler Transform and is suitable for reads ≥70 bp.

---

### Processing Steps

The alignment output is processed using `samtools` as follows:

#### 1. Convert SAM to BAM

```
samtools view -bS
```

#### 2. Remove unmapped reads

```
samtools view -bF 4
```

* `-F 4` removes reads flagged as unmapped

#### 3. Sort BAM by coordinate

```
samtools sort
```

#### 4. Index final BAM

```
samtools index
```

#### 5. Generate mapping statistics

```
samtools flagstat
```

This produces per-sample mapping summaries including:

* Total reads
* Mapped reads
* Properly paired reads
* Singleton reads

---

### Output

| Directory                 | Contents                           |
| ------------------------- | ---------------------------------- |
| `STAGE2/mapped_bam/`      | Final sorted and indexed BAM files |
| `STAGE2/mapping_results/` | Mapping statistics reports         |

Each final BAM file:

* contains only mapped reads
* is coordinate sorted
* is indexed for downstream SNP calling and coverage analysis

---

### Design Rationale

* Intermediate SAM files are avoided to reduce disk usage.
* Unmapped reads are removed to ensure cleaner downstream variant analysis.
* Sorting and indexing are required for most population genomic tools.
* Flagstat summaries allow QC-based sample filtering before SNP calling.

==================================================================


## Stage 3 – Mapping Summary Generation and PCA-based Sample Quality Filtering

Stage 3 serves as a rigorous multivariate quality control (QC) checkpoint. The primary objective is to identify and isolate **technical outliers**—samples that successfully bypassed initial sequence trimming but display anomalous alignment behavior. These anomalies typically stem from library preparation artifacts, DNA degradation, or non-target contamination.

### 1. Overview and Input Requirements

The workflow utilizes high-dimensional mapping statistics to distinguish biological variation from technical noise.

* **Input Files**:
* **Coordinate-sorted BAM files**: Located in `STAGE2/mapped_bam/*.bam`.
* **Alignment Statistics**: Individual text files generated via `samtools flagstat` located in `STAGE2/mapping_results/*.txt`.



---

### 2. Step 1: Mapping Summary Consolidation

A comprehensive summary table, `${PROJECT_NAME}_mapping_summary.csv`, is generated by extracting 14 distinct metrics from each `flagstat` output:

| Category | Extracted Metrics |
| --- | --- |
| **Primary Alignment** | Total reads, Mapped reads, Properly_Paired |
| **Pairing Status** | With_Mate_Mapped, Singletons, Paired_in_Seq |
| **Cross-Chr Mapping** | Mate_Diff_Chr, Mate_Diff_Chr_MapQ5 |
| **Secondary/Filtering** | Secondary, Supplementary, Duplicates |
| **Read Identification** | Read1, Read2 |

---

### 3. Step 2: Derived Quality Metrics

To eliminate bias introduced by variable sequencing depth, the pipeline computes five standardized ratios in R. Analysis proceeds only for samples with complete data across these metrics:

* **`mapping_rate`**: Measures overall alignment success.


* **`properly_paired_rate`**: Reflects library structural integrity.


* **`pairing_efficiency`**: Measures the proportion of mapped reads that are correctly oriented.


* **`singleton_rate`**: Indicates DNA fragmentation or adapter interference.


* **`logTotal`**: Normalizes total read counts for cross-magnitude observation.



---

### 4. Step 3: Principal Component Analysis (PCA) Implementation

Multivariate assessment is performed using `mapping_rate`, `pairing_efficiency`, and `singleton_rate`. These ratios characterize the technical "signature" of a sample independent of absolute read counts.

* **Computation**: The R function `prcomp(center = TRUE, scale. = TRUE)` is utilized.
* **Centering**: Shifts data to a zero mean.
* **Scaling**: Ensures each metric has unit variance, preventing variables with larger numerical ranges from dominating the principal components.


* **Outputs**: Detailed statistics including Eigenvalues (importance of each PC), Variance Explained, and Feature Loadings are recorded in `${PROJECT_NAME}_PCA_statistics.txt`.

---

### 5. Step 4: Outlier Detection using Mahalanobis Distance

The pipeline identifies outliers in the PC1 and PC2 subspace using the **Mahalanobis distance**, which accounts for the covariance structure between metrics.

* **Thresholding**: A cutoff is defined based on the **F-distribution** at a **95% confidence level**.
* **Visualization**: The diagnostic plot `${PROJECT_NAME}_PCA_Mapping_Quality.pdf` visualizes the results:
* **95% Confidence Ellipse**: Defines the boundary for "normal" technical variation.
* **Loading Vectors**: Blue arrows indicate the direction and influence of each metric on the sample distribution.
* **Outlier Labels**: Samples falling outside the ellipse are flagged in red.



---

### 6. Step 5: Decision Logic and Filtering

Upon completion of the statistical detection, the pipeline offers several operational modes:

* **Interactive Mode**: Users manually review the PCA plot and outlier list to decide on sample removal.
* **Automatic Mode**: The system follows a pre-defined strategy to either retain or remove all flagged outliers for high-throughput execution.
* **Filtering Output**: If outliers are removed, a refined list, `${PROJECT_NAME}_after_pca.bamfile`, is generated. This filtered list serves as the exclusive input for Stage 4 (Clone Detection) and Stage 6 (Final SNP Calling), preventing low-quality data from introducing noise into population genetic inferences.

---

### Stage 3 Output Summary

| File | Description |
| --- | --- |
| `_mapping_summary.csv` | Consolidated table of 14 raw statistics and derived ratios |
| `_PCA_statistics.txt` | Detailed summary of eigenvalues, variance explained, and loadings |
| `_PCA_Mapping_Quality.pdf` | Visual diagnostic plot with confidence ellipses and identified outliers |
| `_outliers.txt` | List of samples statistically identified as outliers |
| `_after_pca.bamfile` | The filtered BAM list prepared for downstream analysis |


==================================================================

## Stage 4 – Clone Detection Using IBS-Based Genetic Distance

```markdown
# 階段四：基於 IBS 遺傳距離的複本鑑定 (Clone Identification)

階段四的主要任務是透過計算全基因組範圍內的相同狀態 (Identity By State, IBS) 距離，偵測並處理資料集中的遺傳重複樣本 (Duplicates) 或生物性無性繁殖複本 (Clones)。在群體遺傳學中，維持樣本的遺傳獨立性是所有下游分析的基礎前提。移除重複個體可有效避免等位基因頻率 (Allele Frequency) 估計偏差、連鎖不平衡 (LD) 虛假升高，以及群體結構分析（如 PCA 或 ADMIXTURE）中因樣本過度代表所產生的偏誤。

---

## 1. 複本鑑定的背景與必要性

在進行大規模群體採樣時，非獨立樣本的出現通常源於以下兩個層面，若不加處理，將會嚴重扭曲研究結論：

### 生物性複本 (Biological Clones)

**常見場景：**  
在具備無性生殖能力的生物（如珊瑚、地衣、或是具有根莖繁殖能力的植物）中，於不同地理位置採得的樣本可能具備完全相同的基因型。

**後果：**  
保留這些樣本會人為地降低群體的有效個體數（Effective population size, Ne），導致對群體遺傳多樣性的高估或對近交係數（Fis）的錯誤推論。

### 技術性重複 (Technical Duplicates)

**常見場景：**  
實驗室建庫過程中針對同一個體進行多次文庫製備，或是定序過程中同一樣本被重複分配至不同的通道（Lane）。

**後果：**  
這些重複數據在統計上是高度相關的，若視為獨立個體，會導致 F-statistics (Fst) 運算出現顯著誤差，並在 PCA 分析中形成無意義的緊密群簇，遮蔽真實的空間遺傳分化。

---

## 2. 第一步：ANGSD IBS 矩陣運算與全參數解析

本流程利用 ANGSD (Analysis of Next Generation Sequencing Data) 的機率模型來評估樣本間的兩兩遺傳相似度。對於 MIG-seq 這類中低深度的定序數據，這比傳統的「硬判斷」(Hard-call) 基因型更具科學性。

以下為腳本中所使用的完整參數解析：

### 基本輸入與品質過濾參數

- `-bam "$BAM_LIST"`：指定輸入的 BAM 檔案路徑清單。  
- `-uniqueOnly 1`：僅保留唯一比對（Unique mapping）的序列，剔除可能導致遺傳距離低估的多處比對序列。  
- `-remove_bads 1`：移除被旗標標記為不合格（如 PCR duplicates、failed QC）的 reads。  
- `-minMapQ 20`：設定最低比對品質門檻，排除比對位置具備高度不確定性的序列。  
- `-minQ 30`：設定最低鹼基品質門檻（Q30），確保參與計算的變異位點具備高定序精度。  

### 統計模型與位點篩選參數

- `-GL 1` (Genotype Likelihood)：採用 SAMtools 的基因型似然模型。此模型不強制判定特定基因型，而是保留鹼基品質所隱含的機率資訊，在處理雜合子（Heterozygote）偵測時比純數值計算更穩健。  
- `-P 1`：指定使用的執行緒數量。  
- `-minInd (70%)`：動態計算的樣本門檻，規定位點必須在 70% 的樣本中存在，防止因缺失數據導致樣本間距離被誤判。  
- `-snp_pval 1e-5`：對每個候選位點進行統計檢定，僅保留 P-value 小於 1e-5 的顯著變異位點。  
- `-minMaf 0.05`：剔除次要等位基因頻率（Minor Allele Frequency）低於 5% 的位點，降低定序隨機誤差產生的偽 SNP 影響。  

### 運算輸出與輸出格式參數

- `-doMajorMinor 1`：自動推斷主、次等位基因（Major and Minor alleles）。  
- `-doMaf 1`：估計等位基因頻率。  
- `-doCounts 1`：計算各鹼基的覆蓋數量。  
- `-makeMatrix 1`：指示 ANGSD 生成兩兩樣本間的距離矩陣。  
- `-doIBS 1`：執行基於相同狀態（Identity By State）的相似度計算。  
- `-doCov 1`：計算樣本間的共變異矩陣（Covariance matrix），用於輔助結構分析。  
- `-doGeno 32`：輸出二進位的基因型機率檔案。  
- `-doPost 1`：計算每個位點的後驗機率（Posterior probabilities）。  
- `-doGlf 2`：輸出 Beagle 格式的基因型似然檔案，供後續 R 語言腳本讀取。  
- `-out`：定義輸出檔案的前綴路徑與檔名。  

---

## 3. 第二步：階層式聚類 (Hierarchical Clustering) 的計算機制

R 腳本將 ANGSD 產出的 `.ibsMat` 視為遺傳相似度矩陣，並將其轉換為距離對象（Distance Object）後執行聚類。

### 平均連鎖法 (Average Linkage / UPGMA) 的統計原理

腳本採用 `method = "ave"`（即 UPGMA 演算法）。其核心公式為計算兩個群集 $u$ 與 $v$ 之間所有成員兩兩遺傳距離的算術平均數：

$$
d(u, v) = \frac{1}{|u||v|} \sum_{x \in u} \sum_{y \in v} d(x, y)
$$

**抗噪能力：**  
相較於「單一連鎖法 (Single Linkage)」易受極端樣本（Outliers）影響產生「鏈狀效應」，平均連鎖法能平滑化 MIG-seq 在低深度位點產生的隨機誤差，確保聚類結果反映整體的遺傳背景。

**樹狀高度意義：**  
在該演算法下，分支融合的高度標誌著兩個群集間的平均遺傳距離。較低的高度代表樣本間具備極高的一致性，通常暗示其為同一個體的複本。

---

## 4. 第三步：分支高度間隙偵測 (Branch Height Gap Detection) 的數學邏輯

本流程實作自動化門檻判定，透過分析聚類樹狀圖中的分支高度（Merge Heights）向量 $H$ 演變，尋找技術誤差與生物變異之間的斷層。

### 統計跳躍點偵測演算法

#### 微分處理 (Gap Analysis)

提取聚類節點的高度向量 $H = \{h_1, h_2, ..., h_{n-1}\}$，並計算其一階差分 $G$，即相鄰融合高度間的間隙：

$$
g_i = h_{i+1} - h_i
$$

#### 雜訊基線定義 (Noise Baseline)

演算法假設樹狀圖底部的初始融合（高度最小的前幾個節點）反映的是樣本間的「測序雜訊」或「極微量體細胞突變」。系統擷取 $G$ 向量的前 5 個元素，計算其均值 ($\mu_{noise}$) 與標準差 ($\sigma_{noise}$)。這定義了當前數據集下「無意義變異」的統計分佈。

#### 三倍標準差判定門檻 (3-Sigma Rule)

系統尋找第一個顯著突破雜訊區間的間隙點。門檻 $T$ 定義為：

$$
T = \mu_{noise} + 3 \times \sigma_{noise}
$$

**物理意義：**  
若兩次融合之間的高度跳躍超過了雜訊區間的三倍標準差，代表此融合已跨越了「技術重複」的範疇，進入了「生物獨立個體」的遺傳距離量級。

#### 邊界約束機制 (Safety Fallback)

若數據集的分佈極度集中導致 $\sigma_{noise}$ 過小，系統會啟動回退機制，強制將門檻 $h$ 限制在 $0.05$ 至 $0.2$ 的物理約束區間內，避免過度過濾。

---

## 5. 第四步：結果產出與資料過濾策略

判定完成後，系統執行自動分群與資料清單的導出。

### 代表樣本選取邏輯

在每個被判定為複本的群集 (Cluster) 中，系統僅保留清單中的第一個樣本。這確保了在後續的群體遺傳分析中，每個獨立的「基因型」僅擁有一次權重。

### 路徑繼承與清單更新

系統產出標註後的 `${PROJECT_NAME}_after_clones.bamfile`。在階段五（LD Pruning）與階段六（正式變異位點標定）啟動時，腳本會自動掃描並優先加載此清單，大幅提升了研究結果的重複性。

---

## 階段四產出檔案總結

| 檔案名稱 | 詳細說明與診斷用途 |
|----------|--------------------|
| `${PROJECT_NAME}_identify_clones.r` | 封裝了 IBS 數據讀取、UPGMA 聚類與三倍標準差跳躍偵測邏輯的 R 腳本。 |
| `${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf` | 呈現所有輸入樣本最原始的遺傳親緣關係分布。 |
| `${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf` | 標註自動化門檻線 (紅線) 的診斷圖，呈現統計判定的臨界位置。 |
| `${PROJECT_NAME}_clones_to_review.txt` | 系統建議移除的重複樣本名單，研究者應核對其採集紀錄。 |
| `${PROJECT_NAME}_after_clones.bamfile` | 最終通過去重過濾的高品質樣本路徑清單。 |

---

## 操作建議

完成分析後，應開啟 `Filtered.pdf` 觀察紅線（門檻）是否精確切分了底部的緊密群簇。若紅線位置與生物學常理不符，應參考 R 輸出的 $\mu$ 與 $\sigma$ 數值微調代碼中的倍數參數。
```


---

## Stage 5：連鎖不平衡過濾與位點清單生成 (LD Pruning & Site Map)

本階段的核心邏輯在於從全基因組位點中，篩選出在遺傳上相互獨立的標記。若不進行此步驟，後續的群體遺傳分析（如 PCA 或 ADMIXTURE）會因為高度連鎖的位點（Linked Sites）而被特定基因組區域過度代表，導致分析偏差。

### 為何需要「二次呼叫」？

分析流程採用「初步找位點  計算相關性  篩選獨立位點  正式呼叫」的兩手策略：

1. **初步 SNP 呼叫**：必須先產生全樣本的遺傳變異位點池（Potential SNPs），才能得知哪些位點在同一條染色體上、距離多遠，以及它們之間的關聯程度。
2. **生成 Site Map (`mc1.sites`)**：將初步找到的位點位置（Chromosome 與 Position）提取出來，作為計算 LD 的座標基礎。
3. **計算與修剪 (Pruning)**：利用 `ngsLD` 計算位點間的相關係數 ，再透過 `prune_graph` 剔除彼此相關性過高的位點。
4. **再次 Index Site**：ANGSD 無法直接讀取純文字格式的位點清單。`angsd sites index` 會將修剪後留下的位點轉換成二進位索引（Binary index），使 Stage 6 的正式呼叫能精確且快速地只針對這些高品質、獨立的位點進行計算。

### LD 計算與修剪參數解析

| 工具 | 參數 | 意義與設定理由 |
| --- | --- | --- |
| **ngsLD** | `--max_kb_dist 50` | 限制只計算距離 50 kb 以內的位點對。MIG-seq 片段通常較短，計算超長距離的 LD 缺乏統計意義且耗費運算資源。 |
| **ngsLD** | `--probs 1` | 使用基因型後驗機率（Posterior probabilities）進行計算，這比直接使用 Hard-call 基因型更能容忍定序深度不足帶來的誤差。 |
| **prune_graph** | `dist <= 10000` | 修剪範圍設定在 10 kb 內。在此物理距離內的位點若具備高相關性，通常被視為同一遺傳連鎖群。 |
| **prune_graph** | `r2 >= 0.5` | 決定「關聯強度」的門檻。當兩位點的  超過 0.5 時，`prune_graph` 會從中選擇一個保留，並剔除另一個，以確保標記的獨立性。 |

---

## 補充模組：NCBI Genome Fetcher (Genome 下載模組)

`RefMIG.sh` 內建的 `run_fetch_genome_module` 函式提供自動化的參考基因組檢索與預處理流程，旨在簡化分析前的環境設定。

### 運作邏輯與功能

* **多樣化檢索**：支援物種名稱（Species name）或生物專案編號（BioProject ID）關鍵字搜尋。透過 NCBI EDirect 套件中的 `esearch` 指令即時查詢資料庫。
* **技術指標篩選**：系統會自動抓取並呈現下列 12 項關鍵數據供使用者評估：
* **Scaffold N50**：評估基因組組裝的連續性。
* **Coverage**：該基因組的定序深度。
* **Total Genome Size**：基因組總長度。
* **Assembly Status**：組裝狀態（如 Chromosome、Scaffold 或 Contig 級別）。


* **自動化配置**：
* **下載與解壓縮**：獲取目標基因組的 FTP 位置後，自動進行下載與 `gunzip`。
* **建立索引**：下載完成後立即執行 `bwa index` 與 `samtools faidx`，省去手動操作時間。
* **環境變數持久化**：程式會將該基因組的絕對路徑寫入系統設定檔（`.bashrc` 或 `.zshrc`），並以 `Ref_` 作為前綴命名變數（例如 `Ref_Acropora_millepora`），方便在以後的分析中直接選取。


==================================================================

## Stage 6 – Final SNP Calling Using LD-Pruned Site List

### Purpose

Perform final SNP calling restricted to LD-pruned loci to generate a high-quality, approximately independent SNP dataset in VCF format.

---

### Prerequisite

LD-pruned site list must exist:

```
STAGE5/LDpruned_snp.sites
```

This file defines the genomic positions retained after LD filtering.

---

## Step 1 – Site-Restricted SNP Calling

ANGSD performs SNP calling only at positions specified in the LD-pruned site map.

```
angsd \
  -sites LDpruned_snp.sites \
  -b BAM_LIST \
  -GL 1 \
  -minInd MIN_IND \
  -minMapQ 20 \
  -minQ 25 \
  -sb_pval 1e-5 \
  -Hetbias_pval 1e-5 \
  -skipTriallelic 1 \
  -snp_pval 1e-5 \
  -minMaf 0.05 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doCounts 1 \
  -doGlf 2 \
  -dosnpstat 1 \
  -doPost 1 \
  -doGeno 8 \
  -doBcf 1 \
  -doHWE 1 \
  -ref REF_GENOME \
  -out ${PROJECT_NAME}_snps_final
```

### Parameter Rationale

`-sites`
Restricts analysis to LD-pruned loci.

`-minInd`
Requires site presence in minimum number of individuals.

`-minMapQ 20` and `-minQ 25`
Mapping and base quality filters.

`-skipTriallelic 1`
Restricts to biallelic SNPs.

`-minMaf 0.05`
Removes very rare alleles.

`-snp_pval 1e-5`
Significance threshold for SNP detection.

`-doBcf 1`
Outputs BCF format for downstream conversion.

---

## Step 2 – Convert BCF to VCF

Final SNP dataset is converted to VCF format:

```
bcftools view -O v \
  -o ${PROJECT_NAME}_snps_final.vcf \
  ${PROJECT_NAME}_snps_final.bcf
```

---

## Output

Final SNP dataset:

```
${PROJECT_NAME}_snps_final.vcf
```

This VCF:

* Contains only LD-pruned loci
* Includes only biallelic SNPs
* Applies mapping, base quality, and allele frequency filters
* Is suitable for PCA, ADMIXTURE, STRUCTURE, FST, and other population genetic analyses

---

## Final Summary

The pipeline produces:

* Quality-controlled samples
* Clone-filtered individuals
* LD-pruned SNP sites
* Final high-confidence SNP dataset in VCF format
