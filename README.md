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

### Overview

This stage summarises mapping statistics from `samtools flagstat` outputs and performs multivariate quality assessment using PCA.
The objective is to detect outlier samples that may represent sequencing artefacts, poor library quality, or biological divergence.

The workflow:

1. Collect per-sample mapping statistics
2. Generate a consolidated summary table
3. Perform PCA based on key mapping quality metrics
4. Detect statistical outliers using Mahalanobis distance
5. Produce filtered BAM list for downstream analysis

---

## Input Files

* Sorted BAM files

  ```
  STAGE2/mapped_bam/*.bam
  ```

* Mapping statistics generated from Stage 2

  ```
  STAGE2/mapping_results/*.txt
  ```

Each mapping statistics file originates from:

```
samtools flagstat sample.bam
```

---

## Step 1 – Generate Mapping Summary Table

A consolidated CSV file is produced:

```
${PROJECT_NAME}_mapping_summary.csv
```

The following 14 statistics are extracted from each `flagstat` file:

* Total
* Mapped
* Properly_Paired
* With_Mate_Mapped
* Singletons
* Mate_Diff_Chr
* Mate_Diff_Chr_MapQ5
* Secondary
* Supplementary
* Duplicates
* Paired_in_Seq
* Read1
* Read2

These values are stored in a structured table for quantitative analysis.

---

## Step 2 – Derived Quality Metrics

The following indices are computed per sample:

| Metric                 | Definition               |
| ---------------------- | ------------------------ |
| `mapping_rate`         | Mapped / Total           |
| `properly_paired_rate` | Properly_Paired / Total  |
| `pairing_efficiency`   | Properly_Paired / Mapped |
| `singleton_rate`       | Singletons / Total       |
| `logTotal`             | log10(Total reads)       |

Only samples with non-missing values across selected metrics are retained.

---

## Step 3 – PCA Analysis

PCA is performed using the following three variables:

* mapping_rate
* pairing_efficiency
* singleton_rate

Computation:

```
prcomp(center = TRUE, scale. = TRUE)
```

Outputs:

* Eigenvalues
* Variance explained per component
* Cumulative variance
* Feature loadings

Statistics are saved to:

```
${PROJECT_NAME}_PCA_statistics.txt
```

---

## Step 4 – Outlier Detection

Outliers are identified using Mahalanobis distance on PC1 and PC2.

Cutoff threshold:

95% confidence threshold from F distribution.

Samples exceeding this threshold are flagged as outliers.

Generated files:

* PCA plot (PDF)

  ```
  ${PROJECT_NAME}_PCA_Mapping_Quality.pdf
  ```

* Outlier list

  ```
  ${PROJECT_NAME}_outliers.txt
  ```

The PCA figure contains:

* Sample points
* 95% confidence ellipse
* Variable loading vectors
* Labels for detected outliers

---

## Step 5 – BAM List Filtering

If outliers are removed, a new BAM list file is generated:

```
${PROJECT_NAME}_after_pca.bamfile
```

This filtered BAM list becomes the input for downstream SNP calling.

Users may choose:

* Remove detected outliers
* Keep all samples
* Exit analysis

---

## Output Summary

| File                       | Description              |
| -------------------------- | ------------------------ |
| `_mapping_summary.csv`     | Mapping statistics table |
| `_PCA_statistics.txt`      | PCA eigenvalue summary   |
| `_PCA_Mapping_Quality.pdf` | PCA visualisation        |
| `_outliers.txt`            | Detected outlier samples |
| `_after_pca.bamfile`       | Filtered BAM list        |

---

## Design Considerations

* PCA reduces multiple mapping quality metrics into interpretable axes.
* Mahalanobis distance accounts for covariance structure.
* Using three mapping-derived ratios minimises bias from raw sequencing depth.
* Filtering at this stage prevents low-quality samples from influencing SNP and population genetic inference.

==================================================================

## Stage 4 – Clone Detection Using IBS-Based Genetic Distance

### Purpose

Detect duplicated or clonally identical samples based on genome-wide Identity By State distance calculated from BAM files.
Samples identified as genetic duplicates can be removed prior to downstream population genetic analysis.

---

### Input

BAM list from previous stage

* `${PROJECT_NAME}_bwa_mapped.bamfile`
* or `${PROJECT_NAME}_after_pca.bamfile`

Each line contains one sorted BAM file.

---

## Step 1 – IBS Matrix Calculation with ANGSD

Genome-wide pairwise IBS distances are computed using ANGSD.

```
angsd \
  -bam BAM_LIST \
  -GL 1 \
  -uniqueOnly 1 \
  -remove_bads 1 \
  -minMapQ 20 \
  -minQ 30 \
  -minInd 70%_of_samples \
  -snp_pval 1e-5 \
  -minMaf 0.05 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doCounts 1 \
  -makeMatrix 1 \
  -doIBS 1 \
  -doCov 1 \
  -doGeno 32 \
  -doPost 1 \
  -doGlf 2 \
  -out ${PROJECT_NAME}_clone_identification
```

### Parameter Rationale

`-GL 1`
Uses genotype likelihood model based on SAMtools.

`-uniqueOnly 1`
Excludes multi-mapping reads.

`-remove_bads 1`
Removes reads flagged as low quality.

`-minMapQ 20`
Minimum mapping quality threshold.

`-minQ 30`
Minimum base quality threshold.

`-minInd`
Site must be present in at least 70 percent of samples.

`-minMaf 0.05`
Reduces low frequency noise variants.

`-makeMatrix 1`
Produces pairwise distance matrix.

`-doIBS 1`
Computes IBS similarity between samples.

Primary output:

```
${PROJECT_NAME}_clone_identification.ibsMat
```

This matrix contains pairwise genetic distances.

---

## Step 2 – Hierarchical Clustering

The IBS distance matrix is read into R and converted into a distance object:

```
hc <- hclust(as.dist(ma), method = "ave")
```

Average linkage clustering is applied to group genetically similar samples.

---

## Step 3 – Automatic Threshold Determination

Branch heights from the dendrogram are examined.

* Differences between consecutive merge heights are calculated.
* Early significant increases relative to background noise are used to determine threshold.
* If no early significant jump is detected, a conservative fallback threshold is applied.

This avoids manual threshold specification.

---

## Step 4 – Clone Group Assignment

Clusters are defined using:

```
cutree(hc, h = threshold)
```

If a cluster contains multiple samples:

* They are treated as potential clones.
* One representative per cluster is retained.
* Remaining samples are written to removal list.

---

## Outputs

Raw clustering dendrogram:

```
${PROJECT_NAME}_Clone_Dendrogram_RAW.pdf
```

Threshold-applied dendrogram:

```
${PROJECT_NAME}_Clone_Dendrogram_Filtered.pdf
```

List of samples suggested for removal:

```
${PROJECT_NAME}_clones_to_review.txt
```

Filtered BAM list:

```
${PROJECT_NAME}_after_clones.bamfile
```

---

## Downstream Usage

If clone removal is applied, the filtered BAM list replaces the previous list and is used for subsequent analyses.
==================================================================

## Stage 5 – LD Pruning and Site Map Generation

### Purpose

Generate a linkage disequilibrium filtered SNP site list for downstream population genetic analyses.
This step identifies genome-wide SNPs and removes linked loci to retain approximately independent markers.

---

### Input

* Filtered BAM list from previous stage
* Reference genome

---

## Step 1 – Preliminary SNP Calling

ANGSD is used to identify polymorphic sites across all individuals.

```
angsd \
  -b BAM_LIST \
  -GL 1 \
  -uniqueOnly 1 \
  -remove_bads 1 \
  -minMapQ 30 \
  -baq 1 \
  -setMinDepth 5 \
  -SNP_pval 1e-6 \
  -skipTriallelic 1 \
  -doHWE 1 \
  -Hetbias_pval 0.00001 \
  -minInd MIN_IND \
  -doMajorMinor 1 \
  -doMaf 1 \
  -dosnpstat 1 \
  -doPost 2 \
  -doGeno 32 \
  -doCounts 1 \
  -ref REF_GENOME \
  -out STAGE5/allsnps
```

### Parameter Rationale

`-minMapQ 30`
Higher mapping quality threshold to reduce false positives.

`-setMinDepth 5`
Minimum total depth across individuals.

`-SNP_pval 1e-6`
Strict SNP significance threshold.

`-skipTriallelic 1`
Keeps only biallelic loci.

`-minInd`
Requires site presence in minimum number of individuals.

`-doHWE` and `-Hetbias_pval`
Remove sites with extreme Hardy Weinberg deviations or heterozygosity bias.

Primary outputs include:

* `.mafs.gz`
* `.geno`
* other SNP statistics files

---

## Step 2 – Extract Polymorphic Site Positions

From the MAF file, chromosome and position columns are extracted:

```
chr   position
```

File generated:

```
STAGE5/mc1.sites
```

This file is used for LD computation.

---

## Step 3 – LD Matrix Calculation

Pairwise linkage disequilibrium is computed using ngsLD.

```
ngsLD \
  --geno allsnps.geno \
  --probs 1 \
  --n_ind N_IND \
  --n_sites N_SITES \
  --max_kb_dist 50 \
  --pos mc1.sites \
  --n_threads THREADS \
  --extend_out 1 \
  --out allsnpsites.LD
```

### Key Settings

`--max_kb_dist 50`
Only compare loci within 50 kb window.

`--probs 1`
Use posterior genotype probabilities.

---

## Step 4 – LD Pruning

Sites are pruned based on r2 threshold.

```
prune_graph \
  --weight-field r2 \
  --weight-filter "dist <=10000 && r2 >= 0.5"
```

Criteria:

* Distance ≤ 10 kb
* r2 ≥ 0.5

From each linked cluster, one representative SNP is retained.

---

## Step 5 – Site Map Formatting and Indexing

Pruned positions are converted into ANGSD site format:

```
chr    position
```

Output file:

```
STAGE5/LDpruned_snp.sites
```

Index is created:

```
angsd sites index LDpruned_snp.sites
```

---

## Outputs

* `LDpruned_snp.sites`
  Final list of LD-filtered SNP positions.

* Indexed site file
  Used for downstream SNP analysis, PCA, admixture, and population structure inference.

---

## Summary

* Genome-wide SNP discovery
* LD computation within 50 kb windows
* Removal of strongly linked loci
* Generation of independent SNP site map




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
