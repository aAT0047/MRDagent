# MRD-Agent
[![DOI](https://zenodo.org/badge/920480472.svg)](https://doi.org/10.5281/zenodo.15458496)
## ABSTRACT

### Motivation
Minimal residual disease (MRD) detection is critical for cancer treatment and prognostic evaluation. Next-generation sequencing (NGS)-based variant detection from circulating tumour DNA (ctDNA) represents the state-of-the-art approach for MRD assessment. Crucially, MRD detection results must be stable, as fluctuations in detection performance can severely impact clinical decisions and therapeutic evaluations, yet maintaining stable detection remains challenging. Specifically, ctDNA samples exhibit substantial heterogeneity due to mixtures of normal DNA and multiple tumour subclones, leading to pronounced inter- and intra-sample variability. Consequently, the tools relying on fixed parameter configurations experience severe performance fluctuations across different genomic regions and samples. Although dividing genomic regions into heterogeneous subsets for separate optimisation can theoretically stabilise detection, manual parameter tuning is impractical due to numerous genomic intervals and highly coupled parameters across variant detection and filtering stages. Adaptive optimisation also faces significant challenges, including gradient-free objectives and unknown, dynamically changing constraints. 

### Results
In this study, we propose MRD-Agent, a novel variant detection tool designed specifically for MRD detection. MRD-Agent incorporates an iterative and self-adaptive optimisation framework capable of handling unknown objectives, varying constraints, and highly coupled parameters across stages. A key innovation of MRD-Agent is the integration of a convolutional neural network (CNN)-based meta-model, trained on historical data to enable rapid parameter prediction. This significantly enhances computational efficiency and generalisation performance. Extensive evaluations on simulated and real-world datasets demonstrate MRD-Agent’s superior and stable performance, providing an efficient, reliable solution for MRD detection in clinical and high-throughput research applications.

![Figure 1](https://github.com/aAT0047/MRD-Agent/blob/main/pict/figure1.png)
## Installation

### Prerequisites
The following packages and software need to be installed before running **MRD-Agent**:

#### Python Libraries
Make sure you have the following Python libraries installed:
- `os`
- `csv`
- `math`
- `pysam`
- `numpy`
- `random`
- `argparse`
- `matplotlib`
- `sklearn`
- `scipy`
- `skmultilearn`
- `quickgt`
- `pytorch`

You can install the Python dependencies using `pip`:
```bash
pip install pysam numpy matplotlib sklearn scipy skmultilearn quickgt torch
```
## Data Access

### PACA-CA Project
The dataset used in this project is sourced from the **PACA-CA Project**, which can be accessed via a **DACO-authorised account** at the following URL:
[ICGC ARGO Platform](https://platform.icgc-argo.org/).

### Download Instructions
To download the required data files, you can use the `score-client` command-line tool. Below is an example command for downloading and processing the data:
## Reference Genome
at the following URL:[icgc-argo-workflows](https://github.com/icgc-argo-workflows/argo-reference-files).
The reference genome file `GRCh38_hla_decoy_ebv.fa.gz` is required for data processing. You can download it using the following `wget` command:

```bash
wget https://object.genomeinformatics.org/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa.gz
```
## Data Download Example

To download data from the PACA-CA project, use the `score-client` command-line tool. Below is an example command for downloading and converting the data:

```bash
bin/score-client view \
  --manifest /yourpath/Donor_cramTSVs/DO35116.tsv \
  --reference-file /yourpath/GRCh38_hla_decoy_ebv.fa.gz \
  --bed-query /yourpath/icgccram/bed/DO35116.bed \
  --output-dir cramdata \
  --output-format bam

```
# MRD-Agent workflow
MRD-Agent workflow consists of five key steps：1.Initialize 2. Feature Extractor 3.Adaptive Parameters  4.Meta-model Training

## 1.Initialize 
Split ctDNA Samples & gold VCF files ( **historical data**)

This section describes how to use the **split_genomic_interval.py** script to split BAM and VCF files based on BED regions for each sample, and generate a CSV mapping table.


### Prerequisites

- Python 3  
- samtools  
- bcftools  
- tqdm (`pip install tqdm`)  
- Example directory structure:
  
```bash
python split_genomic_interval.py \
-i /yourpath/ctDNA \
-v /yourpath/VCF \
-b /yourpath/bed \
-B /yourpath/split_bam \
-V /yourpath/split_vcf \
-D /yourpath/split_bed \
-m /yourpath/mapping.csv \
-t 8 \
-r 10
```
| SampleID | RegionIndex | NewFileName  | BEDRegion      | RegionFile                           |
|----------|-------------|--------------|----------------|--------------------------------------|
| sample1  | 1           | sample1_sv_1 | chr1:1000-2000 | /yourpath/split_bed/sample1_sv_1.bed |
| …        | …           | …            | …              | …                                    |

## 2  Feature Extractor
  For each BEDregion, we extracted meta-features X (Supplementary table s3) capturing inter-sample heterogeneity from corresponding BAM files
### 2.1Indel Feature Extractor

This script extracts insertion/deletion (indel) features from BAM files and outputs per-sample CSVs as well as a merged summary.

---

#### Prerequisites

- Python 3.6+  
- [samtools](http://www.htslib.org/download/)  
- [tqdm](https://github.com/tqdm/tqdm) (`pip install tqdm`)  

---
```bash
python extract_indel_features.py \
  -i /your/path/to/bam_files \
  -o /your/path/to/output_dir \
  -w 8 \
  -t 4
```
### 2.2 SNV Feature Extractor

Extract single-nucleotide variant (SNV) features from BAM files per sample and produce a merged summary CSV.

---

#### Prerequisites

- Python 3.6+  
- [samtools](http://www.htslib.org/download/)  
- [bedtools](https://bedtools.readthedocs.io/)  
- [pandas](https://pandas.pydata.org/) (`pip install pandas`)  

---
```bash
python extract_snv_features.py \
  -i /path/to/bam_dir \
  -b /path/to/bed_dir \
  -r /path/to/reference.fa \
  -g /path/to/germline_snp.bed.gz \
  -t 8 \
  -o /path/to/output/all_snv_features.csv
```
 indel & SNV Meta-feature Explanation(X)
| Meta Feature                               | Explanation                                                                                              |
|--------------------------------------------|----------------------------------------------------------------------------------------------------------|
| Reference Allele Count                     | The count of reads that support the reference allele at the variant site.                                |
| Non-reference Allele Count                 | The count of reads that support the non-reference (alternative) allele at the variant site.              |
| Sum of Base Qualities                      | The sum of base quality scores for all reads at the variant site, indicating sequencing confidence.     |
| Tumour Coverage                            | The total number of reads covering the variant site in tumour samples.                                   |
| Mapping Quality                            | The average mapping quality of reads at the variant site, reflecting alignment accuracy.                 |
| Median Read Position                       | The median position of the variant within the sequencing reads, used to assess positional bias.          |
| Homopolymer Rate                           | The proportion of homopolymer regions (repeated bases) near the variant site.                            |
| GC Content                                 | The percentage of guanine (G) and cytosine (C) nucleotides in the surrounding sequence.                  |
| Trinucleotide Sequence Counts              | The frequency of specific trinucleotide sequences near the variant site.                                 |
| Trinucleotide Total Count                  | The total count of trinucleotide sequences in the region surrounding the variant site.                   |
| Distance to Germline SNP                   | The distance between the detected variant and the nearest known germline single nucleotide polymorphism. |
| Repeat Percentage                          | The proportion of repeated sequences in the region around the variant site.                              |
| Short SV Percentage                        | The percentage of short structural variants (SVs) in the region around the variant.                      |
| Middle SV Percentage                       | The percentage of medium-sized structural variants (SVs) in the region around the variant.               |
| Long SV Percentage                         | The percentage of long structural variants (SVs) in the region around the variant.                       |
| Average Read Length                        | The average length of sequencing reads covering the variant site.                                        |
| Small Gap Account                          | The number of small gaps (insertions or deletions) near the variant site.                                |
| RMB (Read Mismatch Bias)                   | The bias in the occurrence of mismatched reads at the variant site.                                      |
| HMDP (High-Quality Mismatch Density %)     | The density of high-quality mismatches in the surrounding sequence.                                      |
| Average Depth                              | The average sequencing depth at the variant site.                                                        |

## 3.Adaptive Parameters in DQN & ADMM

MRD-Agent supports a comprehensive set of **discrete**, **continuous**, and **filtering** parameters for adaptive optimization during the MRD detection process. Below is a categorized list of these parameters.

---

### 3.1 Calling Parameters(y1)
The following discrete parameters control specific thresholds and settings for variant calling:

| Parameter                                           | Explanation                                                                                       | Range      | Type     |
|-----------------------------------------------------|---------------------------------------------------------------------------------------------------|------------|----------|
| base-quality-score-threshold (T)                    | Minimum base quality score required for a base to be considered.                                  | (6, 25)    | Discrete |
| callable-depth (T)                                  | Minimum read depth required to call a site.                                                       | (5, 101)   | Discrete |
| f1r2-median-mq (T)                                  | Median mapping quality of F1R2 reads.                                                             | (30, 71)   | Discrete |
| f1r2-min-bq (T)                                     | Minimum base quality of F1R2 reads.                                                               | (6, 31)    | Discrete |
| max-reads-per-alignment-start (T)                   | Maximum number of reads allowed per alignment start.                                              | (0, 5001)  | Discrete |
| pcr-indel-qual (T)                                  | Quality score threshold for PCR indel filtering.                                                  | [10, 60]   | Discrete |
| pcr-snv-qual (T)                                    | Quality score threshold for PCR SNV filtering.                                                    | [10, 60]   | Discrete |
| assembly-region-padding (T)                         | Padding size added to assembly regions.                                                            | (50, 2001) | Discrete |
| kmer-size (first instance) (T)                      | Size of kmers used in assembly (small regions).                                                   | (5, 25)    | Discrete |
| kmer-size (second instance) (T)                     | Size of kmers used in assembly (large regions).                                                   | (25, 50)   | Discrete |
| max-assembly-region-size (T)                        | Maximum size of assembly regions.                                                                 | (200, 001) | Discrete |
| max-prob-propagation-distance (T)                   | Maximum distance for probability propagation.                                                     | (40, 301)  | Discrete |
| min-assembly-region-size (T)                        | Minimum size of assembly regions.                                                                 | (30, 151)  | Discrete |
| max-unpruned-variants (T)                           | Maximum number of unpruned variants allowed.                                                      | (50, 501)  | Discrete |
| min-dangling-branch-length (T)                      | Minimum length of dangling branches allowed in assembly graph.                                     | (2, 10)    | Discrete |
| phred-scaled-global-read-mismapping-rate (T)        | Phred-scaled global mismapping rate threshold for reads.                                           | (30, 51)   | Discrete |
| pair-hmm-gap-continuation-penalty (T)               | Penalty for gaps in PairHMM calculations.                                                         | (6, 15)    | Discrete |
| mbq (T)                                             | Minimum base quality for read bases to be considered.                                              | (6, 20)    | Discrete |

| Parameter                                                  | Explanation                                                      | Range           | Type       |
|------------------------------------------------------------|------------------------------------------------------------------|-----------------|------------|
| init-lod (X)                                               | Initial LOD threshold for variant calling.                       | (0.5, 3.5)      | Continuous |
| max-af (X)                                                 | Maximum allowed allele frequency for variants.                   | (0.005, 0.05)   | Continuous |
| emit-lod (X)                                               | LOD threshold for emitting potential variants.                   | (1.5, 3.0)      | Continuous |
| active-probability-threshold (X)                           | Probability threshold for activating a region.                   | (0.0005, 0.01)  | Continuous |
| adaptive-pruning-initial-error-rate (X)                    | Initial error rate for adaptive pruning in assembly.             | (0.0005, 0.005) | Continuous |
| pruning-lod-threshold (X)                                  | LOD threshold for pruning assembly graph.                        | (2, 5)          | Continuous |
| flow-probability-threshold (X)                             | Threshold for flow-based probability.                            | (0.002, 0.01)   | Continuous |
| expected-mismatch-rate-for-read-disqualification (X)       | Expected mismatch rate for disqualifying reads.                  | (0.01, 0.05)    | Continuous |
| min-AF (X)                                                 | Minimum allele frequency threshold for variant calling.          | (1e-5, 1e-3)    | Continuous |

---

### 3.2 Filter Parameters(y2)
These parameters are required for **FilterMutectCalls_objectivelast** and control post-calling filtering steps:

| Parameter                                        | Explanation                                                                                               | Range           | Type       |
|--------------------------------------------------|-----------------------------------------------------------------------------------------------------------|-----------------|------------|
| distance_on_haplotype (X)                        | Maximum allowable distance between variants on the same haplotype for consideration.                      | (10, 200)       | Continuous |
| f_score_beta (X)                                 | Weighting factor for precision-recall balance in F1-score calculation.                                    | (0.5, 2.0)      | Continuous |
| false_discovery_rate (X)                         | Threshold for controlling the expected proportion of false discoveries among the identified variants.     | (0.001, 0.2)    | Continuous |
| initial_threshold (X)                            | Initial confidence threshold for variant calling.                                                         | (0.005, 0.2)    | Continuous |
| log_artifact_prior (X)                           | Log-scaled prior probability of artifacts in sequencing data.                                             | (-5.0, -1.5)    | Continuous |
| log_indel_prior (X)                              | Log-scaled prior probability of indels in sequencing data.                                                | (-20.0, -10.0)  | Continuous |
| log_snv_prior (X)                                | Log-scaled prior probability of single-nucleotide variants (SNVs) in sequencing data.                     | (-20.0, -10.0)  | Continuous |
| max_events_in_region (T)                         | Maximum number of variant events allowed in a specific genomic region.                                    | (2, 10)         | Discrete   |
| min_slippage_length (T)                          | Minimum length of homopolymer repeats required to identify PCR slippage events.                           | (4, 12)         | Discrete   |
| pcr_slippage_rate (X)                            | Estimated rate of PCR slippage in sequencing data.                                                        | (0.005, 0.3)    | Continuous |


---



### 3.3 Usage Guide for MRD-Agent

Below is an example usage script for running MRD-Agent using SLURM for parallel task management. This script processes BAM files by splitting them into groups, dynamically assigns tasks to SLURM job arrays, and executes the main MRD-Agent Python program.

```bash
#!/bin/bash
#SBATCH -N 1                       # Request 1 node per task
#SBATCH -n 45                      # Maximum of 45 parallel tasks per job
#SBATCH --ntasks-per-node=45       # Run 45 tasks per node
#SBATCH --partition=9242           # Use partition 9242
#SBATCH --output=%j_%a.out         # Output log file
#SBATCH --error=%j_%a.err          # Error log file
#SBATCH --array=0-10               # Job array index range (11 tasks)

# Configuration paths
FOLDER_PATH="yourpath/split_bam"
OUTPUT_BASE_PATH="yourpath/ICGCCRAM/runbam/"
TASK_OUTPUT_PATH="${OUTPUT_BASE_PATH}group_${SLURM_ARRAY_TASK_ID}/"

# Dynamically create output directory
mkdir -p "${TASK_OUTPUT_PATH}"

# Retrieve all BAM file names
bam_files=($(ls ${FOLDER_PATH}/*.bam | xargs -n 1 basename))
total_files=${#bam_files[@]}

# Grouping logic
total_groups=11  # Total number of groups, matching SLURM_ARRAY_TASK_ID range
group_size=$((total_files / total_groups))
remainder=$((total_files % total_groups))

if [ "$SLURM_ARRAY_TASK_ID" -lt "$remainder" ]; then
    start_idx=$((SLURM_ARRAY_TASK_ID * (group_size + 1)))
    end_idx=$((start_idx + group_size + 1))
else
    start_idx=$((SLURM_ARRAY_TASK_ID * group_size + remainder))
    end_idx=$((start_idx + group_size))
fi

# Extract BAM files for the current group
group_files=("${bam_files[@]:$start_idx:$((end_idx - start_idx))}")

# Save file names to a TXT file (no paths)
bam_files_list="${TASK_OUTPUT_PATH}bam_files_group_${SLURM_ARRAY_TASK_ID}.txt"
printf "%s\n" "${group_files[@]}" > "${bam_files_list}"

# Execute main program, passing the TXT file path
python https://github.com/aAT0047/MRD-Agent/src_code/main.py \
    --bam_files_path "${bam_files_list}" \
    --output_path "${TASK_OUTPUT_PATH}"
```

## 4. Meta-model Training
We viewed parameter recommendation for each genomic region as distinct yet related "tasks". For each task, we extracted meta-features X (Supplementary table s3) capturing inter-sample heterogeneity from corresponding BAM files and combined these with the optimal parameter configurations y, thus forming a meta-dataset. The CNN-based meta-learning model was then trained on this dataset, explicitly learning the common mapping between meta-features and optimal parameters across tasks (Supplementary Methods.3)
```bash
python metamodel.py --features "metamodel/Supplementary Table3.xlsx" --labels "metamodel//Supplementary Table5_sorted.xlsx"
```

### Why an Agent is Needed to Balance Constrained Conditions in Variant Detection
We used MRD-Agent to conduct variant detection. MRD-Agent adopts a DQN framework to facilitate iterative parameter optimisation between the preliminary detection and filtering stages. By dynamically balancing the constraints on false negatives and false positives and integrating stepwise optimisation, MRD-Agent seeks to maximise the overall detection performance.

![Figure 2](https://github.com/aAT0047/MRD-Agent/blob/main/pict/figure4.png)

### Why a Meta-Model is Needed for Rapid Parameter Configuration
As previously described, MRD-Agent’s iterative interactions between the agent and ADMM module result in exponentially increased computational complexity. Moreover, due to the absence of an explicit loss function for gradient calculation, parameter optimisation relies solely on evaluation against a gold-standard reference, significantly restricting its applicability to new ctDNA samples.
, we trained a CNN-based meta-model using historical samples (X: Supplementary Table 3; y:  Supplementary Table 5), We conducted 10-fold cross-validation on 474 samples, using 427 for training and 47 for testing in each fold.


![Figure 5](https://github.com/aAT0047/MRD-Agent/blob/main/pict/figure5.png)

