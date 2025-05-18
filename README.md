# MRD-Agent

## ABSTRACT

### Motivation
Minimal residual disease (MRD) detection holds significant value in the treatment and prognostic evaluation of cancer. However, the underlying circulating tumour DNA (ctDNA) mutations exhibit high heterogeneity both inter- and intra-samples, necessitating customised parameter configurations to achieve stable detection performance. Additionally, small variant detection pipelines commonly encompass at least two key steps: an upper-layer preliminary detection step and a lower-layer filtering step, each involving multiple parameters that are intricately coupled across steps. Overly stringent thresholds in the upper-layer detection step markedly increase the false-negative rate, whereas overly permissive thresholds complicate the subsequent filtering step and elevate false positives. Existing parameter optimisation approaches often fail to incorporate a unified dynamic mechanism to balance false-negative and false-positive rates under real-time feedback conditions, thereby impeding stable ctDNA-based MRD detection.

### Results
To address these challenges, we propose MRD-Agent, which integrates an ADMM (Alternating Direction Method of Multipliers) module into a DQN (Deep Q-Network) framework. This design enables the agent to dynamically adjust multiple constraints across the detection and filtering steps, achieving optimal parameter configurations through iterative optimisation. Furthermore, a meta-learning model is introduced to facilitate rapid parameter recommendations for highly heterogeneous samples. Experiments on both real-world and simulated datasets demonstrate that MRD-Agent achieves superior stability compared with existing software, with a root mean square error (RMSE) of 6.2% and a variance of 0.47% based on F1-measure. These findings indicate that MRD-Agent effectively balances high sensitivity and high specificity for large-scale detection tasks, substantially enhancing the stability of MRD detection.

![Figure 1](https://github.com/aAT0047/MRD-Agent/blob/main/pict/figure1.png)
## Installation

### Prerequisites
The following packages and software need to be installed before running **TMBstable**:

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


## Verifying Installation

After successfully installing the required libraries and tools, verify their versions to ensure compatibility. Use the following commands to check:

```bash
# Verify Python library installations
python -c "import pysam, numpy, matplotlib, sklearn, scipy, skmultilearn, quickgt, torch; print('All libraries are installed correctly')"
```
# Verify external tools
samtools --version
bedtools --version
bcftools --version
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
## Split ctDNA Samples & gold VCF files ( historical data)

This section describes how to use the `split_by_bed.py` script to split BAM and VCF files based on BED regions for each sample, and generate a CSV mapping table.


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

## Indel Feature Extractor

This script extracts insertion/deletion (indel) features from BAM files and outputs per-sample CSVs as well as a merged summary.

---

### Prerequisites

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

## Adaptive Parameters in MRD-Agent

MRD-Agent supports a comprehensive set of **discrete**, **continuous**, and **filtering** parameters for adaptive optimization during the MRD detection process. Below is a categorized list of these parameters.

---

### **Calling Discrete Parameters**
The following discrete parameters control specific thresholds and settings for variant calling:

- `base-quality-score-threshold`
- `callable-depth`
- `f1r2-median-mq`
- `f1r2-min-bq`
- `max-reads-per-alignment-start`
- `pcr-indel-qual`
- `pcr-snv-qual`
- `assembly-region-padding`
- `kmer-size`
- `max-assembly-region-size`
- `max-prob-propagation-distance`
- `min-assembly-region-size`
- `max-unpruned-variants`
- `min-dangling-branch-length`
- `phred-scaled-global-read-mismapping-rate`
- `pair-hmm-gap-continuation-penalty`
- `mbq`

---

### **Calling Continuous Parameters**
The following continuous parameters are used for dynamic threshold adjustments and error handling:

- `init-lod`
- `max-af`
- `emit-lod`
- `active-probability-threshold`
- `adaptive-pruning-initial-error-rate`
- `pruning-lod-threshold`
- `flow-probability-threshold`
- `expected-mismatch-rate-for-read-disqualification`
- `min-AF`

---

### **Filter Parameters**
These parameters are required for **FilterMutectCalls_objectivelast** and control post-calling filtering steps:

- `distance_on_haplotype`
- `f_score_beta`
- `false_discovery_rate`
- `initial_threshold`
- `log_artifact_prior`
- `log_indel_prior`
- `log_snv_prior`
- `max_events_in_region`
- `min_slippage_length`
- `pcr_slippage_rate`

---



## Usage Guide for MRD-Agent

Below is an example usage script for running MRD-Agent using SLURM for parallel task management. This script processes BAM files by splitting them into groups, dynamically assigns tasks to SLURM job arrays, and executes the main MRD-Agent Python program.

---

### SLURM Batch Script Example

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

### Why an Agent is Needed to Balance Constrained Conditions in Variant Detection
We used MRD-Agent to conduct variant detection. MRD-Agent adopts a DQN framework to facilitate iterative parameter optimisation between the preliminary detection and filtering stages. By dynamically balancing the constraints on false negatives and false positives and integrating stepwise optimisation, MRD-Agent seeks to maximise the overall detection performance.

![Figure 2](https://github.com/aAT0047/MRD-Agent/blob/main/pict/figure4.png)

### Why a Meta-Model is Needed for Rapid Parameter Configuration
```bash
python https://github.com/aAT0047/MRD-Agent/metafeature/metafeatureindel.py
python https://github.com/aAT0047/MRD-Agent/metafeature/metafeaturesnv.py
```

![Figure 5](https://github.com/aAT0047/MRD-Agent/blob/main/pict/figure5.png)

