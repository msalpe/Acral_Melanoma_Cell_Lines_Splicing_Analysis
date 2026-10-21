# Acral Melanoma SF3B4 KO Cell Lines – RNA-seq Alternative Splicing Analysis

## Overview

This repository is part of the [Acral Melanoma PDX Models from Latin America](https://github.com/team113sanger/Acral_Melanoma_PDX_models_from_Latin_America) project.

It contains the RNA-seq processing and downstream analysis pipelines used to investigate gene expression and alternative splicing in acral melanoma cell lines.

The analyses focus on two acral melanoma cell lines, **AM007 and AM016**, comparing control conditions (sgControl) with depletion of the spliceosomal factor **SF3B4 (sgSF3B4)**. 
Each condition was performed in triplicate (total: 12 samples)

The main objectives are to:

- Identify changes in alternative splicing following SF3B4 knockout.
- Compare splicing responses between the two cell lines.
- Investigate the predicted functional consequences of differential splicing events.
- Identify differentially expressed genes associated with SF3B4 depletion.

The repository provides scripts for RNA-seq processing using [VAST-TOOLS](https://github.com/vastgroup/vast-tools) and downstream analyses in R, together with a reproducible R package environment managed through `renv`.

## Repository Structure

```text
Acral_Melanoma_Cell_Lines_Splicing_Analysis/
├── processing/
│   ├── 01_BAM_to_FASTQ.sh
│   ├── 02_FASTQC.sh
│   ├── 03_VASTTOOLS_align.sh
│   └── 04_VASTTOOLS_combine.sh
├── downstream_analysis/
│   ├── 00_utils.R
│   ├── 01_Filtering_AS_events.R
│   ├── 02_Diff_spliced_events.R
│   ├── 03_Functional_impact_diff_spliced_events.R
│   └── 04_Diff_Gene_Expression.R
├── renv/
│   ├── activate.R
│   └── settings.json
├── renv.lock
├── .Rprofile
└── .gitignore
```

The repository is organized into two main analysis stages:

- **`processing/`** contains Bash scripts for RNA-seq preprocessing and VAST-TOOLS analysis, executed on the CRG HPC cluster using SLURM.
- **`downstream_analysis/`** contains R scripts for alternative splicing analysis, predicted functional impact annotation, and differential gene expression analysis.

The `renv` configuration records the R package versions required to reproduce the downstream analyses.

## 1. RNA-seq Processing

RNA-seq preprocessing and splicing quantification were performed using Bash scripts designed for execution on an HPC cluster with the SLURM workload manager.

The processing scripts are organized in the following order:

| Script | Description |
|---|---|
| `01_BAM_to_FASTQ.sh` | Converts BAM files transferred from the Wellcome Sanger Institute into FASTQ files for downstream processing with VAST-TOOLS. |
| `02_FASTQC.sh` | Performs sequencing read quality control using FastQC. |
| `03_VASTTOOLS_align.sh` | Processes RNA-seq reads using the VAST-TOOLS `align` module. |
| `04_VASTTOOLS_combine.sh` | Uses the VAST-TOOLS `merge` and `combine` modules to integrate sample results and generate inclusion and gene expression tables. |

### Required configuration files

**FASTQ file list**

The processing workflow uses a `fastq_list.txt` file containing the absolute paths to the input FASTQ files, with one file path per line.

Example:

```text
/path/to/AM007_sgCtrol_R1.fastq.gz
/path/to/AM007_sgCtrol_R2.fastq.gz
/path/to/AM007_sgCtrol_R3.fastq.gz
```

These paths must be adapted to the local computing environment.

**VAST-TOOLS sample merging**

The VAST-TOOLS `merge` module can combine alignment outputs from multiple samples into a new set of output files, allowing biological replicates to be grouped into their corresponding experimental conditions.

The following sample-to-condition mapping was used:

```text
AM007_sgCtrol_R1      AM007_sgCtrol
AM007_sgCtrol_R2      AM007_sgCtrol
AM007_sgCtrol_R3      AM007_sgCtrol
AM007_sgSF3B4_R1      AM007_sgSF3B4
AM007_sgSF3B4_R2      AM007_sgSF3B4
AM007_sgSF3B4_R3      AM007_sgSF3B4
AM016_sgCtrol_R1      AM016_sgCtrol
AM016_sgCtrol_R2      AM016_sgCtrol
AM016_sgCtrol_R3      AM016_sgCtrol
AM016_sgSF3B4_R1      AM016_sgSF3B4
AM016_sgSF3B4_R2      AM016_sgSF3B4
AM016_sgSF3B4_R3      AM016_sgSF3B4
```

For additional information about the `align`, `merge`, and `combine` modules, consult the [VAST-TOOLS documentation](https://github.com/vastgroup/vast-tools).

### Main VAST-TOOLS outputs

The downstream analyses use two principal output files:

- `INCLUSION_LEVELS_FULL-hg38-16-v251.tab`: alternative splicing inclusion levels (PSI) and associated quality information.
- `COUNTS-hg38-16.tab`: gene expression read counts.

These files provide the starting point for the alternative splicing and differential gene expression analyses, respectively.

## 2. Downstream Analysis

Downstream analyses were performed in R. Shared functions are defined in `00_utils.R` and sourced by the individual analysis scripts.

### 00_utils.R – Shared Analysis Functions

This file contains reusable functions for:

- Quality filtering of VAST-TOOLS inclusion values.
- PCA visualization of PSI measurements.
- Filtering alternative splicing events and calculating differential PSI (dPSI).
- Quantification and visualization of differential splicing events.
- Correlation analysis of splicing changes between cell lines.
- Annotation and export of predicted functional impacts.
- Differential gene expression analysis using edgeR.

**Additional methodological information:**

- **PSI quality filtering:** The `replaceNAI()` function is adapted from the original `replaceNAI.R` function available in [SpliceNet](https://github.com/estepi/SpliceNet), with minor modifications for compatibility with VAST-TOOLS v2.5.1 `combine` output. It identifies low-quality PSI measurements and intron retention events affected by intron imbalance, which are subsequently treated as missing values.

- **PSI PCA:** PCA plots were generated for all alternative splicing events and separately for each event type. These analyses correspond to **Extended Data Figure 20a**.

- **dPSI calculation:** Differential PSI values were calculated by subtracting the control PSI from the SF3B4-depleted condition PSI using merged sample values (dPSI = PSI KO − PSI control).

- **Differential splicing:** Quantification of differential splicing events contributed to **Figure 5**.

- **Splicing correlations:** Correlations between the dPSI values of AM007 and AM016 correspond to **Extended Data Figure 20b**.

- **Functional impact annotation:** Predicted protein impact annotations were obtained from [VastDB](https://vastdb.crg.eu/wiki/Downloads#Protein_impact), version 2.3. These results contributed to **Figure 5h** and **Extended Data Figure 21**.

- **Differential gene expression:** Expression analysis was performed using edgeR. These results contributed to **Extended Data Figure 21**.

### 01_Filtering_AS_events.R – Alternative Splicing Quality Control and Filtering

Processes the VAST-TOOLS inclusion table and prepares the data for differential splicing analysis.

**Main steps:**

1. Identify low-quality PSI measurements using the `replaceNAI()` function.
2. Process missing values and remove events without valid inclusion measurements.
3. Perform PCA of PSI values, both globally and separately for different splicing event types.
4. Filter events according to variability between biological replicates (SD < 5).
5. Calculate dPSI values between sgSF3B4 and sgControl within each cell line.
6. Retain events passing the filtering criteria in both AM007 and AM016.

Splicing event categories include cassette exons, retained introns, alternative 5′ splice sites, and alternative 3′ splice sites.

**Outputs:** Filtered PSI and dPSI tables (TSV) and PCA plots (SVG).

### 02_Diff_spliced_events.R – Differential Splicing Analysis

Identifies and characterizes alternative splicing changes following SF3B4 depletion.

**Main steps:**

1. Classify alternative splicing events by event type.
2. Identify differential splicing events using **|dPSI| > 10**.
3. Quantify differential splicing events for each cell line and event category.
4. Calculate Pearson correlations between dPSI values in AM007 and AM016.
5. Generate correlation plots for all events and separately for each splicing event type.
6. Export filtered differential splicing events for each cell line.

**Outputs:** Differential splicing barplots and correlation plots (SVG), and filtered event tables (Excel and RDS).

### 03_Functional_impact_diff_spliced_events.R – Predicted Functional Impact of Alternative Splicing

Integrates differential splicing results with protein impact annotations obtained from **VastDB/VAST-TOOLS (hg38, version 2.3)**.

**Main steps:**

1. Process and correct formatting issues in the protein impact annotation file.
2. Select differentially spliced events using |dPSI| > 10.
3. Classify splicing changes according to increased inclusion or exclusion.
4. Associate differential splicing events with their predicted functional consequences.
5. Simplify functional impact annotations into biologically interpretable categories.

Functional impact annotations are grouped into the following categories:

- Alternative protein isoforms.
- Open reading frame (ORF) disruption.
- Changes affecting untranslated regions (UTRs).
- Noncoding transcript regions.
- Uncertain or unavailable functional impact annotations.

Results are exported both with separate UTR categories and with combined UTR/noncoding annotations.

**Outputs:** Differential splicing events with original and simplified functional impact annotations (Excel and RData).

### 04_Diff_Gene_Expression.R – Differential Gene Expression Analysis

Performs differential gene expression analysis using **edgeR**.

**Main steps:**

1. Extract raw gene expression counts for individual biological replicates.
2. Remove genes with missing count values.
3. Filter lowly expressed genes using `edgeR::filterByExpr()`.
4. Normalize gene expression counts using the Trimmed Mean of M-values (TMM) method.
5. Generate PCA and MDS plots to assess sample relationships.
6. Perform differential expression analysis independently for AM007 and AM016.
7. Fit quasi-likelihood negative binomial generalized linear models and perform statistical testing.
8. Identify differentially expressed genes using **|log2FC| ≥ 1 and FDR ≤ 0.05**.

**Outputs:** Gene expression PCA and MDS plots, library-size barplots (SVG), differential expression tables (Excel and RData), and TMM-normalized log2-CPM values (RData).

## 3. Software and Dependencies

The pipeline uses the following main tools:

**RNA-seq processing**

| Software | Version | Repository |
|---|---|---|
| SLURM | HPC workload manager | [SLURM](https://slurm.schedmd.com/) |
| VAST-TOOLS | 2.5.1 | [VAST-TOOLS GitHub](https://github.com/vastgroup/vast-tools) |
| FastQC | 0.12.1 | [FastQC GitHub](https://github.com/s-andrews/fastqc) |
| SAMtools | 1.20 | [SAMtools GitHub](https://github.com/samtools/samtools) |

**Downstream analysis**

- R
- edgeR
- limma
- dplyr
- tidyr
- stringr
- ggplot2
- Other supporting R packages

R package dependencies and versions are managed using `renv`.

### Restoring the R environment

Clone the repository:

```bash
git clone https://github.com/msalpe/Acral_Melanoma_Cell_Lines_Splicing_Analysis.git
cd Acral_Melanoma_Cell_Lines_Splicing_Analysis
```

Start R from the repository root and restore the package environment:

```r
renv::restore()
```

Check that the environment is consistent:

```r
renv::status()
```

On an HPC cluster, an appropriate R module should be loaded before restoring the environment.

## 4. Input Data and Output Files

The analyses require RNA-seq input files, VAST-TOOLS outputs, and external annotation files used for functional impact assessment.

Large sequencing files, intermediate processing files, and generated results are not included in this repository.

Protein impact annotations must be downloaded separately from the [VastDB Downloads page](https://vastdb.crg.eu/wiki/Downloads#Protein_impact). **Version 2.3 (hg38)** was used for this analysis.

The downstream scripts generate tables and figures for two principal analyses:

- **Alternative splicing (AS):** PSI quality control, differential splicing, and predicted functional impact.
- **Gene expression (GE):** normalized gene expression, sample clustering, and differential expression.

Results can be organized locally into the following directories:

```text
results/
├── AS/
│   ├── Files/
│   └── Graphs/
└── GE/
    ├── Files/
    └── Graphs/
```

Input and output directories must be configured according to the local computing environment before executing the scripts.

## 5. Reproducibility

The repository is intended to facilitate reproducibility of RNA-seq processing and downstream analyses.

To reproduce the workflow:

1. Obtain the appropriate RNA-seq input data and external annotation files.
2. Configure the input, output, and software paths for the computing environment.
3. Execute the scripts in `processing/` in numerical order.
4. Restore the R environment using `renv::restore()`.
5. Execute the scripts in `downstream_analysis/` in numerical order, sourcing `00_utils.R` as required.

The analysis scripts document the principal filtering thresholds, normalization methods, statistical procedures, and output formats.

The `renv.lock` file records the R package versions used for the downstream analyses.

## 6. Contact

For questions, suggestions, or feedback regarding this repository, please contact:

**Mónica Salinas**

- [monicasalinaspena@gmail.com](mailto:monicasalinaspena@gmail.com)
- [monica.salinas@crg.eu](mailto:monica.salinas@crg.eu)
