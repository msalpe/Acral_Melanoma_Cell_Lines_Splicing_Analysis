library(dplyr)

# ======================================================================================
# 04 - DIFFERENTIAL GENE EXPRESSION ANALYSIS USING edgeR
# ======================================================================================
# Input:
#   - VAST-TOOLS combine gene expression count table (COUNTS-hg38-16.tab)
#   - Sample metadata defining cell lines, conditions and biological replicates
#   - Comparison table specifying treatment versus control conditions
#
# Analysis:
#   - Extract raw gene counts from individual biological replicates
#   - Remove genes with missing count values
#   - Filter lowly expressed genes using edgeR::filterByExpr()
#   - Perform TMM normalization of gene expression counts
#   - Generate PCA and MDS plots to assess sample clustering
#   - Perform differential expression analysis independently for AM007 and AM016
#   - Fit edgeR quasi-likelihood negative binomial GLMs and perform statistical tests
#   - Identify differentially expressed genes (abs(log2FC) >= 1, FDR <= 0.05)
#
# Output:
#   - Raw library-size barplot, PCA and MDS plots (SVG)
#   - Differential gene expression results for each comparison (Excel and RData)
#   - Filtered differentially expressed genes (Excel and RData)
#   - TMM-normalized log2-CPM values and analysis metadata (RData)
#
# Dependencies:
#   - VAST-TOOLS combine output (gene count table)
#   - 00_utils.R
#   - dplyr, edgeR, limma, openxlsx and additional packages managed through renv
#
# ======================================================================================
# paths & directories
root_dir <- ""
files_dir <- file.path(root_dir,"GE", "Files")
graphs_dir <- file.path(root_dir,"GE", "Graphs")

dir.create(files_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs_dir, recursive = TRUE, showWarnings = FALSE)

# functions/utils
source(file.path(root_dir,"code_paper/downstream_analyses/00_utils.R"))

# Load read counts table
counts <- read.table(file.path(root_dir,"vast_out/COUNTS-hg38-16.tab"), header = TRUE, sep = "\t") # 41614

# params
log2FC_threshold <- 1
fdr_threshold <- 0.05


# ============================================================================
###### BUILDING NEEDED DATAFRAMES ###### 
# ============================================================================
# Sample_groups df: config to associate replicates to conditions

sample_groups <- data.frame(
  sample    = c("AM007_sgCtrol_R1", "AM007_sgCtrol_R2", "AM007_sgCtrol_R3", 
                "AM007_sgSF3B4_R1", "AM007_sgSF3B4_R2", "AM007_sgSF3B4_R3",
                "AM016_sgCtrol_R1", "AM016_sgCtrol_R2", "AM016_sgCtrol_R3",
                "AM016_sgSF3B4_R1", "AM016_sgSF3B4_R2", "AM016_sgSF3B4_R3"),
  
  condition = c("AM007_sgCtrol","AM007_sgCtrol" ,"AM007_sgCtrol", 
                "AM007_sgSF3B4","AM007_sgSF3B4" ,"AM007_sgSF3B4",
                "AM016_sgCtrol","AM016_sgCtrol" ,"AM016_sgCtrol",
                "AM016_sgSF3B4","AM016_sgSF3B4" ,"AM016_sgSF3B4"),
  stringsAsFactors = FALSE
)

# Comparison_table df with all the comparisons: Cond vs Ctrol
comparison_table <- data.frame(
  condition = c(
    "AM007_sgSF3B4", "AM016_sgSF3B4"),
  control = c(
    "AM007_sgCtrol", "AM016_sgCtrol"),
  stringsAsFactors = FALSE
)


# ============================================================================
# Preparing data to run edgeR 
# ============================================================================

# colnames
colnames(counts)

# remove merge samples. keep only replicates
counts <- counts %>%
  select(1:2, matches("_R[0-9]+$"))

# Remove rows with NA (to avoid issues during normalization or modeling)
counts <- na.omit(counts) # 40188

# Split into gene annotation and count matrix
gene_info <- counts[, 1:2]                  # e.g., ENSEMBL ID and gene name
count_data <- counts[, 3:ncol(counts)]      # raw counts

# Ensure counts are integers (required by edgeR)
count_data <- round(as.matrix(count_data))

# ENSEMBL ID as rownames
rownames(gene_info) <- gene_info[,1]
rownames(count_data) <- gene_info[,1]


# ============================================================================
### Diff gene expression
# ============================================================================
# run edgeR and save outputs
edger_analysis <- run_edger_comparisons(
                                  count_data = count_data,
                                   gene_info = gene_info,
                               sample_groups = sample_groups,
                            comparison_table = comparison_table,
                                   files_dir = files_dir,
                                  graphs_dir = graphs_dir,
                                 output_name = "AM202607_sgSF3B4_vs_sgCtrol_edgeR",
                                       stats = TRUE,
                             logFC_threshold = log2FC_threshold,
                               FDR_threshold = fdr_threshold)


# ===================================================================
# ===================================================================
