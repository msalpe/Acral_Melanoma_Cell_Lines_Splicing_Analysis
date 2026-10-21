
library(dplyr)


# ======================================================================================
# 01 - FILTERING AND QUALITY CONTROL OF ALTERNATIVE SPLICING EVENTS
# ======================================================================================
# Input:
#   - VAST-TOOLS inclusion table (output of vast-tools combine)
#   - Sample metadata (cell line, condition and biological replicates)
#
# Analysis:
#   - Identify low-quality PSI values and IR values (replaceNAI function)
#   - Remove events with missing PSI 
#   - Perform PCA of PSI  for all events and per splicing event type
#   - Filter events based on variability between biological replicates 
#   - Calculate dPSI between sgSF3B4 and sgControl for AM007 and AM016 cell lines
#   - Retain events passing filtering criteria in both cell lines
#
# Output:
#   - Corrected VAST-TOOLS inclusion table
#   - PSI PCA plots (SVG)
#   - Filtered PSI and dPSI tables for each cell line (TSV)
#   - Combined dPSI tables containing events shared between cell lines (TSV)
#
# Dependencies:
#   - 00_utils.R
#   - dplyr and additional R packages managed through renv
#
# ======================================================================================

# paths & directories
root_dir <- ""
files_dir <- file.path(root_dir,"AS", "Files")
graphs_dir <- file.path(root_dir,"AS", "Graphs")

dir.create(files_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(graphs_dir, recursive = TRUE, showWarnings = FALSE)


# functions/utils
source(file.path(root_dir,"00_utils.R"))


# files, parameters and custom colors

# files for replaceNAI.R
INCFile <-file.path(root_dir,"/vast_out/INCLUSION_LEVELS_FULL-hg38-16-v251.tab")
OUTFile <- file.path(root_dir,"/vast_out/INCLUSION_LEVELS_FULL-hg38-16-v251_newNA2_NAI.tab")

# Define annotation colors based on event type
ann_colors <- list(
  EVENT_TYPE = c(
    "Alt3" = "#ef476f",  
    "Alt5" = "#ffd166",  
    "Exon" = "#118ab2", 
    "Intron" = "#06d6a0"
  )
)

# custom color palette 
colors <- c(
  "sgCtrol"    = "gray30",
  "sgSF3B4"    = "#e74c3c")

# params
dataset <- "AM_sgSF3B4"   # just for file names

# metadata
metadata <- data.frame(
  sample    = c("AM007_sgCtrol_R1", "AM007_sgCtrol_R2", "AM007_sgCtrol_R3", 
                "AM007_sgSF3B4_R1", "AM007_sgSF3B4_R2", "AM007_sgSF3B4_R3",
                "AM016_sgCtrol_R1", "AM016_sgCtrol_R2", "AM016_sgCtrol_R3",
                "AM016_sgSF3B4_R1", "AM016_sgSF3B4_R2", "AM016_sgSF3B4_R3"),
  
  condition = c("sgCtrol","sgCtrol" ,"sgCtrol", 
                "sgSF3B4","sgSF3B4" ,"sgSF3B4",
                "sgCtrol","sgCtrol" ,"sgCtrol",
                "sgSF3B4","sgSF3B4" ,"sgSF3B4"),
  
  cell_line = c("AM007", "AM007", "AM007", 
                "AM007", "AM007", "AM007",
                "AM016", "AM016", "AM016",
                "AM016", "AM016", "AM016"),
  stringsAsFactors = FALSE
)

# ======================================================================================
# QC: replaceNAI
# ======================================================================================

replaceNAI(INCFile, OUTFile)


# Inclusion Table (vast-tools output) after applying replaceNAI.R function (NAold & NAnew2 values)
incData <- read.table(OUTFile, stringsAsFactors = FALSE, sep="\t", header=TRUE )  # 721551
colnames(incData)[1] <- "EVENT"   #just change column name, error in replaceNAI.R


# NAold, NAnew2, NAI as "NA" values
incData <- incData %>%
  mutate(across(everything(), ~ ifelse(. %in% c("NAold", "NAnew2", "NAI"), NA, .)))

incData <- incData %>%
  mutate(across(7:ncol(.), as.numeric))

# colnames
colnames(incData)


# number of NAs per sample
colSums(is.na(incData[, 7:ncol(incData)]))

# remove events for which ALL samples are NA

incData <- incData[!apply(incData[, 7:ncol(incData)], 1, function(x) all(is.na(x))), ] #506893
colSums(is.na(incData[, 7:ncol(incData)]))
# AM007_sgCtrol AM007_sgCtrol_R1 AM007_sgCtrol_R2 AM007_sgCtrol_R3    AM007_sgSF3B4 AM007_sgSF3B4_R1 AM007_sgSF3B4_R2 AM007_sgSF3B4_R3    AM016_sgCtrol 
# 49868            76611            95107            73991            50710            85777            82726            80834            37853 
# AM016_sgCtrol_R1 AM016_sgCtrol_R2 AM016_sgCtrol_R3    AM016_sgSF3B4 AM016_sgSF3B4_R1 AM016_sgSF3B4_R2 AM016_sgSF3B4_R3 
# 72194            68649            79299            46775            77136            88815            91257 

# for PCAs, remove merged samples
PSI <- incData %>%
  select(
    1:6, 
    matches("_R[0-9]+$") 
  )
# remove data with NA values for PCA
PSI <- na.omit(PSI) # 374299 



# ======================================================================================
# PCA plots
# ======================================================================================

# list per event type
PSI_per_type <- list(
  all      = PSI,
  exons    = PSI[grep("^HsaEX\\d+", PSI$EVENT), ], 
  introns  = PSI[grep("^HsaINT\\d+", PSI$EVENT), ],               # all introns  
  alt5     = PSI[grep("^HsaALTD\\d+", PSI$EVENT), ],              # all alt5     
  alt3     = PSI[grep("^HsaALTA\\d+", PSI$EVENT), ]               # all alt3     
)


plot_titles <- c("All_events", "Exons", "Introns", "Alt5", "Alt3")


# PCA plot per event type, custom colors by condition and shape by Cell line
save_pca_plots_2variables(
  data_list     = PSI_per_type,
  plot_titles   = plot_titles,
  metadata      = metadata,
  from_col      = 7,  
  output_dir    = graphs_dir, 
  cond1_label   = "Condition",
  cond2_label   = "Cell Line",
  output_svg    = "PCA_Plots_PSI_AM_sgSF3B4_shape_by_cellline.svg",
  svg_width     = 13,
  svg_height    = 6,
  colors_cond1  = colors   
)


##=============================================================
#  FILTERING EVENTS AND CALCULATING dPSI values #
# =============================================================

# using helper function to filter PSI values
# SD between replicates < 5 with no NAs
incData <- na.omit(incData) #370465


results <- lapply(
  c("AM007", "AM016"),
  \(x) filter_PSI_events_calculate_dPSI(incData, cell_line = x)
)

names(results) <- c("AM007", "AM016")

data_AM007 <- results[["AM007"]]
data_AM016 <- results[["AM016"]]


# Just to better compare between cell lines lets select events in both
common_events <- intersect(data_AM007$EVENT, data_AM016$EVENT)

data_AM007 <- data_AM007 %>%
  filter(EVENT %in% common_events) #356613

data_AM016 <- data_AM016 %>%
  filter(EVENT %in% common_events) #356613


# combine df with dPSI values
dpsi_annot <- data_AM007 %>%
  select(1:6, dPSI_AM007) %>%
  inner_join(
    data_AM016 %>% select(EVENT, dPSI_AM016),
    by = "EVENT"
  )

dpsi <- dpsi_annot %>% select(EVENT, dPSI_AM007, dPSI_AM016)


# Save data PSI values
write.table(
  data_AM007, 
  file = file.path(files_dir, paste0("IncTable_PSI_dPSI_values_AM007_n", nrow(data_AM007), ".txt")),
  sep = "\t", 
  quote = FALSE, 
  col.names = TRUE, 
  row.names = FALSE
)

write.table(
  data_AM016, 
  file = file.path(files_dir, paste0("IncTable_PSI_dPSI_values_AM016_n", nrow(data_AM016), ".txt")),
  sep = "\t", 
  quote = FALSE, 
  col.names = TRUE, 
  row.names = FALSE
)

# Save data dPSI values
write.table(
  dpsi_annot, 
  file = file.path(files_dir, paste0("IncTable_dPSI_values_AMcelllines_sgSF3B4_n", nrow(dpsi_annot), ".txt")),
  sep = "\t", 
  quote = FALSE, 
  col.names = TRUE, 
  row.names = FALSE
)

# Save data dPSI values
write.table(
  dpsi, 
  file = file.path(files_dir, paste0("dPSI_values_AMcelllines_sgSF3B4_n", nrow(dpsi), ".txt")),
  sep = "\t", 
  quote = FALSE, 
  col.names = TRUE, 
  row.names = FALSE
)


# ===================================================================
# ===================================================================






