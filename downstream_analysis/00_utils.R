# helper / utils functions to perform downstream analysis of AS and DEG 


## replaceNAI.R function with minor modifications
# NO OUTFile1 (only raplaceNA), now only output is replaceNAI
# Counts N are two (no NAnew2 in new vast-tools version)

# Use dplyr pipe without attaching packages
`%>%` <- dplyr::`%>%`


replaceNAI <- function(INCFile,
                       #OUTFile1,
                       OUTFile) {
  #######################################################################
  #READ INPUT FILE
  INCData <- read.csv(
    INCFile,
    header = T,
    sep = "\t",
    stringsAsFactors = F,
    row.names = 2
  )
  #######################################################################
  #subset PSI values
  onlyValues <- INCData[, 6:ncol(INCData)]
  #subset QUAL values
  Qvals <- grep("\\.Q$", colnames(onlyValues))
  qualValues <- onlyValues[, Qvals]
  psiValues  <- onlyValues[, -Qvals]
  finalValues <- matrix(0, ncol = ncol(psiValues), nrow = nrow(psiValues))
  colnames(finalValues) <- colnames(onlyValues)[-Qvals]
  #######################################################################
  # count "N," for ALL events
  Ns_comma <- t(apply(qualValues, 1, function(x) {
    stringr::str_count(x, "N,")
  }))
  #######################################################################
  #replace
  indexes_N <- which(Ns_comma != 0, arr.in = TRUE)
  indexes_NewN2 <- which(Ns_comma == 2, arr.in = TRUE)
  
  finalValues[indexes_N] <- Ns_comma[indexes_N]
  head(finalValues)
  correctedPsiVal <- psiValues
  
  head(correctedPsiVal)
  correctedPsiVal[indexes_NewN2] <- "NAnew2"
  originalNA <- is.na(psiValues)
  correctedPsiVal[originalNA] <- "NAold"
  fullNAS <- cbind(INCData[, 1:5], correctedPsiVal)
  #####################################################
  #compute and print STATS:
  tVal <- length(originalNA)# 235328422
  oNA <- length(which(originalNA))#95995460 40.8
  new3 <- length(indexes_NewN2) / 2# 827226 28.39
  oNA2 <- length(which(correctedPsiVal == "NAold"))#
  new3b <- length(which(correctedPsiVal == "NAnew2"))#
  m1 <- paste(round(oNA / tVal * 100), "%")
  m3 <- paste(round(new3b / tVal * 100), "%")
  
  message(paste("Num total PSI values:", tVal))
  message(paste("NA original values:",  m1))
  message(paste("NewN2 values:", m3))
  #####################################################
  #PRINT OUTPUT FILE. Monica: No need for OUTFile1 
  #write.table(fullNAS, OUTFile1, sep = "\t", col.names = NA)
  #message(paste("Final corrected table only NewNAs is written in:", OUTFile1))
  #####################################################
  #correct IR:
  message("Number of events by type: ")
  print(table(fullNAS$COMPLEX))
  #IRCqual<-qualValues[fullNAS$COMPLEX=="IR-C" | fullNAS$COMPLEX=="IR-S" ,]
  ii <- which(fullNAS$COMPLEX == "IR")
  IRCqual <- qualValues[ii, ]
  # IRCpsi<-fullNAS[fullNAS$COMPLEX=="IR-C" | fullNAS$COMPLEX=="IR-S" ,6:ncol(fullNAS)]
  IRCpsi <- fullNAS[ii , 6:ncol(fullNAS)]
  #####################################################################
  
  
  pval <- apply(IRCqual, 2, function(x)
  {
    as.numeric(stringr::str_match(matrix(
      unlist(strsplit(x, ",")), ncol = 6, byrow = T
    )[, 5], "(.*)@")[, 2])
  })
  pvalCorrected <- apply(pval, 2, p.adjust)
  pvalCorrected[pvalCorrected < 0.05] <- "NAI"
  ######################################################################
  Nis <- t(apply(pvalCorrected, 1, function(x)
  {
    stringr::str_count(x, "NAI")
  }))
  indexes_NAI <- which(Nis != 0, arr.in = TRUE)
  #####################################################################
  dim(IRCpsi)
  dim(IRCqual)
  rownames(IRCpsi)[1:10]
  rownames(IRCqual)[1:10]
  correctedIRpsi <- IRCpsi
  
  #compute and print STATS IR
  tIR <- dim(IRCpsi)[1] * dim(IRCpsi)[2]
  tNewNAI <- dim(indexes_NAI)[1]
  
  m4 <- paste(round(tNewNAI / tIR * 100), "%")
  message(paste("Num total PSI values:", tIR))
  
  correctedIRpsi[indexes_NAI] <- "NAI"
  df1 <- correctedPsiVal
  
  df2 <- correctedIRpsi
  
  df1[match(rownames(df2) , rownames(df1)),] <- df2
  fullNASIR <- cbind(INCData[, 1:5], df1)
  
  message("pval analyzed (for correction): ", length(pvalCorrected))
  message("pval corrected <0.05: ", m4)
  
  write.table(fullNASIR, OUTFile, sep = "\t", col.names = NA)
  message(paste("Final corrected table NewNAs+NAI is written in:", OUTFile))
  
  
}



plot_PCA_2variables <- function(data,
                                metadata,
                                plot_title = "PCA Plot",
                                columns_to_use = NULL,
                                colors_cond1 = NULL,
                                shapes_cond2 = NULL,
                                cond1_label = "Condition 1",
                                cond2_label = "Condition 2",
                                na_threshold = 0.2,
                                top_variable_fraction = NULL) {
  if (ncol(metadata) < 2) {
    stop("metadata must contain at least two columns: sample and condition1.")
  }
  
  has_condition2 <- ncol(metadata) >= 3
  
  metadata <- metadata %>%
    as.data.frame()
  
  colnames(metadata)[1:2] <- c("Sample", "Condition1")
  
  if (has_condition2) {
    colnames(metadata)[3] <- "Condition2"
  }
  
  metadata <- metadata %>%
    dplyr::mutate(
      Sample = as.character(Sample),
      Condition1 = as.factor(Condition1)
    )
  
  if (has_condition2) {
    metadata <- metadata %>%
      dplyr::mutate(Condition2 = as.factor(Condition2))
  }
  
  if (is.null(columns_to_use)) {
    columns_to_use <- metadata$Sample
  }
  
  if (is.list(columns_to_use)) {
    columns_to_use <- unlist(columns_to_use)
  }
  
  if (is.numeric(columns_to_use)) {
    sample_names <- colnames(data)[columns_to_use]
  } else {
    sample_names <- as.character(columns_to_use)
  }
  
  missing_in_data <- setdiff(sample_names, colnames(data))
  if (length(missing_in_data) > 0) {
    stop("These samples are missing from data: ",
         paste(missing_in_data, collapse = ", "))
  }
  
  missing_in_metadata <- setdiff(sample_names, metadata$Sample)
  if (length(missing_in_metadata) > 0) {
    stop("These samples are missing from metadata: ",
         paste(missing_in_metadata, collapse = ", "))
  }
  
  df_pca <- data %>%
    dplyr::select(dplyr::all_of(sample_names)) %>%
    as.matrix()
  
  rownames(df_pca) <- data[[1]]
  
  x <- t(df_pca)
  
  keep <- colMeans(is.na(x)) <= na_threshold
  x <- x[, keep, drop = FALSE]
  
  for (j in seq_len(ncol(x))) {
    x[is.na(x[, j]), j] <- mean(x[, j], na.rm = TRUE)
  }
  
  vars <- apply(x, 2, var)
  x <- x[, vars > 0, drop = FALSE]
  vars <- vars[vars > 0]
  
  if (!is.null(top_variable_fraction)) {
    if (top_variable_fraction <= 0 || top_variable_fraction > 1) {
      stop("top_variable_fraction must be > 0 and <= 1.")
    }
    
    cutoff <- quantile(vars, probs = 1 - top_variable_fraction, na.rm = TRUE)
    x <- x[, vars >= cutoff, drop = FALSE]
  }
  
  if (ncol(x) < 2) {
    stop("Not enough valid variable events for PCA after filtering.")
  }
  
  pca_result <- prcomp(x, center = TRUE, scale. = FALSE)
  
  variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
  pc1_var <- round(variance_explained[1], 2)
  pc2_var <- round(variance_explained[2], 2)
  
  pca_df <- as.data.frame(pca_result$x) %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::left_join(metadata, by = "Sample")
  
  if (!is.null(colors_cond1) && !is.null(names(colors_cond1))) {
    present_levels <- unique(as.character(pca_df$Condition1))
    ordered_levels <- intersect(names(colors_cond1), present_levels)
    remaining_levels <- setdiff(present_levels, ordered_levels)
    
    pca_df$Condition1 <- factor(
      pca_df$Condition1,
      levels = c(ordered_levels, remaining_levels)
    )
  }
  
  if (has_condition2 && !is.null(shapes_cond2) && !is.null(names(shapes_cond2))) {
    present_levels <- unique(as.character(pca_df$Condition2))
    ordered_levels <- intersect(names(shapes_cond2), present_levels)
    remaining_levels <- setdiff(present_levels, ordered_levels)
    
    pca_df$Condition2 <- factor(
      pca_df$Condition2,
      levels = c(ordered_levels, remaining_levels)
    )
  }
  
  if (has_condition2) {
    pca_plot <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, color = Condition1, shape = Condition2)
    )
  } else {
    pca_plot <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, color = Condition1)
    )
  }
  
  pca_plot <- pca_plot +
    ggplot2::geom_point(size = 4, alpha = 0.85) +
    ggplot2::labs(
      title = plot_title,
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)"),
      color = cond1_label,
      shape = if (has_condition2) cond2_label else NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      legend.key.size = unit(0.4, "cm")
    )
  
  if (!is.null(colors_cond1)) {
    pca_plot <- pca_plot + ggplot2::scale_color_manual(values = colors_cond1)
  }
  
  if (has_condition2 && !is.null(shapes_cond2)) {
    pca_plot <- pca_plot + ggplot2::scale_shape_manual(values = shapes_cond2)
  }
  
  return(pca_plot)
}



save_pca_plots_2variables <- function(data_list,
                                      plot_titles,
                                      metadata,
                                      from_col = 7,
                                      to_col = NULL,
                                      output_svg = "PCA_Plots.svg",
                                      output_dir = getwd(),
                                      colors_cond1 = NULL,
                                      shapes_cond2 = NULL,
                                      cond1_label = "Condition 1",
                                      cond2_label = "Condition 2",
                                      na_threshold = 0.2,
                                      top_variable_fraction = NULL,
                                      svg_width = 13,
                                      svg_height = 6,
                                      ncol_grid = 3,
                                      nrow_grid = 2) {
  if (length(plot_titles) != length(data_list)) {
    stop("Error: The number of plot titles must match the number of dataframes in data_list.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  svg_path <- file.path(output_dir, output_svg)
  pca_plots <- list()
  
  for (i in seq_along(data_list)) {
    
    data <- data_list[[i]]
    plot_title <- plot_titles[i]
    
    current_to_col <- if (is.null(to_col)) ncol(data) else to_col
    columns_to_use <- from_col:current_to_col
    
    pca_plots[[i]] <- plot_PCA_2variables(
      data = data,
      metadata = metadata,
      plot_title = plot_title,
      columns_to_use = columns_to_use,
      colors_cond1 = colors_cond1,
      shapes_cond2 = shapes_cond2,
      cond1_label = cond1_label,
      cond2_label = cond2_label,
      na_threshold = na_threshold,
      top_variable_fraction = top_variable_fraction
    )
  }
  
  svg(svg_path, width = svg_width, height = svg_height)
  
  plot_list <- unname(pca_plots)
  num_plots <- length(plot_list)
  plots_per_page <- ncol_grid * nrow_grid
  
  for (i in seq(1, num_plots, by = plots_per_page)) {
    end_index <- min(i + plots_per_page - 1, num_plots)
    
    gridExtra::grid.arrange(
      grobs = plot_list[i:end_index],
      ncol = ncol_grid,
      nrow = nrow_grid
    )
  }
  
  dev.off()
  
  message("SVG successfully generated at: ", svg_path)
}




# filter splicing events considering SD in replicates per cell line and calculating dPSI values
filter_PSI_events_calculate_dPSI <- function(incData,
                                             cell_line,
                                             target = "SF3B4",
                                             sd_cutoff = 5,
                                             min_non_na = 3) {
  
  merged_control <- paste0(cell_line, "_sgCtrol")
  merged_cond    <- paste0(cell_line, "_sg", target)
  
  control_reps <- grep(
    paste0("^", cell_line, "_sgCtrol_R[0-9]+$"),
    colnames(incData),
    value = TRUE
  )
  
  sample_reps <- grep(
    paste0("^", cell_line, "_sg", target, "_R[0-9]+$"),
    colnames(incData),
    value = TRUE
  )
  
  required_cols <- c(
    merged_control,
    merged_cond,
    control_reps,
    sample_reps
  )
  
  missing_cols <- setdiff(required_cols, colnames(incData))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for ", cell_line, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  data <- incData %>%
    dplyr::select(1:6, dplyr::all_of(required_cols)) %>%
    dplyr::filter(
      rowSums(!is.na(dplyr::across(dplyr::all_of(control_reps)))) >= min_non_na,
      rowSums(!is.na(dplyr::across(dplyr::all_of(sample_reps)))) >= min_non_na
    )
  
  data$SD_Ctrol <- apply(
    data[, control_reps, drop = FALSE],
    1,
    sd,
    na.rm = TRUE
  )
  
  data$SD_cond <- apply(
    data[, sample_reps, drop = FALSE],
    1,
    sd,
    na.rm = TRUE
  )
  
  dpsi_name <- paste0("dPSI_", cell_line)
  
  data %>%
    dplyr::filter(
      SD_Ctrol < sd_cutoff,
      SD_cond < sd_cutoff
    ) %>%
    dplyr::mutate(
      !!dpsi_name := .data[[merged_cond]] - .data[[merged_control]]
    ) %>%
    dplyr::filter(!is.na(.data[[dpsi_name]])) %>%
    dplyr::select(-SD_Ctrol, -SD_cond)
}




plot_diff_spliced_events_barplot <- function(dPSI_data, 
                                             num_cols = 3:ncol(dPSI_data), 
                                             colors_per_event_type = c("Alt3" = "#ef476f",  
                                                                       "Alt5" = "#ffd166",  
                                                                       "Exon" = "#118ab2", 
                                                                       "Intron" = "#06d6a0"),
                                             dPSI_threshold = 10,
                                             output_dir = getwd(),
                                             output_svg = "Barplot_DiffSplicedEvents_per_event_type.svg",
                                             svg_width = 10,
                                             svg_height = 6,
                                             label_size = 5) { # Added customizable text label size parameter
  
  # 1. Resolve column names if passed as numeric indices
  if (is.numeric(num_cols)) {
    num_cols <- names(dPSI_data)[num_cols]
  }
  
  # 2. Clean up event names and reshape data to long format
  dPSI_long <- dPSI_data %>%
    dplyr::mutate(
      EVENT_PREFIX = dplyr::if_else(
        EVENT_TYPE %in% c("Alt3", "Alt5"),
        sub("-[0-9]+/[0-9]+$", "", EVENT),
        EVENT
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(num_cols),
      names_to = "Condition",
      values_to = "Value"
    ) %>%
    dplyr::filter(abs(Value) > dPSI_threshold & !is.na(Value))
  
  # 3. Compute unique event category totals per condition
  summary_dPSI <- dPSI_long %>%
    dplyr::group_by(Condition, EVENT_TYPE) %>% 
    dplyr::summarise(
      TOTAL = dplyr::if_else(         
        dplyr::first(EVENT_TYPE) %in% c("Alt3", "Alt5"),
        dplyr::n_distinct(EVENT_PREFIX),
        dplyr::n()
      ),
      .groups = "drop"
    )
  
  # 4. Enforce structural factor constraints for strict plot alignment
  summary_dPSI$Condition <- factor(summary_dPSI$Condition, levels = num_cols)
  summary_dPSI$EVENT_TYPE <- factor(summary_dPSI$EVENT_TYPE, levels = names(colors_per_event_type))
  
  # 5. Build full file path and check destination folder integrity
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  svg_path <- file.path(output_dir, output_svg)
  
  # 6. Generate and save the multi-layered stacked bar chart
  svg(svg_path, width = svg_width, height = svg_height)
  
  p <- ggplot2::ggplot(summary_dPSI, ggplot2::aes(x = Condition, y = TOTAL, fill = EVENT_TYPE)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    
    # Linked to label_size parameter dynamically
    ggplot2::geom_text(ggplot2::aes(label = TOTAL), position = ggplot2::position_stack(vjust = 0.5), size = label_size, color = "black") + 
    
    ggplot2::scale_fill_manual(values = colors_per_event_type) +
    ggplot2::labs(
      title = paste0("Events with abs(dPSI) > ", dPSI_threshold),
      x = "",
      y = "Total Events",
      fill = "Event Type"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), 
      axis.text.y = ggplot2::element_text(size = 14),  
      axis.title.x = ggplot2::element_text(size = 16, face = "bold"),  
      axis.title.y = ggplot2::element_text(size = 16, face = "bold"),  
      legend.title = ggplot2::element_text(size = 14, face = "bold"),  
      legend.text = ggplot2::element_text(size = 12), 
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)  
    )
  
  print(p)
  dev.off()
  
  message(paste("Barplot showing diff spliced events per event type successfully saved at:", svg_path))
}



# plotting corr plot between dPSI values

make_cor_plot <- function(df, plot_title, x_col, y_col, axis_limits = c(-100, 100)) {
  
  cor_result <- cor.test(
    df[[x_col]],
    df[[y_col]],
    method = "pearson"
  )
  
  annotation <- paste0(
    "r = ", round(unname(cor_result$estimate), 2),
    "\np = ", format.pval(cor_result$p.value, digits = 2),
    "\nn = ", nrow(df)
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_point(alpha = 0.45, size = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = axis_limits[1],
      y = axis_limits[2],
      label = annotation,
      hjust = 0,
      vjust = 1,
      size = 3.5
    ) +
    ggplot2::coord_equal(
      xlim = axis_limits,
      ylim = axis_limits
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "dPSI AM007",
      y = "dPSI AM016"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}



# Filter differential splicing events by dPSI threshold and export results to Excel and RDS.
save_data_diff_events <- function(dPSI_data, 
                                  start_col = 7, 
                                  dPSI_threshold = 10,
                                  output_dir = getwd(),
                                  excel_filename = "dPSI_values_diff_events_per_condition.xlsx",
                                  r_object_filename = "dPSI_values_diff_events_per_condition.rds") {
  # 1. Extract the target column names to iterate over
  target_cols <- names(dPSI_data)[start_col:ncol(dPSI_data)]
  # 2. Build the list of filtered data frames using lapply
  filtered_list <- lapply(target_cols, function(col_name) {
    # Filter rows matching the strict absolute threshold criteria
    # and remove potential NA records in that specific column
    df_filtered <- dplyr::filter(dPSI_data, abs(.data[[col_name]]) > dPSI_threshold & !is.na(.data[[col_name]]))
    # Sort rows in ascending order based on the current column values
    df_filtered <- dplyr::arrange(df_filtered, .data[[col_name]])
    return(df_filtered)
  })
  # 3. Assign names to the list elements matching their corresponding column names
  names(filtered_list) <- target_cols
  # 4. Verify and prepare the destination directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  excel_path <- file.path(output_dir, excel_filename)
  r_obj_path <- file.path(output_dir, r_object_filename)
  # 5. Export as a multi-sheet Excel Workbook using openxlsx
  # openxlsx natively converts a named list into sheets matching the list names
  openxlsx::write.xlsx(filtered_list, file = excel_path, overwrite = TRUE)
  # 6. Save the list as a compressed native R object (.rds)
  saveRDS(filtered_list, file = r_obj_path)
  # 7. Print success messages and return the object to the environment
  message(paste("Excel spreadsheet successfully saved at:", excel_path))
  message(paste("R object (.rds) successfully saved at:", r_obj_path))
  return(filtered_list)
}


# Classify differentially spliced events, annotate their predicted protein impact, and export results.
filter_and_export_impact_diff_events <- function(
    dpsi_df,
    conditions,
    impact_df,
    dpsi_threshold     = 10,
    col_splicing       = "Splicing_change",
    col_event_type     = "EVENT_TYPE",
    col_impact         = "IMPACT",
    env_all_name       = "impact_diff_events_all",
    env_filtered_name  = "impact_diff_events_filtered",
    output_dir         = "Files_dir",
    file_base_name     = "Functional_impact_diff_exons_introns",
    separate_UTR       = "no"
) {
  
  # --- 1. DIRECTORY CONFIGURATION ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  excel_all_path      <- file.path(output_dir, paste0(file_base_name, "_all.xlsx"))
  excel_filtered_path <- file.path(output_dir, paste0(file_base_name, "_filtered.xlsx"))
  rdata_all_path      <- file.path(output_dir, paste0(file_base_name, "_all.RData"))
  rdata_filtered_path <- file.path(output_dir, paste0(file_base_name, "_filtered.RData"))
  
  # --- 2. DEFINE EVENT TYPE ---
  message("Classifying event types...")
  dpsi_annotated <- dpsi_df %>%
    dplyr::mutate(!!rlang::sym(col_event_type) := dplyr::case_when(
      grepl("^HsaEX\\d+", EVENT)  ~ "Exon",
      grepl("^HsaINT\\d+", EVENT) ~ "Intron",
      grepl("^HsaALTD\\d+", EVENT) ~ "Alt5",
      grepl("^HsaALTA\\d+", EVENT) ~ "Alt3"
    )) %>%
    dplyr::select(EVENT, !!rlang::sym(col_event_type), dplyr::everything())
  
  # --- 3. CREATE A LIST OF DF PER CONDITION ---
  message("Creating a list of df with diff spliced events per condition & their predicted impact...")
  all_events_list <- lapply(conditions, function(col_name) {
    
    df_filtered <- dpsi_annotated %>%
      dplyr::mutate(!!rlang::sym(col_splicing) := dplyr::case_when(
        .data[[col_name]] > dpsi_threshold   ~ "Included",
        .data[[col_name]] < -dpsi_threshold  ~ "Excluded"
      )) %>%
      dplyr::filter(abs(.data[[col_name]]) > dpsi_threshold) %>%
      dplyr::arrange(.data[[col_name]]) %>% 
      dplyr::left_join(impact_df, by = "EVENT") 
    
    return(df_filtered)
  })
  names(all_events_list) <- conditions
  
  # --- 4. EXCEL WORKBOOK & RDATA OBJECT FOR ALL (raw impact annotations) ---
  export_excel <- function(data_list, path) {
    wb <- openxlsx::createWorkbook()
    for (sheet_name in names(data_list)) {
      short_name <- substr(sheet_name, 1, 31) 
      openxlsx::addWorksheet(wb, sheetName = short_name)
      openxlsx::writeData(wb, sheet = short_name, x = data_list[[sheet_name]])
    }
    openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
  }
  
  export_rdata <- function(target_object, name_to_save, path) {
    temp_env <- new.env()
    temp_env[[name_to_save]] <- target_object
    save(list = name_to_save, envir = temp_env, file = path)
  }
  
  export_excel(all_events_list, excel_all_path)
  export_rdata(all_events_list, env_all_name, rdata_all_path)
  message("R object and excel workbook with impact annotations saved successfully.")
  
  # --- 5. POLISHING AND FILTERING ---
  message("Simplifying & filtering protein impact definitions ...")
  filtered_events_list <- list()
  
  for (df_name in names(all_events_list)) {
    current_df <- all_events_list[[df_name]]
    
    if (is.null(current_df) || nrow(current_df) == 0) next
    
    # 5a. Apply text polishing on the raw impact column 
    filtered_df <- current_df %>%
      dplyr::mutate(!!rlang::sym(col_impact) := dplyr::case_when(
        grepl("^Alternative protein isoforms", .data[[col_impact]]) ~ "Alt. isoforms",
        grepl("^ORF disruption upon sequence exclusion", .data[[col_impact]]) ~ "ORF disruption upon sequence exclusion",
        grepl("^ORF disruption upon sequence inclusion", .data[[col_impact]]) ~ "ORF disruption upon sequence inclusion",
        grepl("^Protein isoform when splice site is used", .data[[col_impact]]) ~ "Protein isoform when splice site is used",
        TRUE ~ dplyr::coalesce(as.character(.data[[col_impact]]), "NA")
      )) %>%
      
      # 5b. Restrict exclusively to Exons and Introns
      dplyr::filter(.data[[col_event_type]] %in% c("Exon", "Intron")) %>%
      
      # Consider splicing change to associated impact
      dplyr::filter(
        .data[[col_impact]] %in% c("3 UTR", "5 UTR", "NonCoding", "Alt. isoforms", "In the CDS, with uncertain impact", "NA") |
          (.data[[col_splicing]] == "Excluded" & .data[[col_impact]] == "ORF disruption upon sequence exclusion") |
          (.data[[col_splicing]] == "Included" & .data[[col_impact]] == "ORF disruption upon sequence inclusion")
      ) %>%
      
      # 5c. Aggroupate impact categories
      dplyr::mutate(
        !!rlang::sym(col_impact) := stringr::str_replace(.data[[col_impact]], "ORF disruption upon sequence.*", "ORF disruption"),
        
        !!rlang::sym(col_impact) := dplyr::if_else(.data[[col_impact]] %in% c("NA", "In the CDS, with uncertain impact"), 
                                                   "NA, Uncertain", 
                                                   .data[[col_impact]])
      )
    
    # Conditionally combine UTR/NonCoding categories only if separate_UTR == "no"
    if (separate_UTR == "no") {
      filtered_df <- filtered_df %>%
        dplyr::mutate(
          !!rlang::sym(col_impact) := dplyr::if_else(.data[[col_impact]] %in% c("3 UTR", "5 UTR", "NonCoding"), 
                                                     "UTR, NonCoding", 
                                                     .data[[col_impact]])
        )
    }
    
    filtered_events_list[[df_name]] <- filtered_df
  }
  
  # --- 6. EXCEL & RDATA WORKBOOKS FOR FILTERED IMPACT OF EXONS/INTRONS ---
  export_excel(filtered_events_list, excel_filtered_path)
  export_rdata(filtered_events_list, env_filtered_name, rdata_filtered_path)
  
  if (separate_UTR == "no") {
    message("R object and excel workbook with filtered impact annotations:\nORF disruption, UTR/NonCoding, NA/uncertain of diff exons/introns saved successfully")
  } else {
    message("R object and excel workbook with filtered impact annotations:\nORF disruption, 3 UTR, 5 UTR, NonCoding, NA/uncertain of diff exons/introns saved successfully")
  }
  
  # --- 7. list of df into the global environment ---
  assign(env_all_name, all_events_list, envir = .GlobalEnv)
  assign(env_filtered_name, filtered_events_list, envir = .GlobalEnv)
  
  message(paste0("\n>>> END!.\n",
                 "Variables '", env_all_name, "' containing the raw impact annotations and\n'", env_filtered_name,
                 "' containing the filtered ones are now loaded in your active Global Environment."))
}



# Function to run diff gene expression with edgeR
run_edger_comparisons <- function(
    count_data,
    gene_info,
    sample_groups,
    comparison_table,
    files_dir,
    graphs_dir,
    output_name,
    stats = TRUE,
    min_count = 20,
    min_prop = 1,
    logFC_threshold = 1,
    FDR_threshold = 0.05,
    prior_count = 1,
    robust = TRUE,
    save_files = TRUE,
    save_graphs = TRUE
) {
  
  
  ################################
  ### Input checks             ###
  ################################
  
  if (!all(c("sample", "condition") %in% colnames(sample_groups))) {
    stop("sample_groups must contain columns: sample and condition.")
  }
  
  if (!all(c("condition", "control") %in% colnames(comparison_table))) {
    stop("comparison_table must contain columns: condition and control.")
  }
  
  if (!is.logical(stats) || length(stats) != 1) {
    stop("stats must be TRUE or FALSE.")
  }
  
  if (anyDuplicated(sample_groups$sample)) {
    stop("Duplicated samples found in sample_groups.")
  }
  
  if (is.null(rownames(count_data))) {
    stop("count_data must have gene IDs as row names.")
  }
  
  if (is.null(rownames(gene_info))) {
    stop("gene_info must have gene IDs as row names.")
  }
  
  # if (ncol(gene_info) < 2) {
  #   stop("gene_info must contain at least two columns: ID and NAME.")
  # }
  
  missing_metadata <- setdiff(
    colnames(count_data),
    sample_groups$sample
  )
  
  if (length(missing_metadata) > 0) {
    stop(
      "Samples missing from sample_groups: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  
  unknown_samples <- setdiff(
    sample_groups$sample,
    colnames(count_data)
  )
  
  if (length(unknown_samples) > 0) {
    stop(
      "Samples in sample_groups but absent from count_data: ",
      paste(unknown_samples, collapse = ", ")
    )
  }
  
  sample_groups <- sample_groups[
    match(colnames(count_data), sample_groups$sample),
    ,
    drop = FALSE
  ]
  
  requested_conditions <- unique(
    c(
      comparison_table$condition,
      comparison_table$control
    )
  )
  
  missing_conditions <- setdiff(
    requested_conditions,
    unique(sample_groups$condition)
  )
  
  if (length(missing_conditions) > 0) {
    stop(
      "Conditions not found in sample_groups: ",
      paste(missing_conditions, collapse = ", ")
    )
  }
  
  if (!all(rownames(count_data) %in% rownames(gene_info))) {
    stop("Some count_data genes are absent from gene_info.")
  }
  
  gene_info <- gene_info[
    rownames(count_data),
    ,
    drop = FALSE
  ]
  
  dir.create(
    files_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    graphs_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ################################
  ### Helper functions         ###
  ################################
  
  make_sheet_names <- function(x) {
    
    x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
    x <- substr(x, 1, 31)
    
    make.unique(x, sep = "_")
  }
  
  format_threshold <- function(x) {
    
    formatted <- format(
      x,
      scientific = FALSE,
      trim = TRUE
    )
    
    gsub("\\.", "", formatted)
  }
  
  fold_change_threshold <- 2^logFC_threshold
  
  fc_label <- format_threshold(
    fold_change_threshold
  )
  
  fdr_label <- format_threshold(
    FDR_threshold
  )
  
  if (stats) {
    
    filtered_name <- paste0(
      output_name,
      "_fc",
      fc_label,
      "_fdr",
      fdr_label
    )
    
  } else {
    
    filtered_name <- paste0(
      output_name,
      "_fc",
      fc_label
    )
  }
  
  
  ################################
  ### Library-size barplot     ###
  ################################
  
  library_sizes <- data.frame(
    sample = colnames(count_data),
    condition = sample_groups$condition,
    total_counts = colSums(count_data),
    stringsAsFactors = FALSE
  )
  
  if (save_graphs) {
    
    svg(
      filename = file.path(
        graphs_dir,
        "Barplot_gene_counts_per_sample.svg"
      ),
      width = 11,
      height = 7
    )
    
    par(mar = c(10, 5, 3, 1))
    
    barplot(
      height = library_sizes$total_counts,
      names.arg = library_sizes$sample,
      las = 2,
      ylab = "Total raw counts",
      main = "Total gene counts per sample"
    )
    
    dev.off()
  }
  
  
  ################################
  ### Global TMM normalization ###
  ################################
  
  all_groups <- factor(
    sample_groups$condition
  )
  
  dge_all <- edgeR::DGEList(
    counts = count_data,
    genes = gene_info,
    group = all_groups
  )
  
  keep_all <- edgeR::filterByExpr(
    dge_all,
    group = all_groups,
    min.count = min_count,
    min.prop = min_prop
  )
  
  dge_all <- dge_all[
    keep_all,
    ,
    keep.lib.sizes = FALSE
  ]
  
  dge_all <- edgeR::calcNormFactors(
    dge_all,
    method = "TMM"
  )
  
  logCPM_all <- edgeR::cpm(
    dge_all,
    log = TRUE,
    prior.count = prior_count
  )
  
  
  ################################
  ### PCA                      ###
  ################################
  
  pca_model <- prcomp(
    t(logCPM_all),
    center = TRUE,
    scale. = FALSE
  )
  
  variance_explained <- (
    pca_model$sdev^2 /
      sum(pca_model$sdev^2)
  ) * 100
  
  pca <- data.frame(
    sample = rownames(pca_model$x),
    pca_model$x,
    stringsAsFactors = FALSE
  )
  
  pca <- merge(
    pca,
    sample_groups,
    by = "sample",
    sort = FALSE
  )
  
  pca <- pca[
    match(colnames(count_data), pca$sample),
    ,
    drop = FALSE
  ]
  
  pca_df <- data.frame(
    Gene_ID = rownames(logCPM_all),
    logCPM_all,
    check.names = FALSE
  )
  
  if (save_graphs) {
    
    svg(
      filename = file.path(
        graphs_dir,
        "PCA_gene_expression.svg"
      ),
      width = 8,
      height = 7
    )
    
    condition_factor <- factor(
      pca$condition
    )
    
    plot(
      pca$PC1,
      pca$PC2,
      pch = 19,
      col = as.integer(condition_factor),
      xlab = paste0(
        "PC1 (",
        round(variance_explained[1], 1),
        "%)"
      ),
      ylab = paste0(
        "PC2 (",
        round(variance_explained[2], 1),
        "%)"
      ),
      main = "PCA of TMM-normalized log2-CPM"
    )
    
    text(
      pca$PC1,
      pca$PC2,
      labels = pca$sample,
      pos = 3,
      cex = 0.7
    )
    
    legend(
      "topright",
      legend = levels(condition_factor),
      col = seq_along(levels(condition_factor)),
      pch = 19,
      cex = 0.7
    )
    
    dev.off()
  }
  
  
  ################################
  ### edgeR MDS plot           ###
  ################################
  
  if (save_graphs) {
    
    svg(
      filename = file.path(
        graphs_dir,
        "MDS_gene_expression.svg"
      ),
      width = 8,
      height = 7
    )
    
    mds_groups <- factor(
      sample_groups$condition
    )
    
    limma::plotMDS(
      dge_all,
      labels = sample_groups$sample,
      col = as.integer(mds_groups),
      main = "edgeR MDS plot"
    )
    
    legend(
      "topright",
      legend = levels(mds_groups),
      col = seq_along(levels(mds_groups)),
      pch = 16,
      cex = 0.7
    )
    
    dev.off()
  }
  
  
  ################################
  ### Differential expression  ###
  ################################
  
  deg <- list()
  analysis_info <- list()
  
  for (i in seq_len(nrow(comparison_table))) {
    
    treatment <- comparison_table$condition[i]
    control <- comparison_table$control[i]
    
    comparison_name <- paste0(
      treatment,
      "_vs_",
      control
    )
    
    message(
      "Processing comparison ",
      i,
      "/",
      nrow(comparison_table),
      ": ",
      comparison_name
    )
    
    selected <- sample_groups$condition %in%
      c(control, treatment)
    
    selected_samples <- sample_groups$sample[selected]
    selected_conditions <- sample_groups$condition[selected]
    
    counts_subset <- count_data[
      ,
      selected_samples,
      drop = FALSE
    ]
    
    groups_subset <- factor(
      selected_conditions,
      levels = c(control, treatment)
    )
    
    names(groups_subset) <- selected_samples
    
    replicate_counts <- table(
      groups_subset
    )
    
    dge_sub <- edgeR::DGEList(
      counts = counts_subset,
      genes = gene_info[
        rownames(counts_subset),
        ,
        drop = FALSE
      ],
      group = groups_subset
    )
    
    keep <- edgeR::filterByExpr(
      dge_sub,
      group = groups_subset,
      min.count = min_count,
      min.prop = min_prop
    )
    
    dge_sub <- dge_sub[
      keep,
      ,
      keep.lib.sizes = FALSE
    ]
    
    if (nrow(dge_sub) == 0) {
      
      warning(
        "No genes passed filtering for ",
        comparison_name
      )
      
      deg[[comparison_name]] <- data.frame(
        ID = character(),
        NAME = character(),
        logFC = numeric(),
        logCPM = numeric(),
        F = numeric(),
        PValue = numeric(),
        FDR = numeric(),
        stringsAsFactors = FALSE
      )
      
      analysis_info[[comparison_name]] <- list(
        treatment = treatment,
        control = control,
        replicate_counts = replicate_counts,
        analysis_type = "no_genes_passed_filter"
      )
      
      next
    }
    
    dge_sub <- edgeR::calcNormFactors(
      dge_sub,
      method = "TMM"
    )
    
    run_statistics <- (
      stats &&
        all(replicate_counts >= 2)
    )
    
    if (run_statistics) {
      
      design_sub <- model.matrix(
        ~0 + groups_subset
      )
      
      colnames(design_sub) <- c(
        "control",
        "treatment"
      )
      
      dge_sub <- edgeR::estimateDisp(
        dge_sub,
        design_sub,
        robust = robust
      )
      
      fit_sub <- edgeR::glmQLFit(
        dge_sub,
        design_sub,
        robust = robust
      )
      
      test_sub <- edgeR::glmQLFTest(
        fit_sub,
        contrast = c(-1, 1)
      )
      
      result_df <- edgeR::topTags(
        test_sub,
        n = Inf,
        sort.by = "PValue"
      )$table
      
      result_df$ID <- rownames(result_df)
      
      result_df$NAME <- gene_info[
        result_df$ID,
        2
      ]
      
      result_df <- dplyr::select(
        result_df,
        ID,
        NAME,
        logFC,
        logCPM,
        F,
        PValue,
        FDR
      )
      
      rownames(result_df) <- NULL
      
      analysis_type <- "edgeR_QLF"
      
    } else {
      
      if (stats) {
        message(
          "  Insufficient replication: calculating log2FC only."
        )
      } else {
        message(
          "  stats = FALSE: calculating log2FC only."
        )
      }
      
      normalized_logCPM <- edgeR::cpm(
        dge_sub,
        log = TRUE,
        prior.count = prior_count
      )
      
      control_columns <- groups_subset == control
      treatment_columns <- groups_subset == treatment
      
      mean_control <- rowMeans(
        normalized_logCPM[
          ,
          control_columns,
          drop = FALSE
        ]
      )
      
      mean_treatment <- rowMeans(
        normalized_logCPM[
          ,
          treatment_columns,
          drop = FALSE
        ]
      )
      
      gene_ids <- rownames(
        normalized_logCPM
      )
      
      result_df <- data.frame(
        ID = gene_ids,
        NAME = gene_info[
          gene_ids,
          2
        ],
        logFC = mean_treatment - mean_control,
        logCPM = rowMeans(normalized_logCPM),
        F = NA_real_,
        PValue = NA_real_,
        FDR = NA_real_,
        stringsAsFactors = FALSE
      )
      
      result_df <- result_df[
        order(
          abs(result_df$logFC),
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]
      
      rownames(result_df) <- NULL
      
      analysis_type <- "logFC_only"
      design_sub <- NULL
      fit_sub <- NULL
      test_sub <- NULL
    }
    
    deg[[comparison_name]] <- result_df
    
    analysis_info[[comparison_name]] <- list(
      treatment = treatment,
      control = control,
      treatment_samples = names(groups_subset)[
        groups_subset == treatment
      ],
      control_samples = names(groups_subset)[
        groups_subset == control
      ],
      replicate_counts = replicate_counts,
      analysis_type = analysis_type,
      normalization_factors = dge_sub$samples$norm.factors,
      library_sizes = dge_sub$samples$lib.size,
      genes_before_filtering = nrow(counts_subset),
      genes_after_filtering = nrow(dge_sub),
      design = design_sub,
      fit = fit_sub,
      test = test_sub
    )
  }
  
  
  ################################
  ### Preserve comparison order ###
  ################################
  
  comparison_order <- paste0(
    comparison_table$condition,
    "_vs_",
    comparison_table$control
  )
  
  deg <- deg[
    comparison_order
  ]
  
  analysis_info <- analysis_info[
    comparison_order
  ]
  
  
  ################################
  ### Differential gene filter ###
  ################################
  
  deg_diff <- lapply(
    names(deg),
    function(comparison_name) {
      
      df <- deg[[comparison_name]]
      
      if (nrow(df) == 0) {
        return(df)
      }
      
      has_statistics <- (
        analysis_info[[comparison_name]]$analysis_type ==
          "edgeR_QLF"
      )
      
      if (stats && has_statistics) {
        
        df[
          !is.na(df$logFC) &
            !is.na(df$FDR) &
            abs(df$logFC) >= logFC_threshold &
            df$FDR <= FDR_threshold,
          ,
          drop = FALSE
        ]
        
      } else {
        
        df[
          !is.na(df$logFC) &
            abs(df$logFC) >= logFC_threshold,
          ,
          drop = FALSE
        ]
      }
    }
  )
  
  names(deg_diff) <- names(deg)
  
  deg_diff <- deg_diff[
    comparison_order
  ]
  
  
  ################################
  ### Save result files        ###
  ################################
  
  output_files <- NULL
  
  if (save_files) {
    
    deg_rdata_file <- file.path(
      files_dir,
      paste0(output_name, ".RData")
    )
    
    deg_excel_file <- file.path(
      files_dir,
      paste0(output_name, ".xlsx")
    )
    
    diff_rdata_file <- file.path(
      files_dir,
      paste0(filtered_name, ".RData")
    )
    
    diff_excel_file <- file.path(
      files_dir,
      paste0(filtered_name, ".xlsx")
    )
    
    save(
      deg,
      analysis_info,
      pca,
      pca_model,
      pca_df,
      logCPM_all,
      library_sizes,
      file = deg_rdata_file
    )
    
    save(
      deg_diff,
      analysis_info,
      file = diff_rdata_file
    )
    
    deg_excel <- deg
    
    names(deg_excel) <- make_sheet_names(
      names(deg_excel)
    )
    
    deg_diff_excel <- deg_diff
    
    names(deg_diff_excel) <- make_sheet_names(
      names(deg_diff_excel)
    )
    
    openxlsx::write.xlsx(
      deg_excel,
      file = deg_excel_file,
      overwrite = TRUE,
      rowNames = FALSE
    )
    
    openxlsx::write.xlsx(
      deg_diff_excel,
      file = diff_excel_file,
      overwrite = TRUE,
      rowNames = FALSE
    )
    
    output_files <- c(
      DEG_RData = deg_rdata_file,
      DEG_Excel = deg_excel_file,
      filtered_RData = diff_rdata_file,
      filtered_Excel = diff_excel_file
    )
  }
  
  
  ################################
  ### Return results           ###
  ################################
  
  return(
    list(
      deg = deg,
      deg_diff = deg_diff,
      pca = pca,
      pca_model = pca_model,
      logCPM = logCPM_all,
      pca_df = pca_df,
      library_sizes = library_sizes,
      analysis_info = analysis_info,
      output_files = output_files,
      parameters = list(
        stats = stats,
        min_count = min_count,
        min_prop = min_prop,
        logFC_threshold = logFC_threshold,
        fold_change_threshold = fold_change_threshold,
        FDR_threshold = FDR_threshold,
        prior_count = prior_count,
        robust = robust
      )
    )
  )
}