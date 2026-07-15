# change column names

sample_counts_col <- function(p) {
  # draw indices
  idx  <- sample.int(n_feats, size = n_draws, replace = TRUE, prob = p)
  # count occurrences per feature
  tabulate(idx, nbins = n_feats)
}

fix_names <- function(input_cols){
  return(gsub("\\_R1\\|grp1", "", input_cols))
}

# function to remove specific miRNAs from matrix
filter_miRNA_out <- function(mat, remove = NULL) {
  # If 'remove' is NULL or empty, return the matrix unchanged
  if (is.null(remove) || length(remove) == 0) {
    return(mat)
  }
  
  # Otherwise, keep only rows not in 'remove'
  keep_rows <- !rownames(mat) %in% remove
  mat[keep_rows, , drop = FALSE]
}

extract_mouse <- function(x) {
  m <- regexpr("\\d+(?=[SL])", x, perl = TRUE)  # digits immediately before S or L
  out <- ifelse(m > 0, regmatches(x, m), NA_character_)
  as.integer(out)
}

filter_mice <- function(RC_df, column_data, mice_to_drop) {
  # Ensure the column order matches between RC_df and column_data
  common_cols <- intersect(colnames(RC_df), rownames(column_data))
  RC_df <- RC_df[, common_cols, drop = FALSE]
  column_data <- column_data[common_cols, , drop = FALSE]
  
  # Identify columns (samples) to keep
  keep_cols <- rownames(column_data)[!column_data$mouse %in% mice_to_drop]
  
  # Filter both dataframes
  RC_df_filtered <- RC_df[, keep_cols, drop = FALSE]
  column_data_filtered <- column_data[keep_cols, , drop = FALSE]
  
  # Return both filtered dataframes as a list
  return(list(
    RC_df = RC_df_filtered,
    column_data = column_data_filtered
  ))
}



# define groups
define_groups <- function(input_df) {
  groups <- list()
  for (x in c(20, 40, 60, 80)) {
    current_group <- as.character(x)
    print(current_group)
    # columns <- grep(paste0("^",current_group), colnames(input_df))
    columns <- grep(paste0("^",current_group, "(?!.*_2$)"), colnames(input_df), perl = TRUE)
    members <- colnames(input_df)[columns]
    print(members)
    groups[[paste0(current_group, "S")]] <- members
    
  }
  columns <- grep("^[1-5]L(?!.*_.*_)", colnames(input_df), perl = TRUE)
  groups[["0S"]] <- colnames(input_df)[columns]
  columns <- grep("^[1-5]S(?!.*_2$)", colnames(input_df), perl = TRUE)
  groups[["100S"]] <- colnames(input_df)[columns]
  group_list <- names(groups)
  return(list(group_members = groups, group_list = group_list))
}

get_group <- function(entry, group_list) {
  entry <- trimws(entry)
  
  # Special case
  if (entry == "1L_2L_3L_4L_5L") return("0S")
  
  # First: replace '_2' with '_1'
  entry_alt1 <- sub("_2$", "_1", entry)
  for (name in names(group_list)) {
    if (entry_alt1 %in% group_list[[name]]) return(name)
  }
  
  # Second: remove '_2' altogether
  entry_alt2 <- sub("_2$", "", entry)
  if (entry_alt2 != entry_alt1) {
    for (name in names(group_list)) {
      if (entry_alt2 %in% group_list[[name]]) return(name)
    }
  }
  
  return(NA_character_)
}

# function to aggregate per group

geometric_mean_by_group <- function(data, group_list) {
  # Initialize an empty list to store the results
  group_means <- list()
  
  # Iterate over each group in the list
  for (group_name in names(group_list)) {
    # Get the column names for the current group
    cols <- group_list[[group_name]]
    
    # Extract the relevant columns
    group_data <- data[, cols, drop = FALSE]
    
    # Calculate the geometric mean across rows
    geom_mean <- apply(group_data, 1, function(row) {
      exp(mean(log(row + 1), na.rm = TRUE)) # Calculate geometric mean for each row
    })
    
    # Store the result
    group_means[[group_name]] <- geom_mean
  }
  
  # Combine the results into a new dataframe
  result_df <- as.data.frame(group_means)
  colnames(result_df) <- names(group_list)
  desired_order <- c("0S", "20S", "40S", "60S", "80S", "100S")
  result_df <- result_df[, desired_order]
  
  # Return the result
  return(result_df)
}

# check if a row is monotonic

check_monotonicity <- function(x) {
  if (all(diff(x) > 0)) {
    return("increasing")
  } else if (all(diff(x) < 0)) {
    return("decreasing")
  } else {
    return("neither")
  }
}





# perform rlog

perform_rlog <- function(dds) {
  # Ensure the input is a DESeqDataSet object
  if (!inherits(dds, "DESeqDataSet")) {
    stop("Input must be a DESeqDataSet object")
  }
  
  # Define parameter combinations
  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  blinds <- c(TRUE, FALSE)
  
  # Initialize an empty list to store results
  rlog_results <- list()
  
  # Iterate over all combinations
  for (fit_type in fit_types) {
    for (blind in blinds) {
      # Generate a descriptive key name for this combination
      key_name <- paste0("fitType_", fit_type, "_blind_", blind)
      
      # Perform rlog transformation
      rlog_results[[key_name]] <- tryCatch({
        rlog(dds, blind = blind, fitType = fit_type)
      }, error = function(e) {
        message(paste("Error in combination:", key_name, "-", e$message))
        return(NULL)
      })
    }
  }
  
  # Return the list of results
  return(rlog_results)
}

# perform vst
perform_vst <- function(dds) {
  # Ensure the input is a DESeqDataSet object
  if (!inherits(dds, "DESeqDataSet")) {
    stop("Input must be a DESeqDataSet object")
  }
  
  # Define parameter combinations
  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  blinds <- c(TRUE, FALSE)
  
  # Initialize an empty list to store results
  vst_results <- list()
  
  # Iterate over all combinations
  for (fit_type in fit_types) {
    for (blind in blinds) {
      # Generate a descriptive key name for this combination
      key_name <- paste0("fitType_", fit_type, "_blind_", blind)
      
      # Perform VST
      vst_results[[key_name]] <- tryCatch({
        varianceStabilizingTransformation(dds, blind = blind, fitType = fit_type)
      }, error = function(e) {
        message(paste("Error in combination:", key_name, "-", e$message))
        return(NULL)
      })
    }
  }
  
  # Return the list of results
  return(vst_results)
}

# MAP estimator function
# MAP Estimator Function
# Implements: π̂_ij^MAP = w_j * (1/I) + (1-w_j) * π̂_ij^MLE
# where w_j = (I * α_j) / (I * α_j + n_j)

map_estimator <- function(pi_mle, I, alpha_j, n_j) {
  # Calculate the weight w_j
  w_j <- (I * alpha_j) / (I * alpha_j + n_j)
  
  # Calculate the MAP estimate
  pi_map <- w_j * (1/I) + (1 - w_j) * pi_mle
  
  return(pi_map)
}

# # Example usage:
# pi_mle <- 0.01    # MLE estimate (can be a vector of MLEs for multiple miRNAs)
# I <- 500          # Total number of miRNAs measured
# alpha_j <- 1      # Prior parameter for sample j (treat as a tuning parameter)
# n_j <- 1e6         # Total number of reads in sample j
# 
# result <- map_estimator(pi_mle, I, alpha_j, n_j)
# print(result)

generate_monotonic_heatmap_old <- function(
    list_of_transformations,
    groups,
    rpm_grp,
    name_map,
    plot_dir,
    topN = 100,
    sort = "default"
) {
  # 0) setup
  RAW_DF     <- rpm_grp
  average_RPM <- rowMeans(RAW_DF, na.rm = TRUE)
  
  # 1) pick topN by RPM
  top100 <- names(sort(average_RPM, decreasing = TRUE))[seq_len(min(topN, length(average_RPM)))]
  
  # 2) monotonicity for each normalization
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp100 <- geometric_mean_by_group(
      list_of_transformations[[key]], groups
    )[top100, , drop = FALSE]
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = apply(grp100, 1, check_monotonicity),
      stringsAsFactors = FALSE
    )
  })
  heat_df <- dplyr::bind_rows(mono_list)
  
  # 3) add log2FC trend
  idx    <- match(top100, rownames(RAW_DF))
  fc_vec <- log2(RAW_DF[idx, "100S"] / RAW_DF[idx, "0S"])
  fc_df  <- data.frame(
    miRNA          = top100,
    transformation = "log2FC",
    monotonicity   = ifelse(
      fc_vec >  0, "increasing",
      ifelse(fc_vec <  0, "decreasing", "none")
    ),
    stringsAsFactors = FALSE
  )
  heat_df <- dplyr::bind_rows(heat_df, fc_df)
  
  # 4) factor levels
  heat_df <- heat_df %>%
    dplyr::mutate(
      miRNA = factor(miRNA, levels = top100),
      transformation = factor(
        transformation,
        levels = c(names(list_of_transformations), "log2FC")
      )
    )
  
  # 5) build heatmap
  monotone_cols <- c(
    increasing = "#1b9e77",
    decreasing = "#d95f02",
    none       = "#FFFF00"
  )
  gg <- ggplot2::ggplot(
    heat_df,
    aes(transformation, miRNA, fill = monotonicity)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      name   = "Trend",
      values = monotone_cols,
      labels = c(
        increasing = "Increasing",
        decreasing = "Decreasing",
        none       = "Non-monotonic"
      ),
      na.value = "white"
    ) +
    ggplot2::scale_x_discrete(
      labels = c(name_map, log2FC = "log2FC"),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x     = "Normalization method",
      y     = paste0("miRNA (top ", length(top100), " by RPM)"),
      title = "Per-miRNA monotonic trend + log2FC"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 6),
      panel.grid.major = element_line(color = "black")
    )
  
  # 6) save to PDF
  out_dir <- if (length(plot_dir) > 1) plot_dir[1] else plot_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  pdf(file.path(out_dir, "heatmap_wFC.pdf"), width = 8, height = 10)
  print(gg)
  dev.off()
  
  message("done")
  invisible(gg)
}

generate_monotonic_heatmap <- function(
    list_of_transformations,
    groups,
    rpm_grp,
    name_map,
    plot_dir,
    topN = 100,
    sort_by = "default"
) {
  # 0) setup
  RAW_DF     <- rpm_grp
  average_RPM <- rowMeans(RAW_DF, na.rm = TRUE)
  
  # 1) pick topN by RPM
  top100 <- names(sort(average_RPM, decreasing = TRUE))[seq_len(min(topN, length(average_RPM)))]
  
  # 2) monotonicity for each normalization
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp100 <- geometric_mean_by_group(
      list_of_transformations[[key]], groups
    )[top100, , drop = FALSE]
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = apply(grp100, 1, check_monotonicity),
      stringsAsFactors = FALSE
    )
  })
  heat_df <- dplyr::bind_rows(mono_list)
  
  # 3) add log2FC trend
  idx    <- match(top100, rownames(RAW_DF))
  fc_vec <- log2(RAW_DF[idx, "100S"] / RAW_DF[idx, "0S"])
  fc_df  <- data.frame(
    miRNA          = top100,
    transformation = "log2FC",
    monotonicity   = ifelse(
      fc_vec >  0, "increasing",
      ifelse(fc_vec <  0, "decreasing", "none")
    ),
    logFC = fc_vec,
    stringsAsFactors = FALSE
  )
  heat_df <- dplyr::bind_rows(heat_df, fc_df)
  
  #3.1 sort by logFC if specified
  if (identical(sort_by, "logFC")) {
    # order by decreasing logFC
    print("option logfc")
    ordered_miRNA <- fc_df %>%
      dplyr::arrange(dplyr::desc(logFC)) %>%
      dplyr::pull(miRNA)
  } else {
    # default RPM order
    ordered_miRNA <- top100
  }
  
  # 4) factor levels with chosen order
  heat_df <- heat_df %>%
    dplyr::mutate(
      miRNA = factor(miRNA, levels = ordered_miRNA),
      transformation = factor(
        transformation,
        levels = c(names(list_of_transformations), "log2FC")
      )
    )
  
  # 5) build heatmap
  monotone_cols <- c(
    increasing = "#1b9e77",
    decreasing = "#d95f02",
    none       = "#FFFF00"
  )
  gg <- ggplot2::ggplot(
    heat_df,
    aes(transformation, miRNA, fill = monotonicity)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      name   = "Trend",
      values = monotone_cols,
      labels = c(
        increasing = "Increasing",
        decreasing = "Decreasing",
        none       = "Non-monotonic"
      ),
      na.value = "white"
    ) +
    ggplot2::scale_x_discrete(
      labels = c(name_map, log2FC = "log2FC"),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x     = "Normalization method",
      y     = paste0("miRNA (top ", length(top100), " by RPM)"),
      title = "Per-miRNA monotonic trend + log2FC"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 6),
      panel.grid.major = element_line(color = "black")
    )
  
  # 6) save to PDF
  out_dir <- if (length(plot_dir) > 1) plot_dir[1] else plot_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  pdf(file.path(out_dir, "heatmap_FC_sorted.pdf"), width = 8, height = 10)
  print(gg)
  dev.off()
  
  message("done")
  invisible(gg)
  print(ordered_miRNA)
}

generate_correlation_heatmap_old <- function(
    list_of_transformations,
    groups,
    rpm_grp,
    name_map,
    plot_dir,
    ref_miRNA,
    topN = 100
) {
  # 0) Validate that ref_miRNA exists in raw RPM data
  if (!ref_miRNA %in% rownames(rpm_grp)) {
    stop("Reference miRNA '", ref_miRNA, "' not found in rpm_grp.")
  }
  
  # 1) Prepare output directory
  out_dir <- plot_dir[1]
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 2) Select topN miRNAs by raw RPM mean
  avg_rpm <- rowMeans(rpm_grp, na.rm = TRUE)
  miRNAs  <- names(sort(avg_rpm, decreasing = TRUE))[seq_len(min(topN, length(avg_rpm)))]
  
  # 3) Compute Pearson correlation for each transformation
  corr_list <- lapply(names(list_of_transformations), function(key) {
    mat <- list_of_transformations[[key]]
    # Ensure reference miRNA exists in this transformation
    if (!ref_miRNA %in% rownames(mat)) {
      stop("Reference miRNA '", ref_miRNA, "' not found in transformation '", key, "'.")
    }
    # Subset to top miRNAs and all samples
    mat_sub     <- mat[miRNAs, , drop = FALSE]
    ref_vec     <- as.numeric(mat[ref_miRNA, , drop = TRUE])
    # Compute correlations, coercing each profile to numeric
    cors <- apply(mat_sub, 1, function(vals) {
      stats::cor(as.numeric(vals), ref_vec, use = "pairwise.complete.obs")
    })
    data.frame(
      miRNA          = miRNAs,
      transformation = key,
      correlation    = cors,
      stringsAsFactors = FALSE
    )
  })
  corr_df <- dplyr::bind_rows(corr_list)
  
  # 4) Factor levels for plotting
  corr_df <- corr_df %>%
    dplyr::mutate(
      miRNA          = factor(miRNA, levels = rev(miRNAs)),
      transformation = factor(transformation, levels = names(list_of_transformations))
    )
  
  # 5) Build heatmap
  heatmap_plot <- ggplot2::ggplot(
    corr_df,
    ggplot2::aes(x = transformation, y = miRNA, fill = correlation)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      name     = "Pearson R",
      low      = "#d73027",
      mid      = "#ffffbf",
      high     = "#1a9850",
      midpoint = 0,
      na.value = "white"
    ) +
    ggplot2::scale_x_discrete(labels = name_map, expand = c(0, 0)) +
    ggplot2::labs(
      x     = "Normalization method",
      y     = paste0("miRNA (top ", length(miRNAs), " by RPM)"),
      title = paste0("Correlation to ", ref_miRNA)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y      = ggplot2::element_text(size = 6),
      panel.grid.major = ggplot2::element_line(color = "white")
    )
  
  # 6) Save to PDF
  out_file <- file.path(out_dir, paste0("corr_to_", ref_miRNA, ".pdf"))
  ggplot2::ggsave(out_file, heatmap_plot, width = 8, height = 10)
  
  message("Correlation heatmap saved to: ", out_file)
  invisible(heatmap_plot)
}

generate_correlation_heatmap <- function(
    list_of_transformations,
    groups,
    rpm_grp,
    name_map,
    plot_dir,
    ref_miRNA,
    # ref_transform = "RPM_lib",
    ref_transform = "rlog_parametric_blind_TRUE",
    topN = 100
) {
  # 0) Validate reference transformation and miRNA
  if (!ref_transform %in% names(list_of_transformations)) {
    stop("Reference transformation '", ref_transform, "' not found in list_of_transformations.")
  }
  ref_mat <- list_of_transformations[[ref_transform]]
  if (!ref_miRNA %in% rownames(ref_mat)) {
    stop("Reference miRNA '", ref_miRNA, "' not found in transformation '", ref_transform, "'.")
  }
  
  # 1) Prepare output directory
  out_dir <- plot_dir[1]
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 2) Select topN miRNAs by raw RPM mean
  avg_rpm <- rowMeans(rpm_grp, na.rm = TRUE)
  miRNAs  <- names(sort(avg_rpm, decreasing = TRUE))[seq_len(min(topN, length(avg_rpm)))]
  
  # 3) Extract reference profile from chosen transformation
  ref_vec <- as.numeric(ref_mat[ref_miRNA, , drop = TRUE])
  ref_samples <- colnames(ref_mat)
  
  # 4) Compute Pearson correlation for each transformation
  corr_list <- lapply(names(list_of_transformations), function(key) {
    mat <- list_of_transformations[[key]]
    # Subset to top miRNAs and common samples
    mat_sub <- mat[miRNAs, ref_samples, drop = FALSE]
    # Compute correlations
    cors <- apply(mat_sub, 1, function(vals) {
      stats::cor(as.numeric(vals), ref_vec, use = "pairwise.complete.obs")
    })
    data.frame(
      miRNA          = miRNAs,
      transformation = key,
      correlation    = cors,
      stringsAsFactors = FALSE
    )
  })
  corr_df <- dplyr::bind_rows(corr_list)
  
  # 5) Factor levels for plotting
  corr_df <- corr_df %>%
    dplyr::mutate(
      miRNA          = factor(miRNA, levels = rev(miRNAs)),
      transformation = factor(transformation, levels = names(list_of_transformations))
    )
  
  # 6) Build heatmap
  heatmap_plot <- ggplot2::ggplot(
    corr_df,
    ggplot2::aes(x = transformation, y = miRNA, fill = correlation)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      name     = "Pearson R",
      low      = "#d73027",
      mid      = "#ffffbf",
      high     = "#1a9850",
      midpoint = 0,
      na.value = "white"
    ) +
    ggplot2::scale_x_discrete(labels = name_map, expand = c(0, 0)) +
    ggplot2::labs(
      x     = "Normalization method",
      y     = paste0("miRNA (top ", length(miRNAs), " by RPM)"),
      title = paste0("Correlation to ", ref_miRNA, " (", ref_transform, ")")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y      = ggplot2::element_text(size = 6),
      panel.grid.major = ggplot2::element_line(color = "white")
    )
  
  # 7) Save to PDF
  out_file <- file.path(
    out_dir,
    paste0("corr_to_", ref_miRNA, "_", ref_transform, ".pdf")
  )
  ggplot2::ggsave(out_file, heatmap_plot, width = 8, height = 10)
  
  message("Correlation heatmap saved to: ", out_file)
  invisible(heatmap_plot)
}

# "RC"
# [1] "RPM_total"
# [1] "RPM_lib"
# [1] "MAP"
# [1] "rlog_parametric_blind_TRUE"
# [1] "rlog_mean_blind_FALSE"
# [1] "vst_local_blind_TRUE"
# [1] "TMM"









# calculate expected values (average)
# calculate expected values (per mouse)