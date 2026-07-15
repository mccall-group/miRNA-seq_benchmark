library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
conflicts_prefer(base::intersect)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(htmlwidgets)
library(plotly)
library(tidyr)
library(ggnewscale)
library(compositions)
library(zCompositions)
library(tibble)
library(ComplexHeatmap)


source("R/normalization.R")
data_dir <- "./data/"
plot_dir <- "./plots/"
dir.create(plot_dir, recursive = TRUE)
# where to write monotonic miRNAs
monotonicity_path <- "../DE_benchmark/input_data/monotonic.csv"
minRPM <- 50
resample_reads <- FALSE

# parameters to enable easy reanalysis
microRNAs_to_drop <- c()
mice_to_drop <- c()


# name map
name_map <- c(
  "RPM_total"                   = "RPM (total)",
  "RPM_lib"                     = "RPM",
  "TMM"                         = "edgeR - TMM",
  "rlog_parametric_blind_TRUE"  = "DESeq2 - rlog (parametric)",
  "rlog_mean_blind_FALSE"       = "DESeq2 - rlog (mean)",
  "vst_local_blind_TRUE"        = "DESeq2 - VST (local)",
  "RC"                          = "Raw read count",
  "MAP" = "MAP"
)

#### create resampled RC (read count) matrix
if (resample_reads){
  # read original RC matrix
  RC_df <- read.delim("./data/mature_sense_minExpr0_RCadj.mat",
                      check.names = FALSE, row.names = 1)
  
  # drop miRNAs
  RC_df <- filter_miRNA_out(RC_df, microRNAs_to_drop)
  
  # count min amount of total reads
  col_sums <- colSums(RC_df)
  min_am  <- min(col_sums) # 12214826, let's sample 10M
  
  # calculate probabilities
  prob_df <- sweep(RC_df, 
                   MARGIN = 2, 
                   STATS  = colSums(RC_df), 
                   FUN    = "/")
  # sample 10M for each column
  n_feats   <- nrow(prob_df)
  n_draws   <- 1e7  # 10 million
  
  sampled_mat <- sapply(seq_len(ncol(prob_df)), function(i) {
    sample_counts_col(prob_df[[i]])
  })
  
  # Convert to data.frame and restore row/col names
  sampled_df <- as.data.frame(sampled_mat,
                              row.names = rownames(prob_df))
  colnames(sampled_df) <- colnames(prob_df)
  
  RC_df <- sampled_df
  
  
}else{
  # load sRNAbench Read Count matrix
  
  RC_df <- read.delim("./data/mature_sense_minExpr0_RCadj.mat",
                      check.names = FALSE, row.names = 1)
  RC_df <- filter_miRNA_out(RC_df, microRNAs_to_drop)
}

colnames(RC_df) <- fix_names(colnames(RC_df))

# separate groups
groups_object <- define_groups(RC_df)

# some mixture_group are missing
column_data <- data.frame(row.names = colnames(RC_df))
column_data$mixture_group <- sapply(rownames(column_data), 
                                    function(x) get_group(x,groups_object$group_members))

column_data$mouse <- vapply(rownames(column_data), extract_mouse, integer(1))

# filter out mice if any is included in the exclusion list
filtered <- filter_mice(RC_df, column_data, mice_to_drop = mice_to_drop)

RC_df <- filtered$RC_df           # RC_df with columns from mouse 3 removed

# log transform (not included in the final analysis)
logged_df <- log2(RC_df + 1)

# RPM total (RPM reads normalized to total mapped reads)
RPMt_df <- read.delim("./data/mature_sense_minExpr0_RCadj_totalRPM.mat",
                      check.names = FALSE, row.names = 1)

RPMt_df <- filter_miRNA_out(RPMt_df, microRNAs_to_drop)

colnames(RPMt_df) <- fix_names(colnames(RPMt_df))
filtered <- filter_mice(RPMt_df, column_data, mice_to_drop = mice_to_drop)
RPMt_df<- filtered$RC_df 
column_data <- filtered$column_data     # Matching metadata

RPMl_df <- sweep(RC_df, 2, colSums(RC_df), FUN = "/") * 1e6

list_of_transformations <- list(RC = RC_df, 
                                # logged = logged_df, 
                                RPM_total = RPMt_df,
                                RPM_lib = RPMl_df
                                )

# monotonicity examples based on experimental data

if (TRUE) {
  plot_dir <- "./plots/"
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Recreate column_data if missing
  if (!exists("column_data")) {
    groups_object <- define_groups(RC_df)
    column_data <- data.frame(row.names = colnames(RC_df))
    column_data$mixture_group <- sapply(
      rownames(column_data),
      function(x) get_group(x, groups_object$group_members)
    )
  }
  column_data$Sample <- rownames(column_data)
  column_data$mixture_group <- gsub("\\s*\\(.*\\)", "", column_data$mixture_group)
  
  # Colors and group order
  inc_col <- "#1f78b4"   # blue
  dec_col <- "#e31a1c"   # red
  nom_col <- "#878787"   # green
  group_order <- c("0S", "20S", "40S", "60S", "80S", "100S")
  
  # Target miRNAs
  dec_miR <- "mmu-miR-122-5p"
  inc_miR <- "mmu-miR-126a-3p"
  # nomono_miR <- "mmu-miR-101a-3p"
  nomono_miR <- "mmu-miR-21a-5p"
  # nomono_miR <- "mmu-miR-99a-5p"
  # nomono_miR <- ""
  if (!exists("nomono_miR")) {
    stop("Please define nomono_miR <- 'mmu-<id>' for the non-monotonic example.", call. = FALSE)
  }
  
  mirnas_use <- c(dec_miR, inc_miR, nomono_miR)
  mirnas_use <- mirnas_use[mirnas_use %in% rownames(RPMl_df)]
  if (length(mirnas_use) < 3) stop("One of the three miRNAs not found in RPMl_df")
  
  # Legend labels with miRNA names
  dec_lab <- paste0("Monotonic decrease \n(", dec_miR, ")")
  inc_lab <- paste0("Monotonic increase \n(", inc_miR, ")")
  nom_lab <- paste0("No monotonicity \n(", nomono_miR, ")")
  
  # Prepare data
  df_long <- RPMl_df[mirnas_use, , drop = FALSE] %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    left_join(column_data[, c("Sample", "mixture_group")], by = "Sample") %>%
    pivot_longer(cols = all_of(mirnas_use), names_to = "miRNA", values_to = "RPM")
  
  df_summary <- df_long %>%
    group_by(mixture_group, miRNA) %>%
    summarise(
      mean_RPM = mean(RPM, na.rm = TRUE),
      sd_RPM   = sd(RPM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(mixture_group = factor(mixture_group, levels = group_order)) %>%
    group_by(miRNA) %>%
    mutate(
      rel_pct = 100 * mean_RPM / max(mean_RPM, na.rm = TRUE),
      rel_sd  = 100 * sd_RPM  / max(mean_RPM, na.rm = TRUE)
    ) %>%
    ungroup()
  
  
  # Position nudges for three side-by-side bars
  nudge_map <- c(
    setNames(-0.28, inc_miR),
    setNames( 0.00, nomono_miR),
    setNames( 0.28, dec_miR)
  )
  
  df_summary$pos_nudge <- unlist(nudge_map[df_summary$miRNA])
  
  # Map each miRNA to its legend class
  df_summary$legend_class <- dplyr::case_when(
    df_summary$miRNA == dec_miR     ~ dec_lab,
    df_summary$miRNA == inc_miR     ~ inc_lab,
    df_summary$miRNA == nomono_miR  ~ nom_lab,
    TRUE ~ "Other"
  )
  
  legend_levels <- c(inc_lab, nom_lab, dec_lab)
  
  df_summary$legend_class <- factor(df_summary$legend_class, levels = legend_levels)
  
  
  # Plot
  gg <- ggplot(df_summary, aes(x = mixture_group, y = rel_pct)) +
    geom_col(
      aes(fill = legend_class),
      position = position_nudge(x = df_summary$pos_nudge),
      width = 0.28
    ) +
    geom_errorbar(
      aes(
        x = mixture_group,
        ymin = rel_pct - rel_sd,
        ymax = rel_pct + rel_sd
      ),
      position = position_nudge(x = df_summary$pos_nudge),
      width = 0.05,
      linewidth = 0.7,
      color = "black"
    ) +
    # scale_fill_manual(
    #   name = NULL,
    #   values = setNames(c(inc_col, nom_col, dec_col), c(inc_lab, nom_lab, dec_lab))
    # ) +
    scale_fill_manual(
      name = NULL,
      values = setNames(c(inc_col, nom_col, dec_col), c(inc_lab, nom_lab, dec_lab)),
      breaks = legend_levels,
      limits = legend_levels
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    
    labs(
      title = "",
      x = NULL,
      y = "Relative expression"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 14),
      plot.title  = element_text(size = 16, hjust = 0.5),
      legend.position = "bottom"
    )
  
  # pdf(file.path(plot_dir, "example_monotonicity_real_relative_3miR.pdf"), width = 8, height = 6)
  # print(gg)
  # dev.off()
  
  gg <- ggplot(df_summary, aes(x = mixture_group, y = rel_pct)) +
    geom_col(
      aes(fill = legend_class),
      position = position_nudge(x = df_summary$pos_nudge),
      width = 0.28
    ) +
    geom_errorbar(
      aes(
        ymin = rel_pct - rel_sd,
        ymax = rel_pct + rel_sd
      ),
      position = position_nudge(x = df_summary$pos_nudge),
      width = 0.05,
      linewidth = 0.7,
      color = "black"
    ) +
    scale_fill_manual(
      name = NULL,
      values = setNames(c(inc_col, nom_col, dec_col), c(inc_lab, nom_lab, dec_lab))
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    coord_cartesian(ylim = c(0, 110)) +   # <<< key change: allow error bars above 100
    labs(
      x = NULL,
      y = "Relative expression"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      plot.title  = element_text(size = 16, hjust = 0.5),
      legend.position = "bottom",
      legend.text  = element_text(size = 14)
      
    )
  pdf(file.path(plot_dir, "example_monotonicity_bars.pdf"), width = 8, height = 6)
  print(gg)
  dev.off()
  
}

### Normalizations/Transformations

# MAP
if (TRUE){
  # 1. Compute library sizes and counts
  lib_sizes <- colSums(RC_df)
  I <- nrow(RC_df)
  alpha <- 1               # uniform prior concentration
  n_j <- lib_sizes
  alpha_j <- rep(alpha, length(n_j))
  
  # 2. Compute MLE proportions
  pi_mle <- sweep(RC_df, 2, lib_sizes, FUN = "/")
  
  # 3. Apply MAP shrinkage column‐wise
  pi_map <- sapply(seq_len(ncol(pi_mle)), function(j) {
    map_estimator(
      pi_mle = pi_mle[, j],
      I      = I,
      alpha_j= alpha_j[j],
      n_j    = n_j[j]
    )
  })
  # restore dimnames
  rownames(pi_map) <- rownames(pi_mle)
  colnames(pi_map) <- colnames(pi_mle)
}
list_of_transformations$MAP <- as.data.frame(pi_map)

dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                              colData = data.frame(row.names = colnames(RC_df)), 
                              design = ~ 1)
                              
# DESEq transformations
if(TRUE){
# if(FALSE){
  # fitType = c("parametric", "local", "mean", "glmGamPoi")
  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  blind_options <- c(TRUE, FALSE)
  
  # Helper function to handle errors
  safe_transform <- function(transform_func, ...) {
    tryCatch(
      expr = transform_func(...),
      error = function(e) {
        message("Error in transformation: ", conditionMessage(e))
        return(NULL)  # Return NULL if an error occurs
      }
    )
  }
  
  # Add rlog transformations
  for (blind in blind_options) {
    if  (blind){
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ 1)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ mixture_group)
    }
    allowed_combinations <- c(
      "rlog_mean_blind_FALSE",
      "rlog_parametric_blind_TRUE",
      "vst_local_blind_TRUE"
    )
    for (fit in fit_types) {
      transformation_name <- paste0("rlog_", fit, "_blind_", blind)
      # print(paste0("rlog_", fit, "_blind_", blind))
      
      # print(head(result))
      
      if (transformation_name %in% allowed_combinations){
        result <- as.data.frame(assay(safe_transform(rlog, dds, fitType = fit, blind = blind)))
        list_of_transformations[[transformation_name]] <- 2^result  
      }
      
    }
  }
  fit_types <- c("parametric", "local", "mean")
  # Add VST transformations
  for (blind in blind_options) {
    if  (blind){
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ 1)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ mixture_group)
    }
    for (fit in fit_types) {
      # print(paste0("vst_", fit, "_blind_", blind))
      
      # print(head(result))
      transformation_name <- paste0("vst_", fit, "_blind_", blind)
      if (transformation_name %in% allowed_combinations){
        result <- as.data.frame(assay(safe_transform(varianceStabilizingTransformation, dds, fitType = fit, blind = blind)))
      list_of_transformations[[transformation_name]] <- 2^result
      }
    }
  }
}

# edgeR transformation TMM
if(TRUE){dge <- DGEList(counts = RC_df)

# Calculate normalization factors using TMM
dge <- calcNormFactors(dge, method = "TMM")

# To extract the TMM-normalized values (on a per-million scale)
tmm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
list_of_transformations[["TMM"]] <- as.data.frame(tmm_counts)
}

### defining mixture groups
group_obj <- define_groups(RC_df)
groups <- group_obj$group_members
group_list <- group_obj$group_list

# calculate group RPM average
if (TRUE){
  group_RPMs <- list()
  for (group in group_list) {
    # Get the samples corresponding to this group
    group_samples <- groups[[group]]
    
    # Subset the input_df to only the columns corresponding to this group
    group_RPM <- RPMl_df[, group_samples]
    
    # Store the RPM values for this group
    group_RPMs[[group]] <- group_RPM
  }
  
  # Initialize an empty matrix to store the average RPMs
  group_avg_RPM <- matrix(NA, nrow = nrow(RPMl_df), ncol = length(group_list))
  rownames(group_avg_RPM) <- rownames(RPMl_df)
  colnames(group_avg_RPM) <- group_list
  
  for (group in group_list) {
    # Extract the RPM data for this group
    group_RPM <- group_RPMs[[group]]
    
    # Calculate the geometric mean for each miRNA (log space, then exponentiate)
    log_RPM <- log2(group_RPM + 1)
    log_means <- rowMeans(log_RPM)
    group_avg_RPM[, group] <- log_means
  }
  
  group_avg_RPM <- 2 ^ group_avg_RPM
  
  group_extended <- group_avg_RPM
  group_extended <- as.data.frame(group_extended)
  
  # Add the average row value column
  group_extended$average <- rowMeans(group_extended)
  
  # Add the FC (fold change) column between 0S and 100S
  group_extended$FC <- group_extended$`100S` / (group_extended$`0S`) 
  
  group_extended$log2FC <- log2(group_extended$FC)
  group_extended$FC_direction <- ifelse(group_extended$log2FC > 0, "increasing", "decreasing")
  group_extended$abs_log2FC <- abs(group_extended$log2FC)
}

#######
##### evaluate monotonic trend for each method
#######

all_filtered <- list()
all_unfiltered <- list()

if (TRUE){
  final_results <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  final_log2 <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    print(key_name)
    temp_df <- list_of_transformations[[i]]
    # make same row order (miRNAs)
    group_df <- geometric_mean_by_group(temp_df, groups)
    results_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = 0))
    rownames(results_df) <- rownames(group_df)
    results_df$monotonicity <- apply(group_df, 1, check_monotonicity)
    # calculate absolute log2FC
    results_df$FC <- group_df$`100S` / (group_df$`0S`)
    results_df$log2FC <- log2(results_df$FC)
    results_df$putative_direction <- ifelse(results_df$log2FC > 0, "increasing", "decreasing")
    results_df$abs_log2FC <- abs(results_df$log2FC)
    results_df$average <- rowMeans(group_df)
    results_df$average_RPM <- group_extended$average
    # calculate agreement
    results_df$agreement <- ifelse(
      results_df$monotonicity %in% c("increasing", "decreasing"),
      1,
      0
    )
    
    # use RPM to filter so it's equivalent to all transformations
    filtered_df <- results_df[results_df$average_RPM >= minRPM, ]
    filtered_df$transformation <- key_name
    all_filtered[[key_name]] <- filtered_df
    all_unfiltered[[key_name]] <- results_df
    # Sort by average_RPM in descending order
    sorted_df <- filtered_df[order(-filtered_df$average), ]
    top_n_values <- c(5, 10, 20, 30, 40, 50, seq(60, nrow(sorted_df), by = 10))
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results df
      parts <- unlist(strsplit(key_name, "_"))
      tgroup = parts[1]  # Return the second keyword
      final_results <- rbind(final_results, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
    
    sorted_df <- filtered_df[order(-filtered_df$abs_log2FC), ]
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results
      parts <- unlist(strsplit(key_name, "_"))
      tgroup = parts[1]  # Return the second keyword
      
      final_log2 <- rbind(final_log2, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  old_group_extended <- group_extended
  
  methods <- unique(final_results$transformation)
  pal_hex <- unname(grDevices::palette.colors())
  pal_hex <- pal_hex[pal_hex != "#000000"]
  cols_use <- pal_hex[seq_along(methods)]
  colour_map <- setNames(cols_use, methods)
  
  
  ## fix levels to mimic order in legend
  
  order_levels <- final_results %>%
    filter(top_n == max(top_n)) %>%               # pick the last x‐value
    arrange(desc(percentage_agreement)) %>%        # sort so highest % is first
    pull(transformation)                           # extract the transformation names
  
  final_results <- final_results %>%
    mutate(transformation = factor(transformation, levels = order_levels))
  
}

# 100 RPM threshold
rpm_threshold <- 100
rpm_threshold_label    <- "100 RPM"
rpm_threshold_position <- sum(group_extended$average >= rpm_threshold, na.rm = TRUE)

# lineplot % monotonic for top n (expression)
gg <- ggplot(final_results, aes(x = top_n, y = percentage_agreement, 
                                color = transformation) ) +
  geom_line() +
  geom_point() +
  # geom_vline(
  #   xintercept = rpm_threshold_position,
  #   linetype = "dashed",
  #   color = "grey70", 
  #   linewidth = 0.6
  # ) +
  # annotate(
  #   "text",
  #   x = rpm_threshold_position,
  #   y = Inf,
  #   label = rpm_threshold_label,
  #   angle = 15,
  #   vjust = 3,
  #   hjust = 0,
  #   size = 4
  # ) +
  labs(x = "Top n miRNAs (sorted by expression)", y = "Percentage with monotonic trend (%)", 
       title = "Proportion of top-expressed miRNAs displaying a monotonic trend", 
       color = "Normalization method") +
  scale_color_manual(
    name   = "Normalization method",
    values = colour_map,
    labels = name_map,
    breaks = order_levels
  ) +
  theme_minimal(
    base_size = 16
  )

pdf(paste0(plot_dir, "percentage_monotonic_expres_relative.pdf"), width = 8, height = 6)  
print(gg)  
dev.off() 

log2FC_threshold <- 0.5
pure_log2FC <- group_extended %>%
  tibble::rownames_to_column("miRNA") %>%
  filter(average >= minRPM) %>%
  mutate(abs_log2FC_pure = abs(log2(`100S` / `0S`))) %>%
  arrange(desc(abs_log2FC_pure))

log2FC_threshold_pos <- sum(pure_log2FC$abs_log2FC_pure >= log2FC_threshold)


# lineplot % monotonic for top n (absolute log2FC)
gg <- ggplot(final_log2, aes(x = top_n, y = percentage_agreement, 
                             color = transformation) ) +
  geom_line() +
  geom_point() +
  # geom_vline(
  #   xintercept = log2FC_threshold_pos,
  #   linetype = "dashed",
  #   color = "grey70", 
  #   linewidth = 0.6
  # ) +
  # annotate(
  #   "text",
  #   x = log2FC_threshold_pos,
  #   y = Inf,
  #   label = paste0("abs(log[2]*FC) == 0.5"),
  #   parse = TRUE,
  #   angle = 15,
  #   vjust = 3,
  #   hjust = 0,
  #   size = 4
  # ) +
  labs(x = "Top n miRNAs (sorted by  |log2FC|)", y = "Percentage with monotonic trend (%)", 
       title = "Proportion of top DE miRNAs displaying a monotonic trend", 
       color = "normalization method") +
  scale_color_manual(
    name   = "normalization method",
    values = colour_map,
    labels = name_map,
    breaks = order_levels
  ) +
  theme_minimal(
    base_size = 16
  ) +
  coord_cartesian(clip = "off")

pdf(paste0(plot_dir, "percentage_monotonic_log2FC_relative.pdf"), width = 8, height = 6)
print(gg)  
dev.off()

# same plot only using top 100 most expressed miRNAS

if (TRUE){
  
  ## ---- Filter DE-sorted monotonicity plot to top 100 most expressed miRNAs ----
  ## Requires: list_of_transformations, groups, group_extended, name_map, colour_map, order_levels (optional), plot_dir
  
  # 1) Pick the top 100 most expressed miRNAs (robust to column naming / type)
  expr_col <- NULL
  if ("average_RPM" %in% colnames(group_extended)) expr_col <- "average_RPM"
  if (is.null(expr_col) && "average" %in% colnames(group_extended)) expr_col <- "average"
  if (is.null(expr_col)) stop("group_extended must contain a numeric column named 'average_RPM' or 'average'.")
  
  expr_vec <- suppressWarnings(as.numeric(group_extended[[expr_col]]))
  if (all(is.na(expr_vec))) stop(paste0("Column '", expr_col, "' could not be coerced to numeric."))
  
  top100_miR <- rownames(group_extended)[order(-expr_vec)][1:min(100, length(expr_vec))]
  
  # 2) Recompute the monotonicity-vs-topN curve using ONLY those miRNAs
  final_log2_top100 <- data.frame(
    transformation = character(),
    transformation_group = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    key_name <- names(list_of_transformations)[i]
    temp_df  <- list_of_transformations[[i]]
    
    # group geometric means for this transformation
    group_df <- geometric_mean_by_group(temp_df, groups)
    
    # restrict to the top 100 expressed miRNAs
    keep_miR <- intersect(rownames(group_df), top100_miR)
    if (length(keep_miR) == 0) next
    group_df <- group_df[keep_miR, , drop = FALSE]
    
    # compute monotonicity and |log2FC|
    results_df <- data.frame(row.names = rownames(group_df))
    results_df$monotonicity <- apply(group_df, 1, check_monotonicity)
    
    results_df$FC <- group_df$`100S` / group_df$`0S`
    results_df$log2FC <- log2(results_df$FC)
    results_df$abs_log2FC <- abs(results_df$log2FC)
    
    # monotonic in either direction counts as agreement
    results_df$agreement <- ifelse(results_df$monotonicity %in% c("increasing", "decreasing"), 1, 0)
    
    # sort by |log2FC| (top DE)
    sorted_df <- results_df[order(-results_df$abs_log2FC), , drop = FALSE]
    
    # top_n grid (cannot exceed 100 here)
    top_n_values <- c(5, 10, 20, 30, 40, 50, seq(60, nrow(sorted_df), by = 10))
    top_n_values <- unique(pmin(top_n_values, nrow(sorted_df)))
    
    for (top_n in top_n_values) {
      top_subset <- sorted_df[1:top_n, , drop = FALSE]
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      parts  <- unlist(strsplit(key_name, "_"))
      tgroup <- parts[1]
      
      final_log2_top100 <- rbind(final_log2_top100, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 3) Keep the same legend order (if available)
  if (exists("order_levels")) {
    final_log2_top100$transformation <- factor(final_log2_top100$transformation, levels = order_levels)
  }
  
  # 4) Plot (DE-sorted) restricted to top-100 expressed miRNAs
  gg_top100 <- ggplot(final_log2_top100,
                      aes(x = top_n, y = percentage_agreement, color = transformation)) +
    geom_line() +
    geom_point() +
    labs(
      x = "Top n miRNAs (sorted by  |log2FC|)",
      y = "Percentage with monotonic trend (%)",
      title = "Proportion of top DE miRNAs displaying a monotonic trend",
      color = "normalization method"
    ) +
    scale_color_manual(
      name   = "normalization method",
      values = colour_map,
      labels = name_map,
      breaks = if (exists("order_levels")) order_levels else NULL
    ) +
    theme_minimal(
      base_size = 16
    )
  
  pdf(file.path(plot_dir, "percentage_monotonic_log2FC_relative_top100ExprOnly.pdf"), width = 8, height = 6)
  print(gg_top100)
  dev.off()
  
  
  
}

# FC_direction up or down

# RPM_mono_df <- all_filtered[["RPM_lib"]]
RPM_mono_df <- all_unfiltered[["RPM_lib"]]
pvals <- ifelse(RPM_mono_df$monotonicity %in% c("increasing", "decreasing"), 0, 1)
FC_direction <- ifelse(
  RPM_mono_df$monotonicity == "increasing", "up",
  ifelse(RPM_mono_df$monotonicity == "decreasing", "down", "neither")
)

df_monotonicity <- data.frame(
  miRNA      = rownames(RPM_mono_df),
  average    = 0,
  log2FC     = 0,
  pval       = pvals,
  padj       = pvals,
  comparison = "0S_VS_100S",
  DE_method  = "monotonicity",
  parameters = "RPM",
  FC_direction = FC_direction,
  stringsAsFactors = FALSE
)

write.csv(df_monotonicity, file = monotonicity_path, 
          row.names = FALSE, quote = TRUE)


### Upset plot of overlapping monotonicity sets among methods

if (TRUE){
  
  monotonicity_sets <- list()
  
  for (key_name in names(all_filtered)) {
    df <- all_filtered[[key_name]]

    # Increasing miRNAs
    inc <- rownames(df[df$monotonicity == "increasing", , drop = FALSE])
    if (length(inc) > 0) {
      monotonicity_sets[[paste0(key_name, "_up")]] <- inc
    }
    
    # Decreasing miRNAs
    dec <- rownames(df[df$monotonicity == "decreasing", , drop = FALSE])
    if (length(dec) > 0) {
      monotonicity_sets[[paste0(key_name, "_down")]] <- dec
    }
  }
  
  # get current set names
  set_order_vec <- names(monotonicity_sets)
  
  set_order_vec <- gsub("^vst_local_blind_TRUE", "VST", set_order_vec)
  set_order_vec <- gsub("^rlog_mean(_blind_(TRUE|FALSE))?", "Rlog_mean", set_order_vec)
  set_order_vec <- gsub("^rlog_parametric(_blind_(TRUE|FALSE))?", "Rlog_parametric", set_order_vec)
  # convert RPM_lib to just RPM
  set_order_vec <- gsub("^RPM_lib", "RPM", set_order_vec)
  
  
  names(monotonicity_sets) <- set_order_vec
  
  set_order_vec <- c(
    grep("_up$", set_order_vec, value = TRUE),
    grep("_down$", set_order_vec, value = TRUE)
  )
  
  comb_mat <- make_comb_mat(monotonicity_sets)
  
  # order of intersections by size (largest -> smallest)
  comb_order_by_size <- order(comb_size(comb_mat), decreasing = TRUE)
  top_anno <- HeatmapAnnotation(
    "Intersection size" = anno_barplot(
      comb_size(comb_mat),
      border = FALSE,
      gp = gpar(fill = "black"),
      height = unit(2, "cm"),
      axis = TRUE,
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 10),
      numbers_offset = unit(2, "mm")
    )
  )
  ht <- UpSet(
    comb_mat,
    comb_order      = comb_order_by_size,
    set_order       = set_order_vec,
    # top_annotation  = upset_top_annotation(comb_mat),
    top_annotation  = top_anno,
    left_annotation = upset_left_annotation(comb_mat),
    row_names_side  = "right",
    row_names_gp    = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 12)
    
  )
  
  pdf(file.path(plot_dir, "upset_monotonicity.pdf"), width = 7, height = 5)
  draw(ht)
  dev.off()
    }


# aggregated sd pseudo-MA plot
if (TRUE){
  replicates1 <- c("1S", "20_5S_80_5L_1", "40_4S_60_4L_1",  "80_2S_20_2L_1") #"60_3S_40_3L_1",
  replicates2 <- c("1S_2", "20_5S_80_5L_2", "40_4S_60_4L_2", "80_2S_20_2L_2") # "60_3S_40_3L_2"
  
  cum_result <- data.frame(
    A = numeric(),
    sd_M = numeric(),
    transformation = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    if (key_name == "MAP"){
      temp_df <- temp_df * 1000000
    }
    temp_df <- temp_df + 1
    rep1_df <- log2(temp_df[, replicates1])
    rep2_df <- log2(temp_df[, replicates2])
    M_df <- rep1_df - rep2_df
    A_df <- 1/2 * (rep1_df + rep2_df)
    M_vector <- as.vector(as.matrix(M_df))
    A_vector <- as.vector(as.matrix(A_df))
    MA_df <- data.frame(A=A_vector, M=M_vector)
    plot_df <- data.frame(A=A_vector, M=M_vector)
    
    MA_df <- MA_df %>% arrange(A)
    bin_rel_size <- 0.1
    bin_width <- bin_rel_size * diff(range(MA_df$A))
    total_bins <- 1000
    step_size <- 1/total_bins * (max(MA_df$A) - min(MA_df$A))
    # Define the breaks for overlapping bins
    breaks <- seq(min(MA_df$A), max(MA_df$A) - bin_width, by = step_size)
    assign_bins <- function(x) {
      lower_edges <- breaks
      upper_edges <- breaks + bin_width
      which(x >= lower_edges & x < upper_edges)
    }
    MA_df$Bins <- lapply(MA_df$A, assign_bins)
    # Explode the dataframe by bins
    binned_data <- tidyr::unnest(MA_df, cols = Bins)
    # Calculate mean A and mean M for each bin
    result <- binned_data %>%
      group_by(Bins) %>%
      summarise(
        A_mean = mean(A),
        M_mean = stats::sd(M),
        .groups = 'drop'
      )
    # Assign the mean A values as bin centers
    ####
    result$A_mean <- breaks[result$Bins] + bin_width / 2 
    ####
    colnames(result) <- c("bin", "A", "sd_M")
    result$transformation <- key_name
    cum_result <- rbind(cum_result, result)
    
  }
  
  cum_result$transformation <- factor(cum_result$transformation, levels = names(list_of_transformations))
  line_types <- c("solid", "dashed")
  
  line_type_mapping <- rep(line_types, length.out = length(names(list_of_transformations)))
  names(line_type_mapping) <- names(list_of_transformations)
 
  order_levels <- cum_result %>%
    group_by(transformation) %>%
    slice_max(order_by = A, n = 1, with_ties = FALSE) %>%
    arrange(desc(sd_M)) %>%
    pull(transformation)
  
  gg <- ggplot(cum_result, aes(x = A , y = sd_M, color = transformation, linetype = transformation)) +
    # gg <- ggplot(plot_df, aes(x = A, y = M)) +
    geom_line() +
    scale_linetype_manual(values = line_type_mapping,
                          labels = name_map,
                          breaks = order_levels
    ) +
    scale_colour_manual(
      name   = "normalization method",
      values = colour_map,
      breaks = order_levels,
      labels = name_map
    ) +
    scale_linetype_manual(
      name   = "normalization method",
      values = line_type_mapping,
      breaks = order_levels,
      labels = name_map
    ) +
    theme_minimal() + 
    labs(
      title = paste0(""),
      x = "Average expression",
      y = "sd(log2FC)"
    )
  
  pdf(paste0(plot_dir, "sd_mean_plot_raw.pdf"), width = 8, height = 6)
  print(gg)  
  dev.off()
   
}

# pseudo-MA plot for each transformation
if (TRUE){
  
  replicates1 <- c("1S", "20_5S_80_5L_1", "40_4S_60_4L_1", "80_2S_20_2L_1")
  replicates2 <- c("1S_2", "20_5S_80_5L_2", "40_4S_60_4L_2", "80_2S_20_2L_2")
  
  for (method_name in names(list_of_transformations)) {
    
    message("Making aggregated MA plot for: ", method_name)
    
    temp_df <- list_of_transformations[[method_name]]
    
    if (method_name == "MAP") {
      temp_df <- temp_df * 1e6
    }
    
    keep_pairs <- replicates1 %in% colnames(temp_df) &
      replicates2 %in% colnames(temp_df)
    
    common_rep1 <- replicates1[keep_pairs]
    common_rep2 <- replicates2[keep_pairs]
    
    if (length(common_rep1) == 0) {
      warning("No valid replicate pairs found for ", method_name)
      next
    }
    
    temp_df <- temp_df + 1
    
    plot_list <- list()
    
    for (i in seq_along(common_rep1)) {
      
      s1 <- common_rep1[i]
      s2 <- common_rep2[i]
      
      rep1 <- log2(temp_df[, s1])
      rep2 <- log2(temp_df[, s2])
      
      plot_list[[i]] <- data.frame(
        miRNA = rownames(temp_df),
        A = 0.5 * (rep1 + rep2),
        M = rep1 - rep2,
        comparison = paste0(s1, " vs ", s2),
        stringsAsFactors = FALSE
      )
    }
    
    plot_df <- dplyr::bind_rows(plot_list)
    
    plot_df <- plot_df[
      is.finite(plot_df$A) & is.finite(plot_df$M),
      ,
      drop = FALSE
    ]
    
    gg <- ggplot(plot_df, aes(x = A, y = M)) +
      geom_point(alpha = 0.18, size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_minimal(base_size = 14) +
      labs(
        title = method_name,
        x = "A = average log2 expression",
        y = "M = log2 fold-change"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    output_file <- file.path(
      plot_dir,
      paste0("MA_plot_aggregated_", gsub("[^A-Za-z0-9_]+", "_", method_name), ".pdf")
    )
    
    pdf(output_file, width = 7, height = 6)
    print(gg)
    dev.off()
  }
}

# one figure pseudo-MA plots
if (TRUE){
  library(patchwork)
  
  ## Combined aggregated MA plots: 2 columns x 4 rows, A4 proportions
  ## Uses: list_of_transformations, plot_dir
  
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  replicates1 <- c("1S", "20_5S_80_5L_1", "40_4S_60_4L_1", "80_2S_20_2L_1")
  replicates2 <- c("1S_2", "20_5S_80_5L_2", "40_4S_60_4L_2", "80_2S_20_2L_2")
  
  plot_list_all <- list()
  
  for (method_name in names(list_of_transformations)) {
    
    temp_df <- list_of_transformations[[method_name]]
    
    if (method_name == "MAP") {
      temp_df <- temp_df * 1e6
    }
    
    keep_pairs <- replicates1 %in% colnames(temp_df) &
      replicates2 %in% colnames(temp_df)
    
    common_rep1 <- replicates1[keep_pairs]
    common_rep2 <- replicates2[keep_pairs]
    
    if (length(common_rep1) == 0) next
    
    temp_df <- temp_df + 1
    
    plot_list <- list()
    
    for (i in seq_along(common_rep1)) {
      
      s1 <- common_rep1[i]
      s2 <- common_rep2[i]
      
      rep1 <- log2(temp_df[, s1])
      rep2 <- log2(temp_df[, s2])
      
      plot_list[[i]] <- data.frame(
        miRNA = rownames(temp_df),
        A = 0.5 * (rep1 + rep2),
        M = rep1 - rep2,
        comparison = paste0(s1, " vs ", s2),
        stringsAsFactors = FALSE
      )
    }
    
    plot_df <- dplyr::bind_rows(plot_list)
    
    plot_df <- plot_df[
      is.finite(plot_df$A) & is.finite(plot_df$M),
      ,
      drop = FALSE
    ]
    
    plot_list_all[[method_name]] <- ggplot(plot_df, aes(x = A, y = M)) +
      geom_point(alpha = 0.18, size = 0.35) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
      theme_minimal(base_size = 8) +
      labs(
        title = ifelse(method_name %in% names(name_map), name_map[[method_name]], method_name),
        x = "Average log2 expression",
        y = "log2FC"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)
      )
  }
  
  combined_plot <- wrap_plots(plot_list_all, ncol = 2, nrow = 4)
  
  pdf(
    file.path(plot_dir, "MA_plots_aggregated_A4_2col_4row.pdf"),
    width = 8.27,
    height = 11.69
  )
  
  print(combined_plot)
  
  dev.off()
}

#### pseudo-MA plot (common x-axis)
if (TRUE){
  replicates1 <- c("1S", "20_5S_80_5L_1", "40_4S_60_4L_1",  "80_2S_20_2L_1")
  replicates2 <- c("1S_2", "20_5S_80_5L_2", "40_4S_60_4L_2", "80_2S_20_2L_2")
  
  cum_result <- data.frame(
    A = numeric(),
    sd_M = numeric(),
    transformation = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    if (key_name == "MAP") {
      temp_df <- temp_df * 1000000
    }
    temp_df <- temp_df + 1
    
    rep1_df <- log2(temp_df[, replicates1])
    rep2_df <- log2(temp_df[, replicates2])
    M_df <- rep1_df - rep2_df
    A_df <- 0.5 * (rep1_df + rep2_df)
    
    M_vector <- as.vector(as.matrix(M_df))
    A_vector <- as.vector(as.matrix(A_df))
    
    MA_df <- data.frame(A = A_vector, M = M_vector)
    
    # ---- same as before: compute sd(M) in overlapping bins ----
    bin_rel_size <- 0.1
    bin_width <- bin_rel_size * diff(range(MA_df$A))
    total_bins <- 1000
    step_size <- 1 / total_bins * (max(MA_df$A) - min(MA_df$A))
    breaks <- seq(min(MA_df$A), max(MA_df$A) - bin_width, by = step_size)
    
    assign_bins <- function(x) {
      lower_edges <- breaks
      upper_edges <- breaks + bin_width
      which(x >= lower_edges & x < upper_edges)
    }
    
    MA_df$Bins <- lapply(MA_df$A, assign_bins)
    binned_data <- tidyr::unnest(MA_df, cols = Bins)
    
    result <- binned_data %>%
      group_by(Bins) %>%
      summarise(
        A = mean(A),
        sd_M = stats::sd(M),
        .groups = 'drop'
      )
    
    # ---- here's the key fix ----
    # Instead of sorting A before calculating sd, sort only for plotting:
    result <- result[order(result$A), ]
    result$A_order <- seq_len(nrow(result))  # fixed axis 1…N
    
    result$transformation <- key_name
    cum_result <- rbind(cum_result, result)
  }
  
  cum_result$transformation <- factor(cum_result$transformation, levels = names(list_of_transformations))
  line_types <- c("solid", "dashed")
  line_type_mapping <- rep(line_types, length.out = length(names(list_of_transformations)))
  names(line_type_mapping) <- names(list_of_transformations)
  
  order_levels <- cum_result %>%
    group_by(transformation) %>%
    slice_max(order_by = A, n = 1, with_ties = FALSE) %>%
    arrange(desc(sd_M)) %>%
    pull(transformation)
  
  gg <- ggplot(cum_result, aes(x = A_order, y = sd_M, color = transformation, linetype = transformation)) +
    geom_line() +
    scale_linetype_manual(name = "normalization method", values = line_type_mapping, labels = name_map, breaks = order_levels) +
    scale_colour_manual(name = "normalization method", values = colour_map, breaks = order_levels, labels = name_map) +
    theme_minimal(
      base_size = 16
    ) +
    labs(
      title = "",
      x = "Average expression (scaled)",
      y = "sd(log2FC)"
    )
  
  pdf(paste0(plot_dir, "sd_mean_plot_scaled.pdf"), width = 8, height = 6)
  print(gg)
  dev.off()
}

#### pseudo-MA plot using biological replicates within each group
if (TRUE) {
  
  group_samples <- list(
    "0S_100L" = c("1L", "2L", "3L", "4L", "5L"),
    "20S_80L" = c("20_1S_80_1L", "20_2S_80_2L", "20_3S_80_3L", "20_4S_80_4L", "20_5S_80_5L"),
    "40S_60L" = c("40_1S_60_1L", "40_2S_60_2L", "40_3S_60_3L", "40_4S_60_4L", "40_5S_60_5L"),
    "60S_40L" = c("60_1S_40_1L", "60_2S_40_2L", "60_3S_40_3L", "60_4S_40_4L", "60_5S_40_5L"),
    "80S_20L" = c("80_1S_20_1L", "80_2S_20_2L", "80_3S_20_3L", "80_4S_20_4L", "80_5S_20_5L"),
    "100S_0L" = c("1S", "2S", "3S", "4S", "5S")
  )
  
  cum_result <- data.frame(
    A = numeric(),
    sd_M = numeric(),
    transformation = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    
    if (key_name == "MAP") {
      temp_df <- temp_df * 1000000
    }
    temp_df <- temp_df + 1
    
    all_MA <- list()
    
    
    # For each group, compare all biological replicate pairs
    for (group_name in names(group_samples)) {
      
      samples_in_group <- group_samples[[group_name]]
      samples_in_group <- samples_in_group[samples_in_group %in% colnames(temp_df)]
      
      if (length(samples_in_group) < 2) next
      
      pair_mat <- combn(samples_in_group, 2)
      
      for (j in seq_len(ncol(pair_mat))) {
        s1 <- pair_mat[1, j]
        s2 <- pair_mat[2, j]
        
        rep1 <- log2(temp_df[, s1])
        rep2 <- log2(temp_df[, s2])
        
        M_vec <- rep1 - rep2
        A_vec <- 0.5 * (rep1 + rep2)
        
        all_MA[[length(all_MA) + 1]] <- data.frame(
          A = A_vec,
          M = M_vec
        )
      }
    }
    
    MA_df <- do.call(rbind, all_MA)
    
    ##  same binning approach as before
    bin_rel_size <- 0.1
    bin_width <- bin_rel_size * diff(range(MA_df$A))
    total_bins <- 1000
    step_size <- 1 / total_bins * (max(MA_df$A) - min(MA_df$A))
    breaks <- seq(min(MA_df$A), max(MA_df$A) - bin_width, by = step_size)
    
    assign_bins <- function(x) {
      lower_edges <- breaks
      upper_edges <- breaks + bin_width
      which(x >= lower_edges & x < upper_edges)
    }
    
    MA_df$Bins <- lapply(MA_df$A, assign_bins)
    binned_data <- tidyr::unnest(MA_df, cols = Bins)
    
    result <- binned_data %>%
      dplyr::group_by(Bins) %>%
      dplyr::summarise(
        A = mean(A),
        sd_M = stats::sd(M),
        .groups = "drop"
      )
    
    result <- result[order(result$A), ]
    result$A_order <- seq_len(nrow(result))
    
    result$transformation <- key_name
    cum_result <- rbind(cum_result, result)
  }
  
  cum_result$transformation <- factor(
    cum_result$transformation,
    levels = names(list_of_transformations)
  )
  
  line_types <- c("solid", "dashed")
  line_type_mapping <- rep(line_types, length.out = length(names(list_of_transformations)))
  names(line_type_mapping) <- names(list_of_transformations)
  
  order_levels <- cum_result %>%
    dplyr::group_by(transformation) %>%
    dplyr::slice_max(order_by = A, n = 1, with_ties = FALSE) %>%
    dplyr::arrange(desc(sd_M)) %>%
    dplyr::pull(transformation)
  
  gg <- ggplot(cum_result, aes(x = A_order, y = sd_M, color = transformation, linetype = transformation)) +
    geom_line() +
    scale_linetype_manual(
      name = "normalization method",
      values = line_type_mapping,
      labels = name_map,
      breaks = order_levels
    ) +
    scale_colour_manual(
      name = "normalization method",
      values = colour_map,
      breaks = order_levels,
      labels = name_map
    ) +
    theme_minimal(base_size = 16) +
    labs(
      title = "",
      x = "Average expression (scaled)",
      y = "sd(log2FC)"
    )
  
  pdf(paste0(plot_dir, "sd_mean_plot_scaled_biological_replicates.pdf"), width = 8, height = 6)
  print(gg)
  dev.off()
}

#### pseudo-MA plot using 100 random comparisons (same pairs for all methods)
if (TRUE) {
  
  set.seed(123)
  
  
  # Choose the common sample set across all normalization methods
  common_samples <- Reduce(intersect, lapply(list_of_transformations, colnames))
  
  n_random_pairs <- 100
  
  all_possible_pairs <- t(combn(common_samples, 2))
  
  if (nrow(all_possible_pairs) < n_random_pairs) {
    stop("There are fewer than 100 unique sample pairs available.")
  }
  
  selected_idx <- sample(seq_len(nrow(all_possible_pairs)), n_random_pairs, replace = FALSE)
  random_pairs <- all_possible_pairs[selected_idx, , drop = FALSE]
  
  random_pairs_df <- data.frame(
    sample1 = random_pairs[, 1],
    sample2 = random_pairs[, 2],
    stringsAsFactors = FALSE
  )
  
  cum_result <- data.frame(
    A = numeric(),
    sd_M = numeric(),
    transformation = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    
    if (key_name == "MAP") {
      temp_df <- temp_df * 1000000
    }
    
    temp_df <- temp_df + 1
    
    all_MA <- vector("list", nrow(random_pairs_df))
    
    
    # Use the same random sample pairs for each method
    
    for (j in seq_len(nrow(random_pairs_df))) {
      s1 <- random_pairs_df$sample1[j]
      s2 <- random_pairs_df$sample2[j]
      
      rep1 <- log2(temp_df[, s1])
      rep2 <- log2(temp_df[, s2])
      
      M_vec <- rep1 - rep2
      A_vec <- 0.5 * (rep1 + rep2)
      
      all_MA[[j]] <- data.frame(
        A = A_vec,
        M = M_vec
      )
    }
    
    MA_df <- do.call(rbind, all_MA)
    MA_df <- MA_df[is.finite(MA_df$A) & is.finite(MA_df$M), , drop = FALSE]
    
    
    # Same overlapping-bin approach
    
    bin_rel_size <- 0.1
    bin_width <- bin_rel_size * diff(range(MA_df$A))
    total_bins <- 1000
    step_size <- (max(MA_df$A) - min(MA_df$A)) / total_bins
    breaks <- seq(min(MA_df$A), max(MA_df$A) - bin_width, by = step_size)
    
    assign_bins <- function(x) {
      lower_edges <- breaks
      upper_edges <- breaks + bin_width
      which(x >= lower_edges & x < upper_edges)
    }
    
    MA_df$Bins <- lapply(MA_df$A, assign_bins)
    binned_data <- tidyr::unnest(MA_df, cols = Bins)
    
    result <- binned_data %>%
      dplyr::group_by(Bins) %>%
      dplyr::summarise(
        A = mean(A),
        sd_M = stats::sd(M),
        .groups = "drop"
      )
    
    result <- result[order(result$A), ]
    result$A_order <- seq_len(nrow(result))
    result$transformation <- key_name
    
    cum_result <- rbind(cum_result, result)
  }
  
  cum_result$transformation <- factor(
    cum_result$transformation,
    levels = names(list_of_transformations)
  )
  
  line_types <- c("solid", "dashed")
  line_type_mapping <- rep(line_types, length.out = length(names(list_of_transformations)))
  names(line_type_mapping) <- names(list_of_transformations)
  
  order_levels <- cum_result %>%
    dplyr::group_by(transformation) %>%
    dplyr::slice_max(order_by = A, n = 1, with_ties = FALSE) %>%
    dplyr::arrange(desc(sd_M)) %>%
    dplyr::pull(transformation)
  
  gg <- ggplot(cum_result, aes(x = A_order, y = sd_M, color = transformation, linetype = transformation)) +
    geom_line() +
    scale_linetype_manual(
      name = "normalization method",
      values = line_type_mapping,
      labels = name_map,
      breaks = order_levels
    ) +
    scale_colour_manual(
      name = "normalization method",
      values = colour_map,
      breaks = order_levels,
      labels = name_map
    ) +
    theme_minimal(base_size = 16) +
    labs(
      title = "",
      x = "Average expression (scaled)",
      y = "sd(log2FC)"
    )
  
  pdf(paste0(plot_dir, "sd_mean_plot_scaled_random100.pdf"), width = 8, height = 6)
  print(gg)
  dev.off()
  
}

#### dot plot: rlog-mean vs Raw read count (robust version)
if (FALSE) {
  
  ## ------------------------------------------------------------
  ## Select the two matrices (adjust names if needed)
  ## ------------------------------------------------------------
  raw_df  <- list_of_transformations[["RC"]]
  rlog_df <- list_of_transformations[["rlog_mean_blind_FALSE"]]
  
  ## ------------------------------------------------------------
  ## Keep only common miRNAs and samples
  ## ------------------------------------------------------------
  common_rows <- intersect(rownames(raw_df), rownames(rlog_df))
  common_cols <- intersect(colnames(raw_df), colnames(rlog_df))
  
  raw_df  <- raw_df[common_rows, common_cols, drop = FALSE]
  rlog_df <- rlog_df[common_rows, common_cols, drop = FALSE]
  
  ## ensure matrix format (extra safety)
  raw_df  <- as.matrix(raw_df)
  rlog_df <- as.matrix(rlog_df)
  
  ## ------------------------------------------------------------
  ## Convert to long format (safe, no as.table)
  ## ------------------------------------------------------------
  plot_df <- data.frame(
    miRNA = rep(rownames(raw_df), times = ncol(raw_df)),
    sample = rep(colnames(raw_df), each = nrow(raw_df)),
    raw_count = as.vector(raw_df),
    rlog_mean = as.vector(rlog_df),
    stringsAsFactors = FALSE
  )
  
  ## remove invalid values
  plot_df <- plot_df[
    is.finite(plot_df$raw_count) & is.finite(plot_df$rlog_mean),
    ,
    drop = FALSE
  ]
  
  ## ------------------------------------------------------------
  ## Plot (recommended: log-scale on raw counts)
  ## ------------------------------------------------------------
  gg <- ggplot(plot_df, aes(x = log10(raw_count + 1), y = rlog_mean)) +
    geom_point(alpha = 0.15, size = 0.6) +
    theme_minimal(base_size = 16) +
     scale_y_log10() +
    labs(
      title = "",
      x = "log10(raw read count + 1)",
      y = "rlog (mean fit)"
    )
  
  pdf(paste0(plot_dir, "dotplot_rlog_mean_vs_raw_counts.pdf"), width = 7, height = 6)
  print(gg)
  dev.off()
}


#### dot plot: rlog_mean_blind_FALSE vs rlog_parametric_blind_TRUE
if (FALSE) {
  
  ## ------------------------------------------------------------
  ## Set transformation names
  ## ------------------------------------------------------------
  x_name <- "rlog_mean_blind_FALSE"
  y_name <- "rlog_parametric_blind_TRUE"
  
  ## ------------------------------------------------------------
  ## Check availability
  ## ------------------------------------------------------------
  if (!x_name %in% names(list_of_transformations)) {
    stop(paste0(
      "Transformation '", x_name, "' not found.\nAvailable:\n",
      paste(names(list_of_transformations), collapse = ", ")
    ))
  }
  
  if (!y_name %in% names(list_of_transformations)) {
    stop(paste0(
      "Transformation '", y_name, "' not found.\nAvailable:\n",
      paste(names(list_of_transformations), collapse = ", ")
    ))
  }
  
  x_df <- list_of_transformations[[x_name]]
  y_df <- list_of_transformations[[y_name]]
  
  if (is.null(x_df) || is.null(y_df)) {
    stop("One of the selected transformations is NULL.")
  }
  
  ## ------------------------------------------------------------
  ## Ensure matrix format
  ## ------------------------------------------------------------
  x_df <- as.matrix(x_df)
  y_df <- as.matrix(y_df)
  
  ## ------------------------------------------------------------
  ## Match miRNAs and samples
  ## ------------------------------------------------------------
  common_rows <- intersect(rownames(x_df), rownames(y_df))
  common_cols <- intersect(colnames(x_df), colnames(y_df))
  
  if (length(common_rows) == 0) stop("No common miRNAs.")
  if (length(common_cols) == 0) stop("No common samples.")
  
  x_df <- x_df[common_rows, common_cols, drop = FALSE]
  y_df <- y_df[common_rows, common_cols, drop = FALSE]
  
  ##
  
  selected_sample <- "1S"   # or specify explicitly, e.g. "1S"
  selected_sample <- "1L"   # or specify explicitly, e.g. "1S"
  
  x_df <- x_df[, selected_sample, drop = FALSE]
  y_df <- y_df[, selected_sample, drop = FALSE]
  
  
  ## ------------------------------------------------------------
  ## Long format (safe)
  ## ------------------------------------------------------------
  plot_df <- data.frame(
    miRNA = rep(rownames(x_df), times = ncol(x_df)),
    sample = rep(colnames(x_df), each = nrow(x_df)),
    x_val = as.vector(x_df),
    y_val = as.vector(y_df),
    stringsAsFactors = FALSE
  )
  
  plot_df <- plot_df[
    is.finite(plot_df$x_val) & is.finite(plot_df$y_val),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) == 0) {
    stop("No valid points to plot.")
  }
  
  ## ------------------------------------------------------------
  ## Plot
  ## ------------------------------------------------------------
  gg <- ggplot(plot_df, aes(x = x_val, y = y_val)) +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 16) +
     scale_x_log10() + scale_y_log10() +
    labs(
      title = "",
      x = x_name,
      y = y_name
    )
  
  pdf(paste0(plot_dir, "dotplot_", x_name, "_vs_", y_name, ".pdf"), width = 7, height = 6)
  print(gg)
  dev.off()
}


# CVA plot for each transformation (dot per feature)
if (FALSE){
  print("hey")
  replicates1 <- c("1S", "20_5S_80_5L_1", "40_4S_60_4L_1",  "80_2S_20_2L_1") #"60_3S_40_3L_1",
  replicates2 <- c("1S_2", "20_5S_80_5L_2", "40_4S_60_4L_2", "80_2S_20_2L_2") # "60_3S_40_3L_2"
  
  cum_result <- data.frame(
    A = numeric(),
    cv_M = numeric(),     # << changed
    transformation = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    temp_df <- temp_df + 1
    rep1_df <- log2(temp_df[, replicates1])
    rep2_df <- log2(temp_df[, replicates2])
    M_df <- rep1_df - rep2_df
    A_df <- 1/2 * (rep1_df + rep2_df)
    M_vector <- as.vector(as.matrix(M_df))
    A_vector <- as.vector(as.matrix(A_df))
    MA_df <- data.frame(A=A_vector, M=M_vector)
    plot_df <- data.frame(A=A_vector, M=M_vector)
    
    MA_df <- MA_df %>% arrange(A)
    bin_rel_size <- 0.1
    bin_width <- bin_rel_size * diff(range(MA_df$A))
    total_bins <- 1000
    step_size <- 1/total_bins * (max(MA_df$A) - min(MA_df$A))
    breaks <- seq(min(MA_df$A), max(MA_df$A) - bin_width, by = step_size)
    
    assign_bins <- function(x) {
      lower_edges <- breaks
      upper_edges <- breaks + bin_width
      which(x >= lower_edges & x < upper_edges)
    }
    MA_df$Bins <- lapply(MA_df$A, assign_bins)
    binned_data <- tidyr::unnest(MA_df, cols = Bins)
    
    # Calculate mean A and coefficient of variation of M for each bin
    result <- binned_data %>%
      group_by(Bins) %>%
      summarise(
        A_mean = mean(A),
        M_mean = mean(M),
        M_sd   = stats::sd(M),
        cv_M   = ifelse(M_mean == 0, NA, M_sd / abs(M_mean)),  # << changed
        .groups = 'drop'
      )
    
    result$A_mean <- breaks[result$Bins] + bin_width / 2
    result <- result[, c("Bins", "A_mean", "cv_M")]            # << changed
    colnames(result) <- c("bin", "A", "cv_M")                  # << changed
    result$transformation <- key_name
    cum_result <- rbind(cum_result, result)
  }
  
  cum_result$transformation <- factor(cum_result$transformation, levels = names(list_of_transformations))
  line_types <- c("solid", "dashed")
  line_type_mapping <- rep(line_types, length.out = length(names(list_of_transformations)))
  names(line_type_mapping) <- names(list_of_transformations)
  
  order_levels <- cum_result %>%
    group_by(transformation) %>%
    slice_max(order_by = A, n = 1, with_ties = FALSE) %>%
    arrange(desc(cv_M)) %>%             # << changed
    pull(transformation)
  
  gg <- ggplot(cum_result, aes(x = A , y = cv_M, color = transformation, linetype = transformation)) +  # << changed
    geom_line() +
    scale_linetype_manual(values = line_type_mapping,
                          labels = name_map,
                          breaks = order_levels
    ) +
    scale_colour_manual(
      name   = "normalization method",
      values = colour_map,
      breaks = order_levels,
      labels = name_map
    ) +
    scale_linetype_manual(
      name   = "normalization method",
      values = line_type_mapping,
      breaks = order_levels,
      labels = name_map
    ) +
    theme_minimal() + 
    labs(
      title = paste0(""),
      x = "Average expression",
      y = "CV(log2FC)"          # << changed
    )
  
  pdf(paste0(plot_dir, "cv_mean_plot.pdf"), width = 8, height = 6)  # << changed
  print(gg)  
  dev.off()
}



### O/E (observed/expected) values  
if (TRUE){
  # for loop that iterates all transformations
  # calculate average value in pure samples
  # remove anything below RPM threshold (get list of miRNAs over threshold to subsect)
  RPM_threshold <- minRPM
  over50_miRNAs <- rownames(old_group_extended[old_group_extended$average > RPM_threshold, ])
  
  group_lookup <- list()
  for (key in names(groups)) {
    for (value in groups[[key]]) {
      group_lookup[[value]] <- key
    }
  }
  
  OE_ratio_list <- list()
  meme <- 0
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    
    if (key_name == "MAP") {
      temp_df <- list_of_transformations[[i]] * 1000000
    }
    print(key_name)
    # temp_df <- list_of_transformations[[i]]
    temp_df <- temp_df[, colnames(temp_df) %in% names(group_lookup), drop = FALSE]
    group_df <- geometric_mean_by_group(temp_df, groups)
    temp_df <- temp_df[rownames(temp_df) %in% over50_miRNAs, ]
    group_df <- group_df[rownames(group_df) %in% over50_miRNAs, ]
    if (key_name == "MAP") {
      cat("rownames intersect over50_miRNAs:", length(intersect(rownames(temp_df), over50_miRNAs)), "\n")
    }
    
    recalculated_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df)))
    colnames(recalculated_df) <- colnames(temp_df)
    columns_to_drop <- c()
    for (j in seq_along(colnames(temp_df))) {
      # Get the group name and extract the percentages
      col_name <- colnames(temp_df)[j]
      # get group name from group_assignment
      if (col_name %in% names(group_lookup)){
        group_name <- group_lookup[[col_name]]
        # print(col_name)
        # print(group_name)
        percentage <- as.numeric(gsub("S", "", group_name)) / 100
        # print(percentage)
        remaining_percentage <- 1 - percentage
        # Calculate the weighted value for each row using the columns 100S and 0S in summary_df
        recalculated_values <- (percentage * group_df[["100S"]]) + (remaining_percentage * group_df[["0S"]])
        # Assign the recalculated values to the corresponding column in recalculated_df
        recalculated_df[[col_name]] <- recalculated_values
      }else{
        # store columns not in list
        columns_to_drop <- c(columns_to_drop, col_name)
      }
      
    }
    # remove columns_to_drop from both dfs
    temp_df <- temp_df[, !(colnames(temp_df) %in% columns_to_drop), drop = FALSE]
    recalculated_df <- recalculated_df[, !(colnames(recalculated_df) %in% columns_to_drop), drop = FALSE]
    
    if (key_name == "MAP") {
      cat("temp_df dims:", dim(temp_df), "recalc_df dims:", dim(recalculated_df), "\n")
    }
    OE_ratio <- temp_df/recalculated_df
    OE_ratio_list[[key_name]] <- OE_ratio
    
  }
  
  plot_df <- data.frame()
  
  for (i in seq_along(list_of_transformations)){
    temp_OE_ratio <- OE_ratio_list[[i]]
    key_name <- names(list_of_transformations)[i]
    print(key_name)
    print(head(temp_OE_ratio))
    temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
    long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", 
                            variable.name = "Sample", value.name = "ObservedExpectedRatio")
    long_ratio_data$transformation <- key_name
    print(head(long_ratio_data))
    print(head(plot_df))
    plot_df <- rbind(plot_df, long_ratio_data)
  }

  gg <- ggplot(plot_df, aes(x = transformation, y = ObservedExpectedRatio)) +
    geom_boxplot() +  # Boxplot without displaying outliers
    # geom_jitter() +  
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
    labs(x = "transformation", y = "Observed/Expected Ratio") +
    ggtitle("Observed/Expected Ratio") +
    # scale_y_continuous(breaks = seq(0, max_y, by = 1)) +  # Add y-axis ticks at every unit
    # scale_y_continuous() +  # Add y-axis ticks at every unit
    scale_y_log10()+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(plot_dir, "log_OE.pdf"), plot = gg, width = 8, height = 6, device = "pdf")
    
}


### O/E values calculated per mouse
# O/E values predicted per sample using source mouse data


if (TRUE){
  OE_ratio_list <- list()
  relevant_columns <- colnames(recalculated_df)
  group_prefixes <- c("20", "40", "60", "80")
  percentages <- c(20, 40, 60, 80) / 100  # Convert to proportions
  remaining_percentages <- 1 - percentages  # Complementary percentages
  
  
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    print(key_name)
    if (key_name == "MAP") {
      temp_df <- list_of_transformations[[i]] * 1000000
    }
    temp_df <- temp_df[rownames(temp_df) %in% over50_miRNAs, ]
    temp_df <- temp_df[, relevant_columns]
    mouse_recalculated_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = 0))
    rownames(mouse_recalculated_df) <- rownames(temp_df)  # Set row names to match filtered_df
    # recalculated_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df)))
    # colnames(recalculated_df) <- colnames(temp_df)
    print("step 1")
    prefix_columns <- c()
    for (j in seq_along(group_prefixes)) {
      prefix <- group_prefixes[j]
      percentage <- percentages[j]
      remaining_percentage <- remaining_percentages[j]
      print("step 2")
      # Select columns that start with the current prefix
      selected_columns <- grep(paste0("^", prefix), colnames(temp_df), value = TRUE)
      prefix_columns <- c(prefix_columns, selected_columns)
      # Loop through the selected columns and extract the corresponding mouse identifiers
      for (col in selected_columns) {
        # Extract the mouse identifier, which is located between the prefix and the "S" or "L" part
        mouse_id <- sub(paste0("^", prefix, "_(\\d+)[SL]_.*"), "\\1", col)  # Extracts the mouse number
        
        # Define the corresponding "S" and "L" columns using the extracted mouse ID
        S_column <- paste0(mouse_id, "S")
        L_column <- paste0(mouse_id, "L")
        print("step 3")
        # Check if both corresponding columns (e.g., 1S and 1L) exist in filtered_df
        print(S_column)
        print(L_column)
        print(colnames(temp_df))
        if (S_column %in% colnames(temp_df) && L_column %in% colnames(temp_df)) {
          # Calculate the weighted value for each row
          print("step 4")
          recalculated_values <- (percentage * temp_df[[S_column]]) + (remaining_percentage * temp_df[[L_column]])
          print("step 5")
          # Add the recalculated values as a new column in mouse_recalculated_df
          mouse_recalculated_df[[col]] <- recalculated_values
        }
      }
    }
    temp_df <- temp_df[, prefix_columns]
    OE_ratio <- temp_df/mouse_recalculated_df
    OE_ratio_list[[key_name]] <- OE_ratio
  }
  
  # generate plot_df
  
  plot_df <- data.frame()
  
  for (i in seq_along(list_of_transformations)){
    temp_OE_ratio <- OE_ratio_list[[i]]
    key_name <- names(list_of_transformations)[i]
    temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
    long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", 
                            variable.name = "Sample", value.name = "ObservedExpectedRatio")
    long_ratio_data$transformation <- key_name
    # print(head(long_ratio_data))
    plot_df <- rbind(plot_df, long_ratio_data)
  }
  plot_df$transformation <- factor(plot_df$transformation,
                                   levels = names(list_of_transformations))
  gg <- ggplot(plot_df,
               aes(x = transformation,
                   y = ObservedExpectedRatio,
                   fill = transformation)) +
    geom_boxplot(alpha = 0.8) +                # show.legend defaults to TRUE
    guides(alpha = "none") +
    geom_hline(yintercept = 1,
               linetype = "dashed",
               color = "blue") +
    scale_y_log10() +
    scale_x_discrete(
      labels = name_map
    ) +
    
    scale_fill_manual(
      name   = "normalization method",
      values = colour_map,
      labels = name_map
    ) +
    
    theme_minimal(
      base_size = 16
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = NULL,
      y = "Observed/Expected Ratio"
    )
  
  ggsave(paste0(plot_dir, "log_mouse_OE.pdf"), plot = gg, width = 8, height = 6, device = "pdf")
  
}


#### supplementary figures

### heatmaps
if (TRUE){
  rpm_key <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  rpm_key <- "RPM_lib"
  # rpm_key <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  
  # Compute geometric mean by group for RPM
  rpm_grp       <- geometric_mean_by_group(list_of_transformations[[rpm_key]], groups)
  average_RPM   <- rowMeans(rpm_grp)
  
  # Get the top 100 miRNAs by that average
  topN <- 100
  top100 <- names(sort(average_RPM, decreasing = TRUE))[1:topN]
  
  # ADDING FC COLUMN TO HEATMAP
  if (FALSE){
    
    # generate RAW_DF here
    RAW_DF <- rpm_grp  
    print(colnames(RAW_DF))
    
    
    # Compute log2FC from two raw columns (placeholders RAW_DF, COL_A, COL_B)
    # Ensure miRNA names align to top100 order
    fc_vec <- log2(
      RAW_DF[top100, "100S", drop = TRUE] /
        RAW_DF[top100, "0S" , drop = TRUE]
    )
    print(fc_vec)
    
    # Build a data frame for the extra column
    fc_df <- data.frame(
      miRNA          = top100,
      transformation = "log2FC",
      log2FC         = fc_vec,
      stringsAsFactors = FALSE
    )
    
    
    
    # 4. Build a long table of monotonicity for those 100
    mono_list <- lapply(names(list_of_transformations), function(key) {
      grp  <- geometric_mean_by_group(list_of_transformations[[key]], groups)
      grp100 <- grp[top100, , drop = FALSE]                 # only top100 rows
      mono <- apply(grp100, 1, check_monotonicity)
      data.frame(
        miRNA          = rownames(grp100),
        transformation = key,
        monotonicity   = mono,
        stringsAsFactors = FALSE
      )
    })
    heat_df <- bind_rows(mono_list)
    
    # 5. Factor levels to fix the x-axis order
    heat_df <- heat_df %>%
      mutate(
        miRNA          = factor(miRNA, levels = top100),
        transformation = factor(transformation, levels = names(colour_map))
      )
    # Make sure factor levels include the new “log2FC” method
    heat_df <- heat_df %>%
      bind_rows(fc_df %>%
                  mutate(
                    miRNA          = factor(miRNA, levels = top100),
                    transformation = factor(transformation,
                                            levels = c(names(colour_map), "log2FC"))
                  ))
    gg <- ggplot() +
      # 1) existing monotonicity tiles
      geom_tile(
        data = filter(heat_df, transformation != "log2FC"),
        aes(x = transformation, y = miRNA, fill = monotonicity),
        color = "white", linewidth = 0.2
      ) +
      # 2) new log2FC tiles
      geom_tile(
        data = filter(heat_df, transformation == "log2FC"),
        aes(x = transformation, y = miRNA, fill = log2FC),
        color = "white", linewidth = 0.2
      ) +
      # scales for monotonicity
      scale_fill_manual(
        name   = "Trend",
        values = c(
          increasing = "#1b9e77",
          decreasing = "#d95f02",
          none       = "#FFFF00"
        ),
        labels = c(
          increasing = "Increasing",
          decreasing = "Decreasing",
          none       = "Non-monotonic"
        ),
        na.value = "white",
        guide = guide_legend(order = 1)
      ) +
      # continuous scale for the log2FC column only
      scale_fill_gradient2(
        name     = "log₂FC",
        low      = "red",
        mid      = "white",
        high     = "green",
        midpoint = 0,
        guide    = guide_colorbar(order = 2)
      ) +
      scale_x_discrete(
        labels = c(name_map, log2FC = "log₂FC"),
        expand = c(0, 0)
      ) +
      labs(
        x     = "normalization method",
        y     = paste0("miRNA (top ", topN, " by RPM)"),
        title = "Per-miRNA monotonic trend + log₂FC"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        panel.grid   = element_blank()
      )
    
    pdf(paste0(plot_dir, "top100_wFC.pdf"), width = 8, height = 10)
    print(gg)  
    dev.off()
    
    
  }
  
  
  # Build a long table of monotonicity for those 100
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp  <- geometric_mean_by_group(list_of_transformations[[key]], groups)
    grp100 <- grp[top100, , drop = FALSE]                 # only top100 rows
    mono <- apply(grp100, 1, check_monotonicity)
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = mono,
      stringsAsFactors = FALSE
    )
  })
  heat_df <- bind_rows(mono_list)
  
  # Factor levels to fix the x-axis order
  heat_df <- heat_df %>%
    mutate(
      miRNA          = factor(miRNA, levels = top100),
      transformation = factor(transformation, levels = names(colour_map))
    )
  
  # Draw the heatmap
  gg <- ggplot(heat_df,
               aes(x = transformation, y = miRNA, fill = monotonicity)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(
      name   = "trend",
      values = c(increasing = "#1b9e77",
                 decreasing = "#d95f02",
                 none       = "#FFFF00"),
      labels = c(increasing = "Increasing",
                 decreasing = "Decreasing",
                 none       = "Non-monotonic"),
      na.value = "white"
    ) +
    scale_x_discrete(
      labels = name_map,
      expand = c(0, 0)
    ) +
    labs(
      x     = "normalization method",
      y     = paste0("miRNA (top ", topN, " by RPM)"),
      title = "Per-miRNA monotonic trend"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      panel.grid   = element_blank()
    )
  
  pdf(paste0(plot_dir, "top100.pdf"), width = 8, height = 10)
  print(gg)  
  dev.off()
  
  # interactive heatmap
  # Recompute the top 100 miRNAs by RPM
  rpm_key     <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  rpm_grp     <- geometric_mean_by_group(list_of_transformations[[rpm_key]], groups)
  average_RPM <- rowMeans(rpm_grp)
  top100      <- names(sort(average_RPM, decreasing = TRUE))[1:100]
  
  # Build a long df with monotonicity + expression text
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp     <- geometric_mean_by_group(list_of_transformations[[key]], groups)
    grp100  <- grp[top100, , drop = FALSE]
    mono    <- apply(grp100, 1, check_monotonicity)
    # pack all group‐means into one tooltip string per miRNA
    expr_txt <- apply(grp100, 1, function(vals) {
      paste0(names(vals), ": ", round(vals, 2), collapse = "\n")
    })
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = mono,
      expr_text      = expr_txt,
      stringsAsFactors = FALSE
    )
  })
  heat_df <- bind_rows(mono_list) %>%
    mutate(
      miRNA          = factor(miRNA, levels = top100),
      transformation = factor(transformation, levels = names(colour_map)),
      # build the full hover text
      tooltip = paste0(
        "miRNA: ", miRNA, "\n",
        "Method: ", name_map[transformation], "\n",
        "Trend: ", monotonicity, "\n\n",
        expr_text
      )
    )
  
  # Static ggplot with text aesthetic
  p <- ggplot(heat_df,
              aes(x = transformation,
                  y = miRNA,
                  fill = monotonicity,
                  text = tooltip)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_manual(
      name   = "Trend",
      values = c(increasing = "#1b9e77",
                 decreasing = "#d95f02",
                 none       = "#FFFFFF"),
      labels = c(increasing = "Increasing",
                 decreasing = "Decreasing",
                 none       = "Non-monotonic")
    ) +
    scale_x_discrete(labels = name_map, expand = c(0, 0)) +
    labs(x = "normalization method",
         y = "miRNA (top 100 by RPM)",
         title = "Per-miRNA monotonic trend") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      panel.grid   = element_blank()
    )
  
  # Turn it interactive and save as HTML
  interactive_heatmap <- ggplotly(p, tooltip = "text")
  saveWidget(interactive_heatmap,
             file          =  paste0(plot_dir, "interactive_heatmap.html"),
             selfcontained = TRUE)
  
}

# heatmap with FC between pure groups
generate_monotonic_heatmap(
  list_of_transformations,
  groups,
  rpm_grp,
  name_map,
  plot_dir,
  100,
  "logFC"
)

### Heatmap: % monotonic vs log(RPM) expression (rank/quantile bins; common bins for all methods)

### Heatmap: per-miRNA monotonic direction (top 150 by RPM) -- y=miRNA, x=method (colors + method order like generate_monotonic_heatmap)

if (TRUE) {
  
  
  
  #  top 150 miRNAs by RPM
  rpm_df  <- list_of_transformations[["RPM_lib"]]
  rpm_grp <- geometric_mean_by_group(rpm_df, groups)
  if (is.null(rownames(rpm_grp))) rownames(rpm_grp) <- rownames(rpm_df)
  
  average_RPM <- rowMeans(rpm_grp, na.rm = TRUE)
  
  topN <- 150
  top_miR <- names(sort(average_RPM, decreasing = TRUE))[seq_len(min(topN, length(average_RPM)))]
  
  # monotonicity for each normalization
  mono_list <- lapply(names(list_of_transformations), function(key) {
    
    grp_top <- geometric_mean_by_group(list_of_transformations[[key]], groups)[top_miR, , drop = FALSE]
    
    data.frame(
      miRNA          = rownames(grp_top),
      transformation = key,
      monotonicity   = apply(grp_top, 1, check_monotonicity),
      stringsAsFactors = FALSE
    )
  })
  
  heat_df <- dplyr::bind_rows(mono_list)
  
  # Ensure full grid (missing combos -> "none")
  heat_df <- heat_df %>%
    tidyr::complete(
      transformation = names(list_of_transformations),
      miRNA = top_miR,
      fill = list(monotonicity = "none")
    )
  
  # miRNAs sorted by expression, methods as given in list_of_transformations
  heat_df <- heat_df %>%
    dplyr::mutate(
      miRNA = factor(miRNA, levels = top_miR),  # highest expression first (top of plot after reversing)
      transformation = factor(transformation, levels = names(list_of_transformations))
    )
  
  # highest-expressed miRNA at the top (ggplot draws first level at bottom)
  heat_df$miRNA <- factor(heat_df$miRNA, levels = rev(top_miR))
  
  monotone_cols <- c(
    increasing = "#1b9e77",
    decreasing = "#d95f02",
    none       = "#FFFF00"
  )
  
  p <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = transformation, y = miRNA, fill = monotonicity)
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
      labels = name_map,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x     = "normalization method",
      y     = paste0("miRNA (top ", length(top_miR), " by RPM)"),
      title = "Per-miRNA monotonic trend (sorted by expression)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y      = ggplot2::element_text(size = 6),
      panel.grid.major = ggplot2::element_line(color = "black")
    )
  
  ggplot2::ggsave(
    filename = file.path(plot_dir, "heatmap_monotonic_direction_top150_byRPM.pdf"),
    plot = p,
    width = 8,
    height = 14
  )
}


### Heatmap: per-miRNA monotonic direction (top 150 by RPM) y=miRNA, x=method

if (TRUE) {
  
  #  top 150 miRNAs by RPM (geometric mean by group, then average across groups) 
  rpm_df  <- list_of_transformations[["RPM_lib"]]
  rpm_grp <- geometric_mean_by_group(rpm_df, groups)
  if (is.null(rownames(rpm_grp))) rownames(rpm_grp) <- rownames(rpm_df)
  
  rpm_avg <- rowMeans(rpm_grp, na.rm = TRUE)
  
  topN <- 150
  top_miR <- names(sort(rpm_avg, decreasing = TRUE))[seq_len(min(topN, length(rpm_avg)))]
  
  # Method order: reuse your existing order if present
  methods <- names(list_of_transformations)
  if (exists("order_levels")) {
    methods <- as.character(order_levels)
    methods <- methods[methods %in% names(list_of_transformations)]
  }
  
  # monotonic direction per miRNA per method
  mono_list <- lapply(methods, function(meth) {
    
    temp_df <- list_of_transformations[[meth]]
    if (is.null(rownames(temp_df))) stop(paste0("Method ", meth, " has no rownames."))
    
    mirs_use <- intersect(rownames(temp_df), top_miR)
    
    grp <- geometric_mean_by_group(temp_df[mirs_use, , drop = FALSE], groups)
    if (is.null(rownames(grp))) rownames(grp) <- mirs_use
    
    mono <- apply(grp, 1, check_monotonicity)
    if (is.null(names(mono)) || length(names(mono)) != length(mono)) {
      names(mono) <- rownames(grp)
    }
    
    data.frame(
      transformation = meth,
      miRNA          = names(mono),
      trend          = ifelse(mono %in% c("increasing", "decreasing"), mono, "none"),
      stringsAsFactors = FALSE
    )
  })
  
  heat_df <- dplyr::bind_rows(mono_list)
  
  # Complete grid (missing combos -> none)
  heat_df <- heat_df %>%
    tidyr::complete(
      transformation = methods,
      miRNA = top_miR,
      fill = list(trend = "none")
    )
  
  heat_df$miRNA <- factor(heat_df$miRNA, levels = rev(top_miR))  # rev => highest at top
  heat_df$transformation <- factor(heat_df$transformation, levels = methods)
  
  p <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = transformation, y = miRNA, fill = trend)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_x_discrete(labels = name_map, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      name = "Trend",
      values = c(
        increasing = "green",
        decreasing = "red",
        none       = "white"
      ),
      breaks = c("decreasing", "increasing", "none"),
      labels = c("Decreasing", "Increasing", "Non-monotonic")
    ) +
    ggplot2::labs(
      x = "normalization method",
      y = paste0("miRNA (top ", topN, " by RPM)"),
      title = "Per-miRNA monotonic direction (sorted by expression)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6),
      legend.position = "right"
    )
  
  ggplot2::ggsave(
    filename = file.path(plot_dir, "heatmap_monotonic_direction_top150_byRPM.pdf"),
    plot = p,
    width = 10,
    height = 10
  )
}

# individual miRNA plots for supplementary

if (TRUE){
  # build sample → group table from groups list
  .build_sample_to_group <- function(groups) {
    sample_vec <- unlist(groups, use.names = FALSE)
    group_vec  <- rep(names(groups), lengths(groups))
    df <- data.frame(sample = sample_vec, group = group_vec, stringsAsFactors = FALSE)
    df$mouse <- as.numeric(sub(".*?([0-9]+)(S|L).*", "\\1", df$sample))
    df
  }
  
  #  facet plots for selected miRNAs
  plot_selected_miRNAs_faceted <- function(
    mirnas,
    list_of_transformations,
    groups,
    output_pdf      = "./plots/selected_scatter_by_method.pdf",
    ncol_facets     = 4,
    free_y_scales   = TRUE,
    show_summary    = TRUE,
    method_order    = names(list_of_transformations),
    method_labels   = NULL
  ) {
    stopifnot(is.character(mirnas), length(mirnas) > 0)
    
    # 1) map samples to groups and mice
    sample_to_group_df <- .build_sample_to_group(groups)
    
    # 2) long table across all methods but only for requested miRNAs
    long_list <- lapply(names(list_of_transformations), function(method_key) {
      mat <- list_of_transformations[[method_key]]
      # restrict to miRNAs available in this matrix
      present <- intersect(mirnas, rownames(mat))
      if (length(present) == 0) return(NULL)
      
      df <- as.data.frame(mat[present, , drop = FALSE])
      df$miRNA <- rownames(df)
      
      df %>%
        tidyr::pivot_longer(
          cols = -miRNA,
          names_to  = "sample",
          values_to = "value"
        ) %>%
        dplyr::left_join(sample_to_group_df, by = "sample") %>%
        dplyr::filter(!is.na(group)) %>%
        dplyr::mutate(method = method_key)
    })
    
    expr_long <- dplyr::bind_rows(long_list)
    if (nrow(expr_long) == 0) {
      warning("None of the requested miRNAs were found in the provided transformations.")
      return(invisible(NULL))
    }
    
    # 3) factors and labels
    expr_long <- expr_long %>%
      dplyr::mutate(
        miRNA  = factor(miRNA, levels = unique(mirnas)),
        method = factor(method, levels = method_order),
        group  = factor(group, levels = c("0S", "20S", "40S", "60S", "80S", "100S")),
        mouse  = factor(mouse)
      )
    
    # optional pretty labels for facets
    if (!is.null(method_labels)) {
      # ensure labels cover all levels; fallback to raw level if missing
      lvl <- levels(expr_long$method)
      method_lab_vec <- setNames(ifelse(lvl %in% names(method_labels), method_labels[lvl], lvl), lvl)
      lab_fun <- ggplot2::labeller(method = function(x) method_lab_vec[x])
    } else {
      lab_fun <- ggplot2::label_value()
    }
    
    # 4) open device
    pdf(output_pdf, width = 8, height = 6)
    
    # 5) iterate miRNAs and draw
    for (mir in levels(expr_long$miRNA)) {
      df_sub <- dplyr::filter(expr_long, miRNA == mir)
      
      if (show_summary) {
        summ_df <- df_sub %>%
          dplyr::group_by(method, group) %>%
          dplyr::summarise(
            mean_val = mean(value),
            geo_mean = exp(mean(log(value))),
            geo_trim = {
              vals   <- value
              drop_i <- if (length(vals) > 1) which.max(abs(vals - mean(vals))) else integer(0)
              if (length(drop_i)) exp(mean(log(vals[-drop_i]))) else exp(mean(log(vals)))
            },
            .groups = "drop"
          )
      }
      
      p <- ggplot2::ggplot(
        df_sub,
        aes(x = group, y = value, group = mouse, color = mouse)
      ) +
        ggplot2::geom_line(alpha = 0.8, linewidth = 0.7) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        {
          if (show_summary)
            list(
              ggplot2::geom_point(
                data = summ_df, aes(x = group, y = mean_val),
                inherit.aes = FALSE, shape = 45, size = 6, color = "black"
              ),
              ggplot2::geom_point(
                data = summ_df, aes(x = group, y = geo_mean),
                inherit.aes = FALSE, shape = 95, size = 6, color = "darkgreen"
              ),
              ggplot2::geom_point(
                data = summ_df, aes(x = group, y = geo_trim),
                inherit.aes = FALSE, shape = 95, size = 6, color = "darkorange"
              )
            )
          else NULL
        } +
        ggplot2::facet_wrap(
          ~ method,
          ncol   = ncol_facets,
          scales = if (free_y_scales) "free_y" else "fixed",
          labeller = lab_fun
        ) +
        ggplot2::labs(
          title = paste0(mir, " – expression by method"),
          x     = "Mixture group",
          y     = "Expression value"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          strip.text       = ggplot2::element_text(size = 8),
          axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        )
      
      print(p)
    }
    
    dev.off()
    message("Saved: ", output_pdf)
    invisible(output_pdf)
  }
  
  
  plot_selected_miRNAs_faceted_page <- function(
    mirnas,
    list_of_transformations,
    groups,
    output_pdf      = "./plots/selected_scatter_by_method_and_miRNA.pdf",
    show_summary    = TRUE,
    method_order    = names(list_of_transformations),
    method_labels   = NULL
  ) {
    stopifnot(is.character(mirnas), length(mirnas) > 0)
    
    if (!requireNamespace("ggh4x", quietly = TRUE)) {
      stop("Please install ggh4x first: install.packages('ggh4x')")
    }
    
    sample_to_group_df <- .build_sample_to_group(groups)
    
    long_list <- lapply(names(list_of_transformations), function(method_key) {
      mat <- list_of_transformations[[method_key]]
      present <- intersect(mirnas, rownames(mat))
      if (length(present) == 0) return(NULL)
      
      df <- as.data.frame(mat[present, , drop = FALSE])
      df$miRNA <- rownames(df)
      
      df %>%
        tidyr::pivot_longer(
          cols = -miRNA,
          names_to = "sample",
          values_to = "value"
        ) %>%
        dplyr::left_join(sample_to_group_df, by = "sample") %>%
        dplyr::filter(!is.na(group)) %>%
        dplyr::mutate(method = method_key)
    })
    
    expr_long <- dplyr::bind_rows(long_list)
    
    if (nrow(expr_long) == 0) {
      warning("None of the requested miRNAs were found.")
      return(invisible(NULL))
    }
    
    expr_long <- expr_long %>%
      dplyr::mutate(
        miRNA  = factor(miRNA, levels = mirnas),
        method = factor(method, levels = method_order),
        group  = factor(group, levels = c("0S", "20S", "40S", "60S", "80S", "100S")),
        mouse  = factor(mouse)
      ) %>%
      dplyr::filter(!is.na(method), !is.na(miRNA))
    
    if (!is.null(method_labels)) {
      lvl <- levels(expr_long$method)
      method_lab_vec <- setNames(
        ifelse(lvl %in% names(method_labels), method_labels[lvl], lvl),
        lvl
      )
      method_labeller <- function(x) method_lab_vec[x]
    } else {
      method_labeller <- ggplot2::label_value
    }
    
    if (show_summary) {
      summ_df <- expr_long %>%
        dplyr::group_by(method, miRNA, group) %>%
        dplyr::summarise(
          mean_val = mean(value, na.rm = TRUE),
          geo_mean = {
            vals <- value[is.finite(value) & !is.na(value) & value > 0]
            if (length(vals) == 0) NA_real_ else exp(mean(log(vals)))
          },
          geo_trim = {
            vals <- value[is.finite(value) & !is.na(value) & value > 0]
            if (length(vals) == 0) {
              NA_real_
            } else if (length(vals) == 1) {
              exp(mean(log(vals)))
            } else {
              drop_i <- which.max(abs(vals - mean(vals)))
              exp(mean(log(vals[-drop_i])))
            }
          },
          .groups = "drop"
        )
    }
    
    n_methods <- nlevels(expr_long$method)
    n_mirnas  <- length(mirnas)
    
    pdf(
      output_pdf,
      width  = max(3.2 * n_mirnas, 8),
      height = max(2.8 * n_methods, 6)
    )
    
    p <- ggplot2::ggplot(
      expr_long,
      ggplot2::aes(x = group, y = value, group = mouse, color = mouse)
    ) +
      ggplot2::geom_line(alpha = 0.8, linewidth = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      {
        if (show_summary)
          list(
            ggplot2::geom_point(
              data = summ_df,
              ggplot2::aes(x = group, y = mean_val),
              inherit.aes = FALSE,
              shape = 45, size = 6, color = "black"
            ),
            ggplot2::geom_point(
              data = summ_df,
              ggplot2::aes(x = group, y = geo_mean),
              inherit.aes = FALSE,
              shape = 95, size = 6, color = "darkgreen"
            ),
            ggplot2::geom_point(
              data = summ_df,
              ggplot2::aes(x = group, y = geo_trim),
              inherit.aes = FALSE,
              shape = 95, size = 6, color = "darkorange"
            )
          )
        else NULL
      } +
      ggh4x::facet_grid2(
        rows = ggplot2::vars(method),
        cols = ggplot2::vars(miRNA),
        scales = "free_y",
        independent = "y",
        labeller = ggplot2::labeller(
          method = method_labeller,
          miRNA  = ggplot2::label_value
        )
      ) +
      ggplot2::labs(
        title = "Selected miRNAs – expression by method",
        x = "Mixture group",
        y = "Expression value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text.x    = ggplot2::element_text(size = 9, face = "bold"),
        strip.text.y    = ggplot2::element_text(size = 9, face = "bold"),
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid      = ggplot2::element_blank(),
        legend.position = "bottom"
      )
    
    print(p)
    dev.off()
    
    message("Saved: ", output_pdf)
    invisible(output_pdf)
  }
  
  
  
  # choose miRNAs explicitly
  sel <- c(
    "mmu-miR-21a-5p", 
    "mmu-miR-7a-5p",
    "mmu-let-7a-5p",
    "mmu-miR-26a-5p"
    # "mmu-miR-192-5p", 
    # "mmu-miR-99a-5p"
    # "mmu-miR-30c-5p", 
    # "mmu-miR-100-5p", 
    # "mmu-miR-126a-3p", 
    # "mmu-miR-92a-3p", 
    # "mmu-miR-7a-5p", 
    # "mmu-miR-30e-5p", 
    # "mmu-miR-10b-5p", 
    # "mmu-miR-143-3p"
  )
  
  plot_selected_miRNAs_faceted(
    mirnas               = sel,
    list_of_transformations = list_of_transformations,
    groups               = groups,
    output_pdf           = file.path(plot_dir, "selected_scatter_by_method.pdf"),
    ncol_facets          = 3,
    free_y_scales        = TRUE,
    show_summary         = TRUE,
    method_order         = names(list_of_transformations),
    method_labels        = name_map   # your human-readable names
  )
  
  plot_selected_miRNAs_faceted_page(
    mirnas                  = sel,
    list_of_transformations = list_of_transformations,
    groups                  = groups,
    output_pdf              = file.path(plot_dir, "selected_scatter_by_method.pdf"),
    # free_y_scales           = TRUE,
    show_summary            = TRUE,
    method_labels           = name_map
  )
  
  
}

if (TRUE){
  .build_sample_to_group <- function(groups) {
    sample_vec <- unlist(groups, use.names = FALSE)
    group_vec  <- rep(names(groups), lengths(groups))
    
    df <- data.frame(
      sample = sample_vec,
      group  = group_vec,
      stringsAsFactors = FALSE
    )
    
    df$mouse <- as.numeric(sub(".*?([0-9]+)(S|L).*", "\\1", df$sample))
    df
  }
  
  
  plot_selected_miRNAs_2x4 <- function(
    mirnas,
    list_of_transformations,
    groups,
    output_dir     = "./plots",
    method_order   = names(list_of_transformations),
    method_labels  = NULL,
    show_summary   = TRUE,
    page_width     = 14,
    page_height    = 8
  ) {
    stopifnot(is.character(mirnas), length(mirnas) > 0)
    
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Please install patchwork first: install.packages('patchwork')")
    }
    
    sample_to_group_df <- .build_sample_to_group(groups)
    group_levels <- c("0S", "20S", "40S", "60S", "80S", "100S")
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Keep only valid methods, in requested order
    method_order <- method_order[method_order %in% names(list_of_transformations)]
    
    if (length(method_order) == 0) {
      stop("No valid methods found in method_order.")
    }
    
    get_method_label <- function(x) {
      if (!is.null(method_labels) && x %in% names(method_labels)) {
        return(method_labels[[x]])
      }
      x
    }
    
    for (mir in mirnas) {
      
      plot_list <- vector("list", length(method_order))
      names(plot_list) <- method_order
      
      for (i in seq_along(method_order)) {
        method_key <- method_order[i]
        mat <- list_of_transformations[[method_key]]
        
        if (!(mir %in% rownames(mat))) {
          p_empty <- ggplot2::ggplot() +
            ggplot2::annotate(
              "text",
              x = 0.5, y = 0.5,
              label = paste0(get_method_label(method_key), "\n\nmiRNA not found"),
              size = 4
            ) +
            ggplot2::xlim(0, 1) +
            ggplot2::ylim(0, 1) +
            ggplot2::theme_void() +
            ggplot2::ggtitle(get_method_label(method_key)) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10)
            )
          
          plot_list[[i]] <- p_empty
          next
        }
        
        df <- as.data.frame(mat[mir, , drop = FALSE])
        df$miRNA <- rownames(df)
        
        df_long <- df %>%
          tidyr::pivot_longer(
            cols = -miRNA,
            names_to = "sample",
            values_to = "value"
          ) %>%
          dplyr::left_join(sample_to_group_df, by = "sample") %>%
          dplyr::filter(!is.na(group)) %>%
          dplyr::mutate(
            group = factor(group, levels = group_levels),
            mouse = factor(mouse)
          )
        
        if (nrow(df_long) == 0) {
          p_empty <- ggplot2::ggplot() +
            ggplot2::annotate(
              "text",
              x = 0.5, y = 0.5,
              label = paste0(get_method_label(method_key), "\n\nno valid samples"),
              size = 4
            ) +
            ggplot2::xlim(0, 1) +
            ggplot2::ylim(0, 1) +
            ggplot2::theme_void() +
            ggplot2::ggtitle(get_method_label(method_key)) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10)
            )
          
          plot_list[[i]] <- p_empty
          next
        }
        
        if (show_summary) {
          summ_df <- df_long %>%
            dplyr::group_by(group) %>%
            dplyr::summarise(
              mean_val = mean(value, na.rm = TRUE),
              geo_mean = {
                vals <- value[is.finite(value) & !is.na(value) & value > 0]
                if (length(vals) == 0) NA_real_ else exp(mean(log(vals)))
              },
              geo_trim = {
                vals <- value[is.finite(value) & !is.na(value) & value > 0]
                if (length(vals) == 0) {
                  NA_real_
                } else if (length(vals) == 1) {
                  exp(mean(log(vals)))
                } else {
                  drop_i <- which.max(abs(vals - mean(vals)))
                  exp(mean(log(vals[-drop_i])))
                }
              },
              .groups = "drop"
            )
        }
        
        p <- ggplot2::ggplot(
          df_long,
          ggplot2::aes(x = group, y = value, group = mouse, color = mouse)
        ) +
          ggplot2::geom_line(alpha = 0.8, linewidth = 0.7) +
          ggplot2::geom_point(size = 2) +
          ggplot2::scale_color_brewer(palette = "Set1") +
          ggplot2::labs(
            title = get_method_label(method_key),
            x = "Mixture group",
            y = "Expression value"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title       = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
            axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position  = "none"
          )
        
        if (show_summary) {
          p <- p +
            ggplot2::geom_point(
              data = summ_df,
              ggplot2::aes(x = group, y = mean_val),
              inherit.aes = FALSE,
              shape = 45, size = 6, color = "black"
            ) +
            ggplot2::geom_point(
              data = summ_df,
              ggplot2::aes(x = group, y = geo_mean),
              inherit.aes = FALSE,
              shape = 95, size = 6, color = "darkgreen"
            ) +
            ggplot2::geom_point(
              data = summ_df,
              ggplot2::aes(x = group, y = geo_trim),
              inherit.aes = FALSE,
              shape = 95, size = 6, color = "darkorange"
            )
        }
        
        plot_list[[i]] <- p
      }
      
      combined <- patchwork::wrap_plots(plot_list, ncol = 4) +
        patchwork::plot_annotation(
          title = mir,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face  = "bold",
              size  = 14
            )
          )
        )
      
      out_file <- file.path(output_dir, paste0(mir, "_scatter_facet.pdf"))
      
      grDevices::pdf(out_file, width = page_width, height = page_height)
      print(combined)
      grDevices::dev.off()
      
      message("Saved: ", out_file)
    }
    
    invisible(NULL)
  }
  
  plot_selected_miRNAs_2x4(
    mirnas                  = sel,
    list_of_transformations = list_of_transformations,
    groups                  = groups,
    output_dir              = plot_dir,
    method_order            = names(list_of_transformations),
    method_labels           = name_map,
    show_summary            = TRUE
  )
}

if (TRUE){
  # =========================================================
  # Plot top 100 most expressed miRNAs in separate folder
  # =========================================================
  
  top100_plot_dir <- file.path(plot_dir, "top100_most_expressed")
  dir.create(top100_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  top100_miRNAs <- group_extended %>%
    tibble::rownames_to_column("miRNA") %>%
    dplyr::arrange(dplyr::desc(average)) %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::pull(miRNA)
  
  plot_selected_miRNAs_2x4(
    mirnas                  = top100_miRNAs,
    list_of_transformations = list_of_transformations,
    groups                  = groups,
    output_dir              = top100_plot_dir,
    method_order            = names(list_of_transformations),
    method_labels           = name_map,
    show_summary            = TRUE
  )
}