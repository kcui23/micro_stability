suppressPackageStartupMessages({
  library(Maaslin2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

run_maaslin2 <- function(ASV_file, groupings_file, output_file, output_dir, seed = 1234) {
  set.seed(seed)
  print('======seed==========')
  print(seed)

  # Read ASV table
  ASV_table <- read_tsv(ASV_file, comment = "", col_names = TRUE, 
                        skip = ifelse(grepl("Constructed from biom file", readLines(ASV_file, n=1)), 1, 0),
                        show_col_types = FALSE)
  ASV_table <- as.data.frame(ASV_table)
  row.names(ASV_table) <- ASV_table[,1]
  ASV_table <- ASV_table[,-1]
  print('======dim in maaslin2==========')
  print(dim(ASV_table))

  # Read groupings table
  groupings <- read_tsv(groupings_file, col_names = TRUE, show_col_types = FALSE)
  groupings <- as.data.frame(groupings)
  row.names(groupings) <- groupings[,1]

  # Check if the same number of samples are being input
  sample_num <- ncol(ASV_table)
  grouping_num <- nrow(groupings)

  if (sample_num != grouping_num) {
    message("The number of samples in the ASV table and the groupings table are unequal")
    message("Will remove any samples that are not found in either the ASV table or the groupings table")
  }

  # Check if order of samples match up
  if (identical(colnames(ASV_table), rownames(groupings))) {
    message("Groupings and ASV table are in the same order")
  } else {
    rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
    groupings <- groupings[rows_to_keep,, drop=F]
    ASV_table <- ASV_table[, rows_to_keep]
  }

  # Transpose ASV_table for Maaslin2
  ASV_table <- data.frame(t(ASV_table), check.rows = F, check.names = F, stringsAsFactors = F)

  # Run Maaslin2 analysis
  sink(file.path(tempdir(), "maaslin2_output.txt"), append = TRUE)
  fit_data <- Maaslin2(
    ASV_table, 
    groupings, 
    output_dir, 
    transform = "AST",
    fixed_effects = c(colnames(groupings)[2]),
    standardize = FALSE, 
    plot_heatmap = FALSE, 
    plot_scatter = FALSE
  )
  sink()

  # copy the 'all_results.tsv' to output_file
  file.copy(file.path(output_dir, "all_results.tsv"), output_file, overwrite = TRUE)
  
  # Read the all_results.tsv file
  all_results <- read_tsv(file.path(output_dir, "all_results.tsv"), show_col_types = FALSE)
  
  # Rename the 'feature' column to 'asv_name'
  colnames(all_results)[colnames(all_results) == "feature"] <- "asv_name"

  # Ensure 'asv_name' is the first column
  all_results <- all_results[, c("asv_name", setdiff(colnames(all_results), "asv_name"))]

  # Write the updated results to the output file
  write_tsv(all_results, output_file)

  list(
    plot1 = "./plot1.png",
    plot2 = "./plot2.png",
    plot3 = "./plot3.png"
  )
}

visualize_maaslin2 <- function(results_file, output_dir, persistent_temp_dir) {
  # Read Maaslin2 results
  maaslin2_results <- read_tsv(results_file, show_col_types = FALSE)
  
  # Visualization 1: Volcano plot of effect vs. -log10(pvalue)
  maaslin2_results <- maaslin2_results %>%
    mutate(neg_log10_pvalue = -log10(pval),
           point_size = ifelse(qval < 0.05, 3, 1))

  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- maaslin2_results %>%
    dplyr::select(asv_name, coef, pval) %>%
    mutate(method = "Maaslin2", log2FoldChange = coef, pvalue = pval)
  write_tsv(volcano_plot_data, file.path(persistent_temp_dir, "maaslin2_volcano_plot_data.tsv"))
  
  p1 <- ggplot(maaslin2_results, aes(x = coef, y = neg_log10_pvalue)) +
    geom_point(aes(color = qval < 0.05, size=point_size)) +
    scale_color_manual(values = c("gray", "#4C3BCF")) +
    scale_size_continuous(range = c(1, 3), guide = "none") +
    theme_minimal() +
    labs(title = "Volcano plot of Effect Size vs. -log10(pvalue)",
         x = "Effect Size",
         y = "-log10(pvalue)")
  
  ggsave(file.path(output_dir, "maaslin2_plot1.png"), plot = p1)

  # Visualization 2: MA-like scatter plot
  #use 'N.not.0' as a proxy for 'baseMean' and 'coef' for 'log2FoldChange'
  # Check if 'N.not.0' exists, otherwise use 'N.not.zero'
  x_column <- if("N.not.0" %in% colnames(maaslin2_results)) "N.not.0" else "N.not.zero"
  p2 <- ggplot(maaslin2_results, aes(x = .data[[x_column]], y = coef)) +
    geom_point(aes(color = qval < 0.05), alpha = 0.5) +
    scale_color_manual(values = c("gray", "#4C3BCF")) +
    theme_minimal() +
    scale_x_log10() +
    labs(title = "MA-like plot",
         x = "Number of Non-Zero Samples",
         y = "Coefficient (Effect Size)")
  
  ggsave(file.path(output_dir, "maaslin2_plot2.png"), plot = p2)

  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(maaslin2_results, aes(x = pval)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() +
    labs(title = "Histogram of p-value distribution", x = "pvalue", y = "Frequency")

  ggsave(file.path(output_dir, "maaslin2_plot3.png"), plot = p3)

  print('finished maaslin2 visualization')
  
  list(
    plot1 = file.path(output_dir, "maaslin2_plot1.png"),
    plot2 = file.path(output_dir, "maaslin2_plot2.png"),
    plot3 = file.path(output_dir, "maaslin2_plot3.png")
  )
}
