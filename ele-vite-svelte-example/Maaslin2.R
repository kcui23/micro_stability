library(Maaslin2)
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

run_maaslin2 <- function(ASV_file, groupings_file, output_file, output_dir, seed = 1234) {
  set.seed(seed)
  print('======seed==========')
  print(seed)
  print(read_tsv(ASV_file))
  print('======output_file==========')
  print(output_file)

  print('======output_dir==========')
  print(output_dir)

  # Read ASV table
  print("Reading ASV table...")
  ASV_table <- read_tsv(ASV_file, comment = "", col_names = TRUE, skip = ifelse(grepl("Constructed from biom file", readLines(ASV_file, n=1)), 1, 0))
  ASV_table <- as.data.frame(ASV_table)
  row.names(ASV_table) <- ASV_table[,1]
  ASV_table <- ASV_table[,-1]
  print("ASV table read successfully.")

  # Read groupings table
  print("Reading groupings table...")
  groupings <- read_tsv(groupings_file, col_names = TRUE)
  groupings <- as.data.frame(groupings)
  row.names(groupings) <- groupings[,1]
  print("Groupings table read successfully.")

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
    if (identical(colnames(ASV_table), rownames(groupings))) {
      message("Groupings table was re-arranged to be in the same order as the ASV table")
      message("A total of ", sample_num - length(colnames(ASV_table)), " from the ASV_table")
      message("A total of ", grouping_num - length(rownames(groupings)), " from the groupings table")
    } else {
      stop("Unable to match samples between the ASV table and groupings table")
    }
  }

  # Transpose ASV_table for Maaslin2
  print("Transposing ASV table...")
  ASV_table <- data.frame(t(ASV_table), check.rows = F, check.names = F, stringsAsFactors = F)
  print("ASV table transposed successfully.")

  print("'-'-'-'-'-'")
  print(output_file)
  print(output_dir)
  print(seed)
  print("'-'-'-'-'-'")

  # Run Maaslin2 analysis
  print("Running Maaslin2 analysis...")
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
  print("Maaslin2 analysis completed.")

  message("Maaslin2 analysis completed and results saved to ", output_dir)

  # print all files under output_dir
  print("-------------")
  print("Files under output_dir:")
  print(list.files(output_dir))

  # copy the 'all_results.tsv' to output_file
  file.copy(file.path(output_dir, "all_results.tsv"), output_file, overwrite = TRUE)
  print(head(read_tsv(output_file)))

  list(
    plot1 = "./plot1.png",
    plot2 = "./plot2.png",
    plot3 = "./plot3.png"
  )
}

visualize_maaslin2 <- function(results_file, output_dir) {
  # Read Maaslin2 results
  print("Reading Maaslin2 results...")
  maaslin2_results <- read_tsv(results_file)

  print(colnames(maaslin2_results))
  print('-=-=-=-=-=-=-=-==')
  
  # Visualization 1: Volcano plot of effect vs. -log10(pvalue)
  maaslin2_results <- maaslin2_results %>%
    mutate(neg_log10_pvalue = -log10(pval),
           point_size = ifelse(qval < 0.05, 3, 1))
  
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
  p2 <- ggplot(maaslin2_results, aes(x = N.not.0, y = coef)) +
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
  
  list(
    plot1 = file.path(output_dir, "maaslin2_plot1.png"),
    plot2 = file.path(output_dir, "maaslin2_plot2.png"),
    plot3 = file.path(output_dir, "maaslin2_plot3.png")
  )
}
