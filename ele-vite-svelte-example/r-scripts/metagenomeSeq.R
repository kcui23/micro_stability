options(
  warn = -1,            
  verbose = FALSE  
)
suppressPackageStartupMessages({
  library(metagenomeSeq)
  library(tidyverse)
})


run_metagenomeseq <- function(ASV_file, groupings_file, output_file, seed = 1234) {
  suppressWarnings({
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
    print('======dim in metagenomeSeq==========')
    print(dim(ASV_table))

    # Read groupings table
    groupings <- read_tsv(groupings_file, col_names = TRUE, show_col_types = FALSE)
    groupings <- as.data.frame(groupings)
    row.names(groupings) <- groupings[,1]

    # Check for matching samples
    sample_num <- ncol(ASV_table)
    grouping_num <- nrow(groupings)

    if (sample_num != grouping_num) {
      message("The number of samples in the ASV table and the groupings table are unequal")
      message("Removing samples not found in both ASV table and groupings")
    }

    # Align ASV table and metadata
    if (identical(colnames(ASV_table), rownames(groupings))) {
      message("Groupings and ASV table are in the same order")
    } else {
      rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
      groupings <- groupings[rows_to_keep,, drop=F]
      ASV_table <- ASV_table[, rows_to_keep]
    }

    # Create a MRexperiment object
    count_matrix <- as.matrix(ASV_table)
    groupings_meta <- AnnotatedDataFrame(groupings)

    # Create MRexperiment object
    MRexp <- metagenomeSeq::newMRexperiment(counts = count_matrix, phenoData = groupings_meta)
    MRexp_norm <- metagenomeSeq::cumNorm(MRexp)
    comparison_variable <- colnames(groupings)[2]

    # Fit the feature model
    model <- model.matrix(~ groupings[, comparison_variable])
    fit <- metagenomeSeq::fitFeatureModel(MRexp_norm, model)

    # Get differential abundance results
    result <- metagenomeSeq::MRfulltable(fit, number = nrow(fit))
    result <- result %>%
      arrange(pvalues) %>%
      mutate(asv_name = rownames(.)) %>%
      dplyr::select(asv_name, everything())

    write_tsv(as.data.frame(result), output_file)
  })

  return(result)
}



visualize_metagenomeseq <- function(input_file, output_dir, persistent_temp_dir) {
  # Read metagenomeSeq results
  metagenomeSeq_results <- read_tsv(input_file, show_col_types = FALSE)
  
  # Visualization 1: Volcano plot of logFC vs. -log10(pvalue)
  metagenomeSeq_results$log10pvalue <- -log10(metagenomeSeq_results$pvalues)
  metagenomeSeq_results$point_size <- ifelse(metagenomeSeq_results$pvalues < 0.05, 2, 1)

  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- metagenomeSeq_results %>%
    dplyr::select(asv_name, logFC, pvalues) %>%
    mutate(method = "metagenomeSeq", log2FoldChange = logFC, pvalue = pvalues)
  write_tsv(volcano_plot_data, file.path(persistent_temp_dir, "metagenomeseq_volcano_plot_data.tsv"))
  
  p1 <- ggplot(metagenomeSeq_results, aes(x = logFC, y = log10pvalue)) +
    geom_point(aes(color = pvalues < 0.05, size = point_size)) +
    scale_color_manual(values = c("gray", "#4C3BCF")) +
    scale_size_continuous(range = c(1, 3), guide = "none") +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)")
  
  ggsave(file.path(output_dir, "metagenomeseq_plot1.png"), plot = p1)
  
  # Visualization 2: MA plot of baseMean vs. log2 Fold Change
  metagenomeSeq_results$baseMean <- rowMeans(as.matrix(metagenomeSeq_results[, -c(1:4)]))
  
  p2 <- ggplot(metagenomeSeq_results, aes(x = baseMean, y = logFC)) +
    geom_point(aes(color = pvalues < 0.05), alpha = 0.5) +
    scale_color_manual(values = c("gray", "#4C3BCF"), name = "Significant") +
    theme_minimal() +
    scale_x_log10() +
    labs(title = "MA Plot", x = "baseMean", y = "log2 Fold Change")
  
  ggsave(file.path(output_dir, "metagenomeseq_plot2.png"), plot = p2)
  
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(metagenomeSeq_results, aes(x = pvalues)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = "Histogram of p-value Distribution", x = "p-value", y = "Frequency")
  
  ggsave(file.path(output_dir, "metagenomeseq_plot3.png"), plot = p3)

  print('finished metagenomeseq visualization')
  
  # Return the paths to the saved plots
  list(
    plot1 = file.path(output_dir, "metagenomeseq_plot1.png"),
    plot2 = file.path(output_dir, "metagenomeseq_plot2.png"),
    plot3 = file.path(output_dir, "metagenomeseq_plot3.png")
  )
}