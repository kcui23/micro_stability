library(DESeq2)
library(readr)

run_deseq2 <- function(ASV_file, groupings_file, output_file) {
  # Read ASV table
  ASV_table <- read_tsv(ASV_file, comment = "", col_names = TRUE, skip = ifelse(grepl("Constructed from biom file", readLines(ASV_file, n=1)), 1, 0))
  ASV_table <- as.data.frame(ASV_table)
  row.names(ASV_table) <- ASV_table[,1]
  ASV_table <- ASV_table[,-1]

  # Read groupings table
  groupings <- read_tsv(groupings_file, col_names = TRUE)
  groupings <- as.data.frame(groupings)
  row.names(groupings) <- groupings[,1]
  groupings <- groupings[,-1, drop=FALSE]

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

  # Run DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                        colData = groupings,
                                        design = ~ comparison)
  dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")

  res <- DESeq2::results(dds_res, tidy = TRUE, format = "DataFrame")

  rownames(res) <- res$row
  res <- res[,-1]

  # Save results to output file
  write_tsv(as.data.frame(res), output_file)

  message("Results written to ", output_file)
  
  # Return paths to the result files (for example, if there are plots)
  list(
    plot1 = "./plot1.png",
    plot2 = "./plot2.png",
    plot3 = "./plot3.png"
  )
}

visualize_deseq2 <- function(input_file, output_dir) {
  # Read DESeq2 results
  deseq2_results <- read_tsv(input_file)
  
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
  deseq2_results$log10pvalue <- -log10(deseq2_results$pvalue)
  
  p1 <- ggplot(deseq2_results, aes(x = log2FoldChange, y = log10pvalue)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Volcano plot", x = "log2FoldChange", y = "-log10(pvalue)")
  
  ggsave(file.path(output_dir, "deseq2_plot1.png"), plot = p1)
  
  # Visualization 2: MA plot of baseMean vs. log2FoldChange
  p2 <- ggplot(deseq2_results, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_x_log10() +
    labs(title = "MA plot", x = "baseMean", y = "log2FoldChange")
  
  ggsave(file.path(output_dir, "deseq2_plot2.png"), plot = p2)
  
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(deseq2_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = "Histogram of p-value distribution", x = "p-value", y = "Frequency")
  
  ggsave(file.path(output_dir, "deseq2_plot3.png"), plot = p3)
  
  list(
    plot1 = file.path(output_dir, "deseq2_plot1.png"),
    plot2 = file.path(output_dir, "deseq2_plot2.png"),
    plot3 = file.path(output_dir, "deseq2_plot3.png")
  )
}
