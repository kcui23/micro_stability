options(
  warn = -1,            
  message = FALSE,      
  verbose = FALSE,      
  dplyr.summarise.inform = FALSE,  
  readr.show_progress = FALSE      
)
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
})
run_deseq2 <- function(ASV_file, groupings_file, output_file, seed = 1234) {
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
  print('======dim in deseq2==========')
  print(dim(ASV_table))

  # Read groupings table
  groupings <- read_tsv(groupings_file, col_names = TRUE, show_col_types = FALSE)
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
  }

  # Run DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                        colData = groupings,
                                        design = ~ comparison)
  dds_res <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)

  res <- DESeq2::results(dds_res, tidy = TRUE, format = "DataFrame")

  # Rename 'row' column to 'asv_name'
  colnames(res)[colnames(res) == "row"] <- "asv_name"
  
  # Ensure 'asv_name' is the first column
  res <- res[, c("asv_name", setdiff(colnames(res), "asv_name"))]

  # Save results to output file
  write_tsv(as.data.frame(res), output_file)
  

  list(
    plot1 = "./plot1.png",
    plot2 = "./plot2.png",
    plot3 = "./plot3.png"
  )
}

visualize_deseq2 <- function(input_file, output_dir, persistent_temp_dir) {
  # Read DESeq2 results
  deseq2_results <- read_tsv(input_file, show_col_types = FALSE)
  
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
  deseq2_results$log10pvalue <- -log10(deseq2_results$pvalue)
  deseq2_results$point_size <- ifelse(deseq2_results$pvalue < 0.05, 2, 1)
  
  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- deseq2_results %>%
    dplyr::select(asv_name, log2FoldChange, pvalue) %>%
    mutate(method = "DESeq2")
  write_tsv(volcano_plot_data, file.path(persistent_temp_dir, "deseq2_volcano_plot_data.tsv"))
  
  p1 <- ggplot(deseq2_results, aes(x = log2FoldChange, y = log10pvalue)) +
    geom_point(aes(color = pvalue < 0.05, size=point_size)) +
    scale_color_manual(values = c("gray", "#4C3BCF")) +
    scale_size_continuous(range = c(1, 3), guide = "none") +
    theme_minimal() +
    labs(title = "Volcano plot", x = "log2FoldChange", y = "-log10(pvalue)")
  
  ggsave(file.path(output_dir, "deseq2_plot1.png"), plot = p1)
  
  # Visualization 2: MA plot of baseMean vs. log2FoldChange
  p2 <- ggplot(deseq2_results, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05), alpha = 0.5) +
    scale_color_manual(values = c("gray", "#4C3BCF"), name = "Significant") +
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

  print('finished deseq2 visualization')
  
  list(
    plot1 = file.path(output_dir, "deseq2_plot1.png"),
    plot2 = file.path(output_dir, "deseq2_plot2.png"),
    plot3 = file.path(output_dir, "deseq2_plot3.png")
  )
}
