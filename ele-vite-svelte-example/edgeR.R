options(
  warn = -1,            
  verbose = FALSE  
)
suppressPackageStartupMessages({
  library(edgeR)
  library(phyloseq)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Define the function to convert phyloseq object to edgeR DGEList
phyloseq_to_edgeR <- function(physeq, group, method="RLE", ...){
  # Enforce orientation
  if(!taxa_are_rows(physeq)) { physeq <- t(physeq) }
  x <- as(otu_table(physeq), "matrix")
  x <- x + 1  # Add one to protect against overflow, log(0) issues
  
  # Check `group` argument
  if(identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1){
    group <- get_variable(physeq, group)
  }
  
  # Define gene annotations (`genes`) as tax_table
  taxonomy <- tax_table(physeq, errorIfNULL=FALSE)
  if(!is.null(taxonomy)){
    taxonomy <- data.frame(as(taxonomy, "matrix"))
  }
  
  # Now turn into a DGEList
  y <- DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  
  # Calculate the normalization factors
  z <- calcNormFactors(y, method=method)
  
  # Check for division by zero inside `calcNormFactors`
  if(!all(is.finite(z$samples$norm.factors))){
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors, consider changing `method` argument")
  }
  
  # Estimate dispersions
  if (length(unique(group)) > 1 && min(table(group)) > 1) {
    z <- estimateCommonDisp(z)
    z <- estimateTagwiseDisp(z)
  } else {
    z$common.dispersion <- NA
    warning("No replication, setting dispersion to NA.")
  }
  
  return(z)
}

# Define the main function to run edgeR analysis
run_edgeR <- function(ASV_file, groupings_file, output_file, seed = 1234) {
  set.seed(seed)
  print('======seed==========')
  print(seed)
  
  # Read ASV table - suppress messages from read_tsv
  ASV_table <- suppressMessages(read_tsv(ASV_file, comment = "", col_names = TRUE, 
    skip = ifelse(grepl("Constructed from biom file", readLines(ASV_file, n=1)), 1, 0)))
  ASV_table <- as.data.frame(ASV_table)
  row.names(ASV_table) <- ASV_table[,1]
  ASV_table <- ASV_table[,-1]
  print('======dim in edgeR==========')
  print(dim(ASV_table))
  
  # Read groupings table - suppress messages
  groupings <- suppressMessages(read_tsv(groupings_file, col_names = TRUE))
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
  
  OTU <- phyloseq::otu_table(ASV_table, taxa_are_rows = TRUE)
  sampledata <- phyloseq::sample_data(groupings, errorIfNULL = TRUE)
  phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
  
  test <- phyloseq_to_edgeR(physeq = phylo, group=groupings[,2])
  et <- exactTest(test)
  tt <- topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
  res <- tt@.Data[[1]]

  # Add ASV names to the results
  res$asv_name <- rownames(ASV_table)
  res <- res[, c("asv_name", setdiff(colnames(res), "asv_name"))]
  
  write_tsv(as.data.frame(res), output_file)
  
  list(output_file = output_file)
}

# Define the function to visualize edgeR results
visualize_edgeR <- function(input_file, output_dir, persistent_temp_dir) {
  # Read edgeR results - suppress messages
  edgeR_results <- suppressMessages(read_tsv(input_file))
  
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
  edgeR_results <- edgeR_results %>%
    mutate(neg_log10_pvalue = -log10(PValue),
           point_size = ifelse(FDR < 0.05, 3, 1))
  
  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- edgeR_results %>%
    dplyr::select(asv_name, logFC, PValue) %>%
    mutate(method = "edgeR", log2FoldChange = logFC, pvalue = PValue)
  write_tsv(volcano_plot_data, file.path(persistent_temp_dir, "edger_volcano_plot_data.tsv"))
  
  p1 <- ggplot(edgeR_results, aes(x = logFC, y = neg_log10_pvalue)) +
    geom_point(aes(color = FDR < 0.05, size=point_size)) +  # Coloring significant points
    scale_color_manual(values = c("gray", "#4C3BCF")) +  # Red for significant
    scale_size_continuous(range = c(1, 3), guide = "none") +
    theme_minimal() +
    labs(title = "Volcano plot of logFC vs. -log10(pvalue)",
         x = "Log2 Fold Change",
         y = "-log10(pvalue)")
  
  ggsave(file.path(output_dir, "edgeR_plot1.png"), plot = p1)

  # Visualization 2: Scatter plot of logCPM vs. logFC
  p2 <- ggplot(edgeR_results, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = FDR < 0.05), alpha = 0.5) +
    scale_color_manual(values = c("gray", "#4C3BCF"), name = "Significant") +
    theme_minimal() +
    labs(title = "Scatter plot of logCPM vs. logFC", x = "logCPM", y = "logFC")
  
  ggsave(file.path(output_dir, "edgeR_plot2.png"), plot = p2)
  
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(edgeR_results, aes(x = PValue)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() +
    labs(title = "Histogram of p-value distribution", x = "pvalue", y = "Frequency")

  ggsave(file.path(output_dir, "edgeR_plot3.png"), plot = p3)

  print('finished edgeR visualization')
  
  list(
    plot1 = file.path(output_dir, "edger_plot1.png"),
    plot2 = file.path(output_dir, "edger_plot2.png"),
    plot3 = file.path(output_dir, "edger_plot3.png")
  )
}
