library(ALDEx2)
library(tidyr)
library(dplyr)
library(readr)
library(e1071)

run_aldex2 <- function(ASV_file, groupings_file, output_file, seed = 1234) {
  set.seed(seed)

  # Read ASV table
  ASV_table <- read_tsv(ASV_file, comment = "", col_names = TRUE, skip = ifelse(grepl("Constructed from biom file", readLines(ASV_file, n=1)), 1, 0))
  ASV_table <- as.data.frame(ASV_table)
  row.names(ASV_table) <- ASV_table[,1]
  ASV_table <- ASV_table[,-1]

  # Read groupings table
  groupings <- read_tsv(groupings_file, col_names = TRUE)
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

  conditions <- groupings[,2]
  unique_groups <- length(unique(conditions))
  
  if (unique_groups == 2) {
    test_type <- "t"
  } else if (unique_groups > 2) {
    test_type <- "kw"
  }

  zero_proportion <- mean(ASV_table == 0)
  feature_skewness <- skewness(colSums(ASV_table))

  if (zero_proportion > 0.5) {
    if (feature_skewness > 1) {
      denom_type <- "iqlr"
    } else {
      denom_type <- "zero"
    }
  } else {
    if (feature_skewness > 1) {
      denom_type <- "iqlr"
    } else {
      denom_type <- "all"
    }
  }

  # Run ALDEx2 analysis
  results <- aldex(reads = ASV_table, 
                   conditions = conditions, 
                   mc.samples = 128, 
                   test = test_type, 
                   effect = TRUE,
                   include.sample.summary = FALSE, 
                   verbose = TRUE, 
                   denom = denom_type)

  # Add ASV names to the results
  results$asv_name <- rownames(ASV_table)
  results <- results[, c("asv_name", setdiff(colnames(results), "asv_name"))]

  # Save results to output file
  write_tsv(as.data.frame(results), output_file)
  
  list(
    plot1 = "./plot1.png",
    plot2 = "./plot2.png",
    plot3 = "./plot3.png"
  )
}

visualize_aldex2 <- function(input_file, output_dir, persistent_temp_dir) {
  # Read ALDEx2 results
  aldex2_results <- read_tsv(input_file)
  
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
  alog <- aldex2_results %>%
    mutate(neg_log10_pvalue = -log10(we.eBH),
           point_size = ifelse(we.eBH < 0.05, 3, 1))

  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- alog %>%
    dplyr::select(asv_name, effect, we.eBH) %>%
    mutate(method = "ALDEx2", log2FoldChange = effect, pvalue = we.eBH)
  write_tsv(volcano_plot_data, file.path(persistent_temp_dir, "aldex2_volcano_plot_data.tsv"))
  
  p1 <- ggplot(alog, aes(x = effect, y = neg_log10_pvalue)) +
    geom_point(aes(color = we.eBH < 0.05, size=point_size)) +  # Coloring significant points
    scale_color_manual(values = c("gray", "#4C3BCF")) +
    scale_size_continuous(range = c(1, 3), guide = "none") +
    theme_minimal() +
    labs(title = "Volcano plot of Effect Size vs. -log10(pvalue)",
         x = "Log2 Fold Change",
         y = "-log10(pvalue)")
  
  ggsave(file.path(output_dir, "aldex2_plot1.png"), plot = p1)


  # Visualization 2: Scatter plot of diff.btw vs. effect
  p2 <- ggplot(alog, aes(x = rab.all, y = effect)) +
    geom_point(aes(color = we.eBH < 0.05, size=point_size), alpha = 0.5) +
    scale_color_manual(values = c("gray", "#4C3BCF"), name = "Significant") +
    scale_size_continuous(range = c(1, 3), guide = "none") +
    theme_minimal() +
    scale_x_log10() +
    labs(title = "MA-like plot", x = "Mean Abundance", y = "Effect Size")
  
  ggsave(file.path(output_dir, "aldex2_plot2.png"), plot = p2)

# Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(aldex2_results, aes(x = we.eBH)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() +
    labs(title = "Histogram of p-value distribution", x = "pvalue", y = "Frequency")

  ggsave(file.path(output_dir, "aldex2_plot3.png"), plot = p3)
  
  list(
    plot1 = file.path(output_dir, "aldex2_plot1.png"),
    plot2 = file.path(output_dir, "aldex2_plot2.png"),
    plot3 = file.path(output_dir, "aldex2_plot3.png")
  )
}
