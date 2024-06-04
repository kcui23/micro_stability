library(ggplot2)
library(dplyr)
library(readr)

create_overlap_plots <- function(deseq2_file, aldex2_file, output_dir) {
  deseq2_results <- read_tsv(deseq2_file) %>%
    mutate(method = "DESeq2")
  
  aldex2_results <- read_tsv(aldex2_file) %>%
    mutate(method = "ALDEx2", log2FoldChange = effect, pvalue = we.eBH)
  
  combined_results <- bind_rows(deseq2_results, aldex2_results)
  combined_results <- combined_results %>%
    mutate(neg_log10_pvalue = -log10(pvalue),
           significant_0.05 = ifelse(pvalue < 0.05, TRUE, FALSE))
  
  # Overlap volcano plot
  p1 <- ggplot(combined_results, aes(x = log2FoldChange, y = neg_log10_pvalue, color = method, shape = significant_0.05)) +
    geom_point(alpha = 0.5, size=4) +
    scale_color_manual(values = c("brown", "purple")) +
    scale_shape_manual(values = c(1, 16)) +
    theme_minimal() +
    labs(title = "Volcano plot overlap", x = "log2FoldChange", y = "-log10(pvalue)")
  
  ggsave(file.path(output_dir, "overlap_volcano.png"), plot = p1, width = 10, height = 6)
  
  # Overlap p-value distribution
  p2 <- ggplot(combined_results, aes(x = pvalue, fill = method)) +
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity") +
    scale_fill_manual(values = c("brown", "purple")) +
    theme_minimal() +
    labs(title = "P-value distribution overlap", x = "p-value", y = "Frequency")
  
  ggsave(file.path(output_dir, "overlap_pvalue_distribution.png"), plot = p2, width = 10, height = 6)
  
  list(
    overlap_volcano = file.path(output_dir, "overlap_volcano.png"),
    overlap_pvalue_distribution = file.path(output_dir, "overlap_pvalue_distribution.png")
  )
}
