library(ggplot2)
library(dplyr)
library(readr)
library(ggridges)

create_overlap_plots <- function(deseq2_file, aldex2_file, edger_file, maaslin2_file, metagenomeseq_file, output_dir, persistent_temp_dir) {
  # Helper function to safely read and process data
  safe_file_path <- function(...) {
    normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
  }
  safe_read_and_process <- function(file, method_name, fc_col, p_col) {
    tryCatch({
      data <- read_tsv(file, col_types = cols(.default = "c"))
      
      if (!all(c("asv_name", fc_col, p_col) %in% names(data))) {
        stop(paste("Missing required columns in", method_name, "results"))
      }
      
      processed_data <- data %>%
        dplyr::select(asv_name, !!sym(fc_col), !!sym(p_col)) %>%
        mutate(
          log2FoldChange = as.numeric(!!sym(fc_col)),
          pvalue = as.numeric(!!sym(p_col))
        ) %>%
        dplyr::select(asv_name, log2FoldChange, pvalue) %>%
        mutate(method = method_name) %>%
        filter(!is.na(log2FoldChange) & !is.na(pvalue))
      return(processed_data)
    }, error = function(e) {
      message(paste("Error processing", method_name, "results:", e$message))
      if (file.exists(file)) {
        message("File exists. Attempting to read first few lines:")
        tryCatch({
          lines <- readLines(file, n = 5)
          message(paste(lines, collapse = "\n"))
        }, error = function(e) {
          message("Failed to read file contents:", e$message)
        })
      } else {
        message("File does not exist")
      }
      return(NULL)
    })
  }

  # Read and process data
  deseq2_results <- safe_read_and_process(deseq2_file, "DESeq2", "log2FoldChange", "padj")
  aldex2_results <- safe_read_and_process(aldex2_file, "ALDEx2", "effect", "we.eBH")
  edger_results <- safe_read_and_process(edger_file, "edgeR", "logFC", "FDR")
  maaslin2_results <- safe_read_and_process(maaslin2_file, "MaAsLin2", "coef", "qval")
  metagenomeseq_results <- safe_read_and_process(metagenomeseq_file, "metagenomeSeq", "logFC", "pvalues")
  # Combine results, removing any NULL entries
  combined_results <- bind_rows(list(deseq2_results, aldex2_results, edger_results, maaslin2_results, metagenomeseq_results) %>% 
                                  keep(~!is.null(.)))

  if (nrow(combined_results) == 0) {
    stop("No valid data available for plotting. Check the error messages above for details on each method.")
  }

  combined_results <- combined_results %>%
    mutate(neg_log10_pvalue = -log10(pvalue),
           significant = ifelse(pvalue < 0.05, "Significant", "Not Significant"))

  write_tsv(combined_results, safe_file_path(persistent_temp_dir, "overlap_combined_results.tsv"))

  # Overlap volcano plot
  p1 <- ggplot(combined_results, aes(x = log2FoldChange, y = neg_log10_pvalue, color = significant, shape = method)) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c("Significant" = "#0000FF", "Not Significant" = "#808080")) +
    scale_shape_manual(values = c("DESeq2" = 21, "ALDEx2" = 24, "edgeR" = 22, "MaAsLin2" = 25, "metagenomeSeq" = 23)) +
    theme_minimal() +
    labs(title = "Volcano plot overlap", x = "log2FoldChange", y = "-log10(pvalue)") +
    guides(color = guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size = 5))) +
    theme(legend.position = "right")

  ggsave(file.path(output_dir, "overlap_volcano.png"), plot = p1, width = 12, height = 8, dpi = 300)

  # Overlap p-value distribution (ridge plot)
  p2 <- ggplot(combined_results, aes(x = pvalue, y = method, fill = method)) +
    geom_density_ridges(alpha = 0.7, scale = 2) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(title = "P-value distribution overlap", x = "p-value", y = "Method") +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "overlap_pvalue_distribution.png"), plot = p2, width = 12, height = 8, dpi = 300)

  list(
    overlap_volcano = file.path(output_dir, "overlap_volcano.png"),
    overlap_pvalue_distribution = file.path(output_dir, "overlap_pvalue_distribution.png")
  )
}