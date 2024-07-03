library(plumber)
library(readr)
library(ggplot2)
library(base64enc)
library(jsonlite)
library(dplyr)
library(tools)
library(tidyverse)
library(scales)

# Create a persistent temp directory to save the output files
persistent_temp_dir <- tempdir()

#* @apiTitle ASV Data Processing API

#* Upload dataset and select processing method
#* @post /process
#* @param method The processing method to use (deseq2 or aldex2)
#* @parser json
function(req, method) {
  body <- fromJSON(req$postBody)
  seed <- if (!is.null(body$seed)) as.integer(body$seed) else 1234

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  writeLines(as.character(body$asv), temp_asv_file)
  writeLines(as.character(body$groupings), temp_groupings_file)
  
  if (method == "deseq2") {
    source("DESeq2.R")
    result <- run_deseq2(temp_asv_file, temp_groupings_file, output_file, seed)
    plots <- visualize_deseq2(output_file, output_dir)
  } else if (method == "aldex2") {
    source("Aldex2.R")
    result <- run_aldex2(temp_asv_file, temp_groupings_file, output_file, seed)
    plots <- visualize_aldex2(output_file, output_dir)
  } else if (method == "edger") {
    source("edgeR.R")
    result <- run_edgeR(temp_asv_file, temp_groupings_file, output_file, seed)
    plots <- visualize_edgeR(output_file, output_dir)
  } else if (method == "maaslin2") {
    source("Maaslin2.R")
    # This method is different, it requires a directory to store the output files
    result <- run_maaslin2(temp_asv_file, temp_groupings_file, output_file, output_dir, seed)
    plots <- visualize_maaslin2(output_file, output_dir)
  } else {
    stop("Invalid method selected. Please choose 'deseq2', 'aldex2', 'edger', or 'maaslin2'.")
  }

  # Save the output file to the persistent temp directory
  output_file_name <- paste0(method, "_results.tsv")
  permanent_output_file <- file.path(persistent_temp_dir, output_file_name)
  file.copy(output_file, permanent_output_file, overwrite = TRUE)
  
  list(
    plot1 = base64enc::base64encode(plots$plot1),
    plot2 = base64enc::base64encode(plots$plot2),
    plot3 = base64enc::base64encode(plots$plot3)
  )
}

#* Filter ASV data based on threshold
#* @post /filter
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  
  asv_data <- read_tsv(body$asv, col_types = cols(.default = "c"))
  
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  threshold <- as.numeric(body$threshold)
  filtered_data <- asv_data %>% filter_if(is.numeric, any_vars(. > 50 * threshold))
  
  temp_filtered_file <- tempfile(fileext = ".tsv")
  write_tsv(filtered_data, temp_filtered_file)
  filtered_asv_content <- paste(readLines(temp_filtered_file), collapse = "\n")
  
  list(filteredAsv = filtered_asv_content)
}

#* Quick Explore subset of ASV data
#* @post /quick_explore
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  seed <- if (!is.null(body$seed)) as.integer(body$seed) else 1234

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  deseq2_output_file <- tempfile(fileext = ".tsv")
  aldex2_output_file <- tempfile(fileext = ".tsv")
  edgeR_output_file <- tempfile(fileext = ".tsv")
  maaslin2_output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  writeLines(as.character(body$asv), temp_asv_file)
  writeLines(as.character(body$groupings), temp_groupings_file)
  
  asv_data <- read_tsv(temp_asv_file, col_types = cols(.default = "c"))
  set.seed(seed)
  subset_asv_data <- asv_data %>% sample_frac(0.05)
  write_tsv(subset_asv_data, temp_asv_file)

  source("DESeq2.R")
  deseq2_result <- run_deseq2(temp_asv_file, temp_groupings_file, deseq2_output_file, seed)
  deseq2_plots <- visualize_deseq2(deseq2_output_file, output_dir)
  
  source("Aldex2.R")
  aldex2_result <- run_aldex2(temp_asv_file, temp_groupings_file, aldex2_output_file, seed)
  aldex2_plots <- visualize_aldex2(aldex2_output_file, output_dir)

  source("edgeR.R")
  edgeR_result <- run_edgeR(temp_asv_file, temp_groupings_file, edgeR_output_file, seed)
  edgeR_plots <- visualize_edgeR(edgeR_output_file, output_dir)

  source("Maaslin2.R")
  maaslin2_result <- run_maaslin2(temp_asv_file, temp_groupings_file, maaslin2_output_file, output_dir, seed)
  maaslin2_plots <- visualize_maaslin2(maaslin2_output_file, output_dir)
  
  source("overlap_plots.R")
  overlap_plots <- create_overlap_plots(deseq2_output_file, aldex2_output_file, output_dir)

  list(
    deseq2_plot1 = base64enc::base64encode(deseq2_plots$plot1),
    deseq2_plot2 = base64enc::base64encode(deseq2_plots$plot2),
    deseq2_plot3 = base64enc::base64encode(deseq2_plots$plot3),
    aldex2_plot1 = base64enc::base64encode(aldex2_plots$plot1),
    aldex2_plot2 = base64enc::base64encode(aldex2_plots$plot2),
    aldex2_plot3 = base64enc::base64encode(aldex2_plots$plot3),
    edgeR_plot1 = base64enc::base64encode(edgeR_plots$plot1),
    edgeR_plot2 = base64enc::base64encode(edgeR_plots$plot2),
    edgeR_plot3 = base64enc::base64encode(edgeR_plots$plot3),
    maaslin2_plot1 = base64enc::base64encode(maaslin2_plots$plot1),
    maaslin2_plot2 = base64enc::base64encode(maaslin2_plots$plot2),
    maaslin2_plot3 = base64enc::base64encode(maaslin2_plots$plot3),
    overlap_volcano = base64enc::base64encode(overlap_plots$overlap_volcano),
    overlap_pvalue_distribution = base64enc::base64encode(overlap_plots$overlap_pvalue_distribution)
  )
}

#* Download the result file for the selected method
#* @get /download
#* @param method The processing method to use (deseq2, aldex2, edger, maaslin2)
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res, method) {
  output_file_name <- paste0(method, "_results.tsv")
  permanent_output_file <- file.path(persistent_temp_dir, output_file_name)

  print(head(read_tsv(permanent_output_file)))
  
  if (!file.exists(permanent_output_file)) {
    res$status <- 404
    return(list(error = "File not found"))
  }
  
  file_content <- readLines(permanent_output_file)
  res$body <- paste(file_content, collapse = "\n")
  res$headers$`Content-Disposition` <- paste0("attachment; filename=", method, "_results.tsv")
  return(res)
}

#* @apiTitle ASV Data Processing API

# Add a simple test endpoint
#* @get /test
function() {
  return(list(status = "API is running"))
}

#* Check the status of method files
#* @get /check_method_files
function() {
  cat("Checking method files...\n")  # Debug print
  methods <- c("deseq2", "aldex2", "edger", "maaslin2")
  status <- sapply(methods, function(method) {
    file_path <- file.path(persistent_temp_dir, paste0(method, "_results.tsv"))
    exists <- file.exists(file_path)
    cat(sprintf("%s file exists: %s\n", method, exists))  # Debug print
    exists
  })
  cat("Returning status:", jsonlite::toJSON(as.list(status)), "\n")  # Debug print
  return(as.list(status))
}

#* Download a specific method file
#* @get /download_method_file
#* @param method The method to download the file for
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res, method) {
  valid_methods <- c("deseq2", "aldex2", "edger", "maaslin2")
  if (!(method %in% valid_methods)) {
    res$status <- 400
    return(list(error = "Invalid method specified"))
  }

  file_path <- file.path(persistent_temp_dir, paste0(method, "_results.tsv"))
  if (!file.exists(file_path)) {
    res$status <- 404
    return(list(error = "File not found"))
  }

  file_content <- readLines(file_path)
  res$body <- paste(file_content, collapse = "\n")
  res$headers$`Content-Disposition` <- paste0("attachment; filename=", method, "_results.tsv")
  return(res)
}

#* Generate combined results file
#* @post /generate_combined_results
function(req, res) {
  methods <- c("deseq2", "aldex2", "edger", "maaslin2")
  combined_results <- data.frame(feature = character())

  for (method in methods) {
    file_path <- file.path(persistent_temp_dir, paste0(method, "_results.tsv"))
    if (file.exists(file_path)) {
      results <- read_tsv(file_path)
      
      # Determine the significance column based on the method
      sig_col <- switch(method,
                        "deseq2" = "padj",
                        "aldex2" = "we.eBH",
                        "edger" = "FDR",
                        "maaslin2" = "qval")
      
      # Create a new column for this method's significance
      new_col <- paste0(method, "_significant")
      results[[new_col]] <- results[[sig_col]] < 0.05
      
      # If it's the first method, use all features, otherwise join with existing results
      if (nrow(combined_results) == 0) {
        combined_results <- results[, c("asv_name", new_col)]
      } else {
        combined_results <- full_join(combined_results, results[, c("asv_name", new_col)], by = "asv_name")
      }
    }
  }

  # Replace NA with FALSE for missing values
  combined_results[is.na(combined_results)] <- FALSE

  # Write the combined results to a file
  output_file <- file.path(persistent_temp_dir, "combined_results.tsv")
  write_tsv(combined_results, output_file)

  # Return success message
  list(message = "Combined results file generated successfully")
}

#* Download combined results file
#* @get /download_combined_results
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res) {
  output_file <- file.path(persistent_temp_dir, "combined_results.tsv")
  
  if (!file.exists(output_file)) {
    res$status <- 404
    return(list(error = "Combined results file not found"))
  }
  
  res$body <- paste(readLines(output_file), collapse = "\n")
  res$headers$`Content-Disposition` <- "attachment; filename=combined_results.tsv"
  res
}

#* Get list of ASV names
#* @get /get_asv_list
function() {
  combined_results_file <- file.path(persistent_temp_dir, "combined_results.tsv")
  if (!file.exists(combined_results_file)) {
    return(list(error = "Combined results file not found"))
  }
  
  combined_results <- read_tsv(combined_results_file)
  asv_list <- combined_results$asv_name
  
  return(asv_list)
}

#* Get combined results
#* @get /get_combined_results
function() {
  combined_results_file <- file.path(persistent_temp_dir, "combined_results.tsv")
  if (!file.exists(combined_results_file)) {
    return(list(error = "Combined results file not found"))
  }
  
  combined_results <- read_tsv(combined_results_file)
  
  return(combined_results)
}

#* Generate stability plot
#* @get /generate_stability_plot
#* @serializer contentType list(type="image/png")
function(res) {
  # Read the data
  data <- read_tsv(file.path(persistent_temp_dir, "combined_results.tsv"))
  
  # Data processing
  tools <- names(data)[str_detect(names(data), "_significant$")]
  processed_data <- data %>%
    mutate(total_significant = rowSums(dplyr::select(., all_of(tools)))) %>%
    pivot_longer(cols = all_of(tools), names_to = "tool", values_to = "is_significant") %>%
    group_by(total_significant, tool) %>%
    summarise(tool_significant = sum(is_significant), total = n(), .groups = "drop") %>%
    mutate(
      proportion = tool_significant / total,
      tool = str_remove(tool, "_significant") %>% str_to_title()
    )
  
  # Calculate tool totals
  tool_totals <- processed_data %>%
    group_by(tool) %>%
    summarise(total_significant = sum(tool_significant, na.rm = TRUE)) %>%
    arrange(desc(total_significant))
  
  # Create custom color scale
  max_significant <- max(tool_totals$total_significant)
  color_breaks <- seq(0, max_significant, length.out = 6)[-1]
  color_labels <- scales::comma(color_breaks)
  custom_color_scale <- scale_fill_gradientn(
    colors = colorRampPalette(c("#bac6df", "#1b4197"))(5),
    breaks = color_breaks,
    labels = color_labels,
    limits = c(0, max_significant),
    name = "Significant ASVs"
  )
  
  # Create the plot
  p <- ggplot(processed_data, aes(x = factor(total_significant), y = proportion, fill = tool_significant)) +
    geom_col(position = "dodge") +
    facet_wrap(~ factor(tool, levels = tool_totals$tool), ncol = 1, strip.position = "left") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0, 0.5), limits = c(0, 1)) +
    custom_color_scale +
    labs(
      x = "Number of tools finding ASV significant",
      y = "Proportion of ASVs found significant",
      fill = "Significant ASVs"
    ) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.placement = "outside",
      legend.position = "right"
    )
  
  # Save the plot to a file
  output_file <- tempfile(fileext = ".png")
  ggsave(output_file, p, width = 10, height = 8)
  
  # Return the plot as an image
  readBin(output_file, "raw", n = file.info(output_file)$size)
}