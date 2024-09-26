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
persistent_temp_dir <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)

# Global variable to store active WebSocket connections
active_ws_connections <- new.env()

# Helper function to ensure consistent file path handling
safe_file_path <- function(...) {
  normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
}

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
  
  writeLines(as.character(body$asv), temp_asv_file, sep = "\n")
  writeLines(as.character(body$groupings), temp_groupings_file, sep = "\n")
  
  if (tolower(method) == "deseq2") {
    source(safe_file_path("DESeq2.R"))
    result <- run_deseq2(temp_asv_file, temp_groupings_file, output_file, seed)
    plots <- visualize_deseq2(output_file, output_dir)
  } else if (tolower(method) == "aldex2") {
    source(safe_file_path("Aldex2.R"))
    result <- run_aldex2(temp_asv_file, temp_groupings_file, output_file, seed)
    plots <- visualize_aldex2(output_file, output_dir)
  } else if (tolower(method) == "edger") {
    source(safe_file_path("edgeR.R"))
    result <- run_edgeR(temp_asv_file, temp_groupings_file, output_file, seed)
    plots <- visualize_edgeR(output_file, output_dir)
  } else if (tolower(method) == "maaslin2") {
    source(safe_file_path("Maaslin2.R"))
    # This method is different, it requires a directory to store the output files
    result <- run_maaslin2(temp_asv_file, temp_groupings_file, output_file, output_dir, seed)
    plots <- visualize_maaslin2(output_file, output_dir)
  } else {
    stop("Invalid method selected. Please choose 'deseq2', 'aldex2', 'edger', or 'maaslin2'.")
  }

  # Save the output file to the persistent temp directory
  output_file_name <- paste0(tolower(method), "_results.tsv")
  permanent_output_file <- safe_file_path(persistent_temp_dir, output_file_name)
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
  
  writeLines(as.character(body$asv), temp_asv_file, sep = "\n")
  writeLines(as.character(body$groupings), temp_groupings_file, sep = "\n")
  
  asv_data <- read_tsv(temp_asv_file, col_types = cols(.default = "c"))
  set.seed(seed)
  subset_asv_data <- asv_data %>% sample_frac(0.05)
  write_tsv(subset_asv_data, temp_asv_file)

  tryCatch({
    source(safe_file_path("DESeq2.R"))
    deseq2_result <- run_deseq2(temp_asv_file, temp_groupings_file, deseq2_output_file, seed)
    deseq2_plots <- visualize_deseq2(deseq2_output_file, output_dir)
    
    source(safe_file_path("Aldex2.R"))
    aldex2_result <- run_aldex2(temp_asv_file, temp_groupings_file, aldex2_output_file, seed)
    aldex2_plots <- visualize_aldex2(aldex2_output_file, output_dir)

    source(safe_file_path("edgeR.R"))
    edgeR_result <- run_edgeR(temp_asv_file, temp_groupings_file, edgeR_output_file, seed)
    edgeR_plots <- visualize_edgeR(edgeR_output_file, output_dir)

    source(safe_file_path("Maaslin2.R"))
    maaslin2_result <- run_maaslin2(temp_asv_file, temp_groupings_file, maaslin2_output_file, output_dir, seed)
    maaslin2_plots <- visualize_maaslin2(maaslin2_output_file, output_dir)
    
    source(safe_file_path("overlap_plots.R"))
    overlap_plots <- create_overlap_plots(deseq2_output_file, aldex2_output_file, edgeR_output_file, maaslin2_output_file, output_dir, persistent_temp_dir)

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
  }, error = function(e) {
    message("Error in quick_explore: ", e$message)
    list(error = paste("Error in quick_explore:", e$message))
  })
}

#* Download the result file for the selected method
#* @get /download
#* @param method The processing method to use (deseq2, aldex2, edger, maaslin2)
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res, method) {
  output_file_name <- paste0(tolower(method), "_results.tsv")
  permanent_output_file <- safe_file_path(persistent_temp_dir, output_file_name)

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
  methods <- c("deseq2", "aldex2", "edger", "maaslin2")
  status <- sapply(methods, function(method) {
    file_path <- safe_file_path(persistent_temp_dir, paste0(method, "_results.tsv"))
    file.exists(file_path)
  })
  return(as.list(status))
}

#* Download a specific method file
#* @get /download_method_file
#* @param method The method to download the file for
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res, method) {
  valid_methods <- c("deseq2", "aldex2", "edger", "maaslin2")
  if (!(tolower(method) %in% valid_methods)) {
    res$status <- 400
    return(list(error = "Invalid method specified"))
  }

  file_path <- safe_file_path(persistent_temp_dir, paste0(tolower(method), "_results.tsv"))
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
    file_path <- safe_file_path(persistent_temp_dir, paste0(method, "_results.tsv"))
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
  output_file <- safe_file_path(persistent_temp_dir, "combined_results.tsv")
  write_tsv(combined_results, output_file)

  # Return success message
  list(message = "Combined results file generated successfully")
}

#* Download combined results file
#* @get /download_combined_results
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res) {
  output_file <- safe_file_path(persistent_temp_dir, "combined_results.tsv")
  
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
  combined_results_file <- safe_file_path(persistent_temp_dir, "combined_results.tsv")
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
  combined_results_file <- safe_file_path(persistent_temp_dir, "combined_results.tsv")
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
  data <- read_tsv(safe_file_path(persistent_temp_dir, "combined_results.tsv"))
  
  # Data processing
  tools <- names(data)[stringr::str_detect(names(data), "_significant$")]
  processed_data <- data %>%
    mutate(total_significant = rowSums(dplyr::select(., all_of(tools)))) %>%
    pivot_longer(cols = all_of(tools), names_to = "tool", values_to = "is_significant") %>%
    group_by(total_significant, tool) %>%
    summarise(tool_significant = sum(is_significant), total = n(), .groups = "drop") %>%
    mutate(
      proportion = tool_significant / total,
      tool = stringr::str_remove(tool, "_significant") %>% stringr::str_to_title()
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
    colors = colorRampPalette(c("#bac6df", "#1b4197"))(100),
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

#* @websocket /ws
function(ws) {
  # Generate a unique ID for this connection
  ws_id <- uuid::UUIDgenerate()
  
  # Store the WebSocket connection
  active_ws_connections[[ws_id]] <- ws
  
  # Clean up the connection when it's closed
  ws$onClose(function() {
    rm(list = ws_id, envir = active_ws_connections)
  })
}

# Function to send a message to a specific WebSocket
send_ws_message <- function(ws_id, message) {
  if (exists(ws_id, envir = active_ws_connections)) {
    ws <- active_ws_connections[[ws_id]]
    ws$send(jsonlite::toJSON(message))
  }
}

#* Perform shuffled analysis
#* @post /shuffled_analysis
#* @param iterations The number of iterations to perform
#* @parser json
function(req, iterations = 10, ws_id) {
  body <- fromJSON(req$postBody)
  seed <- if (!is.null(body$seed)) as.integer(body$seed) else 1234
  set.seed(seed)

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  writeLines(as.character(body$asv), temp_asv_file, sep = "\n")
  writeLines(as.character(body$groupings), temp_groupings_file, sep = "\n")

  asv_data <- read_tsv(temp_asv_file)
  groupings_data <- read_tsv(temp_groupings_file)

  methods <- c("deseq2", "aldex2", "edger", "maaslin2")
  results <- list()

  for (i in 1:iterations) {
    set.seed(seed + i)
    # Send progress update to the client
    send_ws_message(ws_id, list(progress = i, total = iterations))

    # Shuffle the groupings using the second column
    shuffled_groupings <- groupings_data %>%
      mutate(!!names(groupings_data)[2] := sample(!!sym(names(groupings_data)[2])))

    print("=====================Shuffled groupings:=====================")
    print(shuffled_groupings)

    
    # Write shuffled groupings to a temporary file
    temp_shuffled_groupings <- tempfile(fileext = ".tsv")
    write_tsv(shuffled_groupings, temp_shuffled_groupings)

    for (method in methods) {
      output_file <- tempfile(fileext = ".tsv")
      
      if (method == "deseq2") {
        source("DESeq2.R")
        run_deseq2(temp_asv_file, temp_shuffled_groupings, output_file, seed)
      } else if (method == "aldex2") {
        source("Aldex2.R")
        run_aldex2(temp_asv_file, temp_shuffled_groupings, output_file, seed)
      } else if (method == "edger") {
        source("edgeR.R")
        run_edgeR(temp_asv_file, temp_shuffled_groupings, output_file, seed)
      } else if (method == "maaslin2") {
        source("Maaslin2.R")
        run_maaslin2(temp_asv_file, temp_shuffled_groupings, output_file, tempdir(), seed)
      }

      # Read results and store significant ASVs
      result <- read_tsv(output_file)
      significant_asvs <- get_significant_asvs(result, method)
      
      if (is.null(results[[method]])) {
        results[[method]] <- list()
      }
      results[[method]][[i]] <- significant_asvs
    }
  }

  # Process results
  processed_results <- process_shuffled_results(results, iterations)

  # Generate plot
  plot <- generate_stability_plot(processed_results)

  list(
    results = processed_results,
    plot = base64enc::base64encode(plot)
  )
}

# Helper function to get significant ASVs based on method
get_significant_asvs <- function(result, method) {
  if (method == "deseq2") {
    return(result$asv_name[result$padj < 0.05])
  } else if (method == "aldex2") {
    return(result$asv_name[result$we.eBH < 0.05])
  } else if (method == "edger") {
    return(result$asv_name[result$FDR < 0.05])
  } else if (method == "maaslin2") {
    return(result$asv_name[result$qval < 0.05])
  }
}

# Helper function to process shuffled results
process_shuffled_results <- function(results, iterations) {
  processed <- list()
  for (method in names(results)) {
    all_asvs <- unique(unlist(results[[method]]))
    if (length(all_asvs) > 0) {
      asv_counts <- sapply(all_asvs, function(asv) {
        sum(sapply(results[[method]], function(iter) asv %in% iter))
      })
      processed[[method]] <- data.frame(
        asv = names(asv_counts),
        count = as.numeric(asv_counts),
        percentage = as.numeric(asv_counts) / iterations * 100
      )
    } else {
      # If no ASVs found, create an empty data frame with the correct structure
      processed[[method]] <- data.frame(
        asv = character(0),
        count = numeric(0),
        percentage = numeric(0)
      )
    }
  }
  return(processed)
}

# Helper function to generate stability plot
generate_stability_plot <- function(processed_results) {
  plot_data <- bind_rows(lapply(names(processed_results), function(method) {
    if (nrow(processed_results[[method]]) > 0) {
      data.frame(
        method = method,
        percentage = processed_results[[method]]$percentage,
        count = processed_results[[method]]$count
      )
    } else {
      # If no data for this method, create a single row with count 0
      data.frame(
        method = method,
        percentage = 0,
        count = 0
      )
    }
  }))

  # Check if we have any data to plot
  if (nrow(plot_data) == 0) {
    # If no data at all, create a blank plot with a message
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant ASVs found in any method") +
      theme_void()
  } else {
    # Create the plot as before
    p <- ggplot(plot_data, aes(x = percentage, y = count)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ method, scales = "free_y") +
      labs(x = "Percentage of iterations with significant ASV", 
           y = "Number of ASVs", 
           title = "Stability of significant ASVs across shuffled analyses") +
      theme_minimal()
  }

  # Save plot to a temporary file
  temp_plot <- tempfile(fileext = ".png")
  ggsave(temp_plot, p, width = 10, height = 8)
  
  return(temp_plot)
}

#* Download a specific plot image
#* @get /download_image
#* @param method The method of the plot (e.g., deseq2, aldex2, edger, maaslin2, overlap)
#* @param plot The specific plot to download (e.g., plot1, plot2, plot3, volcano, pvalue_distribution)
#* @serializer contentType list(type="image/png")
function(req, res, method, plot) {
  tryCatch({
    valid_methods <- c("deseq2", "aldex2", "edger", "maaslin2", "overlap")
    if (!(tolower(method) %in% valid_methods)) {
      res$status <- 400
      return(list(error = "Invalid method specified"))
    }

    # Construct the file path
    if (method == "overlap") {
      file_name <- paste0("overlap_", plot, ".png")
    } else {
      file_name <- paste0(method, "_", plot, ".png")
    }
    file_path <- safe_file_path(persistent_temp_dir, file_name)

    if (!file.exists(file_path)) {
      res$status <- 404
      return(list(error = paste("Image not found:", file_path)))
    }

    # Read the image file
    image_data <- readBin(file_path, "raw", n = file.info(file_path)$size)

    # Set the appropriate headers
    res$setHeader("Content-Type", "image/png")
    res$setHeader("Content-Disposition", paste0('attachment; filename="', file_name, '"'))

    # Return the image data
    image_data
  }, error = function(e) {
    res$status <- 500
    list(error = paste("Server error:", e$message))
  })
}

#* Provide the combined results TSV file
#* @get /overlap_combined_results_tsv
#* @serializer contentType list(type="text/tab-separated-values")
function(req, res) {
  file_path <- safe_file_path(persistent_temp_dir, "overlap_combined_results.tsv")
  
  if (!file.exists(file_path)) {
    res$status <- 404
    return(list(error = "Overlap combined results file not found"))
  }

  file_content <- readLines(file_path)
  print("overlap_combined_results_tsv")
  print(file_content)
  
  res$body <- file_content
  res$headers$`Content-Disposition` <- "attachment; filename=overlap_combined_results.tsv"
  res
}

#* Download selected points as a .tsv file
#* @post /download_selected_points
#* @serializer contentType list(type = "text/tab-separated-values")
function(req, res) {
  body <- fromJSON(req$postBody)
  
  temp_file <- tempfile(fileext = ".tsv")
  write_tsv(as.data.frame(body$points), temp_file)
  
  res$body <- paste(readLines(temp_file), collapse = "\n")
  res$headers$`Content-Disposition` <- "attachment; filename=selected_points.tsv"
  res
}

#* Update leaf data with random data points
#* @post /update_leaf_data
#* @serializer json
#* @response 200 list(message="JSON file successfully updated")
#* @response 500 list(error="Error updating JSON file") 
function(req, res) {
  json_file_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json"
  
  tryCatch({
    json_data <- fromJSON(json_file_path)
    
    for (leaf in names(json_data)) {
      json_data[[leaf]]$data_point <- generate_random_pair()
    }
    
    write_json(json_data, json_file_path, pretty = TRUE, auto_unbox = TRUE)
    
    res$status <- 200
    return(list(message = "JSON file successfully updated"))
  }, error = function(e) {
    
    res$status <- 500
    return(list(error = paste("Error updating the JSON file:", e$message)))
  })
}

generate_random_pair <- function() {
    return(runif(2, min = 0, max = 10))
  }

#* Calculate stability metric
#* @post /calculate_stability_metric
#* @serializer json
#* @response 200 list(message="Stability metric calculated successfully")
#* @response 500 list(error="Error calculating stability metric")
function(req, res) {
  tryCatch({
    body <- fromJSON(req$postBody)
    destroy <- body$destroy
    print("=====destroy=====")
    print(destroy)
    leaf_json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json"
    tree_json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json"
    leaf_data <- fromJSON(leaf_json_path, simplifyVector = TRUE)
    tree_data <- fromJSON(tree_json_path, simplifyVector = FALSE)
    if (destroy) {
      print("=====destroy is true=====")
      for (leaf in names(leaf_data)) {
        leaf_data[[leaf]]$data_point <- c(0, 0)
      }
      write_json(leaf_data, leaf_json_path, pretty = TRUE, auto_unbox = TRUE)
      res$status <- 200
      return(list(message = "Bye."))
    }
    asv <- body$asv
    groupings <- body$groupings
    method <- body$method
    missing_methods <- body$missing_methods

    for (leaf in names(leaf_data)) {
      tmp_path <- find_path_from_id(leaf, tree_data)
      if (tmp_path[5] %in% missing_methods) {
        leaf_data[[leaf]]$data_point <- c(0, 0)
      }
    }

    stability_metric <- list()
    for (leaf in names(leaf_data)) {
      path <- find_path_from_id(leaf, tree_data)
      if (path[5] == method) {
        leaf_data <- calculate_stability_metric(asv, groupings, path, leaf, leaf_data, leaf_json_path, test = TRUE)
      }
    }

    write_json(leaf_data, leaf_json_path, pretty = TRUE, auto_unbox = TRUE)

    res$status <- 200
    return(list(message = "Stability metric calculated successfully", stability_metric = stability_metric))
  }, error = function(e) {
    res$status <- 500
    return(list(error = paste("Error calculating stability metric:", e$message)))
  })
}

find_path_from_id <- function(id, tree = tree_data) {
    recurse <- function(node, target_id, current_path) {
        new_path <- c(current_path, node$name)
        
        if (!is.null(node$id) && node$id == target_id) {
        return(new_path)
        }
        
        if (!is.null(node$children)) {
        for (child in node$children) {
            result <- recurse(child, target_id, new_path)
            if (!is.null(result)) {
            return(result)
            }
        }
        }
        
        return(NULL)
    }
    path <- recurse(tree, id, c())
    if (is.null(path)) {
        stop(paste("ID", id, "not found in the tree data."))
    }
    return(path)
}

calculate_stability_metric <- function(asv, groupings, path, leaf, json_data = NULL, json_file_path = NULL, test = FALSE) {
  if (test) {
    Sys.sleep(0.1)
    new_data_point <- generate_random_pair()
    json_data[[leaf]]$data_point <- new_data_point
    return(json_data)
  }
  else {
    asv_df <- read_tsv(asv)
    groupings_df <- read_tsv(groupings)
    shuffled_groupings <- groupings_df %>%
      mutate(!!names(groupings_df)[2] := sample(!!sym(names(groupings_df)[2])))

    temp_asv_file <- tempfile(fileext = ".tsv")
    temp_shuffled_groupings <- tempfile(fileext = ".tsv")
    write_tsv(asv_df, temp_asv_file)
    write_tsv(shuffled_groupings, temp_shuffled_groupings)

    output_file <- tempfile(fileext = ".tsv")

    if (method == "deseq2") {
      source(safe_file_path("DESeq2.R"))
      run_deseq2(temp_asv_file, temp_shuffled_groupings, output_file, seed = 1234)
    } else if (method == "aldex2") {
      source(safe_file_path("Aldex2.R"))
      run_aldex2(temp_asv_file, temp_shuffled_groupings, output_file, seed = 1234)
    } else if (method == "edger") {
      source(safe_file_path("edgeR.R"))
      run_edgeR(temp_asv_file, temp_shuffled_groupings, output_file, seed = 1234)
    } else if (method == "maaslin2") {
      source(safe_file_path("Maaslin2.R"))
      run_maaslin2(temp_asv_file, temp_shuffled_groupings, output_file, output_dir = tempfile(), seed = 1234)
    } else {
      stop("Invalid method specified")
    }

    result <- read_tsv(output_file)
    significant_asvs <- get_significant_asvs(result, method)

    stability_metric <- list(
      significant_asvs = significant_asvs,
      count = length(significant_asvs)
    )

    return(stability_metric)
  }
}