suppressPackageStartupMessages({
  library(plumber)
  library(readr)
  library(ggplot2)
  library(base64enc)
  library(jsonlite)
  library(dplyr)
  library(tools)
  library(tidyverse)
  library(scales)
  library(impute)
  library(metagenomeSeq)
  library(edgeR)
  library(uwot)
})

# Create a persistent temp directory to save the output files
persistent_temp_dir <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
# Directory to store uploaded files
uploaded_files_dir <- file.path(persistent_temp_dir, "uploaded_files")
dir.create(uploaded_files_dir, showWarnings = FALSE, recursive = TRUE)
code_dir <- file.path(persistent_temp_dir, "code")
dir.create(code_dir, showWarnings = FALSE, recursive = TRUE)


asv_file_path <- file.path(uploaded_files_dir, "asv_data.tsv")
groupings_file_path <- file.path(uploaded_files_dir, "groupings_data.tsv")

#* Generate all R code files
#* @post /generate_all_r_code
#* @serializer json
function(req, res) {
  message(("Generating all R code files"))
  # check if methods_sig_vectors.rds already exists, if so, delete it
  if (file.exists(file.path(persistent_temp_dir, "methods_sig_vectors.rds"))) {
    file.remove(file.path(persistent_temp_dir, "methods_sig_vectors.rds"))
  }
  tryCatch({
    source(safe_file_path("generate_r_code.R"))
    
    params <- list(
      "Low Abundance Filtering" = 5,
      "Prevalence Filtering" = 10,
      "Variance Filtering" = 0,
      "Pseudocount Addition" = 1,
      "knn" = 5,
      "knn_bound" = 50,
      "Random Seed" = 1234
    )
    tryCatch({
      generate_all_r_code(params, code_dir)
      message("Number of files in code_dir: ", crayon::green$bold(length(list.files(code_dir))))
      # zip the code_dir with quiet flag
      zip_file <- file.path(persistent_temp_dir, "code.zip")
      files_to_zip <- list.files(code_dir, full.names = TRUE)
      utils::zip(zip_file, files_to_zip, flags="-q")
      
      message("Code directory zipped to: ", zip_file)
      local_path <- file.path("~/Downloads", basename(zip_file))
      file.copy(zip_file, local_path, overwrite = TRUE)
      message(paste("Downloaded file to:", local_path))
    }, error = function(e) {
      print(paste("Error generating R code:", e$message))
    })

    # generate the path data
    leaf_json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json"
    tree_json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json"
    leaf_data <- fromJSON(leaf_json_path, simplifyVector = TRUE)
    tree_data <- fromJSON(tree_json_path, simplifyVector = FALSE)

    path_df <- data.frame(
      leaf_id = character(),
      path = I(list()),
      stringsAsFactors = FALSE
    )

    for (leaf in names(leaf_data)) {
      tmp_path <- find_path_from_id(leaf, tree_data)
      path_df <- rbind(path_df, data.frame(
        leaf_id = leaf,
        path = I(list(tmp_path)),
        stringsAsFactors = FALSE
      ))
    }

    message("path data.frame saved.")
    saveRDS(path_df, file.path(persistent_temp_dir, "path_df.rds"))

    
    res$status <- 200
    return(list(message = "All R code files generated successfully", directory = code_dir))
  }, error = function(e) {
    res$status <- 500
    return(list(error = paste("Error generating R code files:", e$message)))
  })
}

#* Download R code file
#* @post /download_code
#* @param type The type of code to download
#* @serializer contentType list(type="text/plain")
function(req, res) {
  body <- fromJSON(req$postBody)
  selectedOperations <- body$selectedOperations
  params <- body$params
  source(safe_file_path("generate_r_code.R"))
  code <- generate_r_code(selectedOperations, params)

  res$body <- code
  res$headers$`Content-Disposition` <- paste0("attachment; filename=r_code.R")
  res
}

#* Store uploaded files
#* @post /store_files
#* @param asv The ASV file content
#* @param groupings The groupings file content
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  
  writeLines(as.character(body$asv), asv_file_path, sep = "\n")
  writeLines(as.character(body$groupings), groupings_file_path, sep = "\n")

  list(
    message = "Files stored successfully",
    asv_file = asv_file_path,
    groupings_file = groupings_file_path
  )
}

# Global variable to store active WebSocket connections
active_ws_connections <- new.env()

# Helper function to ensure consistent file path handling
safe_file_path <- function(...) {
  normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
}

#* Preview ASV and grouping data
#* @get /preview_data
function() {
  asv_data <- read_tsv(asv_file_path, show_col_types = FALSE)
  groupings_data <- read_tsv(groupings_file_path, show_col_types = FALSE)
  
  asv_preview <- asv_data %>%
    dplyr::select(1:min(6, ncol(.))) %>%
    dplyr::slice(1:min(5, nrow(.)))
  groupings_preview <- groupings_data %>%
    dplyr::slice(1:min(5, nrow(.)))

  asv_preview <- asv_preview %>%
    mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.))))
  groupings_preview <- groupings_preview %>%
    mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.))))
  
  response <- list(
    asv_preview = asv_preview,
    asv_dimensions = dim(asv_data),
    groupings_preview = groupings_preview,
    groupings_dimensions = dim(groupings_data)
  )
  
  return(response)
}

#* Get overlap combined results
#* @get /get_overlap_combined_results
#* @serializer contentType list(type="text/tab-separated-values")
function(req, res, specific_interact, selectedMethod) {
  if (specific_interact) {
    file_path <- safe_file_path(persistent_temp_dir, paste0(selectedMethod, "_volcano_plot_data.tsv"))
  } else {
    file_path <- safe_file_path(persistent_temp_dir, "overlap_combined_results.tsv")
  }
  if (!file.exists(file_path)) {
    res$status <- 404
    return(list(error = "Overlap combined results file not found"))
  }
  file_content <- readLines(file_path)
  
  res$body <- paste(file_content, collapse = "\n")
  res$headers$`Content-Disposition` <- "attachment; filename=overlap_combined_results.tsv"
  res
}

#* Upload dataset and select processing method
#* @post /process
#* @param method The processing method to use (deseq2 or aldex2)
#* @parser json
function(req, method) {
  body <- fromJSON(req$postBody)
  seed <- if (!is.null(body$seed)) as.integer(body$seed) else 1234

  output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  if (tolower(method) == "deseq2") {
    source(safe_file_path("DESeq2.R"))
    result <- run_deseq2(asv_file_path, groupings_file_path, output_file, seed)
    plots <- visualize_deseq2(output_file, output_dir, persistent_temp_dir)
  } else if (tolower(method) == "aldex2") {
    source(safe_file_path("Aldex2.R"))
    result <- run_aldex2(asv_file_path, groupings_file_path, output_file, seed)
    plots <- visualize_aldex2(output_file, output_dir, persistent_temp_dir)
  } else if (tolower(method) == "edger") {
    source(safe_file_path("edgeR.R"))
    result <- run_edgeR(asv_file_path, groupings_file_path, output_file, seed)
    plots <- visualize_edgeR(output_file, output_dir, persistent_temp_dir)
  } else if (tolower(method) == "maaslin2") {
    source(safe_file_path("Maaslin2.R"))
    result <- run_maaslin2(asv_file_path, groupings_file_path, output_file, output_dir, seed)
    plots <- visualize_maaslin2(output_file, output_dir, persistent_temp_dir)
  } else if (tolower(method) == "metagenomeseq") {
    source(safe_file_path("metagenomeSeq.R"))
    result <- run_metagenomeseq(asv_file_path, groupings_file_path, output_file, seed)
    plots <- visualize_metagenomeseq(output_file, output_dir, persistent_temp_dir)
  } else {
    stop("Invalid method selected. Please choose 'deseq2', 'aldex2', 'edger', 'maaslin2', or 'metagenomeseq'.")
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
  
  asv_data <- read_tsv(asv_file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  filter_method <- body$filter_method
  threshold <- as.numeric(body$threshold)
  
  low_abundance_filter <- function(data, threshold) {
    data %>% mutate_if(is.numeric, ~ifelse(. < threshold, 0, .))
  }

  prevalence_filter <- function(data, threshold) {
    total_samples <- ncol(data) - 1
    data %>%
      mutate(prevalence = rowSums(select(., -1) > 0) / total_samples) %>%
      filter(prevalence >= threshold) %>%
      dplyr::select(-prevalence)
  }

  variance_filter <- function(data, threshold) {
    data %>%
      mutate(variance = apply(select(., -1), 1, var)) %>%
      filter(variance > threshold) %>%
      dplyr::select(-variance)
  }

  # Use switch to call the appropriate function
  filtered_data <- switch(filter_method,
    "Low Abundance Filtering" = low_abundance_filter(asv_data, threshold),
    "Prevalence Filtering" = prevalence_filter(asv_data, threshold),
    "Variance Filtering" = variance_filter(asv_data, threshold),
    stop("Invalid filter method")
  )
  
  write_tsv(filtered_data, asv_file_path)
  # add some log messages about the method and threshold used
  message("Filtered ASV data with ", filter_method, " with a threshold of ", threshold)

  list(success = TRUE)
}

#* Zero handling for ASV data
#* @post /zero_handling
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  
  asv_data <- read_tsv(asv_file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  zero_method <- body$zero_handling_method
  param_value <- as.numeric(body$param_value)
  
  pseudocount_addition <- function(data, param_value) {
    data %>% mutate(across(-1, ~. + param_value))
  }

  knn_imputation <- function(data, param_value) {
    message("KNN Imputation with k = ", param_value[1], " and bound = ", param_value[2])
    tmp_data <- data %>% select(-asv)
    zero_percentages <- rowMeans(tmp_data == 0) * 100
    asvs_to_impute <- zero_percentages < param_value[2]
    
    if (any(asvs_to_impute)) {
      data_to_impute <- tmp_data[asvs_to_impute, ]
      data_to_impute <- data_to_impute %>%
        mutate(across(where(is.numeric), ~ifelse(. == 0, NA, .)))
      tryCatch({
        imputed_data <- impute.knn(as.matrix(data_to_impute), k = param_value[1], rowmax = 0.5)$data
        tmp_data[asvs_to_impute, ] <- imputed_data
      }, error = function(e) {
        warning("Error in k-NN imputation: ", e$message)
        # Fallback to mean imputation for rows with high missing values
        high_missing <- rowMeans(is.na(data_to_impute)) > 0.5
        data_to_impute[high_missing, ] <- apply(data_to_impute[high_missing, ], 1, function(x) {
          x[is.na(x)] <- mean(x, na.rm = TRUE)
          return(x)
        })
        tmp_data[asvs_to_impute, ] <- data_to_impute
      })
    }
    tmp_data[is.na(tmp_data)] <- 0
    tmp_data <- tmp_data %>% 
      mutate(asv = data$asv) %>%
      dplyr::select(asv, everything())
    
    return(tmp_data)
  }

  # Use switch to call the appropriate function
  zero_handled_data <- switch(zero_method,
    "Pseudocount Addition" = pseudocount_addition(asv_data, param_value),
    "k-NN Imputation" = knn_imputation(asv_data, param_value),
    stop("Invalid zero handling method")
  )
  
  write_tsv(zero_handled_data, asv_file_path)
  message("Applied zero handling with ", zero_method, " with a parameter of ", param_value)

  list(success = TRUE)
}

#* Generate a simple ggplot2 image
#* @get /zero_distribution_plot
#* @serializer contentType list(type="image/png")
function() {
  asv_data <- read_tsv(asv_file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  zero_percentages <- asv_data %>%
    rowwise() %>%
    mutate(ZeroPercentage = mean(c_across(-1) == 0) * 100) %>%
    ungroup() %>%
    dplyr::select(1, ZeroPercentage)
  
  p <- ggplot(zero_percentages, aes(x = ZeroPercentage)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(x = "Percentage of Zeros",
         y = "Count") +
    theme_minimal()

  temp_file <- tempfile(fileext = ".png")
  ggsave(temp_file, p, width = 6, height = 2, dpi = 300)
  readBin(temp_file, "raw", n = file.info(temp_file)$size)
}

#* Normalize ASV data
#* @post /normalization
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  
  asv_data <- read_tsv(asv_file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  norm_method <- body$norm_method

  tss_normalization <- function(data) {
    data <- data %>% mutate(across(-1, ~./sum(.)))
    return(data)
  }
  
  css_normalization <- function(data) {
    tmp_data <- data %>% select(-asv)
    MR_obj <- newMRexperiment(tmp_data)
    MR_obj_normalized <- cumNorm(MR_obj)
    normalized_counts <- MRcounts(MR_obj_normalized, norm = TRUE)
    normalized_counts <- as.data.frame(normalized_counts) %>% 
      mutate(asv = data$asv) %>% 
      dplyr::select(asv, everything())
    return(normalized_counts)
  }
  
  tmm_normalization <- function(data) {
    tmp_data <- data %>% select(-asv)
    dge <- DGEList(counts = tmp_data)
    dge <- calcNormFactors(dge, method = "TMM")
    tmm_normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
    tmm_normalized_counts <- as.data.frame(tmm_normalized_counts) %>% 
      mutate(asv = data$asv) %>% 
      dplyr::select(asv, everything())
    return(tmm_normalized_counts)
  }
  
  clr_normalization <- function(data) {
    tmp_data <- data %>% select(-asv)
    count_matrix <- tmp_data + 1
    geom_mean <- function(x) {
      exp(mean(log(x)))
    }
    clr_manual <- t(apply(count_matrix, 1, function(row) {
      log(row) - log(geom_mean(row))
    }))
    clr_manual <- as.data.frame(clr_manual) %>% 
      mutate(asv = data$asv) %>% 
      dplyr::select(asv, everything())
    return(clr_manual)
  }

  normalized_data <- switch(norm_method,
    "TSS" = tss_normalization(asv_data),
    "CSS" = css_normalization(asv_data),
    "TMM" = tmm_normalization(asv_data),
    "CLR" = clr_normalization(asv_data),
    stop("Invalid normalization method")
  )
  
  write_tsv(normalized_data, asv_file_path)
  message("Applied normalization with ", norm_method)

  list(success = TRUE)
}

#* Normalize ASV data
#* @post /transformation
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  
  asv_data <- read_tsv(asv_file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  trans_method <- body$trans_method

  log_transformation <- function(data) {
    data <- data %>% mutate(across(-1, ~log(.)))
    return(data)
  }
  
  logit_transformation <- function(data) {
    logit <- function(p) {
      log(p / (1 - p))
    }
    data <- data %>% mutate(across(-1, ~logit(.)))
    return(data)
  }
  
  ast_transformation <- function(data) {
    data <- data %>% mutate(across(-1, ~asin(sqrt(.))))
    return(data)
  }

  transformed_data <- switch(trans_method,
    "Log" = log_transformation(asv_data),
    "Logit" = logit_transformation(asv_data),
    "AST" = ast_transformation(asv_data),
    stop("Invalid transformation method")
  )
  
  write_tsv(transformed_data, asv_file_path)
  message("Applied transformation with ", trans_method)

  list(success = TRUE)
}

#* Quick Explore subset of ASV data
#* @post /quick_explore
#* @parser json
function(req) {
  body <- fromJSON(req$postBody)
  seed <- if (!is.null(body$seed)) as.integer(body$seed) else 1234

  temp_asv_file <- tempfile(fileext = ".tsv")
  deseq2_output_file <- tempfile(fileext = ".tsv")
  aldex2_output_file <- tempfile(fileext = ".tsv")
  edgeR_output_file <- tempfile(fileext = ".tsv")
  maaslin2_output_file <- tempfile(fileext = ".tsv")
  metagenomeseq_output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  asv_data <- read_tsv(asv_file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  
  set.seed(seed)
  subset_asv_data <- asv_data %>% sample_frac(0.05)
  write_tsv(subset_asv_data, temp_asv_file)

  tryCatch({
    source(safe_file_path("DESeq2.R"))
    deseq2_result <- run_deseq2(temp_asv_file, groupings_file_path, deseq2_output_file, seed)
    deseq2_plots <- visualize_deseq2(deseq2_output_file, output_dir, persistent_temp_dir)
    
    source(safe_file_path("Aldex2.R"))
    aldex2_result <- run_aldex2(temp_asv_file, groupings_file_path, aldex2_output_file, seed)
    aldex2_plots <- visualize_aldex2(aldex2_output_file, output_dir, persistent_temp_dir)

    source(safe_file_path("edgeR.R"))
    edgeR_result <- run_edgeR(temp_asv_file, groupings_file_path, edgeR_output_file, seed)
    edgeR_plots <- visualize_edgeR(edgeR_output_file, output_dir, persistent_temp_dir)

    source(safe_file_path("Maaslin2.R"))
    maaslin2_result <- run_maaslin2(temp_asv_file, groupings_file_path, maaslin2_output_file, output_dir, seed)
    maaslin2_plots <- visualize_maaslin2(maaslin2_output_file, output_dir, persistent_temp_dir)

    source(safe_file_path("metagenomeSeq.R"))
    metagenomeseq_result <- run_metagenomeseq(temp_asv_file, groupings_file_path, metagenomeseq_output_file, seed)
    metagenomeseq_plots <- visualize_metagenomeseq(metagenomeseq_output_file, output_dir, persistent_temp_dir)

    source(safe_file_path("overlap_plots.R"))
    overlap_plots <- create_overlap_plots(deseq2_output_file, aldex2_output_file, edgeR_output_file, maaslin2_output_file, metagenomeseq_output_file, output_dir, persistent_temp_dir)

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
      metagenomeseq_plot1 = base64enc::base64encode(metagenomeseq_plots$plot1),
      metagenomeseq_plot2 = base64enc::base64encode(metagenomeseq_plots$plot2),
      metagenomeseq_plot3 = base64enc::base64encode(metagenomeseq_plots$plot3),
      overlap_volcano = base64enc::base64encode(overlap_plots$overlap_volcano),
      overlap_pvalue_distribution = base64enc::base64encode(overlap_plots$overlap_pvalue_distribution)
    )
  }, error = function(e) {
    message("Error in quick_explore: ", e$message)
    print("Debug: Error details")
    print(e)
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
  methods <- c("deseq2", "aldex2", "edger", "maaslin2", "metagenomeseq")
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
  valid_methods <- c("deseq2", "aldex2", "edger", "maaslin2", "metagenomeseq")
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
  methods <- c("deseq2", "aldex2", "edger", "maaslin2", "metagenomeseq")
  combined_results <- data.frame(feature = character())

  for (method in methods) {
    file_path <- safe_file_path(persistent_temp_dir, paste0(method, "_results.tsv"))
    if (file.exists(file_path)) {
      results <- read_tsv(file_path, show_col_types = FALSE)
      
      # Determine the significance column based on the method
      sig_col <- switch(method,
                        "deseq2" = "padj",
                        "aldex2" = "we.eBH",
                        "edger" = "FDR",
                        "maaslin2" = "qval",
                        "metagenomeseq" = "pvalues")
      
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
  
  combined_results <- read_tsv(combined_results_file, show_col_types = FALSE)
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
  
  combined_results <- read_tsv(combined_results_file, show_col_types = FALSE)
  
  return(combined_results)
}

#* Generate stability plot
#* @get /generate_stability_plot
#* @serializer contentType list(type="image/png")
function(res) {
  # Read the data
  data <- read_tsv(safe_file_path(persistent_temp_dir, "combined_results.tsv"), show_col_types = FALSE)
  
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

  asv_data <- read_tsv(temp_asv_file, show_col_types = FALSE)
  groupings_data <- read_tsv(temp_groupings_file, show_col_types = FALSE)

  methods <- c("deseq2", "aldex2", "edger", "maaslin2", "metagenomeseq")
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
      } else if (method == "metagenomeseq") {
        source("metagenomeSeq.R")
        run_metagenomeseq(temp_asv_file, temp_shuffled_groupings, output_file, tempdir(), seed)
      }

      # Read results and store significant ASVs
      result <- read_tsv(output_file, show_col_types = FALSE)
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
  if (!("asv_name" %in% colnames(result))) {
    warning("No asv_name column found in results")
    return(character(0))
  }

  significant <- character(0)
  
  if (method == "deseq2" && "padj" %in% colnames(result)) {
    significant <- result$asv_name[!is.na(result$padj) & result$padj < 0.05]
  } else if (method == "aldex2" && "we.eBH" %in% colnames(result)) {
    significant <- result$asv_name[!is.na(result$we.eBH) & result$we.eBH < 0.05]
  } else if (method == "edger" && "FDR" %in% colnames(result)) {
    significant <- result$asv_name[!is.na(result$FDR) & result$FDR < 0.05]
  } else if (method == "maaslin2" && "qval" %in% colnames(result)) {
    significant <- result$asv_name[!is.na(result$qval) & result$qval < 0.05]
  } else if (method == "metagenomeseq" && "pvalues" %in% colnames(result)) {
    significant <- result$asv_name[!is.na(result$pvalues) & result$pvalues < 0.05]
  } else {
    warning(paste("Required p-value column not found for method:", method))
    return(character(0))
  }
  
  return(unique(significant[!is.na(significant)]))
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
    valid_methods <- c("deseq2", "aldex2", "edger", "maaslin2", "metagenomeseq", "overlap")
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

    message("Starting stability metric calculation...")

    leaf_json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json"
    tree_json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json"
    leaf_data <- fromJSON(leaf_json_path, simplifyVector = TRUE) # simplifyVector = TRUE -> data.frame
    tree_data <- fromJSON(tree_json_path, simplifyVector = FALSE) # simplifyVector = FALSE -> list like [[1]] $key [1] value [[2]]...
    message("After reading JSON files...")
    if (destroy) {
      print("=====destroy is true=====")
      for (leaf in names(leaf_data)) {
        leaf_data[[leaf]]$data_point <- c(0, 0)
      }
      write_json(leaf_data, leaf_json_path, pretty = TRUE, auto_unbox = TRUE)
      res$status <- 200
      return(list(message = "Bye."))
    }
    method <- body$method
    missing_methods <- body$missing_methods

    message("Reading path_df.rds...")
    path_df <- readRDS(file.path(persistent_temp_dir, "path_df.rds"))
    all_asvs <- read_tsv(asv_file_path, show_col_types = FALSE)$asv

    message("Filtering paths...")
    filtered_df <- path_df[sapply(path_df$path, function(x) x[6] == method), ]
    connected_paths <- sapply(filtered_df$path, function(x) paste(x, collapse = "_"))
    names(connected_paths) <- filtered_df$leaf_id # has length of all paths that have deseq2 in their path
    print("length of connected_paths:")
    print(length(connected_paths))

    message("Checking if methods_sig_vectors.rds exists...")
    if (file.exists(file.path(persistent_temp_dir, "methods_sig_vectors.rds"))) {
      methods_sig_vectors <- read_rds(file.path(persistent_temp_dir, "methods_sig_vectors.rds"))
    } else {
      message("methods_sig_vectors.rds does not exist, creating new data frame...")
      methods_sig_vectors <- data.frame(
        leaf_id = character(),
        is_significant = I(list()),
        stringsAsFactors = FALSE
      )
    }

    bad_paths <- c()
    current_method_count <- 0

    for (i in 1:length(connected_paths)) {
      path <- connected_paths[[i]]
      id <- names(connected_paths[i])
      matching_files <- list.files(code_dir, pattern = paste0("^", gsub("_", "_", path)), full.names = TRUE)
      
      if (length(matching_files) == 1) {
        r_path <- matching_files[1]
        print("r_path:")
        print(r_path)
        # Download the R script file to user's Downloads directory
        download_path <- file.path("~/Downloads/cal_codes", basename(r_path))
        file.copy(r_path, download_path, overwrite = TRUE)
      } else if (length(matching_files) > 1) {
        message("length(matching_files) > 1")
        r_path <- matching_files[1]
      } else {
        message("length(matching_files) == 0")
        r_path <- NULL
      }
      tryCatch({
        env <- new.env()
        env$ASV_file_path <- asv_file_path
        env$groupings_file_path <- groupings_file_path
        suppressWarnings(suppressMessages(source(r_path, local = env)))
        message("source(r_path) done")

        # save the result data
        tmp_file_name <- basename(r_path) # e.g. "Raw data_Low Abundance Filtering_No Zero-Handling_No Normalization_No Transformation_aldex2_FDR_leaf-528.R"
        tmp_file_name <- substr(tmp_file_name, 1, nchar(tmp_file_name) - 2) # remove .R
        tmp_file_name <- paste0(tmp_file_name, ".rds")
        saveRDS(env$result_data, file.path(persistent_temp_dir, tmp_file_name))
        file.copy(file.path(persistent_temp_dir, tmp_file_name),
                  file.path("~/Downloads/cal_codes", tmp_file_name),
                  overwrite = TRUE)

        sig_asvs <- get_significant_asvs(env$result_data, method)
        is_significant <- as.numeric(all_asvs %in% sig_asvs) # convert to 0-1 vector, length = length(all_asvs)
        methods_sig_vectors <- rbind(methods_sig_vectors, data.frame(
          leaf_id = id,
          is_significant = I(list(is_significant)),
            stringsAsFactors = FALSE
          ))
        current_method_count <- current_method_count + 1
        message(crayon::green$bold(paste("Current method ", method, " count: ", current_method_count)))
        # if (current_method_count > 2) {
        #   break
        # }
      }, error = function(e) {
        bad_paths <- c(bad_paths, id) # just in case
        print(e$message)
      })
    }

    saveRDS(methods_sig_vectors, file.path(persistent_temp_dir, "methods_sig_vectors.rds"))
    file.copy(file.path(persistent_temp_dir, "methods_sig_vectors.rds"),
              file.path("~/Downloads", "methods_sig_vectors.rds"),
              overwrite = TRUE)

    sig_matrix <- do.call(rbind, methods_sig_vectors$is_significant)
    n_neighbors <- as.integer(dim(sig_matrix)[1]/3)+1
    umap_result <- umap(sig_matrix,
                      n_neighbors = n_neighbors,     # Due to large data size, use more neighbors
                      min_dist = 0.1,
                      metric = "hamming",
                      n_components = 2,
                      n_epochs = 200,       # Due to large data size, increase training epochs
                      init = "spectral",    # For large datasets, use spectral initialization
                      verbose = FALSE)       # Show progress

    plot_data <- data.frame(
      x = umap_result[,1],
      y = umap_result[,2],
      leaf_id = methods_sig_vectors$leaf_id
    )

    # Update leaf_data with UMAP coordinates for each leaf_id
    for (i in 1:nrow(plot_data)) {
      leaf_id <- plot_data$leaf_id[i]
      leaf_data[[leaf_id]]$data_point <- c(
        plot_data$x[i],
        plot_data$y[i]
      )
    }

    write_json(leaf_data, leaf_json_path, pretty = TRUE, auto_unbox = TRUE)

    res$status <- 200
    return(list(message = "Stability metric calculated successfully"))
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

#* Check path and visualize results
#* @post /check_path_and_visualize
#* @parser json
function(req, res) {
  tryCatch({
    body <- fromJSON(req$postBody)
    path_array <- body$path
    # print the list of files all ends with png
    png_files <- list.files(persistent_temp_dir, pattern = "\\.png$", full.names = TRUE)
    message("PNG files found:")
    for (file in png_files) {
      message(basename(file))
    }

    # Join path elements with underscores
    path_string <- paste(path_array, collapse = "_")
    # Look for matching RDS files in persistent_temp_dir
    pattern <- paste0("^", path_string)
    matching_files <- list.files(persistent_temp_dir, pattern = pattern, full.names = TRUE)
    
    if (length(matching_files) == 0) {
      res$status <- 404
      return(list(error = "No matching RDS file found"))
    }
    
    # Create output directory for visualizations
    output_dir <- tempfile()
    dir.create(output_dir)
    
    # Source the appropriate R file
    method <- switch(path_array[6],
      "deseq2" = "DESeq2",
      "aldex2" = "ALDEx2",
      "edger" = "edgeR",
      "maaslin2" = "MaAsLin2",
      "metagenomeseq" = "metagenomeSeq"
    )
    source(safe_file_path(paste0(method, ".R")))

    # Load RDS data
    result_data <- readRDS(matching_files[1])
    print("result_data:")
    print(head(result_data))
    
    # Call visualization function based on method
    plots <- switch(method,
      "DESeq2" = visualize_deseq2(result_data, output_dir, persistent_temp_dir),
      "ALDEx2" = visualize_aldex2(result_data, output_dir, persistent_temp_dir),
      "edgeR" = visualize_edgeR(result_data, output_dir, persistent_temp_dir),
      "MaAsLin2" = visualize_maaslin2(result_data, output_dir, persistent_temp_dir),
      "metagenomeSeq" = visualize_metagenomeseq(result_data, output_dir, persistent_temp_dir),
      stop("Invalid method")
    )
    
    # Return both result data and plots
    list(
      result_data = result_data,
      plot1 = base64enc::base64encode(plots$plot1),
      plot2 = base64enc::base64encode(plots$plot2),
      plot3 = base64enc::base64encode(plots$plot3)
    )
  }, error = function(e) {
    res$status <- 500
    list(error = paste("Error processing request:", e$message))
  })
}