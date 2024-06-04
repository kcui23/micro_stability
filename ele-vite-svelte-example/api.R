library(plumber)
library(readr)
library(ggplot2)
library(base64enc)
library(jsonlite)
library(dplyr)

#* @apiTitle ASV Data Processing API

#* Upload dataset and select processing method
#* @post /process
#* @param method The processing method to use (deseq2 or aldex2)
#* @parser json
function(req, method) {
  body <- fromJSON(req$postBody)

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  writeLines(as.character(body$asv), temp_asv_file)
  writeLines(as.character(body$groupings), temp_groupings_file)
  
  if (method == "deseq2") {
    source("DESeq2.R")
    result <- run_deseq2(temp_asv_file, temp_groupings_file, output_file)
    plots <- visualize_deseq2(output_file, output_dir)
  } else if (method == "aldex2") {
    source("Aldex2.R")
    result <- run_aldex2(temp_asv_file, temp_groupings_file, output_file)
    plots <- visualize_aldex2(output_file, output_dir)
  } else {
    stop("Invalid method selected. Please choose 'deseq2' or 'aldex2'.")
  }

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

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  deseq2_output_file <- tempfile(fileext = ".tsv")
  aldex2_output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  writeLines(as.character(body$asv), temp_asv_file)
  writeLines(as.character(body$groupings), temp_groupings_file)
  
  asv_data <- read_tsv(temp_asv_file, col_types = cols(.default = "c"))
  subset_asv_data <- asv_data %>% sample_frac(0.05)
  write_tsv(subset_asv_data, temp_asv_file)

  source("DESeq2.R")
  deseq2_result <- run_deseq2(temp_asv_file, temp_groupings_file, deseq2_output_file)
  deseq2_plots <- visualize_deseq2(deseq2_output_file, output_dir)
  
  source("Aldex2.R")
  aldex2_result <- run_aldex2(temp_asv_file, temp_groupings_file, aldex2_output_file)
  aldex2_plots <- visualize_aldex2(aldex2_output_file, output_dir)
  
  source("overlap_plots.R")
  overlap_plots <- create_overlap_plots(deseq2_output_file, aldex2_output_file, output_dir)

  list(
    deseq2_plot1 = base64enc::base64encode(deseq2_plots$plot1),
    deseq2_plot2 = base64enc::base64encode(deseq2_plots$plot2),
    deseq2_plot3 = base64enc::base64encode(deseq2_plots$plot3),
    aldex2_plot1 = base64enc::base64encode(aldex2_plots$plot1),
    aldex2_plot2 = base64enc::base64encode(aldex2_plots$plot2),
    aldex2_plot3 = base64enc::base64encode(aldex2_plots$plot3),
    overlap_volcano = base64enc::base64encode(overlap_plots$overlap_volcano),
    overlap_pvalue_distribution = base64enc::base64encode(overlap_plots$overlap_pvalue_distribution)
  )
}
