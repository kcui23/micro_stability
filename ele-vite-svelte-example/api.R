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
  # Extract the ASV and groupings content from the JSON body
  body <- fromJSON(req$postBody)

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  # Write ASV and groupings content to temporary files
  writeLines(as.character(body$asv), temp_asv_file)
  writeLines(as.character(body$groupings), temp_groupings_file)
  
  if (method == "deseq2") {
    # Execute the DESeq2 script
    source("DESeq2.R")
    result <- run_deseq2(temp_asv_file, temp_groupings_file, output_file)
    plots <- visualize_deseq2(output_file, output_dir)
  } else if (method == "aldex2") {
    # Execute the Aldex2 script
    source("Aldex2.R")
    result <- run_aldex2(temp_asv_file, temp_groupings_file, output_file)
    plots <- visualize_aldex2(output_file, output_dir)
  } else {
    stop("Invalid method selected. Please choose 'deseq2' or 'aldex2'.")
  }

  # Return the plots as base64 encoded strings
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
  # Extract the ASV content and threshold from the JSON body
  body <- fromJSON(req$postBody)
  
  # Read the ASV data from the request
  asv_data <- read_tsv(body$asv, col_types = cols(.default = "c"))
  
  # Convert columns to numeric where applicable (excluding the first column if it's row names)
  asv_data <- asv_data %>% mutate(across(-1, as.numeric))
  
  # Apply the threshold filter (e.g., filter out rows where all values are below the threshold)
  threshold <- as.numeric(body$threshold)
  filtered_data <- asv_data %>% filter_if(is.numeric, any_vars(. > 50 * threshold))
  print("=======diff after filter=======")
  print((nrow(asv_data) - nrow(filtered_data)) / nrow(asv_data) * 100)
  
  # Convert the filtered data back to a TSV string
  temp_filtered_file <- tempfile(fileext = ".tsv")
  write_tsv(filtered_data, temp_filtered_file)
  filtered_asv_content <- paste(readLines(temp_filtered_file), collapse = "\n")
  
  
  # Return the filtered content
  list(filteredAsv = filtered_asv_content)
}

