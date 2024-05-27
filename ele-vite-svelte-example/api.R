library(plumber)
library(readr)
library(ggplot2)
library(base64enc)
library(jsonlite)

#* @apiTitle ASV Data Processing API

#* Upload dataset and select processing method
#* @post /process
#* @param method The processing method to use (deseq2 or aldex2)
#* @parser text
function(req, method) {
  # Extract the ASV and groupings content from the JSON body
  body <- fromJSON(req$postBody)

  temp_asv_file <- tempfile(fileext = ".tsv")
  temp_groupings_file <- tempfile(fileext = ".tsv")
  output_file <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()

  dir.create(output_dir)
  
  # Assuming req$postBody contains the ASV content and groupings content separately
  writeLines(body$asv, temp_asv_file)
  writeLines(body$groupings, temp_groupings_file)
  
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
