library(plumber)
library(readr)
library(base64enc)

#* @apiTitle ASV Data Processing API

#* Upload dataset and select processing method
#* @post /process
#* @param method The processing method to use (deseq2 or aldex2)
#* @parser text
function(req, method) {
  temp_file <- tempfile(fileext = ".tsv")
  writeLines(req$postBody, temp_file)
  
  if (method == "deseq2") {
    # Execute the DESeq2 script
    source("a.R")
    result <- a(temp_file)
  } else if (method == "aldex2") {
    # Execute the Aldex2 script
    source("b.R")
    result <- b(temp_file)
  } else {
    stop("Invalid method selected. Please choose 'deseq2' or 'aldex2'.")
  }

  # Assuming result contains paths to output plots
  list(
    plot1 = base64enc::base64encode(result$plot1),
    plot2 = base64enc::base64encode(result$plot2),
    plot3 = base64enc::base64encode(result$plot3)
  )
}
