assign(".lib.loc", Sys.getenv("R_LIB_PATHS"), envir = environment(.libPaths))

options(repos = c(CRAN = "https://cloud.r-project.org"))



install_cran_packages <- function(
    cran_pkgs, library_path, type, decompress,
    remove_dirs = c(
      "help", "doc", "tests", "html",
      "include", "unitTests",
      file.path("libs", "*dSYM")
      )) {
  installed <- list.files(library_path) # check installed packages

  cran_to_install <- sort(setdiff(
    unique(unlist(
      c(cran_pkgs,
        tools::package_dependencies(cran_pkgs,
          recursive = TRUE,
          which = c("Depends", "Imports", "LinkingTo")
      ))
    )),
    installed
  ))
  
  if (!length(cran_to_install)) {
    message("No packages to install")
  } else {
    td <- tempdir()
    downloaded <- download.packages(cran_to_install, destdir = td, type = type)
    apply(downloaded, 1, function(x) decompress(x[2], exdir = library_path))
    unlink(downloaded[, 2])
  }
  
  z <- lapply(
    list.dirs(library_path, full.names = TRUE, recursive = FALSE),
    function(x) {
      unlink(file.path(x, remove_dirs), force = TRUE, recursive = TRUE)
    }
  )
  invisible(NULL)
}

install_bioc_packages <- function(packages) {
  # Print installed packages in blue using cat with ANSI color codes
  cat('\033[34m')  # Start blue text
  cat("Currently installed packages:\n")
  cat(paste(sort(installed.packages()[,"Package"]), collapse="\n"))
  cat('\033[0m\n') # Reset color
  library(BiocManager)
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages) > 0) {
    options(repos = BiocManager::repositories())
    phyloseq_idx <- which(missing_packages == "phyloseq")
    if (length(phyloseq_idx) > 0) {
      BiocManager::install("phyloseq", ask = FALSE, update = FALSE, type = "source")
      missing_packages <- missing_packages[-phyloseq_idx]
    }

    if (length(missing_packages) > 0) {
      BiocManager::install(missing_packages, ask = FALSE, update = FALSE, type = "binary")
    }
  }
}

cran_packages <- c(
  "BiocManager",
  "plumber",
  "readr",
  "ggplot2",
  "base64enc",
  "jsonlite",
  "dplyr",
  "tidyr",
  "tidyverse",
  "scales",
  "uwot",
  "httr"
)

bioc_packages <- c(
  "phyloseq",
  "impute",
  "metagenomeSeq",
  "edgeR",
  "Maaslin2",
  "DESeq2",
  "ALDEx2"
)

if (dir.exists("r-mac")) {
  install_cran_packages(
    cran_pkgs = cran_packages, library_path = Sys.getenv("R_LIB_PATHS"),
    type = "mac.binary.big-sur-arm64", decompress = untar
  )
}
install_bioc_packages(bioc_packages)


library(plumber)
plumber::plumb(file='api.R')$run(host='0.0.0.0', port=8000)