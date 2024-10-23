library(readr)
library(tidyr)
library(dplyr)
library(readr)
library(e1071)
library(phyloseq)
library(ggplot2)
threshold <- 0
param_value <- c(5, 50)
seed <- 1234
method <- 'metagenomeseq'

data <- read_tsv(ASV_file_path, comment = '', col_names = TRUE, skip = ifelse(grepl('Constructed from biom file', readLines(ASV_file_path, n=1)), 1, 0))
groupings <- read_tsv(groupings_file_path, col_names = TRUE)

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
        warning('Error in k-NN imputation: ', e$message)
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
      select(asv, everything())
    data <- tmp_data


set.seed(seed)
    ASV_table <- as.data.frame(data)
    row.names(ASV_table) <- ASV_table[,1]
    ASV_table <- ASV_table[,-1]
    groupings <- as.data.frame(groupings)
    row.names(groupings) <- groupings[,1]

  # Check if the same number of samples are being input
  sample_num <- ncol(ASV_table)
  grouping_num <- nrow(groupings)

  if (sample_num != grouping_num) {
    message('The number of samples in the ASV table and the groupings table are unequal')
    message('Will remove any samples that are not found in either the ASV table or the groupings table')
  }

  # Check if order of samples match up
  if (identical(colnames(ASV_table), rownames(groupings))) {
    message('Groupings and ASV table are in the same order')
  } else {
    rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
    groupings <- groupings[rows_to_keep,, drop=F]
    ASV_table <- ASV_table[, rows_to_keep]
  }
library(metagenomeSeq)
# Create a MRexperiment object
  count_matrix <- as.matrix(ASV_table)
  groupings_meta <- AnnotatedDataFrame(groupings)

  # Create MRexperiment object
  MRexp <- metagenomeSeq::newMRexperiment(counts = count_matrix, phenoData = groupings_meta)
  MRexp_norm <- metagenomeSeq::cumNorm(MRexp)
  comparison_variable <- colnames(groupings)[2]

  # Fit the feature model
  model <- model.matrix(~ groupings[, comparison_variable])
  fit <- metagenomeSeq::fitFeatureModel(MRexp_norm, model)

  # Get differential abundance results
  result <- metagenomeSeq::MRfulltable(fit, number = nrow(fit))
  result <- result %>%
    arrange(pvalues) %>%
    mutate(asv_name = rownames(.)) %>%
    dplyr::select(asv_name, everything())
result_data <- as.data.frame(result)
