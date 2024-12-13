library(tidyr)
library(jsonlite)

start_lib <- "library(readr)
library(tidyr)
library(dplyr)
library(readr)
library(e1071)
library(phyloseq)
library(ggplot2)"

start_read_files <- "data <- read_tsv(ASV_file_path, comment = '', col_names = TRUE, skip = ifelse(grepl('Constructed from biom file', readLines(ASV_file_path, n=1)), 1, 0))
groupings <- read_tsv(groupings_file_path, col_names = TRUE)"
filter_low_abundance <- "data %>% mutate_if(is.numeric, ~ifelse(. < threshold, 0, .))"
filter_prevalence <- "total_samples <- ncol(data) - 1
    data %>%
      mutate(prevalence = rowSums(select(., -1) > 0) / total_samples) %>%
      filter(prevalence >= threshold) %>%
      select(-prevalence)"
filter_variance <- "data %>%
      mutate(variance = apply(select(., -1), 1, var)) %>%
      filter(variance > threshold) %>%
      select(-variance)"
filter_no <- ""

zero_pseudocount <- "data %>% mutate(across(-1, ~. + param_value))"
zero_knn <- "library(edgeR)
library(impute)
tmp_data <- data %>% dplyr::select(-asv)
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
      dplyr::select(asv, everything())
    data <- tmp_data"
zero_no <- ""

norm_tss <- "data <- data %>% mutate(across(-1, ~./sum(.)))"
norm_css <- "tmp_data <- data %>% dplyr::select(-asv)
    MR_obj <- newMRexperiment(tmp_data)
    MR_obj_normalized <- cumNorm(MR_obj)
    normalized_counts <- MRcounts(MR_obj_normalized, norm = TRUE)
    normalized_counts <- as.data.frame(normalized_counts) %>% 
      mutate(asv = data$asv) %>% 
      dplyr::select(asv, everything())
    data <- normalized_counts"
norm_tmm <- "tmp_data <- data %>% dplyr::select(-asv)
    dge <- DGEList(counts = tmp_data)
    dge <- calcNormFactors(dge, method = 'TMM')
    tmm_normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
    tmm_normalized_counts <- as.data.frame(tmm_normalized_counts) %>% 
      mutate(asv = data$asv) %>% 
      dplyr::select(asv, everything())
    data <- tmm_normalized_counts"
norm_clr <- "tmp_data <- data %>% dplyr::select(-asv)
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
    data <- clr_manual"
norm_no <- ""

trans_log <- "data <- data %>% mutate(across(-1, ~log(.)))"
trans_logit <- "logit <- function(p) {
      log(p / (1 - p))
    }
    data <- data %>% mutate(across(-1, ~logit(.)))"
trans_ast <- "data <- data %>% mutate(across(-1, ~asin(sqrt(.))))"
trans_no <- ""

pre_method <- "set.seed(seed)
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
  }"

method_deseq2 <- "library(DESeq2)

  groupings <- groupings[,-1, drop=FALSE]
  # Run DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                        colData = groupings,
                                        design = ~ comparison)
  dds_res <- DESeq2::DESeq(dds, sfType = 'poscounts')
  res <- DESeq2::results(dds_res, tidy = TRUE, format = 'DataFrame')

  # Rename 'row' column to 'asv_name'
  colnames(res)[colnames(res) == 'row'] <- 'asv_name'
  # Ensure 'asv_name' is the first column
  res <- res[, c('asv_name', setdiff(colnames(res), 'asv_name'))]
  result_data <- as.data.frame(res)
  "
vis_deseq2 <- "deseq2_results <- result_data
  
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
deseq2_results$log10pvalue <- -log10(deseq2_results$pvalue)
deseq2_results$point_size <- ifelse(deseq2_results$pvalue < 0.05, 2, 1)
  
# Save data for drawing p1 to the persistent directory
  volcano_plot_data <- deseq2_results %>%
    dplyr::select(asv_name, log2FoldChange, pvalue) %>%
    mutate(method = 'DESeq2')

p1 <- ggplot(deseq2_results, aes(x = log2FoldChange, y = log10pvalue)) +
    geom_point(aes(color = pvalue < 0.05, size=point_size)) +
    scale_color_manual(values = c('gray', '#4C3BCF')) +
    scale_size_continuous(range = c(1, 3), guide = 'none') +
    theme_minimal() +
    labs(title = 'Volcano plot', x = 'log2FoldChange', y = '-log10(pvalue)')

# Visualization 2: MA plot of baseMean vs. log2FoldChange
p2 <- ggplot(deseq2_results, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = pvalue < 0.05), alpha = 0.5) +
    scale_color_manual(values = c('gray', '#4C3BCF'), name = 'Significant') +
    theme_minimal() +
    scale_x_log10() +
    labs(title = 'MA plot', x = 'baseMean', y = 'log2FoldChange')
    
# Visualization 3: Histogram of p-value distribution
p3 <- ggplot(deseq2_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = 'Histogram of p-value distribution', x = 'p-value', y = 'Frequency')"

method_aldex2 <-"library(ALDEx2)

  conditions <- groupings[,2]
  unique_groups <- length(unique(conditions))
  
  if (unique_groups == 2) {
    test_type <- 't'
  } else if (unique_groups > 2) {
    test_type <- 'kw'
  }

  zero_proportion <- mean(ASV_table == 0)
  feature_skewness <- skewness(colSums(ASV_table))

  if (zero_proportion > 0.5) {
    if (feature_skewness > 1) {
      denom_type <- 'iqlr'
    } else {
      denom_type <- 'zero'
    }
  } else {
    if (feature_skewness > 1) {
      denom_type <- 'iqlr'
    } else {
      denom_type <- 'all'
    }
  }

  # Run ALDEx2 analysis
  results <- aldex(reads = ASV_table, 
                   conditions = conditions, 
                   mc.samples = 128, 
                   test = test_type, 
                   effect = TRUE,
                   include.sample.summary = FALSE, 
                   verbose = TRUE, 
                   denom = denom_type)

  # Add ASV names to the results
  results$asv_name <- rownames(ASV_table)
  results <- results[, c('asv_name', setdiff(colnames(results), 'asv_name'))]
  result_data <- as.data.frame(results)"
vis_aldex2 <- "aldex2_results <- result_data
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
  alog <- aldex2_results %>%
    mutate(neg_log10_pvalue = -log10(we.eBH),
           point_size = ifelse(we.eBH < 0.05, 3, 1))
  
  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- alog %>%
    dplyr::select(asv_name, effect, we.eBH) %>%
    mutate(method = 'ALDEx2', log2FoldChange = effect, pvalue = we.eBH)
  
  p1 <- ggplot(alog, aes(x = effect, y = neg_log10_pvalue)) +
    geom_point(aes(color = we.eBH < 0.05, size=point_size)) +  # Coloring significant points
    scale_color_manual(values = c('gray', '#4C3BCF')) +
    scale_size_continuous(range = c(1, 3), guide = 'none') +
    theme_minimal() +
    labs(title = 'Volcano plot of Effect Size vs. -log10(pvalue)',
         x = 'Log2 Fold Change',
         y = '-log10(pvalue)')

  # Visualization 2: Scatter plot of diff.btw vs. effect
  p2 <- ggplot(alog, aes(x = rab.all, y = effect)) +
    geom_point(aes(color = we.eBH < 0.05, size=point_size), alpha = 0.5) +
    scale_color_manual(values = c('gray', '#4C3BCF'), name = 'Significant') +
    scale_size_continuous(range = c(1, 3), guide = 'none') +
    theme_minimal() +
    scale_x_log10() +
    labs(title = 'MA-like plot', x = 'Mean Abundance', y = 'Effect Size')
    
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(aldex2_results, aes(x = we.eBH)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() +
    labs(title = 'Histogram of p-value distribution', x = 'pvalue', y = 'Frequency')"

method_edger <- "library(edgeR)
phyloseq_to_edgeR <- function(physeq, group, method='RLE', ...){
  require('edgeR')
  require('phyloseq')
  if(!taxa_are_rows(physeq)) { physeq <- t(physeq) }
  x <- as(otu_table(physeq), 'matrix')
  x <- x + 1  # Add one to protect against overflow, log(0) issues
  if(identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1){
    group <- get_variable(physeq, group)
  }
  taxonomy <- tax_table(physeq, errorIfNULL=FALSE)
  if(!is.null(taxonomy)){
    taxonomy <- data.frame(as(taxonomy, 'matrix'))
  }
  y <- DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  z <- calcNormFactors(y, method=method)
  
  if(!all(is.finite(z$samples$norm.factors))){
    stop('Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors, consider changing `method` argument')
  }
  
  # Estimate dispersions
  if (length(unique(group)) > 1 && min(table(group)) > 1) {
    z <- estimateCommonDisp(z)
    z <- estimateTagwiseDisp(z)
  } else {
    z$common.dispersion <- NA
    warning('No replication, setting dispersion to NA.')
  }
  
  return(z)
}
OTU <- phyloseq::otu_table(ASV_table, taxa_are_rows = TRUE)
  sampledata <- phyloseq::sample_data(groupings, errorIfNULL = TRUE)
  phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
  
  test <- phyloseq_to_edgeR(physeq = phylo, group=groupings[,2])
  et <- exactTest(test)
  tt <- topTags(et, n=nrow(test$table), adjust.method='fdr', sort.by='PValue')
  res <- tt@.Data[[1]]

  # Add ASV names to the results
  res$asv_name <- rownames(ASV_table)
  res <- res[, c('asv_name', setdiff(colnames(res), 'asv_name'))]
  result_data <- as.data.frame(res)"
vis_edger <- "edgeR_results <- result_data
  # Visualization 1: Volcano plot of log2FoldChange vs. -log10(pvalue)
  edgeR_results <- edgeR_results %>%
    mutate(neg_log10_pvalue = -log10(PValue),
           point_size = ifelse(FDR < 0.05, 3, 1))

  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- edgeR_results %>%
    dplyr::select(asv_name, logFC, PValue) %>%
    mutate(method = 'edgeR', log2FoldChange = logFC, pvalue = PValue)
  
  p1 <- ggplot(edgeR_results, aes(x = logFC, y = neg_log10_pvalue)) +
    geom_point(aes(color = FDR < 0.05, size=point_size)) +  # Coloring significant points
    scale_color_manual(values = c('gray', '#4C3BCF')) +  # Red for significant
    scale_size_continuous(range = c(1, 3), guide = 'none') +
    theme_minimal() +
    labs(title = 'Volcano plot of logFC vs. -log10(pvalue)',
         x = 'Log2 Fold Change',
         y = '-log10(pvalue)')

  # Visualization 2: Scatter plot of logCPM vs. logFC
  p2 <- ggplot(edgeR_results, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = FDR < 0.05), alpha = 0.5) +
    scale_color_manual(values = c('gray', '#4C3BCF'), name = 'Significant') +
    theme_minimal() +
    labs(title = 'Scatter plot of logCPM vs. logFC',
         x = 'logCPM',
         y = 'logFC')
         
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(edgeR_results, aes(x = PValue)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() +
    labs(title = 'Histogram of p-value distribution', x = 'pvalue', y = 'Frequency')"

method_maaslin2 <- "library(Maaslin2)
ASV_table <- data.frame(t(ASV_table), check.rows = F, check.names = F, stringsAsFactors = F)
temp_dir <- tempfile()
dir.create(temp_dir)
# Run Maaslin2 analysis
fit_data <- Maaslin2(
  ASV_table, 
  groupings, 
  output = temp_dir,
  transform = 'AST',
  fixed_effects = c(colnames(groupings)[2]),
  standardize = FALSE, 
  plot_heatmap = FALSE, 
  plot_scatter = FALSE
)

all_results <- fit_data$results
colnames(all_results)[colnames(all_results) == 'feature'] <- 'asv_name'
result_data <- all_results[, c('asv_name', setdiff(colnames(all_results), 'asv_name'))]"
vis_maaslin2 <- "maaslin2_results <- result_data
# Visualization 1: Volcano plot of effect vs. -log10(pvalue)
  maaslin2_results <- maaslin2_results %>%
    mutate(neg_log10_pvalue = -log10(pval),
           point_size = ifelse(qval < 0.05, 3, 1))
  
  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- maaslin2_results %>%
    dplyr::select(asv_name, coef, pval) %>%
    mutate(method = 'Maaslin2', log2FoldChange = coef, pvalue = pval)
  
  p1 <- ggplot(maaslin2_results, aes(x = coef, y = neg_log10_pvalue)) +
    geom_point(aes(color = qval < 0.05, size=point_size)) +
    scale_color_manual(values = c('gray', '#4C3BCF')) +
    scale_size_continuous(range = c(1, 3), guide = 'none') +
    theme_minimal() +
    labs(title = 'Volcano plot of Effect Size vs. -log10(pvalue)',
         x = 'Effect Size',
         y = '-log10(pvalue)')

  # Visualization 2: MA-like scatter plot
  #use 'N.not.0' as a proxy for 'baseMean' and 'coef' for 'log2FoldChange'
  # Check if 'N.not.0' exists, otherwise use 'N.not.zero'
  x_column <- if('N.not.0' %in% colnames(maaslin2_results)) 'N.not.0' else 'N.not.zero'
  p2 <- ggplot(maaslin2_results, aes(x = .data[[x_column]], y = coef)) +
    geom_point(aes(color = qval < 0.05), alpha = 0.5) +
    scale_color_manual(values = c('gray', '#4C3BCF')) +
    theme_minimal() +
    scale_x_log10() +
    labs(title = 'MA-like plot',
         x = 'Number of Non-Zero Samples',
         y = 'Coefficient (Effect Size)')
         
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(maaslin2_results, aes(x = pval)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() +
    labs(title = 'Histogram of p-value distribution', x = 'pvalue', y = 'Frequency')"

method_metagenomeseq <- "library(metagenomeSeq)
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
    result_data <- as.data.frame(result)"
vis_metagenomeseq <- "metagenomeSeq_results <- result_data
# Visualization 1: Volcano plot of logFC vs. -log10(pvalue)
  metagenomeSeq_results$log10pvalue <- -log10(metagenomeSeq_results$pvalues)
  metagenomeSeq_results$point_size <- ifelse(metagenomeSeq_results$pvalues < 0.05, 2, 1)

  # Save data for drawing p1 to the persistent directory
  volcano_plot_data <- metagenomeSeq_results %>%
    dplyr::select(asv_name, logFC, pvalues) %>%
    mutate(method = 'metagenomeSeq', log2FoldChange = logFC, pvalue = pvalues)
  
  p1 <- ggplot(metagenomeSeq_results, aes(x = logFC, y = log10pvalue)) +
    geom_point(aes(color = pvalues < 0.05, size = point_size)) +
    scale_color_manual(values = c('gray', '#4C3BCF')) +
    scale_size_continuous(range = c(1, 3), guide = 'none') +
    theme_minimal() +
    labs(title = 'Volcano Plot', x = 'log2 Fold Change', y = '-log10(p-value)')
  
  # Visualization 2: MA plot of baseMean vs. log2 Fold Change
  metagenomeSeq_results$baseMean <- rowMeans(as.matrix(metagenomeSeq_results[, -c(1:4)]))  # Assuming baseMean is not directly in the file
  
  p2 <- ggplot(metagenomeSeq_results, aes(x = baseMean, y = logFC)) +
    geom_point(aes(color = pvalues < 0.05), alpha = 0.5) +
    scale_color_manual(values = c('gray', '#4C3BCF'), name = 'Significant') +
    theme_minimal() +
    scale_x_log10() +
    labs(title = 'MA Plot', x = 'baseMean', y = 'log2 Fold Change')
    
  # Visualization 3: Histogram of p-value distribution
  p3 <- ggplot(metagenomeSeq_results, aes(x = pvalues)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = 'Histogram of p-value Distribution', x = 'p-value', y = 'Frequency')"
  
end <- "result_plots <- list(p1, p2, p3)

write_tsv(result_data, paste0(method, '_results.tsv'))
for (i in 1:length(result_plots)) {
  ggsave(paste0(method, '_plot', i, '.png'), plot = result_plots[[i]])
}"


json_path <- "/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json"
data <- fromJSON(json_path, simplifyVector = FALSE)

paths_list <- list()
ids_list <- character()

traverse <- function(node, current_path) {
    new_path <- c(current_path, node$name)
    if (!is.null(node$children)) {
        for (child in node$children) {
        traverse(child, new_path)
        }
    } else {
        paths_list[[length(paths_list) + 1]] <<- new_path
        ids_list <<- c(ids_list, node$id)
    }
}
traverse(data, c())

generate_all_r_code <- function(params, code_dir) {
    parameters <- unlist(params)
    for (i in 1:length(paths_list)) {
        cstr <- paste(paths_list[[i]][-1], collapse = "_")
        code <- generate_r_code_for_each_combination(cstr, parameters)
        writeLines(code, file.path(code_dir, paste0("Raw data_", cstr, "_", ids_list[i], ".R")))
    }
}

generate_r_code <- function(selectedOperations, params) {
    processed <- process_frontend_inputs(selectedOperations, params)
    operations <- processed$operations
    parameters <- processed$parameters

    cartesian_strings_tidyr <- generate_cartesian_product(operations)
    code <- generate_r_code_for_each_combination(cartesian_strings_tidyr[1], parameters, complete_code = TRUE)
    return(code)
}

process_frontend_inputs <- function(selectedOperations, params) {
    operations <- lapply(selectedOperations, function(op) {
        if (length(op) == 0) {
        NULL
        } else {
        unlist(op)
        }
    })
    operations <- operations[!sapply(operations, is.null)]
    parameters <- unlist(params)

    return(list(operations = operations, parameters = parameters))
}

generate_cartesian_product <- function(operations) {
    filtering <- operations$Filtering
    zero_handling <- operations$`Zero-Handling`
    normalization <- operations$Normalization
    transformation <- operations$Transformation
    model_perturbation <- operations$`Model Perturbation`

    lists <- list(filtering, zero_handling, normalization, transformation, model_perturbation)
    cartesian_product <- expand.grid(lists)
    cartesian_strings_tidyr <- apply(cartesian_product, 1, paste, collapse = "_")
    return(cartesian_strings_tidyr)
}

generate_r_code_for_each_combination <- function(cartesian_string, parameters, complete_code = FALSE) {
    split_string <- strsplit(cartesian_string, "_")[[1]]
    filtering <- split_string[1]
    zero_handling <- split_string[2]
    normalization <- split_string[3]
    transformation <- split_string[4]
    model_perturbation <- split_string[5]

    threshold <- switch(filtering,
        'Low Abundance Filtering' = parameters[['Low Abundance Filtering']],
        'Prevalence Filtering' = parameters[['Prevalence Filtering']],
        'Variance Filtering' = parameters[['Variance Filtering']],
        'No Filtering' = 0
    )
    param_value <- switch(zero_handling,
        'Pseudocount Addition' = parameters[['Pseudocount Addition']],
        'k-NN Imputation' = c(parameters[['knn']], parameters[['knn_bound']]),
        'No Zero-Handling' = 0
    )
    seed <- parameters[['Random Seed']]

    # set up those parameters
    newline <- "\n"
    param_set_up <- paste("threshold <- ", threshold, newline,
                          "param_value <- ", 
                          ifelse(length(param_value) > 1, 
                                paste("c(", paste(param_value, collapse = ", "), ")", sep = ""),
                                param_value),
                          newline,
                          "seed <- ", seed, newline,
                          "method <- '", model_perturbation, "'", newline, sep = "")

    filter_code <- switch(filtering,
        'Low Abundance Filtering' = filter_low_abundance,
        'Prevalence Filtering' = filter_prevalence,
        'Variance Filtering' = filter_variance,
        'No Filtering' = filter_no
    )
    zero_code <- switch(zero_handling,
        'Pseudocount Addition' = zero_pseudocount,
        'k-NN Imputation' = zero_knn,
        'No Zero-Handling' = zero_no
    )
    normalization_code <- switch(normalization,
        'TSS' = norm_tss,
        'CSS' = norm_css,
        'TMM' = norm_tmm,
        'CLR' = norm_clr,
        'No Normalization' = norm_no
    )
    trans_code <- switch(transformation,
        'Log' = trans_log,
        'Logit' = trans_logit,
        'AST' = trans_ast,
        'No Transformation' = trans_no
    )
    model_code <- switch(model_perturbation,
        'deseq2' = method_deseq2,
        'aldex2' = method_aldex2,
        'edger' = method_edger,
        'maaslin2' = method_maaslin2,
        'metagenomeseq' = method_metagenomeseq
    )
    vis_code <- switch(model_perturbation,
        'deseq2' = vis_deseq2,
        'aldex2' = vis_aldex2,
        'edger' = vis_edger,
        'maaslin2' = vis_maaslin2,
        'metagenomeseq' = vis_metagenomeseq
    )

    if (complete_code) {
      code <- paste(start_lib, param_set_up, start_read_files, 
                    filter_code, zero_code, normalization_code, trans_code, 
                    pre_method, model_code, vis_code, end, sep = '\n')
    } else {
      code <- paste(start_lib, param_set_up, start_read_files, 
                    filter_code, zero_code, normalization_code, trans_code, 
                    pre_method, model_code, vis_code, sep = '\n')
    }
    return(code)
}