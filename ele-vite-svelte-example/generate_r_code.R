library(tidyr)

generate_all_r_code <- function(selectedOperations, params, code_dir) {
    processed <- process_frontend_inputs(selectedOperations, params)
    operations <- processed$operations
    parameters <- processed$parameters

    cartesian_strings_tidyr <- generate_cartesian_product(operations)
    for (i in 1:length(cartesian_strings_tidyr)) {
        code <- generate_r_code_for_each_combination(cartesian_strings_tidyr[i], parameters)
        writeLines(code, file.path(code_dir, paste0(cartesian_strings_tidyr[i], ".R")))
    }
}

generate_r_code <- function(selectedOperations, params) {
    processed <- process_frontend_inputs(selectedOperations, params)
    operations <- processed$operations
    parameters <- processed$parameters

    cartesian_strings_tidyr <- generate_cartesian_product(operations)
    code <- generate_r_code_for_each_combination(cartesian_strings_tidyr[1], parameters)
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
    stability_metric <- operations$`Stability Metric`

    lists <- list(filtering, zero_handling, normalization, transformation, model_perturbation, stability_metric)
    cartesian_product <- expand.grid(lists)
    cartesian_strings_tidyr <- apply(cartesian_product, 1, paste, collapse = "_")
    return(cartesian_strings_tidyr)
}

generate_r_code_for_each_combination <- function(cartesian_string, parameters) {
    return(paste0("Code for ", cartesian_string))
}