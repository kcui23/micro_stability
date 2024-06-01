library(ALDEx2)
library(tidyr)
library(dplyr)
library(readr)

run_aldex2 <- function(ASV_file, groupings_file, output_file) {
  # Read ASV table
  ASV_table <- read_tsv(ASV_file, comment = "", col_names = TRUE, skip = ifelse(grepl("Constructed from biom file", readLines(ASV_file, n=1)), 1, 0))
  ASV_table <- as.data.frame(ASV_table)
  row.names(ASV_table) <- ASV_table[,1]
  ASV_table <- ASV_table[,-1]

  # Read groupings table
  groupings <- read_tsv(groupings_file, col_names = TRUE)
  groupings <- as.data.frame(groupings)
  row.names(groupings) <- groupings[,1]

  # Check if the same number of samples are being input
  sample_num <- ncol(ASV_table)
  grouping_num <- nrow(groupings)

  if (sample_num != grouping_num) {
    message("The number of samples in the ASV table and the groupings table are unequal")
    message("Will remove any samples that are not found in either the ASV table or the groupings table")
  }

  # Check if order of samples match up
  if (identical(colnames(ASV_table), rownames(groupings))) {
    message("Groupings and ASV table are in the same order")
  } else {
    rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
    groupings <- groupings[rows_to_keep,, drop=F]
    ASV_table <- ASV_table[, rows_to_keep]
    if (identical(colnames(ASV_table), rownames(groupings))) {
      message("Groupings table was re-arranged to be in the same order as the ASV table")
      message("A total of ", sample_num - length(colnames(ASV_table)), " from the ASV_table")
      message("A total of ", grouping_num - length(rownames(groupings)), " from the groupings table")
    } else {
      stop("Unable to match samples between the ASV table and groupings table")
    }
  }

  # Run ALDEx2 analysis
  results <- aldex(reads = ASV_table, 
                   conditions = groupings[,2], 
                   mc.samples = 128, 
                   test = "t", 
                   effect = TRUE,
                   include.sample.summary = FALSE, 
                   verbose = TRUE, 
                   denom = "all")

  # Save results to output file
  write_tsv(as.data.frame(results), output_file)

  message("Results table saved to ", output_file)
  
  # Return paths to the result files (for example, if there are plots)
  list(
    plot1 = "./plot1.png",
    plot2 = "./plot2.png",
    plot3 = "./plot3.png"
  )
}

visualize_aldex2 <- function(input_file, output_dir) {
  # Read ALDEx2 results
  aldex2_results <- read_tsv(input_file)
  
  # Visualization 1: Boxplot of rab.all values
  p1 <- ggplot(aldex2_results, aes(x = "All", y = rab.all)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Boxplot of rab.all", x = "All", y = "rab.all")

  ggsave(file.path(output_dir, "aldex2_plot1.png"), plot = p1)

  # Visualization 2: Scatter plot of diff.btw vs. effect
  p2 <- ggplot(aldex2_results, aes(x = diff.btw, y = effect)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Scatter plot of diff.btw vs. effect", x = "diff.btw", y = "effect")
  
  ggsave(file.path(output_dir, "aldex2_plot2.png"), plot = p2)

# Visualization 3: Density plot of rab.win.Rhizosphere and rab.win.Soil
  print(colnames(aldex2_results))
  aldex2_long <- aldex2_results %>% 
    dplyr::select(rab.win.Rhizosphere, rab.win.Soil) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Condition", values_to = "Value")

  p3 <- ggplot(aldex2_long, aes(x = Value, fill = Condition)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Density plot of rab.win.Rhizosphere and rab.win.Soil", x = "Value", y = "Density")

  ggsave(file.path(output_dir, "aldex2_plot3.png"), plot = p3)
  
  list(
    plot1 = file.path(output_dir, "aldex2_plot1.png"),
    plot2 = file.path(output_dir, "aldex2_plot2.png"),
    plot3 = file.path(output_dir, "aldex2_plot3.png")
  )
}
