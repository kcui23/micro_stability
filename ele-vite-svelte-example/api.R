library(plumber)
library(ggplot2)
library(readr)
library(base64enc)

#* @apiTitle Iris Dataset Visualization API

#* Upload dataset and get visualizations
#* @post /upload
#* @parser text
function(req) {
  temp_file <- tempfile(fileext = ".csv")
  writeLines(req$postBody, temp_file)
  
  # Read the CSV file
  iris_df <- tryCatch({
    read_csv(temp_file)
  }, error = function(e) {
    stop("Failed to read CSV file. Ensure it is a valid CSV with appropriate columns.")
  })
  
  # Check if required columns exist
  required_cols <- c("Species", "Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
  if (!all(required_cols %in% names(iris_df))) {
    stop("CSV file must contain the columns: 'Species', 'Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'")
  }

  # plots
  p1 <- ggplot(iris_df, aes(x = Species, y = Sepal.Length, fill = Species)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Sepal Length by Species", x = "Species", y = "Sepal Length")
  
  p2 <- ggplot(iris_df, aes(x = Petal.Length, fill = Species)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Petal Length Density by Species", x = "Petal Length", y = "Density")
  
  p3 <- ggplot(iris_df, aes(x = Petal.Length, y = Petal.Width, color = Species)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "Petal Length vs Petal Width", x = "Petal Length", y = "Petal Width")

  # Save the plots to temporary files
  ggsave("plot1.png", plot = p1, device = "png")
  ggsave("plot2.png", plot = p2, device = "png")
  ggsave("plot3.png", plot = p3, device = "png")

  # Return the plots as base64 encoded strings
  list(
    plot1 = base64enc::base64encode("plot1.png"),
    plot2 = base64enc::base64encode("plot2.png"),
    plot3 = base64enc::base64encode("plot3.png")
  )
}