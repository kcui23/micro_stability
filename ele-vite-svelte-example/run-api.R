library(plumber)
pr <- plumber::plumb("api.R")
pr$run(port = 8000)