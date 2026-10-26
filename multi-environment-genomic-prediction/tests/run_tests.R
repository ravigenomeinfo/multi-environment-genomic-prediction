# Run from the repository root:  Rscript tests/run_tests.R
library(testthat)
source("R/functions.R")
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
