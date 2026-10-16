args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
testthat::test_local(root, reporter = "summary")
