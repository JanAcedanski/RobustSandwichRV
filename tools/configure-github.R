args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1L]])) {
  stop("Usage: Rscript tools/configure-github.R GITHUB_OWNER")
}

owner <- args[[1L]]
files <- c("README.Rmd", "README.md")
for (path in files[file.exists(files)]) {
  text <- readLines(path, warn = FALSE, encoding = "UTF-8")
  text <- gsub("OWNER/RobustSandwichRV", paste0(owner, "/RobustSandwichRV"),
               text, fixed = TRUE)
  writeLines(text, path, useBytes = TRUE)
}

description <- readLines("DESCRIPTION", warn = FALSE, encoding = "UTF-8")
description <- description[!grepl("^(URL|BugReports):", description)]
description <- c(
  description,
  paste0("URL: https://github.com/", owner, "/RobustSandwichRV"),
  paste0("BugReports: https://github.com/", owner,
         "/RobustSandwichRV/issues")
)
writeLines(description, "DESCRIPTION", useBytes = TRUE)

pkgdown <- readLines("_pkgdown.yml", warn = FALSE, encoding = "UTF-8")
pkgdown <- pkgdown[!grepl("^url:", pkgdown)]
pkgdown <- c(paste0("url: https://", owner,
                    ".github.io/RobustSandwichRV/"), pkgdown)
writeLines(pkgdown, "_pkgdown.yml", useBytes = TRUE)

message("Configured GitHub links for ", owner, "/RobustSandwichRV")
