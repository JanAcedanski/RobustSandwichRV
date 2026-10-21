args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1L]])) {
  stop("Usage: Rscript tools/configure-github.R GITHUB_OWNER")
}

owner <- args[[1L]]
files <- c("README.Rmd", "README.md")
for (path in files[file.exists(files)]) {
  text <- readLines(path, warn = FALSE, encoding = "UTF-8")
  text <- gsub("OWNER/RVwithSandwichSE", paste0(owner, "/RVwithSandwichSE"),
               text, fixed = TRUE)
  writeLines(text, path, useBytes = TRUE)
}

description <- readLines("DESCRIPTION", warn = FALSE, encoding = "UTF-8")
drop <- logical(length(description))
in_link_field <- FALSE
for (i in seq_along(description)) {
  if (!grepl("^[[:space:]]", description[[i]])) {
    in_link_field <- grepl("^(URL|BugReports):", description[[i]])
  }
  drop[[i]] <- in_link_field
}
description <- description[!drop]
description <- c(
  description,
  paste0("URL: https://github.com/", owner, "/RVwithSandwichSE,"),
  paste0("    https://", owner, ".github.io/RVwithSandwichSE/,"),
  "    https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7639842",
  paste0("BugReports: https://github.com/", owner,
         "/RVwithSandwichSE/issues")
)
writeLines(description, "DESCRIPTION", useBytes = TRUE)

pkgdown <- readLines("_pkgdown.yml", warn = FALSE, encoding = "UTF-8")
pkgdown <- pkgdown[!grepl("^url:", pkgdown)]
pkgdown <- c(paste0("url: https://", owner,
                    ".github.io/RVwithSandwichSE/"), pkgdown)
writeLines(pkgdown, "_pkgdown.yml", useBytes = TRUE)

message("Configured GitHub links for ", owner, "/RVwithSandwichSE")
