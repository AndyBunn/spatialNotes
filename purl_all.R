# purl_all.R
# Run this from the project root to extract R scripts from all chapter .qmds.
# Outputs go to rScripts/, with data paths adjusted for the new location.

library(knitr)

# make sure we're in the right place
if (!file.exists("index.qmd")) {
  stop("Run this script from the project root (where index.qmd lives).")
}

dir.create("rScripts", showWarnings = FALSE)

# grab all chapter qmds - skip files with no real code
skip <- c("index.qmd", "00Welcome.qmd", "08PackagesUsed.qmd")
qmds <- setdiff(list.files(pattern = "\\.qmd$"), skip)

for (qmd in qmds) {
  stem   <- tools::file_path_sans_ext(qmd)
  out    <- file.path("rScripts", paste0(stem, ".R"))
  
  message("Purling: ", qmd, " -> ", out)
  knitr::purl(qmd, output = out, documentation = 1)
  
  # fix data paths: "data/ -> "../data/ so the scripts run from rScripts/
  lines <- readLines(out)
  lines <- gsub('"data/', '"../data/', lines, fixed = TRUE)
  writeLines(lines, out)
}

message("Done. Scripts are in rScripts/")
