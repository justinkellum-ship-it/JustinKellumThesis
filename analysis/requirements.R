# requirements.R ---------------------------------------------------------------
# Installs the packages used by the thesis analysis (run once):
#   Rscript R/requirements.R
# Versions used for the thesis are recorded in results/sessionInfo.txt.
# ------------------------------------------------------------------------------
pkgs <- c(
  # data handling
  "dplyr", "tidyr", "tibble", "readr", "purrr", "stringr", "lubridate", "jsonlite",
  # mixed models and their summaries
  "lme4", "lmerTest", "performance", "car",
  # geospatial
  "sf", "terra",
  # figures
  "ggplot2", "patchwork", "scales", "cowplot",
  # walkthrough document
  "rmarkdown", "knitr", "kableExtra"
)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org") else message("all packages present")
