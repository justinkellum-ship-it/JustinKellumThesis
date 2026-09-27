#!/usr/bin/env Rscript
# run_all.R --------------------------------------------------------------------
# Regenerates every analysis table, result and figure of the thesis from the
# raw logger record and the public geospatial inputs, in order:
#
#   00_prepare_inputs.R    binning, network median, night selection (seconds)
#   01_site_predictors.R   C-CAP / FEMA / SVF / shoreline descriptors (~7 min)
#   02_models.R            all mixed models and checks (~1 min)
#   03_figures.R           the 19 figures (~10 min; the whole-island C-CAP
#                          mosaics are aggregated once and cached)
#   06_elevation_check.R   the terrain check: site elevation added to the models,
#                          plus its figure (~1 min)
#   07_sensitivity_checks.R  site checks: models without HNL08, with 3 m instead of 5 m
#                          for buildings without a recorded height (and without the sites
#                          dominated by them), and with the sky view factor of positions
#                          inside footprints recomputed (~7 min)
#   08_appendix_figures.R  the land-cover and terrain maps and the sky view factor
#                          diagram of Appendix B (~1 min)
#
# Usage:  Rscript analysis/run_all.R          (from the project folder)
# Requires R >= 4.3 with the packages listed in analysis/requirements.R.
# ------------------------------------------------------------------------------
if (!isTRUE(l10n_info()[["UTF-8"]])) suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
Sys.setenv(PROJ_NETWORK = "OFF")
here <- {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath("analysis")
}
Sys.setenv(THESIS_ROOT = normalizePath(file.path(here, "..")))
t_all <- Sys.time()
for (s in c("00_prepare_inputs.R", "01_site_predictors.R", "02_models.R", "03_figures.R", "06_elevation_check.R",
            "07_sensitivity_checks.R", "08_appendix_figures.R")) {
  cat("\n=== ", s, "\n", sep = "")
  t0 <- Sys.time()
  status <- system2("Rscript", c(shQuote(file.path(here, s))), env = c(paste0("THESIS_ROOT=", Sys.getenv("THESIS_ROOT")), "PROJ_NETWORK=OFF"))
  if (status != 0) stop("script ", s, " failed with status ", status)
  cat(sprintf("    done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat(sprintf("\nall done in %.1f min; results in results/\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("THESIS_ROOT"), "results", "sessionInfo.txt"))
