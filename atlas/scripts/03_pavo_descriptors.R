# Montgomerie (2006) colour descriptors computed with the R package pavo (Maia et al. 2013, 2019).
#
# Usage:
#     Rscript 03_pavo_descriptors.R <data folder>
#
# Reads every cleaned_spectra/*.csv (I_deg, C_deg, 380 ... 780; R in % of the white reference, 1 nm),
# applies procspec(fixneg = "zero") and summary.rspec over 400-700 nm, and writes
# <data folder>/pavo_descriptors.csv with B1, B2, B3, H1, S6, S8 for every geometry.
#
#   B1  total brightness   sum of R
#   B2  mean brightness    mean of R
#   B3  intensity          max R
#   H1  peak wavelength    wavelength of max R
#   S6  contrast           max R - min R
#   S8  chroma             (max R - min R) / B2
#
# Requires pavo >= 2.9 (tested with pavo 2.11.0, R 4.5.3).

suppressPackageStartupMessages(library(pavo))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript 03_pavo_descriptors.R <data folder>")
data_dir <- args[1]

files <- sort(list.files(file.path(data_dir, "cleaned_spectra"), pattern = "\\.csv$", full.names = TRUE))
out <- list()
for (f in files) {
  dataset <- sub("\\.csv$", "", basename(f))
  w <- read.csv(f, check.names = FALSE)
  geometry <- paste0("I", w$I_deg, "_C", w$C_deg)
  refl <- t(as.matrix(w[, -(1:2)]))                 # rows = wavelengths, columns = geometries
  colnames(refl) <- geometry
  # as.rspec warns that some values are negative; procspec(fixneg = "zero") sets them to zero right after
  spec <- suppressMessages(suppressWarnings(as.rspec(data.frame(wl = as.numeric(colnames(w)[-(1:2)]), refl, check.names = FALSE), lim = c(380, 780))))
  spec <- suppressMessages(procspec(spec, fixneg = "zero"))
  s <- summary(spec, subset = c("B1", "B2", "B3", "H1", "S6", "S8"), lim = c(400, 700))
  out[[dataset]] <- data.frame(
    species_id = sub("__previous_sweep$", "", dataset),
    variant = ifelse(grepl("__previous_sweep$", dataset), "previous_sweep", "main"),
    I_deg = w$I_deg, C_deg = w$C_deg, s, row.names = NULL
  )
  cat(sprintf("%-28s %4d geometries\n", dataset, nrow(w)))
}
res <- do.call(rbind, out)
write.csv(res, file.path(data_dir, "pavo_descriptors.csv"), row.names = FALSE)
cat("written", nrow(res), "rows to", file.path(data_dir, "pavo_descriptors.csv"), "\n")
