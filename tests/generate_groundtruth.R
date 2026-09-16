## generate_groundtruth.R
##
## Generates reference decodings that TestImportAgilent.m compares the MATLAB
## importers against, using the CRAN chromConverter package as an independent
## reference. Works on any folder of Agilent ".D" data directories.
##
## For every *.D folder found (recursively) under the dataset root, it decodes:
##   *.uv  -> uv_<key>.csv  : rt column + one column per wavelength (raw)
##   *.ch  -> ch_<key>.csv  : rt, raw                                (raw)
##   *.MS  -> ms_<key>.csv  : rt, tic  (+ an intensity aggregate)    (raw)
## where <key> = sanitized "<.D folder name>__<signal file name>".
##
## Files whose version chromConverter cannot read are skipped (a warning is
## printed); the MATLAB test then simply has no reference for them.
##
## Usage (from a shell, with R + chromConverter installed):
##   Rscript generate_groundtruth.R [datasetRoot] [outputDir]
## Defaults: datasetRoot = ../datasets (relative to this script),
##           outputDir   = <system temp>/ca_agilent_truth
## The MATLAB test reads outputDir from the same default (or $AGILENT_TEST_TRUTH).

## If chromConverter was installed into a personal library that is not on the
## default search path, add the usual Windows location (harmless if absent).
userlib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R",
                     "win-library", paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = "."))
if (dir.exists(userlib)) .libPaths(c(userlib, .libPaths()))
suppressPackageStartupMessages(library(chromConverter))

## ---- argument / path resolution ------------------------------------------
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
}
args <- commandArgs(trailingOnly = TRUE)
script_dir  <- get_script_dir()
datasetRoot <- if (length(args) >= 1) args[1] else file.path(script_dir, "..", "datasets")
outdir      <- if (length(args) >= 2) args[2] else file.path(dirname(tempdir()), "ca_agilent_truth")

if (!dir.exists(datasetRoot)) {
  stop("Dataset root not found: ", normalizePath(datasetRoot, mustWork = FALSE),
       "\nPass it as the first argument, e.g. Rscript generate_groundtruth.R C:/my/data")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sanitize <- function(s) gsub("[^A-Za-z0-9]", "_", s)

## Find every .D directory under the dataset root.
all_dirs <- list.dirs(datasetRoot, recursive = TRUE)
d_folders <- all_dirs[grepl("\\.[Dd]$", all_dirs)]

ms_summary <- data.frame(key = character(), intensity_sum = numeric(),
                         intensity_max = numeric(), stringsAsFactors = FALSE)
n_written <- 0

for (dfolder in d_folders) {
  dkey <- sanitize(basename(dfolder))
  files <- list.files(dfolder, full.names = TRUE)

  for (p in files) {
    ext <- tolower(tools::file_ext(p))
    key <- paste0(dkey, "__", sanitize(basename(p)))
    tryCatch({
      if (ext == "uv") {
        x <- read_chemstation_uv(p, format_out = "matrix", scale = FALSE)
        out <- cbind(rt = as.numeric(rownames(x)), x)   # keep wavelength colnames
        write.csv(out, file.path(outdir, paste0("uv_", key, ".csv")), row.names = FALSE)
        n_written <- n_written + 1
      } else if (ext == "ch") {
        x <- read_chemstation_ch(p, format_out = "matrix", scale = FALSE)
        out <- data.frame(rt = as.numeric(rownames(x)), raw = as.numeric(x[, 1]))
        write.csv(out, file.path(outdir, paste0("ch_", key, ".csv")), row.names = FALSE)
        n_written <- n_written + 1
      } else if (ext == "ms") {
        x <- read_chemstation_ms(p)                      # list(MS1, BPC, TIC)
        tic <- as.data.frame(x$TIC); colnames(tic) <- c("rt", "tic")
        write.csv(tic, file.path(outdir, paste0("ms_", key, ".csv")), row.names = FALSE)
        ms_summary <- rbind(ms_summary, data.frame(
          key = key,
          intensity_sum = sum(x$MS1[, "intensity"]),
          intensity_max = max(x$MS1[, "intensity"])))
        n_written <- n_written + 1
      }
    }, error = function(e) {
      message("skipped ", basename(p), " in ", basename(dfolder), ": ", conditionMessage(e))
    })
  }
}
write.csv(ms_summary, file.path(outdir, "ms_summary.csv"), row.names = FALSE)

cat("wrote", n_written, "reference CSVs (+ ms_summary.csv) to", normalizePath(outdir), "\n")
