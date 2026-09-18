
# Rscript 22q_compile_PFN_loadings_gam_chunks_into_hdf5.R \
# --indir "/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_22q_vert_gam_results_nonstd_sex_avg_logTC_subclass/" \
# --out-h5 "/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_22q_vert_gam_results_nonstd_sex_avg_logTC_subclass.h5" \
# --total-len 235355 \
# --compression-level 4 \
# --chunk-rows 5000
#
# Compile per-chunk GAM outputs (.rds) into a single HDF5 file.
# Writes incrementally so we do NOT keep full pred matrices in RAM.

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
})

option_list <- list(
  make_option("--indir", type="character",
              help="Directory containing gam_chunk_*.rds files"),
  make_option("--out-h5", type="character",
              help="Output HDF5 path"),
  make_option("--total-len", type="integer", default=235355,
              help="Total number of masked features (default: 235355)"),
  make_option("--pattern", type="character",
              default="^gam_chunk_\\d+_feat_\\d+_\\d+\\.rds$",
              help="Regex pattern to match chunk files"),
  make_option("--overwrite", action="store_true", default=FALSE,
              help="If set, overwrite existing output HDF5"),
  make_option("--compression-level", type="integer", default=4,
              help="Gzip compression level (0-9). Default 4."),
  make_option("--chunk-rows", type="integer", default=5000,
              help="Row chunk size for pred matrices in HDF5. Default 5000."),
  make_option("--strict", action="store_true", default=FALSE,
              help="If set, error on overlaps or missing coverage")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$indir) || is.null(opt$`out-h5`)) {
  stop("Must supply --indir and --out-h5")
}

INDIR <- opt$indir
OUT_H5 <- opt$`out-h5`
TOTAL_LEN <- as.integer(opt$`total-len`)
PATTERN <- opt$pattern
OVERWRITE <- isTRUE(opt$overwrite)
GZIP_LEVEL <- as.integer(opt$`compression-level`)
CHUNK_ROWS <- as.integer(opt$`chunk-rows`)
STRICT <- isTRUE(opt$strict)

if (GZIP_LEVEL < 0 || GZIP_LEVEL > 9) stop("--compression-level must be 0..9")
if (CHUNK_ROWS < 1) stop("--chunk-rows must be >= 1")

# --------------------------
# Discover chunk files
# --------------------------
files <- list.files(INDIR, pattern=PATTERN, full.names=TRUE)
if (length(files) == 0) stop(sprintf("No chunk files found in %s matching %s", INDIR, PATTERN))

# Sort by chunk_id parsed from filename
get_chunk_id <- function(f) {
  m <- regexec("gam_chunk_([0-9]+)_feat_", basename(f))
  g <- regmatches(basename(f), m)[[1]]
  if (length(g) < 2) return(NA_integer_)
  as.integer(g[2])
}
chunk_ids <- vapply(files, get_chunk_id, integer(1))
ord <- order(chunk_ids)
files <- files[ord]
chunk_ids <- chunk_ids[ord]

cat(sprintf("Found %d chunk files under %s\n", length(files), INDIR))

# --------------------------
# Load first file to get schema
# --------------------------
first <- readRDS(files[1])

required_first <- c("pred",
                    "age_ref", "euler_ref", "fd_ref", "logTC_ref", "isMale_ref",
                    "chunk_id", "start", "end",
                    "features_idx", "net_id", "vertex_id", "feature_col",
                    "beta_stat_22q", "se_stat_22q", "p_stat_22q", "edf_age", "p_age")
missing_first <- setdiff(required_first, names(first))
if (length(missing_first) > 0) {
  stop(sprintf("First chunk missing required fields: %s", paste(missing_first, collapse=", ")))
}

# --------------------------
# Prepare output HDF5
# --------------------------
if (file.exists(OUT_H5)) {
  if (!OVERWRITE) {
    stop(sprintf("Output HDF5 already exists: %s (use --overwrite to replace)", OUT_H5))
  }
  file.remove(OUT_H5)
}

dir.create(dirname(OUT_H5), recursive=TRUE, showWarnings=FALSE)

h5createFile(OUT_H5)

# Helper to create 1D datasets
create_vec <- function(name, dtype) {
  h5createDataset(
    file = OUT_H5,
    dataset = name,
    dims = c(TOTAL_LEN),
    storage.mode = dtype,
    chunk = c(min(TOTAL_LEN, 100000)),
    level = GZIP_LEVEL
  )
}

# Helper to create 2D prediction datasets
create_pred <- function(name) {
  # dims: TOTAL_LEN x 2
  # chunk: CHUNK_ROWS x 2
  h5createDataset(
    file = OUT_H5,
    dataset = name,
    dims = c(TOTAL_LEN, 2),
    storage.mode = "double",    # store as double; you could use "single" if desired
    chunk = c(min(TOTAL_LEN, CHUNK_ROWS), 2),
    level = GZIP_LEVEL
  )
}

# Create groups
h5createGroup(OUT_H5, "meta")
h5createGroup(OUT_H5, "pred")

# Meta datasets (small)
h5write(c("stat_22q0", "stat_22q1"), OUT_H5, "meta/pred_colnames")
h5write(as.numeric(first$age_ref), OUT_H5, "meta/age_ref")
h5write(as.numeric(first$euler_ref), OUT_H5, "meta/euler_ref")
h5write(as.numeric(first$fd_ref), OUT_H5, "meta/fd_ref")
h5write(as.numeric(first$logTC_ref), OUT_H5, "meta/logTC_ref")
h5write(as.numeric(first$isMale_ref), OUT_H5, "meta/isMale_ref")
h5write(as.integer(TOTAL_LEN), OUT_H5, "meta/total_len")

# Create identifier + stats datasets
create_vec("feature_col", "integer")
create_vec("features_idx", "integer")   # 0-based
create_vec("net_id", "integer")         # 1..17
create_vec("vertex_id", "integer")      # 1..59412

create_vec("beta_stat_22q", "double")
create_vec("se_stat_22q",   "double")
create_vec("p_stat_22q",    "double")
create_vec("edf_age",    "double")
create_vec("p_age",      "double")

# Create prediction dataset
create_pred("pred/pred")


# Coverage tracker
filled <- rep(FALSE, TOTAL_LEN)


# --------------------------
# Write each chunk into HDF5
# --------------------------
for (k in seq_along(files)) {
  f <- files[k]
  obj <- readRDS(f)
  
  start <- as.integer(obj$start)
  end   <- as.integer(obj$end)
  if (start < 1 || end > TOTAL_LEN || start > end) {
    stop(sprintf("Bad start/end in %s: %d..%d (TOTAL_LEN=%d)", f, start, end, TOTAL_LEN))
  }
  
  idx <- start:end
  p <- length(idx)
  
  if (any(filled[idx])) {
    msg <- sprintf("Overlap detected for indices %d..%d in %s", start, end, basename(f))
    if (STRICT) stop(msg) else warning(msg)
  }
  
  # --- write identifiers ---
  h5write(as.integer(obj$feature_col),  OUT_H5, "feature_col",  index=list(idx))
  h5write(as.integer(obj$features_idx), OUT_H5, "features_idx", index=list(idx))
  h5write(as.integer(obj$net_id),       OUT_H5, "net_id",       index=list(idx))
  h5write(as.integer(obj$vertex_id),    OUT_H5, "vertex_id",    index=list(idx))
  
  # --- write stats ---
  h5write(as.numeric(obj$beta_stat_22q), OUT_H5, "beta_stat_22q", index=list(idx))
  h5write(as.numeric(obj$se_stat_22q),   OUT_H5, "se_stat_22q",   index=list(idx))
  h5write(as.numeric(obj$p_stat_22q),    OUT_H5, "p_stat_22q",    index=list(idx))
  h5write(as.numeric(obj$edf_age),    OUT_H5, "edf_age",    index=list(idx))
  h5write(as.numeric(obj$p_age),      OUT_H5, "p_age",      index=list(idx))
  
  # --- write predictions ---
  mat <- obj$pred
  # h5write expects index=list(rows, cols)
  if (!is.matrix(mat) || nrow(mat) != p || ncol(mat) != 2) {
    stop(sprintf("Pred matrix mismatch in %s: got %s expected %dx2",
                 basename(f),
                 paste(dim(mat), collapse="x"),
                 p))
  }
  h5write(mat, OUT_H5, "pred/pred", index=list(idx, 1:2))
  
  filled[idx] <- TRUE
  
  if (k %% 10 == 0 || k == length(files)) {
    cat(sprintf("Wrote %d / %d chunks\n", k, length(files)))
  }
}

# --------------------------
# Final coverage check
# --------------------------
n_missing <- sum(!filled)
if (n_missing > 0) {
  msg <- sprintf("WARNING: %d / %d features never filled (missing chunk files?)", n_missing, TOTAL_LEN)
  if (STRICT) stop(msg) else warning(msg)
}

# Save coverage summary
h5write(as.integer(which(!filled)), OUT_H5, "meta/missing_feature_rows")

cat(sprintf("Done. Wrote HDF5: %s\n", OUT_H5))
cat(sprintf("Missing feature rows: %d\n", n_missing))

# close all
H5close()
