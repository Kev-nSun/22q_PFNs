#!/usr/bin/env Rscript
# This script runs vertex-level gams of the impact of 22q on PFN loading weighted surface areas of each vertex
# There are 235,355 non all-zero loadings across all PFN loadings (59,412 x 17)
# These models are chunked into groups of 1000
# This script also gives model-fit predictions based on 22q status

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(mgcv)
  library(data.table)
})

# ---------------------------
# CLI args
# ---------------------------
option_list <- list(
  make_option("--h5", type="character", help="Path to compiled HDF5 file (contains dataset Y and subject id)"),
  make_option("--csv", type="character", help="Path to 22q_PNC_matched_all.csv"),
  make_option("--outdir", type="character", help="Output directory"),
  make_option("--chunk-id", type="integer", help="0-based chunk index"),
  make_option("--chunk-size", type="integer", default=1000, help="Number of features per chunk")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- Coerce + validate scalar CLI args ----
chunk_id   <- as.integer(opt$`chunk-id`)
chunk_size <- as.integer(opt$`chunk-size`)

check_scalar <- function(x, name) {
  if (length(x) != 1 || is.na(x)) {
    stop(sprintf("%s must be a single value. Got: %s", name, paste(x, collapse=",")))
  }
}

check_scalar(chunk_id,   "--chunk-id")
check_scalar(chunk_size, "--chunk-size")

if (is.null(opt$h5) || is.null(opt$csv) ||
    is.null(opt$outdir) || is.null(opt$`chunk-id`)) {
  stop("Must supply --h5, --csv, --outdir, --chunk-id")
}

# ---------------------------
# Read in covariates
# ---------------------------
covar <- fread(opt$csv, colClasses = list(character = "subject"))

required_cols <- c(
  "subject",
  "stat_22q",
  "sex",
  "euler",
  "median_fd_combined",
  "age",
  "subclass"
)

missing_cols <- setdiff(required_cols, names(covar))
if (length(missing_cols) > 0) {
  stop(paste0(
    "22q_PNC_matched_all.csv missing required columns: ",
    paste(missing_cols, collapse=", ")
  ))
}

# Function for digits-only enforcement for csv list
clean_digits_only <- function(x) {
  s <- trimws(as.character(x))
  s <- sub("\\.0$", "", s)
  if (!grepl("^[0-9]+$", s)) {
    stop(paste0("Bad subject id: ", x, " -> ", s))
  }
  return(s)
}

covar[, subj_digits := vapply(subject, clean_digits_only, character(1))]

# ---------------------------
# Read H5 subject order
# ---------------------------
h5_ids <- as.character(h5read(opt$h5, "subject"))

# Convert "sub-00997818717" -> "997818717"
h5_to_digits <- function(s) {
  s2 <- sub("^sub-", "", s)
  s2 <- sub("^0+", "", s2)
  if (nchar(s2) == 0) s2 <- "0"
  return(s2)
}

h5_digits <- vapply(h5_ids, h5_to_digits, character(1))

# Align all_phenos rows to H5 rows
covar_index <- match(h5_digits, covar$subj_digits)

if (anyNA(covar_index)) {
  n_missing <- sum(is.na(covar_index))
  missing_examples <- head(h5_ids[is.na(covar_index)], 10)
  stop(paste0(
    "Failed to match ", n_missing, " H5 subjects to csv covariates. Examples: ",
    paste(missing_examples, collapse=", ")
  ))
}

all_phenos <- covar[covar_index]

# ---------------------------
# Precompute variables used in every model
# ---------------------------
all_phenos[, stat_22q := factor(stat_22q, levels=c(0,1))]
all_phenos[, isMale := as.integer(sex == "Male")]
all_phenos[, subclass := factor(subclass)]

subclass_ref <- levels(all_phenos$subclass)[1]

# Covariate reference values
isMale_ref <- mean(all_phenos$isMale, na.rm=TRUE)
age_ref    <- median(all_phenos$age, na.rm=TRUE)
euler_ref  <- median(all_phenos$euler, na.rm=TRUE)
fd_ref     <- median(all_phenos$median_fd_combined, na.rm=TRUE)

# Function to construct reference data for prediction for each 22q status
make_newdata <- function() {
  status <- c(0, 1)
  data.frame(
    stat_22q = factor(status, levels = c(0, 1)),
    isMale = rep(isMale_ref, length(status)),
    age = rep(age_ref, length(status)),
    euler = rep(euler_ref, length(status)),
    median_fd_combined = rep(fd_ref, length(status)),
    subclass = factor(
      rep(subclass_ref, length(status)),
      levels = levels(all_phenos$subclass)
    )
  )
}

# ---------------------------
# Determine chunk feature range
# ---------------------------
total_len <- as.integer(h5read(opt$h5, "total_len"))
chunk_id <- as.integer(chunk_id)
chunk_size <- as.integer(chunk_size)

start <- as.integer(chunk_id * chunk_size + 1L)
end   <- as.integer(min(total_len, (chunk_id + 1L) * chunk_size))
feat_n <- end - start + 1L

if (feat_n <= 0) stop("Chunk has no features.")

idx_cols <- seq.int(start, end)

# ---------------------------
# Read feature identifiers for this chunk
# ---------------------------
features_idx_all <- h5read(opt$h5, "features_idx")
net_id_all       <- h5read(opt$h5, "net_id")
vertex_id_all    <- h5read(opt$h5, "vertex_id")

features_idx_chunk <- features_idx_all[idx_cols]
net_id_chunk       <- net_id_all[idx_cols]
vertex_id_chunk    <- vertex_id_all[idx_cols]

# ---------------------------
# Read outcome block
# ---------------------------
n_h5 <- length(h5_ids)
stopifnot(nrow(all_phenos) == n_h5)

message("chunk_id = ", chunk_id)
message("chunk_size = ", chunk_size)
message("total_len = ", total_len)
message("start = ", start, " (", class(start), ")")
message("end = ", end, " (", class(end), ")")
message("n_subs = ", nrow(all_phenos))

h5info <- h5ls(opt$h5)
Yrow <- h5info[h5info$name=="Y", ]
message("H5 Y dims: ", paste(Yrow$dim, collapse=" "))

# READ IN VERTEX OUTCOMES "Y"
Y_block <- h5read(opt$h5, "Y", index = list(idx_cols, seq_len(n_h5)))
Y_block <- t(Y_block)

stopifnot(dim(Y_block)[1] == n_h5)
stopifnot(dim(Y_block)[2] == length(idx_cols))

# ---------------------------
# Fit models feature-by-feature
# ---------------------------
n <- nrow(all_phenos)

beta_stat_22q <- rep(NA_real_, feat_n)
se_stat_22q   <- rep(NA_real_, feat_n)
p_stat_22q    <- rep(NA_real_, feat_n)

edf_age <- rep(NA_real_, feat_n)
p_age   <- rep(NA_real_, feat_n)

# Predictions mat
pred_mat <- matrix(NA_real_, nrow=feat_n, ncol=2)

# Covariate model frame
cov_template <- as.data.frame(
  all_phenos[, .(
    stat_22q,
    isMale,
    euler,
    median_fd_combined,
    age,
    subclass
  )]
)

# Check number of subjects
n_cov_template <- nrow(cov_template)
if (n_cov_template < 200) stop("Too few subjects.")

# Fit function for one feature using covariate model frame
fit_one <- function(y) {
  dat2 <- cov_template
  dat2$Y <- y
  
  mod <- bam(
    Y ~ stat_22q +
      isMale +
      euler +
      median_fd_combined +
      s(age, k=4) +
      s(subclass, bs="re"),
    data = dat2,
    method = "fREML",
    discrete = TRUE
  )
  
  return(mod)
}

# Loop through Y_block
for (feat in seq_len(feat_n)) {
  y <- as.numeric(Y_block[, feat])
  mod <- tryCatch(fit_one(y), error=function(e) NULL)
  
  if (is.null(mod)) next
  
  # Parametric term for stat_22q
  mod_sum <- summary(mod)
  
  ptab <- mod_sum$p.table
  if ("stat_22q1" %in% rownames(ptab)) {
    beta_stat_22q[feat] <- ptab["stat_22q1", "Estimate"]
    se_stat_22q[feat]   <- ptab["stat_22q1", "Std. Error"]
    p_stat_22q[feat]    <- ptab["stat_22q1", "Pr(>|t|)"]
  }
  
  # Smooth term for age
  stab <- mod_sum$s.table
  sname <- rownames(stab)
  idx_s <- which(grepl("^s\\(age\\)", sname))
  
  if (length(idx_s) == 1) {
    edf_age[feat] <- stab[idx_s, "edf"]
    p_age[feat]   <- stab[idx_s, "p-value"]
  }
  
  # Predictions at reference covariate values, excluding subclass random effect
  nd <- make_newdata()
  
  pred_mat[feat, ] <- predict(
    mod,
    newdata = nd,
    type = "response",
    exclude = "s(subclass)"
  )
  
  if (feat %% 100 == 0) {
    message(sprintf("  fit %d / %d", feat, feat_n))
  }
}

# ---------------------------
# Save outputs
# ---------------------------
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

out_file <- file.path(
  opt$outdir,
  sprintf("gam_chunk_%04d_feat_%d_%d.rds", chunk_id, start, end)
)

out <- list(
  chunk_id = chunk_id,
  start = start,
  end = end,
  chunk_size = chunk_size,
  
  # --- feature identifiers ---
  feature_col = idx_cols,
  features_idx = features_idx_chunk,
  net_id = net_id_chunk,
  vertex_id = vertex_id_chunk,
  
  # --- prediction info ---
  isMale_ref = isMale_ref,
  age_ref = age_ref,
  euler_ref = euler_ref,
  fd_ref = fd_ref,
  
  # --- model outputs ---
  beta_stat_22q = beta_stat_22q,
  se_stat_22q = se_stat_22q,
  p_stat_22q = p_stat_22q,
  edf_age = edf_age,
  p_age = p_age,
  pred = pred_mat
)

saveRDS(out, out_file, compress = "xz")
message(sprintf("Saved: %s", out_file))