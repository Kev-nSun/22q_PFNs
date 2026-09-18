suppressPackageStartupMessages({
  library(optparse)
  library(R.matlab)
  library(ciftiTools)
})

## ------------------------------------------------------------
## 1. Parse command-line options
## ------------------------------------------------------------
option_list <- list(
  make_option(
    c("--subject_list"), type="character", help="Text file with one subject ID per line"
  ),
  make_option(
    c("--in_dir"), type="character", help="Directory containing sub*/FN.mat folders"
  ),
  make_option(
    c("--out_dir"), type="character", help="Output directory root (per subject folders created)"
  ),
  make_option(
    c("--array_id"), type="integer", help="Array ID in slurm array (based on index in subject_list)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check required arguments
required <- c("subject_list", "in_dir", "out_dir", "array_id")
missing <- required[sapply(opt[required], function(x) is.null(x) || is.na(x))]

if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

## ------------------------------------------------------------
## 2. Load inputs
## ------------------------------------------------------------

# read subject list
subjects <- readLines(opt$subject_list)
subjects <- subjects[nzchar(subjects)]  # remove blank lines

ciftiTools.setOption("wb_path", "/cbica/projects/bbl_22q/software/workbench")

cortex_template <- read_cifti('/cbica/projects/bbl_22q/analysis/allometry/inputs/hardparcel_group.dscalar.nii',
                              brainstructures = c("left", "right"))

# Expected hemisphere vertex counts (from your original script)
n_vert_lh <- 29696
n_vert_rh <- 29716
n_vert_total <- n_vert_lh + n_vert_rh

## ------------------------------------------------------------
## 3. Subject processing function
## ------------------------------------------------------------
process_subject <- function(subj_id) {
  message("---- Processing ", subj_id, " ----")
  
  fn_mat_path <- file.path(opt$in_dir, subj_id, "FN.mat")
  
  if (!file.exists(fn_mat_path)) {
    warning("FN.mat not found for subject ", subj_id, ": ", fn_mat_path)
    return(invisible(NULL))
  }
  
  out_dir <- file.path(opt$out_dir, subj_id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Read FN.mat
  readin <- readMat(fn_mat_path)
  pfn_loadings <- as.matrix(readin$FN)
  
  if (nrow(pfn_loadings) != n_vert_total) {
    warning("Vertex count mismatch for ", subj_id, 
         ": found ", nrow(pfn_loadings), " expected ", n_vert_total)
    return(invisible(NULL))
  }
  
  # Normalize loadings to sum to 1 per vertex
  vertex_sums <- rowSums(pfn_loadings)
  pfn_loadings_norm <- pfn_loadings / vertex_sums
  pfn_loadings_norm[!is.finite(pfn_loadings_norm)] <- 0
  
  # Save pfn_loadings_norm matrix
  write.csv(
    pfn_loadings_norm,
    file = file.path(out_dir, paste0(subj_id, "_PFN_loadings_normed.csv")),
    row.names = FALSE
  )
  
  # Write out per-network CIFTIs
  n_nets <- ncol(pfn_loadings_norm)
  
  for (net in seq_len(n_nets)) {
    lh_vals <- matrix(pfn_loadings_norm[1:n_vert_lh, net], n_vert_lh, 1)
    rh_vals <- matrix(pfn_loadings_norm[(n_vert_lh+1):n_vert_total, net], n_vert_rh, 1)
    
    Soft_map <- cortex_template
    Soft_map$data$cortex_left  <- lh_vals
    Soft_map$data$cortex_right <- rh_vals
    
    outfile <- file.path(
      out_dir,
      sprintf("PFN%02d_soft_parcel_normed", net)
    )
    
    write_cifti(Soft_map, outfile)
  }
  
  message("---- Finished ", subj_id, " ----")
}

## ------------------------------------------------------------
## 4. Loop over subjects (or one if job array)
## ------------------------------------------------------------
array_id <- opt$array_id

# Safety: ensure array_id is in range
if (array_id < 1 || array_id > length(subjects)) {
  stop("array_id out of range: ", array_id,
       " (must be between 1 and ", length(subjects), ")")
}

# Select just this subject
subjects <- subjects[array_id]

for (s in subjects) {
  process_subject(s)
}