library(ciftiTools)
library(R.matlab)

ciftiTools.setOption("wb_path", "/workbench")

# =========================================================================
# SETTINGS
# =========================================================================

# -------------------------------------------------------------------------
# Input directories
# -------------------------------------------------------------------------
# PFN set A: TCA PFN effects
PFN_A_results_dir <- paste0(
  "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/",
  "22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/091426_TCA_vert_normed_loadings_FDR_pval_beta_logTC/"
)

# PFN set B: 22q PFN effects
PFN_B_results_dir <- paste0(
  "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/",
  "22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/091426_22q_logTC_vert_normed_loadings_FDR_pval_beta_subclass/"
)

# Output directory expected by the downstream MATLAB visualization script
VisualizeFolder <- paste0(
  "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/",
  "22q_Project/Analyses/Results/TCA_22q_PFN_conjunction/"
)

dir.create(
  VisualizeFolder,
  showWarnings = FALSE,
  recursive = TRUE
)

# -------------------------------------------------------------------------
# Surface files for clustering
# -------------------------------------------------------------------------
left_surf <- paste0(
  "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/",
  "22q_Project/Analyses/Inputs/tpl-fsLR_den-32k_hemi-L_midthickness.surf.gii"
)

right_surf <- paste0(
  "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/",
  "22q_Project/Analyses/Inputs/tpl-fsLR_den-32k_hemi-R_midthickness.surf.gii"
)

# Template CIFTI used only for temporary binary maps required by Workbench
PFNs_hardparcel <- read_cifti(
  paste0(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/",
    "22q_Project/hardparcel_group.dscalar.nii"
  )
)

# -------------------------------------------------------------------------
# Cortical dimensions / analysis settings
# -------------------------------------------------------------------------
n_lh <- 29696
n_rh <- 29716
n_total <- n_lh + n_rh
n_PFN <- 17
fdr_threshold <- 0.05
cluster_area_mm2 <- 50

# Winner-take-all is based ONLY on the 22q PFN effect map (PFN_B).
# Within positive conjunctions, the largest positive 22q beta wins.
# Within negative conjunctions, the most negative 22q beta wins.
# The winning TCA PFN (PFN_A) is retained as the paired identity.
wta_metric <- "22q_beta"

# -------------------------------------------------------------------------
# Workbench command
# -------------------------------------------------------------------------
wb_command <- file.path("/workbench", "wb_command")

if (!file.exists(wb_command)) {
  wb_command <- Sys.which("wb_command")
}

if (length(wb_command) == 0 || wb_command == "") {
  stop("Could not locate wb_command.")
}

# =========================================================================
# LOAD ALL PFN EFFECT MAPS
# =========================================================================

load_pfn_maps <- function(results_dir, filename_prefix, label) {
  
  beta_matrix <- matrix(
    NA_real_,
    nrow = n_total,
    ncol = n_PFN
  )
  
  fdr_matrix <- matrix(
    NA_real_,
    nrow = n_total,
    ncol = n_PFN
  )
  
  colnames(beta_matrix) <- paste0("PFN", 1:n_PFN)
  colnames(fdr_matrix) <- paste0("PFN", 1:n_PFN)
  
  for (PFN in 1:n_PFN) {
    
    infile <- file.path(
      results_dir,
      sprintf(filename_prefix, PFN)
    )
    
    if (!file.exists(infile)) {
      stop(
        paste0(
          label, " PFN ", PFN,
          " input file does not exist: ", infile
        )
      )
    }
    
    dat <- read.csv(infile)
    
    if (nrow(dat) < n_total) {
      stop(
        paste0(
          label, " PFN ", PFN,
          " has fewer than ", n_total, " vertices."
        )
      )
    }
    
    beta_matrix[, PFN] <- as.numeric(dat[1:n_total, 2])
    fdr_matrix[, PFN] <- as.numeric(dat[1:n_total, 3])
  }
  
  cat(
    label,
    ": loaded ", n_PFN, " PFN maps (",
    n_total, " cortical vertices each)\n",
    sep = ""
  )
  
  list(
    beta = beta_matrix,
    fdr = fdr_matrix
  )
}

PFN_A <- load_pfn_maps(
  results_dir = PFN_A_results_dir,
  filename_prefix = "logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",
  label = "PFN_A_TCA"
)

PFN_B <- load_pfn_maps(
  results_dir = PFN_B_results_dir,
  filename_prefix = "stat_22q_gam_vert_loadings_PFN_%02d_beta_fdr.csv",
  label = "PFN_B_22q"
)

# Basic dimension checks
stopifnot(
  nrow(PFN_A$beta) == n_total,
  nrow(PFN_B$beta) == n_total,
  ncol(PFN_A$beta) == n_PFN,
  ncol(PFN_B$beta) == n_PFN
)

# =========================================================================
# TEMPORARY WORKING DIRECTORY
# =========================================================================

temp_cluster_dir <- tempfile(
  pattern = "TCA_22q_PFN_conjunction_50mm_"
)

dir.create(
  temp_cluster_dir,
  recursive = TRUE
)

# Temporary Workbench files are removed at the end of the script.
# NOTE: Do not call the cleanup function here: doing so would delete
# temp_cluster_dir before the first PFN pair is processed.
on.exit_cleanup <- function() {
  if (exists("temp_cluster_dir") && dir.exists(temp_cluster_dir)) {
    unlink(
      temp_cluster_dir,
      recursive = TRUE,
      force = TRUE
    )
  }
}

# =========================================================================
# FUNCTION: PFN-PAIR CONJUNCTION -> CLUSTER FILTER -> WINNER-TAKE-ALL
# =========================================================================
#
# For each directional analysis:
#   1. Loop through matched PFN pairs only: PFN 1 x PFN 1, ..., PFN 17 x PFN 17.
#   2. Define significance for both matched PFNs using FDR < 0.05 and the
#      requested beta directions.
#   3. Intersect the two PFN masks: this is the matched PFN conjunction.
#   4. Apply the 50 mm2 surface-area cluster threshold to the conjunction.
#   5. Preserve each matched PFN pair that independently survives clustering.
#   6. At vertices where multiple matched PFN pairs survive, winner-take-all
#      is applied across PFN1 through PFN17, using the 22q effect.
#   7. Save only the winning 22q PFN ID for L/R.
#
# Important difference from the original scripts:
#   There is NO SA mask anywhere in this analysis.
#   Conjunction is strictly PFN_A x PFN_B.
#
# Winner metric: the signed 22q PFN beta (PFN_B$beta). This means that
# positive analyses choose the largest positive 22q beta and negative analyses
# choose the most negative 22q beta.
# =========================================================================

cluster_pfn_pairs <- function(
    PFN_A,
    PFN_B,
    A_direction = c("positive", "negative"),
    B_direction = c("positive", "negative"),
    sign_name,
    n_lh = 29696,
    n_rh = 29716,
    cluster_area_mm2 = 50,
    fdr_threshold = 0.05,
    wta_metric = "22q_beta"
) {
  
  A_direction <- match.arg(A_direction)
  B_direction <- match.arg(B_direction)
  n_total <- n_lh + n_rh
  n_PFN <- ncol(PFN_A$beta)
  
  cat(
    "\n========================================\n",
    "Processing: ", sign_name, "\n",
    "PFN_A direction: ", A_direction, "\n",
    "PFN_B direction: ", B_direction, "\n",
    "Cluster threshold: ", cluster_area_mm2, " mm2\n",
    "Number of matched PFN pairs: ", n_PFN, "\n",
    "WTA metric: ", wta_metric, "\n",
    "========================================\n",
    sep = ""
  )
  
  # -----------------------------------------------------------------------
  # Vertex-wise significance masks for each source set
  # -----------------------------------------------------------------------
  
  if (A_direction == "positive") {
    A_sig <-
      !is.na(PFN_A$beta) &
      !is.na(PFN_A$fdr) &
      PFN_A$fdr < fdr_threshold &
      PFN_A$beta > 0
  } else {
    A_sig <-
      !is.na(PFN_A$beta) &
      !is.na(PFN_A$fdr) &
      PFN_A$fdr < fdr_threshold &
      PFN_A$beta < 0
  }
  
  if (B_direction == "positive") {
    B_sig <-
      !is.na(PFN_B$beta) &
      !is.na(PFN_B$fdr) &
      PFN_B$fdr < fdr_threshold &
      PFN_B$beta > 0
  } else {
    B_sig <-
      !is.na(PFN_B$beta) &
      !is.na(PFN_B$fdr) &
      PFN_B$fdr < fdr_threshold &
      PFN_B$beta < 0
  }
  
  # -----------------------------------------------------------------------
  # Storage across all PFN pairs
  # -----------------------------------------------------------------------
  
  # Only matched PFNs are compared: A PFN1 x B PFN1, ..., A PFN17 x B PFN17.
  n_pairs <- n_PFN
  
  pair_A <- rep(NA_integer_, n_pairs)
  pair_B <- rep(NA_integer_, n_pairs)
  pair_label <- character(n_pairs)
  
  unclustered_membership <- matrix(
    FALSE,
    nrow = n_total,
    ncol = n_pairs
  )
  
  clustered_membership <- matrix(
    FALSE,
    nrow = n_total,
    ncol = n_pairs
  )
  
  # Strength used for WTA. The metric is stored for each pair/vertex only
  # after cluster filtering, so nonsurviving pairs can never win.
  clustered_wta_matrix <- matrix(
    NA_real_,
    nrow = n_total,
    ncol = n_pairs
  )
  
  col_idx <- 0L
  
  # -----------------------------------------------------------------------
  # PFN-A x PFN-B conjunctions and clustering
  # -----------------------------------------------------------------------
  
  for (PFN in 1:n_PFN) {
    
    A_PFN <- PFN
    B_PFN <- PFN
    
    col_idx <- col_idx + 1L
    pair_A[col_idx] <- A_PFN
    pair_B[col_idx] <- B_PFN
    pair_label[col_idx] <- paste0(
      "PFN", PFN, "__PFN", PFN
    )
    
    cat(
      "Pair ", col_idx, "/", n_pairs,
      " (A=PFN", A_PFN,
      ", B=PFN", B_PFN, ")\n",
      sep = ""
    )
    
    # PFN-pair conjunction: BOTH source PFNs must be significant at the
    # same cortical vertex in the requested directions.
    conjunction_mask <-
      A_sig[, A_PFN] &
      B_sig[, B_PFN]
    
    unclustered_membership[, col_idx] <- conjunction_mask
    
    n_unclust <- sum(conjunction_mask, na.rm = TRUE)
    
    cat(
      "  Unclustered conjunction vertices = ",
      n_unclust,
      "\n",
      sep = ""
    )
    
    if (n_unclust == 0) {
      next
    }
    
    # -------------------------------------------------------------------
    # Temporary binary conjunction CIFTI for this PFN pair
    # -------------------------------------------------------------------
    
    binary_vec <- as.numeric(conjunction_mask)
    
    binary_cifti <- PFNs_hardparcel
    
    binary_cifti$data$cortex_left <- matrix(
      binary_vec[1:n_lh],
      ncol = 1
    )
    
    binary_cifti$data$cortex_right <- matrix(
      binary_vec[(n_lh + 1):n_total],
      ncol = 1
    )
    
    binary_outbase <- file.path(
      temp_cluster_dir,
      paste0(sign_name, "_", pair_label[col_idx], "_binary")
    )
    
    write_cifti(
      binary_cifti,
      binary_outbase
    )
    
    binary_infile <- paste0(
      binary_outbase,
      ".dscalar.nii"
    )
    
    # -------------------------------------------------------------------
    # Apply the surface-area cluster threshold to this PFN pair
    # -------------------------------------------------------------------
    
    cluster_file <- file.path(
      temp_cluster_dir,
      paste0(sign_name, "_", pair_label[col_idx], "_cluster.dscalar.nii")
    )
    
    cmd <- paste0(
      '"', wb_command, '" -cifti-find-clusters ',
      '"', binary_infile, '" ',
      '1e-12 ', cluster_area_mm2, ' ',
      '0 0 ',
      'COLUMN ',
      '"', cluster_file, '" ',
      '-left-surface ',
      '"', left_surf, '" ',
      '-right-surface ',
      '"', right_surf, '"'
    )
    
    status <- system(cmd)
    
    if (status != 0) {
      stop(
        paste(
          "wb_command clustering failed for",
          sign_name,
          "A PFN",
          A_PFN,
          "x B PFN",
          B_PFN
        )
      )
    }
    
    # -------------------------------------------------------------------
    # Read surviving clusters
    # -------------------------------------------------------------------
    
    cluster_cifti <- read_cifti(cluster_file)
    
    keep_lh <-
      !is.na(cluster_cifti$data$cortex_left) &
      cluster_cifti$data$cortex_left > 0
    
    keep_rh <-
      !is.na(cluster_cifti$data$cortex_right) &
      cluster_cifti$data$cortex_right > 0
    
    keep <- c(
      as.vector(keep_lh),
      as.vector(keep_rh)
    )
    
    stopifnot(length(keep) == n_total)
    
    clustered_membership[, col_idx] <- keep
    
    # -------------------------------------------------------------------
    # WTA strength for this surviving pair
    # -------------------------------------------------------------------
    # WTA is based ONLY on the 22q PFN effect (PFN_B).
    # Because all surviving pairs in a directional analysis have the same
    # 22q sign, the signed beta directly gives the desired winner rule:
    #   positive  -> largest 22q beta wins
    #   negative  -> most negative 22q beta wins
    
    if (wta_metric == "22q_beta") {
      
      wta_strength <- PFN_B$beta[, B_PFN]
      
    } else {
      
      stop(
        paste(
          "Unsupported wta_metric:",
          wta_metric,
          ". Use '22q_beta'."
        )
      )
    }
    
    clustered_wta_matrix[keep, col_idx] <- wta_strength[keep]
    
    cat(
      "  Vertices surviving ",
      cluster_area_mm2,
      " mm2 = ",
      sum(keep),
      "\n",
      sep = ""
    )
    
    # Remove temporary PFN-pair CIFTIs immediately
    unlink(
      c(binary_infile, cluster_file),
      force = TRUE
    )
  }
  
  # -----------------------------------------------------------------------
  # Union across all PFN pairs
  # -----------------------------------------------------------------------
  
  union_unclustered <-
    rowSums(unclustered_membership) > 0
  
  union_clustered <-
    rowSums(clustered_membership) > 0
  
  n_union_unclustered <- sum(union_unclustered)
  n_union_clustered <- sum(union_clustered)
  
  pct_union_of_cortex_unclustered <-
    100 * n_union_unclustered / n_total
  
  pct_union_of_cortex_clustered <-
    100 * n_union_clustered / n_total
  
  # -----------------------------------------------------------------------
  # WINNER-TAKE-ALL ACROSS 22q PFNs AFTER PAIR-SPECIFIC CLUSTER FILTERING
  # -----------------------------------------------------------------------
  #
  # The TCA PFNs determine whether a given 22q PFN has at least one
  # cluster-surviving TCA x 22q conjunction at each vertex. Once that
  # eligibility step is complete, WTA is performed ONLY across 22q PFNs.
  #
  # Thus, there is one output identity per vertex:
  #   winner_22q_PFN_ID
  #
  # Positive 22q conjunctions: largest positive 22q beta wins.
  # Negative 22q conjunctions: most negative 22q beta wins.
  # Exact beta ties across different 22q PFNs favor the lower PFN number.
  # -----------------------------------------------------------------------
  
  # Because only matched PFN pairs are tested, each 22q PFN corresponds
  # to exactly one conjunction pair (same-numbered TCA PFN).
  eligible_22q <- clustered_membership
  
  winner_22q_PFN_ID <- rep(0L, n_total)
  winner_22q_beta <- rep(0, n_total)
  
  surviving_vertices <- which(rowSums(eligible_22q) > 0)
  
  for (vert in surviving_vertices) {
    
    surviving_22q <- which(eligible_22q[vert, ])
    vals <- PFN_B$beta[vert, surviving_22q]
    
    if (B_direction == "positive") {
      winner_local <- which.max(vals)
    } else {
      winner_local <- which.min(vals)
    }
    
    winner_B <- surviving_22q[winner_local]
    
    winner_22q_PFN_ID[vert] <- winner_B
    winner_22q_beta[vert] <- PFN_B$beta[vert, winner_B]
  }
  
  stopifnot(
    sum(winner_22q_PFN_ID > 0) == n_union_clustered
  )
  
  # -----------------------------------------------------------------------
  # Per-pair summary counts
  # -----------------------------------------------------------------------
  
  n_unclust_by_pair <- colSums(unclustered_membership)
  n_clustered_by_pair <- colSums(clustered_membership)
  
  n_winner_by_22q_PFN <- tabulate(
    winner_22q_PFN_ID[winner_22q_PFN_ID > 0],
    nbins = n_PFN
  )
  
  n_multi_pair_50mm <- sum(
    rowSums(clustered_membership) > 1
  )
  
  max_pairs_at_vertex_50mm <- if (n_union_clustered > 0) {
    max(rowSums(clustered_membership))
  } else {
    0
  }
  
  # -----------------------------------------------------------------------
  # Write ONLY the MATLAB-required winner files: L/R for this direction
  # -----------------------------------------------------------------------
  
  # Keep the MAT output simple: one PFN-ID vector per hemisphere.
  # The value is the winning 22q PFN (1:17); 0 means no surviving conjunction.
  winner_lh_fields <- list(
    winner_22q_PFN_ID = as.numeric(winner_22q_PFN_ID[1:n_lh])
  )
  
  winner_rh_fields <- list(
    winner_22q_PFN_ID = as.numeric(winner_22q_PFN_ID[(n_lh + 1):n_total])
  )
  
  map_base <- paste0(
    "TCAx22q_",
    sign_name,
    "_22q_PFN_ID_clust_50mm"
  )
  
  writeMat(
    file.path(
      VisualizeFolder,
      paste0(map_base, "_L.mat")
    ),
    x = winner_lh_fields
  )
  
  writeMat(
    file.path(
      VisualizeFolder,
      paste0(map_base, "_R.mat")
    ),
    x = winner_rh_fields
  )
  
  # -----------------------------------------------------------------------
  # Summary table: one row per matched PFN pair
  # -----------------------------------------------------------------------
  
  data.frame(
    conjunction = sign_name,
    PFN_A_direction = A_direction,
    PFN_B_direction = B_direction,
    PFN_A = pair_A,
    PFN_B = pair_B,
    pair_ID = seq_len(n_pairs),
    n_union_unclustered = n_union_unclustered,
    pct_union_of_cortex_unclustered = pct_union_of_cortex_unclustered,
    n_union_clustered_50mm = n_union_clustered,
    pct_union_of_cortex_clustered_50mm = pct_union_of_cortex_clustered,
    n_vertices_multiPFNpair_50mm = n_multi_pair_50mm,
    max_PFNpairs_at_vertex_50mm = max_pairs_at_vertex_50mm,
    n_vertices_unclustered_pair = as.integer(n_unclust_by_pair),
    n_vertices_50mm_pair = as.integer(n_clustered_by_pair),
    n_vertices_50mm_winner_22q_PFN = as.integer(n_winner_by_22q_PFN[pair_B])
  )
}

# =========================================================================
# RUN FOUR PFN-PAIR CONJUNCTION ANALYSES
# =========================================================================

# 1. Convergent positive:
#    PFN_A positive + PFN_B positive
con_pos_50mm <- cluster_pfn_pairs(
  PFN_A = PFN_A,
  PFN_B = PFN_B,
  A_direction = "positive",
  B_direction = "positive",
  sign_name = "con_pos_pos",
  cluster_area_mm2 = cluster_area_mm2,
  fdr_threshold = fdr_threshold,
  wta_metric = wta_metric
)

# 2. Convergent negative:
#    PFN_A negative + PFN_B negative
con_neg_50mm <- cluster_pfn_pairs(
  PFN_A = PFN_A,
  PFN_B = PFN_B,
  A_direction = "negative",
  B_direction = "negative",
  sign_name = "con_neg_neg",
  cluster_area_mm2 = cluster_area_mm2,
  fdr_threshold = fdr_threshold,
  wta_metric = wta_metric
)

# 3. Divergent:
#    PFN_A positive + PFN_B negative
div_Apos_Bneg_50mm <- cluster_pfn_pairs(
  PFN_A = PFN_A,
  PFN_B = PFN_B,
  A_direction = "positive",
  B_direction = "negative",
  sign_name = "div_pos_neg",
  cluster_area_mm2 = cluster_area_mm2,
  fdr_threshold = fdr_threshold,
  wta_metric = wta_metric
)

# 4. Divergent:
#    PFN_A negative + PFN_B positive
div_Aneg_Bpos_50mm <- cluster_pfn_pairs(
  PFN_A = PFN_A,
  PFN_B = PFN_B,
  A_direction = "negative",
  B_direction = "positive",
  sign_name = "div_neg_pos",
  cluster_area_mm2 = cluster_area_mm2,
  fdr_threshold = fdr_threshold,
  wta_metric = wta_metric
)

# =========================================================================
# WRITE SINGLE SUMMARY CSV
# =========================================================================

summary_df <- rbind(
  con_pos_50mm,
  con_neg_50mm,
  div_Apos_Bneg_50mm,
  div_Aneg_Bpos_50mm
)

summary_df$pct_union_of_cortex_unclustered <- round(
  summary_df$pct_union_of_cortex_unclustered,
  3
)

summary_df$pct_union_of_cortex_clustered_50mm <- round(
  summary_df$pct_union_of_cortex_clustered_50mm,
  3
)

write.csv(
  summary_df,
  file.path(
    VisualizeFolder,
    "TCA_22q_PFN_pairwise_conjunction_50mm_summary.csv"
  ),
  row.names = FALSE
)

# =========================================================================
# CLEAN UP TEMPORARY CIFTI FILES
# =========================================================================

on.exit_cleanup()


# =========================================================================
# FINAL CONSOLE SUMMARY
# =========================================================================

cat(
  "\n========================================\n",
  "50 mm2 PFN-PAIR CONJUNCTION ANALYSIS COMPLETE\n",
  "========================================\n",
  "\nOutputs written to:\n",
  VisualizeFolder,
  "\n\nCreated:\n",
  "  8 MATLAB 22q-PFN winner files (L/R for 4 conjunctions)\n",
  "  1 summary CSV (17 matched PFN rows per conjunction)\n",
  "\nWinner metric: ",
  wta_metric,
  "\n",
  sep = ""
)
