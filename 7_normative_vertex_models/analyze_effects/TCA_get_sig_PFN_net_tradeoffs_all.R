library(readr)

dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/NORMED/Chord_diag/"
dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

indir <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/091426_TCA_vert_normed_loadings_FDR_pval_beta_logTC/"

PFN_names <- paste0("PFN", 1:17)
all_pfns <- 1:17

#---------------------------
# Pre-load all PFN data
#---------------------------
pfn_results <- vector("list", 17)

for (PFN in all_pfns) {
  pfn_results[[PFN]] <- read.csv(
    sprintf(
      "%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",
      indir,
      PFN
    )
  )
}

#---------------------------
# Output matrices / vectors
#---------------------------

# Mean target beta within directional significant reference mask
mean_beta_mat <- matrix(
  NA_real_,
  nrow = 17,
  ncol = 17,
  dimnames = list(PFN_names, PFN_names)
)

# Number of non-missing target betas within reference mask
n_overlap_mat <- matrix(
  NA_integer_,
  nrow = 17,
  ncol = 17,
  dimnames = list(PFN_names, PFN_names)
)

# Number of target vertices that are FDR significant
# within the reference mask
n_target_sig_mat <- matrix(
  NA_integer_,
  nrow = 17,
  ncol = 17,
  dimnames = list(PFN_names, PFN_names)
)

# Number of target vertices that are:
#   1. FDR significant, AND
#   2. opposite in sign to the reference PFN
n_target_sig_opposite_mat <- matrix(
  NA_integer_,
  nrow = 17,
  ncol = 17,
  dimnames = list(PFN_names, PFN_names)
)

# Proportion of significant target vertices that are
# opposite in sign to the reference PFN
prop_target_sig_opposite_mat <- matrix(
  NA_real_,
  nrow = 17,
  ncol = 17,
  dimnames = list(PFN_names, PFN_names)
)

# Reference-network summaries
n_vert_vec <- rep(NA_integer_, 17)
direction_vec <- rep(NA_character_, 17)
mean_sig_beta_vec <- rep(NA_real_, 17)

#---------------------------
# Loop over all 17 reference PFNs
#---------------------------
for (ref_pfn in all_pfns) {
  
  ref_dat <- pfn_results[[ref_pfn]]
  
  #-------------------------------------------------------
  # Step 1:
  # Find all significant vertices in the reference PFN
  #-------------------------------------------------------
  sig_mask <- !is.na(ref_dat[, 3]) &
    ref_dat[, 3] < 0.05 &
    !is.na(ref_dat[, 2])
  
  n_sig <- sum(sig_mask)
  
  if (n_sig > 0) {
    
    #-------------------------------------------------------
    # Step 2:
    # Determine overall reference-network direction
    # using mean beta among ALL significant vertices
    #-------------------------------------------------------
    mean_sig_beta <- mean(
      ref_dat[sig_mask, 2],
      na.rm = TRUE
    )
    
    mean_sig_beta_vec[ref_pfn] <- mean_sig_beta
    
    #-------------------------------------------------------
    # Step 3:
    # Build directional significant reference mask
    #-------------------------------------------------------
    if (mean_sig_beta > 0) {
      
      direction_vec[ref_pfn] <- "positive"
      
      mask <- sig_mask &
        ref_dat[, 2] > 0
      
    } else if (mean_sig_beta < 0) {
      
      direction_vec[ref_pfn] <- "negative"
      
      mask <- sig_mask &
        ref_dat[, 2] < 0
      
    } else {
      
      direction_vec[ref_pfn] <- "zero"
      
      # No meaningful directional mask if mean beta == 0
      mask <- rep(FALSE, nrow(ref_dat))
    }
    
    # Number of vertices in the final directional
    # significant reference mask
    n_vert_vec[ref_pfn] <- sum(mask)
    
    #-------------------------------------------------------
    # Step 4:
    # Examine all 17 target PFNs within that mask
    #-------------------------------------------------------
    for (other_pfn in all_pfns) {
      
      other_dat <- pfn_results[[other_pfn]]
      
      #-----------------------------------------------------
      # All usable target vertices within reference mask
      #-----------------------------------------------------
      overlap_mask <- mask &
        !is.na(other_dat[, 2])
      
      n_overlap <- sum(overlap_mask)
      
      n_overlap_mat[ref_pfn, other_pfn] <- n_overlap
      
      # Mean target beta across ALL usable vertices
      # in the reference mask
      if (n_overlap > 0) {
        
        mean_beta_mat[ref_pfn, other_pfn] <- mean(
          other_dat[overlap_mask, 2],
          na.rm = TRUE
        )
      }
      
      #-----------------------------------------------------
      # Target-significant vertices within reference mask
      #-----------------------------------------------------
      target_sig_mask <- mask &
        !is.na(other_dat[, 2]) &
        !is.na(other_dat[, 3]) &
        other_dat[, 3] < 0.05
      
      n_target_sig <- sum(target_sig_mask)
      
      n_target_sig_mat[ref_pfn, other_pfn] <- n_target_sig
      
      #-----------------------------------------------------
      # Significant target vertices in OPPOSITE direction
      #-----------------------------------------------------
      if (direction_vec[ref_pfn] == "positive") {
        
        target_sig_opposite_mask <- target_sig_mask &
          other_dat[, 2] < 0
        
      } else if (direction_vec[ref_pfn] == "negative") {
        
        target_sig_opposite_mask <- target_sig_mask &
          other_dat[, 2] > 0
        
      } else {
        
        target_sig_opposite_mask <- rep(
          FALSE,
          nrow(other_dat)
        )
      }
      
      n_target_sig_opposite <- sum(
        target_sig_opposite_mask
      )
      
      n_target_sig_opposite_mat[
        ref_pfn,
        other_pfn
      ] <- n_target_sig_opposite
      
      #-----------------------------------------------------
      # Proportion of significant target vertices that
      # show the opposite direction
      #-----------------------------------------------------
      if (n_target_sig > 0) {
        
        prop_target_sig_opposite_mat[
          ref_pfn,
          other_pfn
        ] <- n_target_sig_opposite / n_target_sig
      }
    }
    
  } else {
    
    #-------------------------------------------------------
    # No significant vertices in reference PFN
    #-------------------------------------------------------
    direction_vec[ref_pfn] <- "none"
    mean_sig_beta_vec[ref_pfn] <- NA_real_
    n_vert_vec[ref_pfn] <- 0
    
    n_overlap_mat[ref_pfn, ] <- 0
    n_target_sig_mat[ref_pfn, ] <- 0
    n_target_sig_opposite_mat[ref_pfn, ] <- 0
  }
}

#===========================================================
# SAVE OUTPUTS
#===========================================================

#---------------------------
# 1. Mean target betas
#---------------------------
mean_beta_df <- data.frame(
  Reference_PFN = PFN_names,
  Direction = direction_vec,
  Mean_Beta_Among_All_Sig_Vertices = mean_sig_beta_vec,
  N_Ref_Directional_Sig_Vertices = n_vert_vec,
  mean_beta_mat,
  check.names = FALSE
)

write_csv(
  mean_beta_df,
  file.path(
    dirout,
    "logTC_tradeoff_mean_betas_all17.csv"
  )
)

#---------------------------
# 2. Available vertex counts
#---------------------------
overlap_df <- data.frame(
  Reference_PFN = PFN_names,
  Direction = direction_vec,
  N_Ref_Directional_Sig_Vertices = n_vert_vec,
  n_overlap_mat,
  check.names = FALSE
)

write_csv(
  overlap_df,
  file.path(
    dirout,
    "logTC_tradeoff_overlap_counts_all17.csv"
  )
)

#---------------------------
# 3. N target-significant vertices
#---------------------------
target_sig_df <- data.frame(
  Reference_PFN = PFN_names,
  Direction = direction_vec,
  N_Ref_Directional_Sig_Vertices = n_vert_vec,
  n_target_sig_mat,
  check.names = FALSE
)

write_csv(
  target_sig_df,
  file.path(
    dirout,
    "logTC_tradeoff_target_sig_counts_all17.csv"
  )
)

#---------------------------
# 4. N significant opposite-direction target vertices
#---------------------------
target_sig_opposite_df <- data.frame(
  Reference_PFN = PFN_names,
  Direction = direction_vec,
  N_Ref_Directional_Sig_Vertices = n_vert_vec,
  n_target_sig_opposite_mat,
  check.names = FALSE
)

write_csv(
  target_sig_opposite_df,
  file.path(
    dirout,
    "logTC_tradeoff_target_sig_opposite_counts_all17.csv"
  )
)

#---------------------------
# 5. Proportion of significant target vertices
#    that are opposite direction
#---------------------------
target_sig_opposite_prop_df <- data.frame(
  Reference_PFN = PFN_names,
  Direction = direction_vec,
  N_Ref_Directional_Sig_Vertices = n_vert_vec,
  prop_target_sig_opposite_mat,
  check.names = FALSE
)

write_csv(
  target_sig_opposite_prop_df,
  file.path(
    dirout,
    "logTC_tradeoff_target_sig_opposite_proportions_all17.csv"
  )
)

#===========================================================
# SELECT TRADEOFF NETWORKS WITHIN 90% OF STRONGEST
#
# Rule:
#   1. Exclude the reference PFN itself
#   2. Require >= 25% overlap with the reference mask
#      (target significance NOT required)
#   3. Require target mean beta to be opposite in sign
#      to the reference direction
#   4. Find the strongest tradeoff target
#   5. Retain all targets with |mean beta| >= 90% of
#      the strongest tradeoff magnitude
#===========================================================

min_overlap_prop <- 0.25
relative_tradeoff_threshold <- 0.90

# Store results here
selected_tradeoff_list <- list()

result_counter <- 1

for (ref_pfn in all_pfns) {
  
  #---------------------------------------------------------
  # Skip networks without a usable directional reference mask
  #---------------------------------------------------------
  if (
    is.na(n_vert_vec[ref_pfn]) ||
    n_vert_vec[ref_pfn] == 0 ||
    !direction_vec[ref_pfn] %in% c("positive", "negative")
  ) {
    next
  }
  
  #---------------------------------------------------------
  # Proportion of reference-mask vertices represented
  # in each target PFN
  #---------------------------------------------------------
  overlap_prop <- n_overlap_mat[ref_pfn, ] /
    n_vert_vec[ref_pfn]
  
  #---------------------------------------------------------
  # Initial candidate target PFNs
  #
  # Requirements:
  #   - not the reference PFN itself
  #   - >= 25% overlap
  #   - non-missing mean beta
  #---------------------------------------------------------
  eligible <- all_pfns[
    all_pfns != ref_pfn &
      !is.na(overlap_prop) &
      overlap_prop >= min_overlap_prop &
      !is.na(mean_beta_mat[ref_pfn, ])
  ]
  
  if (length(eligible) == 0) {
    next
  }
  
  #---------------------------------------------------------
  # Require target mean beta to be opposite in direction
  #---------------------------------------------------------
  if (direction_vec[ref_pfn] == "negative") {
    
    # Negative reference -> positive target
    eligible <- eligible[
      mean_beta_mat[ref_pfn, eligible] > 0
    ]
    
  } else if (direction_vec[ref_pfn] == "positive") {
    
    # Positive reference -> negative target
    eligible <- eligible[
      mean_beta_mat[ref_pfn, eligible] < 0
    ]
  }
  
  if (length(eligible) == 0) {
    next
  }
  
  #---------------------------------------------------------
  # Tradeoff strength = absolute target mean beta
  #---------------------------------------------------------
  tradeoff_strength <- abs(
    mean_beta_mat[ref_pfn, eligible]
  )
  
  # Strongest tradeoff magnitude
  strongest_tradeoff <- max(
    tradeoff_strength,
    na.rm = TRUE
  )
  
  #---------------------------------------------------------
  # Keep all targets within 90% of strongest
  #---------------------------------------------------------
  keep <- tradeoff_strength >=
    relative_tradeoff_threshold * strongest_tradeoff
  
  selected_pfns <- eligible[keep]
  
  #---------------------------------------------------------
  # Sort selected networks by tradeoff strength
  # strongest first
  #---------------------------------------------------------
  selected_pfns <- selected_pfns[
    order(
      abs(mean_beta_mat[ref_pfn, selected_pfns]),
      decreasing = TRUE
    )
  ]
  
  #---------------------------------------------------------
  # Store each selected tradeoff as one row
  #---------------------------------------------------------
  for (rank_i in seq_along(selected_pfns)) {
    
    target_pfn <- selected_pfns[rank_i]
    
    target_beta <- mean_beta_mat[
      ref_pfn,
      target_pfn
    ]
    
    relative_strength <- abs(target_beta) /
      strongest_tradeoff
    
    selected_tradeoff_list[[result_counter]] <-
      data.frame(
        Reference_PFN = PFN_names[ref_pfn],
        
        Reference_Direction =
          direction_vec[ref_pfn],
        
        Mean_Beta_Among_All_Sig_Ref_Vertices =
          mean_sig_beta_vec[ref_pfn],
        
        N_Ref_Directional_Sig_Vertices =
          n_vert_vec[ref_pfn],
        
        Selected_Tradeoff_PFN =
          PFN_names[target_pfn],
        
        Selected_Tradeoff_Mean_Beta =
          target_beta,
        
        Tradeoff_Strength =
          abs(target_beta),
        
        Relative_To_Strongest =
          relative_strength,
        
        Tradeoff_Rank =
          rank_i,
        
        N_Overlap =
          n_overlap_mat[
            ref_pfn,
            target_pfn
          ],
        
        Overlap_Proportion =
          overlap_prop[target_pfn],
        
        N_Target_Sig =
          n_target_sig_mat[
            ref_pfn,
            target_pfn
          ],
        
        N_Target_Sig_Opposite =
          n_target_sig_opposite_mat[
            ref_pfn,
            target_pfn
          ],
        
        Prop_Target_Sig_Opposite =
          prop_target_sig_opposite_mat[
            ref_pfn,
            target_pfn
          ],
        
        stringsAsFactors = FALSE
      )
    
    result_counter <- result_counter + 1
  }
}

#===========================================================
# COMBINE RESULTS
#===========================================================

if (length(selected_tradeoff_list) > 0) {
  
  selected_tradeoff_df <- do.call(
    rbind,
    selected_tradeoff_list
  )
  
} else {
  
  selected_tradeoff_df <- data.frame()
}

#===========================================================
# SAVE SELECTED TRADEOFF NETWORKS
#===========================================================

write_csv(
  selected_tradeoff_df,
  file.path(
    dirout,
    "logTC_tradeoff_selected_networks_within90pct_all17.csv"
  )
)