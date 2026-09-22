
library(readr)

#TCA
dirin <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/091426_22q_logTC_vert_normed_loadings_FDR_pval_beta_subclass/"
dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_TCA/"
dir.create(dirout)
cov_mod <- "TCA"

#no TCA
# dirin <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/091426_22q_no_TCA_vert_normed_loadings_FDR_pval_beta_subclass/"
# dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_no_TCA/"
# dir.create(dirout)
# cov_mod <- "no_TCA"

var_verts <- read.csv("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Variable_vertices_PFN_masks/PFN_var_verts.csv")

#COUNT NUMBER OF SIG VERTICES, FDR THRESHOLDED ACROSS ALL NETWORKS:
sig_mat <- matrix(NA_real_, 17, 9)

for (PFN in 1:17) {
  
  TCA_PFN_gam_results <- read.csv(
    sprintf(
      "%s/stat_22q_gam_vert_loadings_PFN_%02d_beta_fdr.csv",
      dirin,
      PFN
    )
  )
  
  beta  <- TCA_PFN_gam_results[, 2]
  p_fdr <- TCA_PFN_gam_results[, 3]
  
  # Vertices that have valid modeled values
  valid <- !is.na(beta) & !is.na(p_fdr)
  
  # FDR-significant modeled vertices
  sig <- valid & p_fdr < 0.05
  
  # Denominator based on external variable vertex mask
  n_var <- var_verts[PFN,2]
  
  # Counts
  sig_verts <- sum(sig)
  pos_sig_verts <- sum(sig & beta > 0, na.rm = TRUE)
  neg_sig_verts <- sum(sig & beta < 0, na.rm = TRUE)
  
  # Mean beta among all FDR-significant vertices
  mean_sig_beta <- if (sig_verts > 0) {
    mean(beta[sig], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  sig_mat[PFN, 1] <- n_var
  sig_mat[PFN, 2] <- sig_verts
  sig_mat[PFN, 3] <- sig_verts / n_var
  sig_mat[PFN, 4] <- pos_sig_verts
  sig_mat[PFN, 5] <- neg_sig_verts
  sig_mat[PFN, 6] <- pos_sig_verts / n_var
  sig_mat[PFN, 7] <- neg_sig_verts / n_var
  sig_mat[PFN, 8] <- (pos_sig_verts - neg_sig_verts) / n_var
  sig_mat[PFN, 9] <- mean_sig_beta
}

sig_mat <- data.frame(sig_mat)

colnames(sig_mat) <- c(
  "Variable_Verts",
  "FDR_Sig_Verts",
  "Sig_Vert_Prop",
  "Sig_Pos",
  "Sig_Neg",
  "Pos_Prop",
  "Neg_Prop",
  "Directional_Diff",
  "Mean_Sig_Beta"
)

PFN_names <- c("DM-B", "SM-foot", "FP-A", "SM-face", "DA-post", "VS-peri", "Salience", "Language", "CingOp", "VS-cent", "SM-lh", "DM-A", "SM-rh", "DA-ant", "ParMem", "Auditory", "FP-B")

sig_22q_summary <- data.frame(PFN = PFN_names, sig_mat)
write_csv(sig_22q_summary,paste0(dirout,"/22q_gam_vert_normed_loadings_FDR_beta_summary_",cov_mod,".csv"))
