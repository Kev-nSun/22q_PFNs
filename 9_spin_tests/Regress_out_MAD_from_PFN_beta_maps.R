library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

#Read in maps
PFNs_hardparcel<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')
Scaling_max_map<-read_cifti('../../Results/TCA_PFN_loadings_gam_results/NORMED/Maps/Unthresholded/ABS_MAX_logTC_PFN_loading_gam_weights_LOG1.5_unthres.dscalar.nii', mwall_values = NULL)
del_22q_max_map<-read_cifti('../../Results/22q_PFN_loadings_gam_results_subclass/Normed_TCA/Maps/Unthresholded/ABS_MAX_22q_PFN_loading_gam_weights_unthres.dscalar.nii', mwall_values = NULL)
MAD_avg_map<-read_cifti('../../../PNC_data/MAD_map/MAD_avg.dscalar.nii', mwall_values = NULL)

# Save out effect sizes
Scaling_betas <- c(Scaling_max_map$data$cortex_left,Scaling_max_map$data$cortex_right)
del_22q_betas <- c(del_22q_max_map$data$cortex_left,del_22q_max_map$data$cortex_right)
MAD_avgs <- c(MAD_avg_map$data$cortex_left,MAD_avg_map$data$cortex_right)

hemi_ids <- factor(
  c(
    rep("LH", 32492),
    rep("RH", 32492)
  ),
  levels = c("LH", "RH")
)

vertex_vals <- data.frame(Scaling_betas = Scaling_betas,del_22q_betas = del_22q_betas,MAD_avgs = MAD_avgs,hemi_ids = hemi_ids)

# exclude any vertices with zero variability
vertex_vals$MAD_avgs[vertex_vals$MAD_avgs == 0] <- NA_real_

# --- Scaling_resid ---
scaling_lm <- lm(Scaling_betas ~ MAD_avgs, data = vertex_vals, na.action = na.exclude) #residual model
Scaling_betas_resid <- residuals(scaling_lm)
scaling_lm_hemi <- lm(Scaling_betas ~ MAD_avgs + hemi_ids, data = vertex_vals, na.action = na.exclude) #residual model with hemisphere covariate
Scaling_betas_resid_hemi <- residuals(scaling_lm_hemi)
scaling_hemi_corr <- cor(Scaling_betas_resid, Scaling_betas_resid_hemi, use = "complete.obs") #check correlation between 2 models, r = 0.9955

Scaling_resid_vert_lh <- matrix(0,32492,1)
Scaling_resid_vert_rh <- matrix(0,32492,1)

#L hemi
for (vert in c(1:32492)) {
  Scaling_resid_vert_lh[vert,] <- Scaling_betas_resid[vert]
}
#R hemi
for (vert in c(1:32492)) {
  Scaling_resid_vert_rh[vert,] <- Scaling_betas_resid[vert+32492]
}

#define map variable
Scaling_MAD_resid_map <- Scaling_max_map

#assign network values for each hemi into cifti 
Scaling_MAD_resid_map$data$cortex_left <- Scaling_resid_vert_lh
Scaling_MAD_resid_map$data$cortex_right <- Scaling_resid_vert_rh

outfile <- paste0('../../Results/Spin_Tests/MAD_&_resid_maps/Scaling_MAD_resid_normed')
write_cifti(Scaling_MAD_resid_map,outfile) # save out cifti

# --- del_22q_resid ---
del_22q_lm <- lm(del_22q_betas ~ MAD_avgs, data = vertex_vals, na.action = na.exclude)
del_22q_betas_resid <- residuals(del_22q_lm)
del_22q_lm_hemi <- lm(del_22q_betas ~ MAD_avgs + hemi_ids, data = vertex_vals, na.action = na.exclude) #residual model with hemisphere covariate
del_22q_betas_resid_hemi <- residuals(del_22q_lm_hemi)
del_22q_hemi_corr <- cor(del_22q_betas_resid, del_22q_betas_resid_hemi, use = "complete.obs") #check correlation between 2 models, r = 1.0

del_22q_resid_vert_lh <- matrix(0,32492,1)
del_22q_resid_vert_rh <- matrix(0,32492,1)

#L hemi
for (vert in c(1:32492)) {
  del_22q_resid_vert_lh[vert,] <- del_22q_betas_resid[vert]
}
#R hemi
for (vert in c(1:32492)) {
  del_22q_resid_vert_rh[vert,] <- del_22q_betas_resid[vert+32492]
}

#define map variable
del_22q_MAD_resid_map <- del_22q_max_map

#assign network values for each hemi into cifti 
del_22q_MAD_resid_map$data$cortex_left <- del_22q_resid_vert_lh
del_22q_MAD_resid_map$data$cortex_right <- del_22q_resid_vert_rh

outfile <- paste0('../../Results/Spin_Tests/MAD_&_resid_maps/del_22q_MAD_resid_normed')
write_cifti(del_22q_MAD_resid_map,outfile) # save out cifti