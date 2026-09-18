library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

Hard_parcel_map<-read_cifti('../../../PNC_data/PNC_group_atlas/PNC_group_hard_parcellation.dscalar.nii') #hard parcel cifti used as template
PNC_PFN_scaling_effects <- read.csv("../../Results/PNC_PFN_Age/PNC_PFNs_log_log_summary_results_bonf_REML_netnames.csv")
PNC_group_scaling_effects <- read.csv("../../Results/PNC_Group_Age/PNC_group_log_log_summary_results_bonf_REML_netnames.csv")
Scaling_effect_diff <- read.csv("../../Results/PNC_Delta_models/PNC_Group_vs_PFN_boot_results.csv")
Func_rep_scaling_effects <- read.csv("../../Results/PNC_PFN_Func_Temp/normed_temp_func_area_summary_results.csv")

#Get network IDs
vertex_ID_lh <- Hard_parcel_map$data$cortex_left
vertex_ID_rh <- Hard_parcel_map$data$cortex_right


# --- PNC_PFN_scaling ---
PFN_effect_for_vert <- PNC_PFN_scaling_effects$beta - 1 #sets betas relative to 0 for easier visualization
#PFN_effect_for_vert <- PNC_PFN_scaling_effects$beta
PFN_effect_vert_lh <- matrix(0,29696,1)
PFN_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      # && PNC_PFN_scaling_effects$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    PFN_effect_vert_lh[vert,] <- PFN_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      # && PNC_PFN_scaling_effects$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    PFN_effect_vert_rh[vert,] <- PFN_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_PFN_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_PFN_scaling_effect_map$data$cortex_left <- PFN_effect_vert_lh
PNC_PFN_scaling_effect_map$data$cortex_right <- PFN_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_PFN_scaling_effect_map_relative_0_unthres') # change depending on sig or not
write_cifti(PNC_PFN_scaling_effect_map,outfile) # save out cifti


#--- PNC_group_scaling ---
group_effect_for_vert <- PNC_group_scaling_effects$beta - 1 #sets betas relative to 0 for easier visualization
#group_effect_for_vert <- PNC_group_scaling_effects$beta
group_effect_vert_lh <- matrix(0,29696,1)
group_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      # && PNC_group_scaling_effects$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    group_effect_vert_lh[vert,] <- group_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      # && PNC_group_scaling_effects$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    group_effect_vert_rh[vert,] <- group_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_group_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_group_scaling_effect_map$data$cortex_left <- group_effect_vert_lh
PNC_group_scaling_effect_map$data$cortex_right <- group_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_group_scaling_effect_map_relative_0_unthres') # change depending on sig or not
write_cifti(PNC_group_scaling_effect_map,outfile) # save out cifti



#---Diff_group_scaling---
Diff_effect_for_vert <- Scaling_effect_diff$delta_beta
Diff_effect_vert_lh <- matrix(0,29696,1)
Diff_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      # && Scaling_effect_diff$p_delta_beta_bonf[Net] < 0.05) { #can comment out to get unthresholded
    Diff_effect_vert_lh[vert,] <- Diff_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      # && Scaling_effect_diff$p_delta_beta_bonf[Net] < 0.05) { #can comment out to get unthresholded
    Diff_effect_vert_rh[vert,] <- Diff_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_Diff_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_Diff_scaling_effect_map$data$cortex_left <- Diff_effect_vert_lh
PNC_Diff_scaling_effect_map$data$cortex_right <- Diff_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_Diff_scaling_effect_map_unthres') # change depending on sig or not
write_cifti(PNC_Diff_scaling_effect_map,outfile) # save out cifti


#---Func_temp_scaling---
Func_rep_effect_for_vert <- Func_rep_scaling_effects$beta*log(1.5)
Func_rep_effect_vert_lh <- matrix(0,29696,1)
Func_rep_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      #&& Func_rep_scaling_effects$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    Func_rep_effect_vert_lh[vert,] <- Func_rep_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      #&& Func_rep_scaling_effects$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    Func_rep_effect_vert_rh[vert,] <- Func_rep_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_Func_rep_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_Func_rep_scaling_effect_map$data$cortex_left <- Func_rep_effect_vert_lh
PNC_Func_rep_scaling_effect_map$data$cortex_right <- Func_rep_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_Func_temp_scaling_effect_map_unthres')  # change depending on sig or not
write_cifti(PNC_Func_rep_scaling_effect_map,outfile) # save out cifti


