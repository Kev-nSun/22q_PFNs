library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

Hard_parcel_map<-read_cifti('../../../PNC_data/PNC_group_atlas/PNC_group_hard_parcellation.dscalar.nii') #hard parcel cifti used as template
PNC_PFN_scaling_effects <- read.csv("../../Results/PNC_PFNs_net_areas_log-log_age_cov/PNC_PFNs_log_log_summary_results.csv")
PNC_group_scaling_effects <- read.csv("../../Results/PNC_group_net_areas_log-log_age_cov/PNC_group_log_log_summary_results.csv")
Scaling_effect_diff <- read.csv("../../Results/Area_stats/PNC_Group_vs_PFN_scaling_beta_diff_results.csv")

#Get network IDs
vertex_ID_lh <- Hard_parcel_map$data$cortex_left
vertex_ID_rh <- Hard_parcel_map$data$cortex_right


#PNC_PFN_scaling
PFN_effect_for_vert <- PNC_PFN_scaling_effects$beta - 1 #sets betas relative to 0 for easier visualization
#PFN_effect_for_vert <- PNC_PFN_scaling_effects$beta
PFN_effect_vert_lh <- matrix(0,29696,1)
PFN_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0 
      && PNC_PFN_scaling_effects$pvals_fdr[Net] < 0.05) {
    PFN_effect_vert_lh[vert,] <- PFN_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0
      && PNC_PFN_scaling_effects$pvals_fdr[Net] < 0.05) {
    PFN_effect_vert_rh[vert,] <- PFN_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_PFN_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_PFN_scaling_effect_map$data$cortex_left <- PFN_effect_vert_lh
PNC_PFN_scaling_effect_map$data$cortex_right <- PFN_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_PFN_scaling_effect_map_relative_0_sig_only') # change depending on rel 0 or not
write_cifti(PNC_PFN_scaling_effect_map,outfile) # save out cifti


#PNC_group_scaling
group_effect_for_vert <- PNC_group_scaling_effects$beta - 1 #sets betas relative to 0 for easier visualization
#group_effect_for_vert <- PNC_group_scaling_effects$beta
group_effect_vert_lh <- matrix(0,29696,1)
group_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0 
      && PNC_group_scaling_effects$pvals_fdr[Net] < 0.05) {
    group_effect_vert_lh[vert,] <- group_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0
      && PNC_group_scaling_effects$pvals_fdr[Net] < 0.05) {
    group_effect_vert_rh[vert,] <- group_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_group_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_group_scaling_effect_map$data$cortex_left <- group_effect_vert_lh
PNC_group_scaling_effect_map$data$cortex_right <- group_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_group_scaling_effect_map_relative_0_sig_only') # change depending on rel 0 or not
write_cifti(PNC_group_scaling_effect_map,outfile) # save out cifti



#Diff_group_scaling
Diff_effect_for_vert <- Scaling_effect_diff$beta
Diff_effect_vert_lh <- matrix(0,29696,1)
Diff_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0 
      && Scaling_effect_diff$pvals_fdr[Net] < 0.05) {
    Diff_effect_vert_lh[vert,] <- Diff_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0
      && Scaling_effect_diff$pvals_fdr[Net] < 0.05) {
    Diff_effect_vert_rh[vert,] <- Diff_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_Diff_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_Diff_scaling_effect_map$data$cortex_left <- Diff_effect_vert_lh
PNC_Diff_scaling_effect_map$data$cortex_right <- Diff_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/PNC_Diff_scaling_effect_map_sig_only') # change depending on rel 0 or not
write_cifti(PNC_Diff_scaling_effect_map,outfile) # save out cifti

