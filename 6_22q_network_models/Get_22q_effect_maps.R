library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

Hard_parcel_map<-read_cifti('../../../PNC_data/PNC_group_atlas/PNC_group_hard_parcellation.dscalar.nii') #hard parcel cifti used as template
del22q_PFN_TCA_cont <- read.csv("../../Results/22q_PFN_Age/subclass/22q_PFN_TCA_euler_age_summary_results_bonf_REML_subclass.csv")
del22q_PFN_TCA_uncont <- read.csv("../../Results/22q_PFN_Age/No_TCA/subclass/22q_PFN_no_TCA_summary_results_bonf_REML_subclass.csv")
del22q_group_TCA_cont <- read.csv("../../Results/22q_Group_Age/subclass/22q_Group_TCA_euler_age_summary_results_bonf_REML_netnames_subclass.csv")
del22q_group_TCA_uncont <- read.csv("../../Results/22q_Group_Age/No_TCA/subclass/22q_Group_no_TCA_age_summary_results_bonf_REML_netnames_subclass.csv")
del22q_func_size_TCA_cont <- read.csv("../../Results/22q_PFN_Func_Temp/subclass/log_log_22q_temp_func_area_summary_results_subclass.csv")
del22q_func_size_TCA_uncont <- read.csv("../../Results/22q_PFN_Func_Temp/subclass/log_log_22q_temp_func_area_no_TCA_summary_results_subclass.csv")

#Get network IDs
vertex_ID_lh <- Hard_parcel_map$data$cortex_left
vertex_ID_rh <- Hard_parcel_map$data$cortex_right

# --- TCA_cont_effect ---
del_TCA_effect_for_vert <- exp(del22q_PFN_TCA_cont$beta) - 1 #exp betas and set relative to 1 for interpretation
del_TCA_effect_vert_lh <- matrix(0,29696,1)
del_TCA_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      #&& del22q_PFN_TCA_cont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_effect_vert_lh[vert,] <- del_TCA_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      #&& del22q_PFN_TCA_cont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_effect_vert_rh[vert,] <- del_TCA_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_PFN_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_PFN_scaling_effect_map$data$cortex_left <- del_TCA_effect_vert_lh
PNC_PFN_scaling_effect_map$data$cortex_right <- del_TCA_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/22q_effect_on_PFNs_TCA_unthres_subclass') # change depending on sig or not
write_cifti(PNC_PFN_scaling_effect_map,outfile) # save out cifti


#--- TCA_uncont_effect ---
del_TCA_uncont_effect_for_vert <- exp(del22q_PFN_TCA_uncont$beta) - 1 #exp betas and set relative to 1 for interpretation
del_TCA_uncont_effect_vert_lh <- matrix(0,29696,1)
del_TCA_uncont_effect_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      #&& del22q_PFN_TCA_uncont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_uncont_effect_vert_lh[vert,] <- del_TCA_uncont_effect_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      #&& del22q_PFN_TCA_uncont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_uncont_effect_vert_rh[vert,] <- del_TCA_uncont_effect_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_group_scaling_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_group_scaling_effect_map$data$cortex_left <- del_TCA_uncont_effect_vert_lh
PNC_group_scaling_effect_map$data$cortex_right <- del_TCA_uncont_effect_vert_rh

outfile <- paste0('../../Results/Network_level_maps/22q_effect_on_PFNs_TCA_uncont_unthres_subclass') # change depending on sig or not
write_cifti(PNC_group_scaling_effect_map,outfile) # save out cifti



# --- Func_temp TCA_cont_effect ---
del_TCA_func_size_for_vert <- del22q_func_size_TCA_cont$beta
del_TCA_func_size_vert_lh <- matrix(0,29696,1)
del_TCA_func_size_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      #&& del22q_func_size_TCA_cont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_func_size_vert_lh[vert,] <- del_TCA_func_size_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      #&& del22q_func_size_TCA_cont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_func_size_vert_rh[vert,] <- del_TCA_func_size_for_vert[Net]
  }
}

#define map variable using hard parcel as template
del_TCA_func_size_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
del_TCA_func_size_effect_map$data$cortex_left <- del_TCA_func_size_vert_lh
del_TCA_func_size_effect_map$data$cortex_right <- del_TCA_func_size_vert_rh

outfile <- paste0('../../Results/Network_level_maps/22q_func_temp_TCA_unthres_subclass') # change depending on sig or not
write_cifti(del_TCA_func_size_effect_map,outfile) # save out cifti


#--- TCA_uncont_effect ---
del_TCA_uncont_func_size_for_vert <- del22q_func_size_TCA_uncont$beta
del_TCA_uncont_func_size_vert_lh <- matrix(0,29696,1)
del_TCA_uncont_func_size_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
      #&& del22q_func_size_TCA_uncont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_uncont_func_size_vert_lh[vert,] <- del_TCA_uncont_func_size_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
      #&& del22q_func_size_TCA_uncont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_uncont_func_size_vert_rh[vert,] <- del_TCA_uncont_func_size_for_vert[Net]
  }
}

#define map variable using hard parcel as template
del_no_TCA_func_size_effect_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
del_no_TCA_func_size_effect_map$data$cortex_left <- del_TCA_uncont_func_size_vert_lh
del_no_TCA_func_size_effect_map$data$cortex_right <- del_TCA_uncont_func_size_vert_rh

outfile <- paste0('../../Results/Network_level_maps/22q_func_temp_TCA_uncont_unthres_subclass') # change depending on sig or not
write_cifti(del_no_TCA_func_size_effect_map,outfile) # save out cifti



# --- TCA_cont_group ---
del_TCA_group_for_vert <- exp(del22q_group_TCA_cont$beta) - 1 #exp betas and set relative to 1 for interpretation
del_TCA_group_vert_lh <- matrix(0,29696,1)
del_TCA_group_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
    #&& del22q_group_TCA_cont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_group_vert_lh[vert,] <- del_TCA_group_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
    #&& del22q_group_TCA_cont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_group_vert_rh[vert,] <- del_TCA_group_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_PFN_scaling_group_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_PFN_scaling_group_map$data$cortex_left <- del_TCA_group_vert_lh
PNC_PFN_scaling_group_map$data$cortex_right <- del_TCA_group_vert_rh

outfile <- paste0('../../Results/Network_level_maps/22q_group_TCA_unthres_subclass') # change depending on sig or not
write_cifti(PNC_PFN_scaling_group_map,outfile) # save out cifti


#--- TCA_uncont_group ---
del_TCA_uncont_group_for_vert <- exp(del22q_group_TCA_uncont$beta) - 1 #exp betas and set relative to 1 for interpretation
del_TCA_uncont_group_vert_lh <- matrix(0,29696,1)
del_TCA_uncont_group_vert_rh <- matrix(0,29716,1)

#L hemi
for (vert in c(1:29696)) {
  Net <- vertex_ID_lh[vert]
  if (Net != 0) {
    #&& del22q_group_TCA_uncont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_uncont_group_vert_lh[vert,] <- del_TCA_uncont_group_for_vert[Net]
  }
}
#R hemi
for (vert in c(1:29716)) {
  Net <- vertex_ID_rh[vert]
  if (Net != 0) {
    #&& del22q_group_TCA_uncont$pvals_bonf[Net] < 0.05) { #can comment out to get unthresholded
    del_TCA_uncont_group_vert_rh[vert,] <- del_TCA_uncont_group_for_vert[Net]
  }
}

#define map variable using hard parcel as template
PNC_group_scaling_group_map <- Hard_parcel_map

#assign network values for each hemi into cifti 
PNC_group_scaling_group_map$data$cortex_left <- del_TCA_uncont_group_vert_lh
PNC_group_scaling_group_map$data$cortex_right <- del_TCA_uncont_group_vert_rh

outfile <- paste0('../../Results/Network_level_maps/22q_group_TCA_uncont_unthres_subclass') # change depending on sig or not
write_cifti(PNC_group_scaling_group_map,outfile) # save out cifti