# Read in 1st and 99th percentile frames to calculate difference between maps for visualization

library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

#Output dir
map_dir <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/Normed_delta_maps/"
dir.create(map_dir)

#Read in frames
perc_1st <- read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/091426_VisualizeFolder/Frame_001/Frame_001_mwall.dlabel.nii')
perc_99th <- read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/091426_VisualizeFolder/Frame_100/Frame_100_mwall.dlabel.nii')

perc_1st_lh <- perc_1st$data$cortex_left
perc_1st_rh <- perc_1st$data$cortex_right
perc_99th_lh <- perc_99th$data$cortex_left
perc_99th_rh <- perc_99th$data$cortex_right

# mask of PFN topography stable to scaling
stable_mask_lh <-  perc_1st_lh == perc_99th_lh
stable_mask_rh <-  perc_1st_rh == perc_99th_rh

#Determine network topography predicted to be stable to scaling
stable_lh <- perc_1st_lh
stable_rh <- perc_1st_rh

stable_lh[!stable_mask_lh] <- 0
stable_rh[!stable_mask_rh] <- 0

#Determine network topography predicted to associated with 1st percentile of TCA
delta_1st_lh <- perc_1st_lh
delta_1st_rh <- perc_1st_rh

delta_1st_lh[stable_mask_lh] <- 0
delta_1st_rh[stable_mask_rh] <- 0

#Determine network topography predicted to associated with 99th percentile of TCA
delta_99th_lh <- perc_99th_lh
delta_99th_rh <- perc_99th_rh

delta_99th_lh[stable_mask_lh] <- 0
delta_99th_rh[stable_mask_rh] <- 0

#Save left + right hem to mat files, converted to gifti--> cifti in matlab "Labeling_PFNs_into_cifti.m"

# For L hemi, ensure that both the list and its elements are named
stable_names_list <- list()
stable_names_list$stable_lh <- stable_lh[,1]
writeMat(paste0(map_dir,"stable_scaling_pred_topo_L.mat"), x = stable_names_list)
# For R hemi, ensure that both the list and its elements are named
stable_names_list <- list()
stable_names_list$stable_rh <- stable_rh[,1]
writeMat(paste0(map_dir,"stable_scaling_pred_topo_R.mat"), x = stable_names_list)

# For L hemi, ensure that both the list and its elements are named
delta_1st_list <- list()
delta_1st_list$delta_1st_lh <- delta_1st_lh[,1]
writeMat(paste0(map_dir,"delta_1st_scaling_pred_topo_L.mat"), x = delta_1st_list)
# For R hemi, ensure that both the list and its elements are named
delta_1st_list <- list()
delta_1st_list$delta_1st_rh <- delta_1st_rh[,1]
writeMat(paste0(map_dir,"delta_1st_scaling_pred_topo_R.mat"), x = delta_1st_list)

# For L hemi, ensure that both the list and its elements are named
delta_99th_list <- list()
delta_99th_list$delta_99th_lh <- delta_99th_lh[,1]
writeMat(paste0(map_dir,"delta_99th_scaling_pred_topo_L.mat"), x = delta_99th_list)
# For R hemi, ensure that both the list and its elements are named
delta_99th_list <- list()
delta_99th_list$delta_99th_rh <- delta_99th_rh[,1]
writeMat(paste0(map_dir,"delta_99th_scaling_pred_topo_R.mat"), x = delta_99th_list)

