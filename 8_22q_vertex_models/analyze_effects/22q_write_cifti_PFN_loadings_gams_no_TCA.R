
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/091426_22q_no_TCA_vert_normed_loadings_FDR_pval_beta_subclass"

dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_no_TCA/Maps/"
dir.create(dirout)

abs_dirout <- sprintf("%s/ABS_maps/",dirout)
dir.create(abs_dirout)

bin_dirout <- sprintf("%s/Binary_cluster_maps/",dirout)
dir.create(bin_dirout)

clust_dirout <- sprintf("%s/Cluster_maps/",dirout)
dir.create(clust_dirout)

unthres_dirout <- sprintf("%s/Unthresholded/",dirout)
dir.create(unthres_dirout)

PFNs_hardparcel<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')

#Surface files needed for clustering
L_surf_file <- 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Inputs/tpl-fsLR_den-32k_hemi-L_midthickness.surf.gii'
R_surf_file <- 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Inputs/tpl-fsLR_den-32k_hemi-R_midthickness.surf.gii'

#GET ALL ZERO MASK FOR VISUALIZATION
PFNs_all_zero_mask<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/PNC_PFN_all_zero_mask.dscalar.nii')
all_zero_lh <- PFNs_all_zero_mask$data$cortex_left
all_zero_rh <- PFNs_all_zero_mask$data$cortex_right

#PFN BETA MAPS FOR 22Q GAM MODELS, FDR THRESHOLDED ACROSS ALL NETWORKS
for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  stat_22q_PFN_gam_results <- read.csv(sprintf("%s/stat_22q_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    if (is.na(stat_22q_PFN_gam_results[i,3])) {
      next
    }
    if (stat_22q_PFN_gam_results[i,3] < 0.05) { #threshold by fdr p-vals
      vertex_weights_lh[i,1] <- stat_22q_PFN_gam_results[i,2]
    }
  }
  
  #assign weights to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    if (is.na(stat_22q_PFN_gam_results[i+29696,3])) {
      next
    }
    if (stat_22q_PFN_gam_results[i+29696,3] < 0.05) { #threshold by fdr p-vals
      vertex_weights_rh[i,1] <- stat_22q_PFN_gam_results[i+29696,2]
    }
  }
  
  #define map variable using hard parcel cifti file
  Weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Weight_map$data$cortex_left <- vertex_weights_lh
  Weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Weight_map_med <- move_from_mwall(Weight_map)
  
  #define output file for map
  outfile <- paste0(dirout,'PFN',PFN,'_22q_PFN_loading_gam_weights_fdr_thresholded_no_TCA')
  
  #save out cifti files
  write_cifti(Weight_map_med,outfile)
  
  #save out ABS beta cifti files (MED WALL NOT ADDED BACK IN)
  Weight_map_abs <- PFNs_hardparcel
  abs_vertex_weights_lh <- abs(vertex_weights_lh)
  abs_vertex_weights_rh <- abs(vertex_weights_rh)
  Weight_map_abs$data$cortex_left <- abs_vertex_weights_lh
  Weight_map_abs$data$cortex_right <- abs_vertex_weights_rh
  outfile_abs <- paste0(abs_dirout,'/PFN',PFN,'_22q_PFN_loading_gam_weights_fdr_thresholded_ABS_no_TCA')
  write_cifti(Weight_map_abs,outfile_abs)
  
  #find clusters using 50 mm2 surface area threshold in cifti, based on abs value beta maps
  infile <- outfile_abs
  clust_file <- paste0(bin_dirout,'/PFN',PFN,'_50mm_cluster_binary_22q_no_TCA.dscalar.nii')
  cmd <- paste0('wb_command -cifti-find-clusters ',infile,'.dscalar.nii ',1e-12,' ',50,' ',0,' ',0,' COLUMN ',clust_file,' -left-surface ',L_surf_file,' -right-surface ',R_surf_file)
  system(cmd)
  
  #read in clustered vertices and relabel them with vertex gam weight
  clust_cifti<-read_cifti(clust_file) #cluster binary cifti
  #L hemi
  vertex_weight_clust_lh <- matrix(0,29696)
  clusts_lh <- clust_cifti$data$cortex_left
  for (i in c(1:29696)) {
    if (clusts_lh[i] > 0){ #if vertex passed cluster threshold
      vertex_weight_clust_lh[i] <- stat_22q_PFN_gam_results[i,2] #if in cluster, set vertex to beta weight
    }
  }
  #R hemi
  vertex_weight_clust_rh <- matrix(0,29716)
  clusts_rh <- clust_cifti$data$cortex_right
  for (i in c(1:29716)) {
    if (clusts_rh[i] > 0){ #if vertex passed cluster threshold
      vertex_weight_clust_rh[i,1] <- stat_22q_PFN_gam_results[i+29696,2] #if in cluster, set vertex to beta weight
    }
  }
  
  #define clustered map variable using weights_cifti file
  Weight_map_clustered <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Weight_map_clustered$data$cortex_left <- vertex_weight_clust_lh
  Weight_map_clustered$data$cortex_right <- vertex_weight_clust_rh
  
  #add back in medial wall
  Weight_map_clustered_med <- move_to_mwall(Weight_map_clustered)
  
  #define output file for averaged map
  clust_outfile <- paste0(clust_dirout,'/PFN',PFN,'_22q_PFN_loading_gam_weights_fdr_thresholded_50mm_CLUST_no_TCA')
  
  #save out cifti files
  write_cifti(Weight_map_clustered_med,clust_outfile)
}

#PFN BETA MAPS FOR 22Q GAM MODELS, UNTHRESHOLDED
for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  stat_22q_PFN_gam_results <- read.csv(sprintf("%s/stat_22q_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_weights_lh[i,1] <- stat_22q_PFN_gam_results[i,2]
  }
  
  #assign weights to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_weights_rh[i,1] <- stat_22q_PFN_gam_results[i+29696,2]
  }
  
  #define map variable using hard parcel cifti file
  Weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Weight_map$data$cortex_left <- vertex_weights_lh
  Weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Weight_map_med <- move_from_mwall(Weight_map)
  
  #define output file for averaged map
  outfile <- paste0(unthres_dirout,'PFN',PFN,'_22q_PFN_loading_gam_weights_unthresholded_no_TCA')
  
  #save out cifti files
  write_cifti(Weight_map_med,outfile)
}


#MAX ABS PFN BETA MAPS FOR 22q GAM MODELS, UNTHRESHOLDED
vertex_weights_lh <- matrix(0,29696,1)
vertex_weights_rh <- matrix(0,29716,1)

for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  stat_22q_PFN_gam_results <- read.csv(sprintf("%s/stat_22q_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  for (i in c(1:29696)) {
    if (is.na(stat_22q_PFN_gam_results[i,2])) {
      next
    }
    if (abs(stat_22q_PFN_gam_results[i,2]) > vertex_weights_lh[i,1]) { # if abs weight greater than current weight
      vertex_weights_lh[i,1] <- abs(stat_22q_PFN_gam_results[i,2])
    }
  }
  
  #assign weights to right hemi
  for (i in c(1:29716)) {
    if (is.na(stat_22q_PFN_gam_results[i+29696,2])) {
      next
    }
    if (abs(stat_22q_PFN_gam_results[i+29696,2]) > vertex_weights_rh[i,1]) { # if abs weight greater than current weight
      vertex_weights_rh[i,1] <- abs(stat_22q_PFN_gam_results[i+29696,2])
    }
  }
}

#Mask by low signal regions
vertex_weights_lh[all_zero_lh == 1] <- NA
vertex_weights_rh[all_zero_rh == 1] <- NA

#define map variable using hard parcel cifti file
ABS_weight_map <- PFNs_hardparcel

#assign weight values into hard parcel as template
ABS_weight_map$data$cortex_left <- vertex_weights_lh
ABS_weight_map$data$cortex_right <- vertex_weights_rh

#add back in medial wall
ABS_weight_map_med <- move_from_mwall(ABS_weight_map)

#define output file for averaged map
outfile <- paste0(unthres_dirout,'ABS_MAX_22q_PFN_loading_gam_weights_unthres_no_TCA')

#save out cifti file
write_cifti(ABS_weight_map_med,outfile)
