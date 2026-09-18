
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/091426_TCA_vert_normed_loadings_FDR_pval_beta_logTC/"

#Define out map directory and subdirs
dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/NORMED/Maps/"
dir.create(dirout)

#absolute value dirout
abs_dirout <- sprintf("%s/ABS_maps/",dirout)
dir.create(abs_dirout)

#binary dirout
bin_dirout <- sprintf("%s/Binary_cluster_maps/",dirout)
dir.create(bin_dirout)

#cluster dirout
clust_dirout <- sprintf("%s/Cluster_maps/",dirout)
dir.create(clust_dirout)

#unthres max map dirout
unthres_dirout <- sprintf("%s/Unthresholded/",dirout)
dir.create(unthres_dirout)

#Hard parcel map used as template
PFNs_hardparcel<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')

#Surface files needed for clustering
L_surf_file <- 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Inputs/tpl-fsLR_den-32k_hemi-L_midthickness.surf.gii'
R_surf_file <- 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Inputs/tpl-fsLR_den-32k_hemi-L_midthickness.surf.gii'

#GET ALL ZERO MASK FOR VISUALIZATION
PFNs_all_zero_mask<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/PNC_PFN_all_zero_mask.dscalar.nii')
all_zero_lh <- PFNs_all_zero_mask$data$cortex_left
all_zero_rh <- PFNs_all_zero_mask$data$cortex_right

#PFN BETA MAPS FOR TCA GAM MODELS, FDR THRESHOLDED ACROSS ALL NETWORKS
for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  TCA_PFN_gam_results <- read.csv(sprintf("%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    if (is.na(TCA_PFN_gam_results[i,3])) {
      next
    }
    if (TCA_PFN_gam_results[i,3] < 0.05) { #threshold by fdr p-vals
      vertex_weights_lh[i,1] <- TCA_PFN_gam_results[i,2]*log(1.5) #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
    }
  }
  
  #assign weights to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    if (is.na(TCA_PFN_gam_results[i+29696,3])) {
      next
    }
    if (TCA_PFN_gam_results[i+29696,3] < 0.05) { #threshold by fdr p-vals
      vertex_weights_rh[i,1] <- TCA_PFN_gam_results[i+29696,2]*log(1.5) #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
    }
  }
  
  #define map variable using hard parcel cifti file
  Weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Weight_map$data$cortex_left <- vertex_weights_lh
  Weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Weight_map_med <- move_from_mwall(Weight_map)
  
  #define output file for averaged map
  outfile <- paste0(dirout,'PFN',PFN,'_logTC_PFN_loading_gam_weights_LOG1.5_fdr_thresholded')
  
  #save out cifti files
  write_cifti(Weight_map_med,outfile)
  
  #save out ABS beta cifti files (MED WALL NOT ADDED BACK IN)
  Weight_map_abs <- PFNs_hardparcel
  abs_vertex_weights_lh <- abs(vertex_weights_lh)
  abs_vertex_weights_rh <- abs(vertex_weights_rh)
  Weight_map_abs$data$cortex_left <- abs_vertex_weights_lh
  Weight_map_abs$data$cortex_right <- abs_vertex_weights_rh
  outfile_abs <- paste0(abs_dirout,'/PFN',PFN,'_logTC_PFN_loading_gam_weights_fdr_thresholded_ABS')
  write_cifti(Weight_map_abs,outfile_abs)

  #find clusters using 50 mm2 surface area threshold in cifti, based on abs value beta maps
  infile <- outfile_abs
  clust_file <- paste0(bin_dirout,'/PFN',PFN,'_50mm_cluster_binary_TCA.dscalar.nii')
  cmd <- paste0('wb_command -cifti-find-clusters ',infile,'.dscalar.nii ',1e-12,' ',50,' ',0,' ',0,' COLUMN ',clust_file,' -left-surface C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii -right-surface C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii')
  system(cmd)

  #read in clustered vertices and relabel them with vertex gam weight
  clust_cifti<-read_cifti(clust_file) #cluster binary cifti
  #L hemi
  vertex_weight_clust_lh <- matrix(0,29696)
  clusts_lh <- clust_cifti$data$cortex_left
  for (i in c(1:29696)) {
    if (clusts_lh[i] > 0){ #if vertex passed cluster threshold, set vertex to beta weight
      vertex_weight_clust_lh[i] <- TCA_PFN_gam_results[i,2]*log(1.5) #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
    }
  }
  #R hemi
  vertex_weight_clust_rh <- matrix(0,29716)
  clusts_rh <- clust_cifti$data$cortex_right
  for (i in c(1:29716)) {
    if (clusts_rh[i] > 0){ #if vertex passed cluster threshold, set vertex to beta weight
      vertex_weight_clust_rh[i,1] <- TCA_PFN_gam_results[i+29696,2]*log(1.5) #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
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
  clust_outfile <- paste0(clust_dirout,'/PFN',PFN,'_logTC_PFN_loading_gam_weights_LOG1.5_fdr_thresholded_50mm_CLUST')

  #save out cifti files
  write_cifti(Weight_map_clustered_med,clust_outfile)
}

#PFN BETA MAPS FOR TCA GAM MODELS, UNTHRESHOLDED
for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  TCA_PFN_gam_results <- read.csv(sprintf("%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_weights_lh[i,1] <- TCA_PFN_gam_results[i,2]*log(1.5) #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
  }
  
  #assign weights to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_weights_rh[i,1] <- TCA_PFN_gam_results[i+29696,2]*log(1.5) #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
  }
  
  #define map variable using hard parcel cifti file
  Weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Weight_map$data$cortex_left <- vertex_weights_lh
  Weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Weight_map_med <- move_from_mwall(Weight_map)
  
  #define output file for averaged map
  outfile <- paste0(unthres_dirout,'PFN',PFN,'_logTC_PFN_loading_gam_weights_LOG1.5_unthresholded')
  
  #save out cifti files
  write_cifti(Weight_map_med,outfile)
}

#MAX ABS PFN BETA MAPS FOR TCA GAM MODELS, UNTHRESHOLDED
vertex_weights_lh <- matrix(0,29696,1)
vertex_weights_rh <- matrix(0,29716,1)

for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  TCA_PFN_gam_results <- read.csv(sprintf("%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  for (i in c(1:29696)) {
    if (is.na(TCA_PFN_gam_results[i,2])) {
      next
    }
    if (abs(TCA_PFN_gam_results[i,2]) > vertex_weights_lh[i,1]) { # if abs weight greater than current weight
      vertex_weights_lh[i,1] <- abs(TCA_PFN_gam_results[i,2])
    }
  }
  
  #assign weights to right hemi
  for (i in c(1:29716)) {
    if (is.na(TCA_PFN_gam_results[i+29696,2])) {
      next
    }
    if (abs(TCA_PFN_gam_results[i+29696,2]) > vertex_weights_rh[i,1]) { # if abs weight greater than current weight
      vertex_weights_rh[i,1] <- abs(TCA_PFN_gam_results[i+29696,2])
    }
  }
}

  #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
  vertex_weights_lh[,1] <- vertex_weights_lh[,1]*log(1.5)
  vertex_weights_rh[,1] <- vertex_weights_rh[,1]*log(1.5)
  
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
  outfile <- paste0(unthres_dirout,'ABS_MAX_logTC_PFN_loading_gam_weights_LOG1.5_unthres')
  
  #save out cifti file
  write_cifti(ABS_weight_map_med,outfile)

#MAX POS PFN BETA MAPS FOR TCA GAM MODELS, UNTHRESHOLDED
vertex_weights_lh <- matrix(0,29696,1)
vertex_weights_rh <- matrix(0,29716,1)

for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  TCA_PFN_gam_results <- read.csv(sprintf("%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  for (i in c(1:29696)) {
    if (is.na(TCA_PFN_gam_results[i,2])) {
      next
    }
    if (TCA_PFN_gam_results[i,2] > 0 && TCA_PFN_gam_results[i,2] > vertex_weights_lh[i,1]) { # if weight is pos and greater than current weight
      vertex_weights_lh[i,1] <- TCA_PFN_gam_results[i,2]
    }
  }
  
  #assign weights to right hemi
  for (i in c(1:29716)) {
    if (is.na(TCA_PFN_gam_results[i+29696,2])) {
      next
    }
    if (TCA_PFN_gam_results[i+29696,2] > 0 && TCA_PFN_gam_results[i+29696,2] > vertex_weights_rh[i,1]) { # if weight is pos and greater than current weight
      vertex_weights_rh[i,1] <- TCA_PFN_gam_results[i+29696,2]
    }
  }
}

  #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
  vertex_weights_lh[,1] <- vertex_weights_lh[,1]*log(1.5)
  vertex_weights_rh[,1] <- vertex_weights_rh[,1]*log(1.5)
  
  #Mask by low signal regions
  vertex_weights_lh[all_zero_lh == 1] <- NA
  vertex_weights_rh[all_zero_rh == 1] <- NA
  
  #define map variable using hard parcel cifti file
  Pos_weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Pos_weight_map$data$cortex_left <- vertex_weights_lh
  Pos_weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Pos_weight_map_med <- move_from_mwall(Pos_weight_map)
  
  #define output file for averaged map
  outfile <- paste0(unthres_dirout,'POS_MAX_logTC_PFN_loading_gam_weights_LOG1.5_unthres')
  
  #save out cifti file
  write_cifti(Pos_weight_map_med,outfile)


#MAX NEG PFN BETA MAPS FOR TCA GAM MODELS, UNTHRESHOLDED
vertex_weights_lh <- matrix(0,29696,1)
vertex_weights_rh <- matrix(0,29716,1)

for (PFN in c(1:17))
{
  #read in weights and p-vals of GAMs for each PFN
  TCA_PFN_gam_results <- read.csv(sprintf("%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  #assign weights to left hemi
  for (i in c(1:29696)) {
    if (is.na(TCA_PFN_gam_results[i,2])) {
      next
    }
    if (TCA_PFN_gam_results[i,2] < 0 && TCA_PFN_gam_results[i,2] < vertex_weights_lh[i,1]) { # if weight is neg and less than current weight
      vertex_weights_lh[i,1] <- TCA_PFN_gam_results[i,2]
    }
  }
  
  #assign weights to right hemi
  for (i in c(1:29716)) {
    if (is.na(TCA_PFN_gam_results[i+29696,2])) {
      next
    }
    if (TCA_PFN_gam_results[i+29696,2] < 0 && TCA_PFN_gam_results[i+29696,2] < vertex_weights_rh[i,1]) { # if weight is neg and less than current weight
      vertex_weights_rh[i,1] <- TCA_PFN_gam_results[i+29696,2]
    }
  }
}

  #TRANSFORM BETA FOR INTERPRETABILITY (log(1.5) EQUIVALENT TO 1.5x TCA)
  vertex_weights_lh[,1] <- vertex_weights_lh[,1]*log(1.5)
  vertex_weights_rh[,1] <- vertex_weights_rh[,1]*log(1.5)
  
  #Mask by low signal regions
  vertex_weights_lh[all_zero_lh == 1] <- NA
  vertex_weights_rh[all_zero_rh == 1] <- NA
  
  #define map variable using hard parcel cifti file
  Neg_weight_map <- PFNs_hardparcel
  
  #assign weight values into hard parcel as template
  Neg_weight_map$data$cortex_left <- vertex_weights_lh
  Neg_weight_map$data$cortex_right <- vertex_weights_rh
  
  #add back in medial wall
  Neg_weight_map_med <- move_from_mwall(Neg_weight_map)
  
  #define output file for averaged map
  outfile <- paste0(unthres_dirout,'NEG_MAX_logTC_PFN_loading_gam_weights_LOG1.5_unthres')
  
  #save out cifti file
  write_cifti(Neg_weight_map_med,outfile)