
library(R.matlab)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

# This script creates data for 2 maps of PFN IDS- one corresponding to positive beta clusters of PFN loading ~ 22q, and one for negative beta clusters

library(R.matlab)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

# This script creates data for 2 maps of PFN IDS- one corresponding to positive beta clusters of PFN loading ~ 22q, and one for negative beta clusters

dirin <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/091426_22q_no_TCA_vert_normed_loadings_FDR_pval_beta_subclass/"
dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_no_TCA/Maps/"
dir.create(dirout)

bin_dirout <- sprintf("%s/Binary_cluster_maps/",dirout)

clust_PFN_ID_dir <- sprintf("%s/PFN_ID_clusts/",dirout)
dir.create(clust_PFN_ID_dir)

POS_vertex_weight_clust_lh <- matrix(0,29696)
NEG_vertex_weight_clust_lh <- matrix(0,29696)
POS_vertex_weight_clust_rh <- matrix(0,29716)
NEG_vertex_weight_clust_rh <- matrix(0,29716)

for (PFN in c(1:17)) {
  #read in weights and p-vals of GAMs for each PFN
  stat_22q_PFN_gam_results <- read.csv(sprintf("%s/stat_22q_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  clust_file <- paste0(bin_dirout,'/PFN',PFN,'_50mm_cluster_binary_22q_no_TCA.dscalar.nii')
  
  #read in clustered vertices
  clust_cifti<-read_cifti(clust_file) #cluster binary cifti
  
  #L hemi
  clusts_lh <- clust_cifti$data$cortex_left
  for (i in c(1:29696)) {
    if (clusts_lh[i] > 0) { # if vertex passed cluster threshold
      if (stat_22q_PFN_gam_results[i,2] > 0) { # if vertex beta is POS, add to POS cluster map
        POS_vertex_weight_clust_lh[i] <- PFN
      } else if (stat_22q_PFN_gam_results[i,2] < 0) { # if vertex beta is NEG, add to NEG cluster map
        NEG_vertex_weight_clust_lh[i] <- PFN
      }
    }
  }
  #R hemi
  clusts_rh <- clust_cifti$data$cortex_right
  for (i in c(1:29716)) {
    if (clusts_rh[i] > 0){ # if vertex passed cluster threshold
      if (stat_22q_PFN_gam_results[i+29696,2] > 0) { # if vertex beta is POS, add to POS cluster map
        POS_vertex_weight_clust_rh[i,1] <- PFN
      } else if (stat_22q_PFN_gam_results[i+29696,2] < 0) { # if vertex beta is NEG, add to NEG cluster map
        NEG_vertex_weight_clust_rh[i] <- PFN
      }
    }
  }
}

#For POS clusters, save clustered left + right hem to mat files, converted to gifti--> cifti in matlab "Labeling_PFNs_into_cifti.m"

# For L hemi, ensure that both the list and its elements are named
pos_names_list <- list()
pos_names_list$POS_vertex_weight_clust_lh <- POS_vertex_weight_clust_lh
writeMat(paste0(clust_PFN_ID_dir,"22q_POS_PFN_ID_clust_50mm_no_TCA_L.mat"), x = pos_names_list)
# For R hemi, ensure that both the list and its elements are named
pos_names_list <- list()
pos_names_list$POS_vertex_weight_clust_rh <- POS_vertex_weight_clust_rh
writeMat(paste0(clust_PFN_ID_dir,"22q_POS_PFN_ID_clust_50mm_no_TCA_R.mat"), x = pos_names_list)

#For NEG clusters, save clustered left + right hem to mat files, converted to gifti--> cifti in matlab "Labeling_PFNs_into_cifti.m"
# For L hemi, ensure that both the list and its elements are named
neg_names_list <- list()
neg_names_list$NEG_vertex_weight_clust_lh <- NEG_vertex_weight_clust_lh
writeMat(paste0(clust_PFN_ID_dir,"22q_NEG_PFN_ID_clust_50mm_no_TCA_L.mat"), x = neg_names_list)
# For R hemi, ensure that both the list and its elements are named
neg_names_list <- list()
neg_names_list$NEG_vertex_weight_clust_rh <- NEG_vertex_weight_clust_rh
writeMat(paste0(clust_PFN_ID_dir,"22q_NEG_PFN_ID_clust_50mm_no_TCA_R.mat"), x = neg_names_list)