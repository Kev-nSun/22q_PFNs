
library(R.matlab)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

# This script creates data for 2 maps of PFN IDS- one corresponding to positive beta clusters of PFN loading ~ TCA, and one for negative beta clusters (all effects)

dirin <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/091426_TCA_vert_normed_loadings_FDR_pval_beta_logTC/"
bin_dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/NORMED/Maps/Binary_cluster_maps/"
clust_PFN_ID_dir <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/NORMED/Maps/PFN_ID_clusts/"
dir.create(clust_PFN_ID_dir)

POS_vertex_weight_clust_lh <- matrix(0,29696,2)
NEG_vertex_weight_clust_lh <- matrix(0,29696,2)
POS_vertex_weight_clust_rh <- matrix(0,29716,2)
NEG_vertex_weight_clust_rh <- matrix(0,29716,2)

for (PFN in c(1:17)) {
  #read in weights and p-vals of GAMs for each PFN
  TCA_PFN_gam_results <- read.csv(sprintf("%s/logTC_gam_vert_loadings_PFN_%02d_beta_fdr.csv",dirin,PFN))
  
  clust_file <- paste0(bin_dirout,'/PFN',PFN,'_50mm_cluster_binary_TCA.dscalar.nii')
  
  #read in clustered vertices
  clust_cifti<-read_cifti(clust_file) #cluster binary cifti
  
  #L hemi
  clusts_lh <- clust_cifti$data$cortex_left
  for (i in c(1:29696)) {
    if (clusts_lh[i] > 0) { # if vertex passed cluster threshold
      if (TCA_PFN_gam_results[i,2] > 0 && TCA_PFN_gam_results[i,2] > POS_vertex_weight_clust_lh[i,2]) { # if vertex beta is POS & greater than current beta, add to POS cluster map
        POS_vertex_weight_clust_lh[i,1] <- PFN
        POS_vertex_weight_clust_lh[i,2] <- TCA_PFN_gam_results[i,2]
      } else if (TCA_PFN_gam_results[i,2] < 0 && TCA_PFN_gam_results[i,2] < NEG_vertex_weight_clust_lh[i,2]) { # if vertex beta is NEG & more neg than current beta, add to NEG cluster map
        NEG_vertex_weight_clust_lh[i,1] <- PFN
        NEG_vertex_weight_clust_lh[i,2] <- TCA_PFN_gam_results[i,2]
      }
    }
  }
  #R hemi
  clusts_rh <- clust_cifti$data$cortex_right
  for (i in c(1:29716)) {
    if (clusts_rh[i] > 0){ # if vertex passed cluster threshold
      if (TCA_PFN_gam_results[i+29696,2] > 0 && TCA_PFN_gam_results[i+29696,2] > POS_vertex_weight_clust_rh[i,2]) { # if vertex beta is POS & greater than current beta, add to POS cluster map
        POS_vertex_weight_clust_rh[i,1] <- PFN
        POS_vertex_weight_clust_rh[i,2] <- TCA_PFN_gam_results[i+29696,2]
      } else if (TCA_PFN_gam_results[i+29696,2] < 0 && TCA_PFN_gam_results[i+29696,2] < NEG_vertex_weight_clust_rh[i,2]) { # if vertex beta is NEG & more neg than current beta, add to NEG cluster map
        NEG_vertex_weight_clust_rh[i,1] <- PFN
        NEG_vertex_weight_clust_rh[i,2] <- TCA_PFN_gam_results[i+29696,2]
      }
    }
  }
}

#For POS clusters, save clustered left + right hem to mat files, converted to gifti--> cifti in matlab "Labeling_PFNs_into_cifti.m"

# For L hemi, ensure that both the list and its elements are named
pos_names_list <- list()
pos_names_list$POS_vertex_weight_clust_lh <- POS_vertex_weight_clust_lh[,1]
writeMat(paste0(clust_PFN_ID_dir,"logTC_POS_PFN_ID_clust_50mm_L.mat"), x = pos_names_list)
# For R hemi, ensure that both the list and its elements are named
pos_names_list <- list()
pos_names_list$POS_vertex_weight_clust_rh <- POS_vertex_weight_clust_rh[,1]
writeMat(paste0(clust_PFN_ID_dir,"logTC_POS_PFN_ID_clust_50mm_R.mat"), x = pos_names_list)

#For NEG clusters, save clustered left + right hem to mat files, converted to gifti--> cifti in matlab "Labeling_PFNs_into_cifti.m"
# For L hemi, ensure that both the list and its elements are named
neg_names_list <- list()
neg_names_list$NEG_vertex_weight_clust_lh <- NEG_vertex_weight_clust_lh[,1]
writeMat(paste0(clust_PFN_ID_dir,"logTC_NEG_PFN_ID_clust_50mm_L.mat"), x = neg_names_list)
# For R hemi, ensure that both the list and its elements are named
neg_names_list <- list()
neg_names_list$NEG_vertex_weight_clust_rh <- NEG_vertex_weight_clust_rh[,1]
writeMat(paste0(clust_PFN_ID_dir,"logTC_NEG_PFN_ID_clust_50mm_R.mat"), x = neg_names_list)