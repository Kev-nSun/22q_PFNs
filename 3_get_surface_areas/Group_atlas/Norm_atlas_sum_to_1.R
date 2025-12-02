library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)
library(dplyr)
library(ggplot2)
library(rhdf5)
library(hdf5r)

#read in group atlas to look into if all loadings across networks sum to 1 for each vertex
readin <- readMat('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PFNs/22q_PFNs/Personalized_FN/sub-011706/FN.mat')
PNC_group_atlas <- readin$initV
vertex_sums <- data.frame(rowSums(PNC_group_atlas))

#Divide loadings by the sum of loadings across networks for each vertex
PNC_group_atlas_normed <- PNC_group_atlas/vertex_sums$rowSums.PNC_group_atlas.
vertex_sums_normed <- data.frame(rowSums(PNC_group_atlas_normed))

write.csv(PNC_group_atlas_normed,"C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/PNC_group_atlas_normed/PNC_group_atlas_normed.csv")

#Replace Na with 0
PNC_group_atlas_normed[is.na(PNC_group_atlas_normed)] <- 0

PFNs_hardparcel<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii') #hard parcel cifti used as template

for (net in c(1:17)) {
  
  #assign network values to left hemi
  vertex_weights_lh <- matrix(0,29696,1)
  for (i in c(1:29696)) {
    vertex_weights_lh[i,1] <- PNC_group_atlas_normed[i,net]
  }
  
  #assign network values to right hemi
  vertex_weights_rh <- matrix(0,29716,1)
  for (i in c(1:29716)) {
    vertex_weights_rh[i,1] <- PNC_group_atlas_normed[i+29696,net]
  }
  
  #define map variable using hard parcel cifti file
  Soft_map <- PFNs_hardparcel
  
  #assign network values for each hemi into cifti 
  Soft_map$data$cortex_left <- vertex_weights_lh
  Soft_map$data$cortex_right <- vertex_weights_rh
  
  outfile <- paste0('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/PNC_group_atlas_normed/PFN', net, '_soft_parcel_normed')
  write_cifti(Soft_map,outfile) # save out cifti for each PFN
}

