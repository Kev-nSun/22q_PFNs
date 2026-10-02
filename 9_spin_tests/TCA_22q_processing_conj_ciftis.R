# Fixes the Matlab output files (*.dlabel.nii) which are missing the medial wall mask values
# Based on hardparcel_group.dscalar.nii as template

library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

#medial wall mask template is hard parcel file
mwall_template <- read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')

maps_list <- list("con_pos_pos","con_neg_neg","div_pos_neg","div_neg_pos")
for (i in seq_along(maps_list)) {
  map_name <- paste0('TCAx22q_',maps_list[[i]],"_22q_PFN_ID_clust_50mm")
  dir <- paste0('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_22q_PFN_conjunction/')
  map <- read_cifti(paste0(dir,'/',map_name,'_AtlasLabel.dlabel.nii'))
  
  #add in mwall mask values based on template
  map$meta$cortex$medial_wall_mask <- mwall_template$meta$cortex$medial_wall_mask
  
  #save out new cifti file
  outfile <- paste0(dir,"/",map_name,'_mwall')
  write_cifti(map,outfile)
}



