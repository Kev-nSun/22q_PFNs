# Fixes the Matlab output files (*.dlabel.nii) which are missing the medial wall mask values
# Based on hardparcel_group.dscalar.nii as template

library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

#medial wall mask template is hard parcel file
mwall_template <- read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')
dir <- paste0('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_results/NORMED/Maps/PFN_ID_clusts/')

for (i in c(1:2)) {
  if (i == 1) {clust = "POS"}
  else if (i == 2) {clust = "NEG"}
  
  map_name <- paste0('logTC_',clust,"_PFN_ID_clust_main_effects_50mm")
  #map_name <- paste0('logTC_',clust,"_PFN_ID_clust_50mm")
  map <- read_cifti(paste0(dir,'/',map_name,'_AtlasLabel.dlabel.nii'))
  
  #add in mwall mask values based on template
  map$meta$cortex$medial_wall_mask <- mwall_template$meta$cortex$medial_wall_mask
  
  #save out new cifti file
  outfile <- paste0(dir,"/",map_name,'_mwall')
  write_cifti(map,outfile)
}



