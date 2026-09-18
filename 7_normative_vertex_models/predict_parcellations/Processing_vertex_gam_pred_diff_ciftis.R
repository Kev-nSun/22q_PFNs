# Fixes the Matlab output files (*.dlabel.nii) which are missing the medial wall mask values
# Based on hardparcel_group.dscalar.nii as template

library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')

#medial wall mask template is hard parcel file
mwall_template <- read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')
dir <- paste0(sprintf('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/Normed_delta_maps'))

for (delta in c("stable","delta_1st","delta_99th")) {
  #read in map that needs medial wall mask values
  map_name <- paste0(delta,'_scaling_pred_topo')
  map <- read_cifti(paste0(dir,'/',map_name,'_AtlasLabel.dlabel.nii'))
  
  #add in mwall mask values based on template
  map$meta$cortex$medial_wall_mask <- mwall_template$meta$cortex$medial_wall_mask
  
  #save out new cifti file
  outfile <- paste0(dir,"/",map_name,'_mwall')
  write_cifti(map,outfile)
}



