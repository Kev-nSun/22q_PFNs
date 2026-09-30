
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirout <- "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/MAD_map/"
dir.create(dirout)

PFN_template<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/hardparcel_group.dscalar.nii')
readin <- readMat("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/MAD_map/PFN_loading_normed_MAD.mat")
MAD_values <- readin$MAD
MAD_row <- t(MAD_values)

#GET ALL ZERO MASK FOR VISUALIZATION
PFNs_all_zero_mask<-read_cifti('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/PNC_PFN_all_zero_mask.dscalar.nii')
all_zero_lh <- PFNs_all_zero_mask$data$cortex_left
all_zero_rh <- PFNs_all_zero_mask$data$cortex_right

#Individual PFN maps of MAD
for (PFN in c(1:17))
{
  #assign MAD to left hemi
  vertex_lh <- matrix(0,29696)
  for (vert in c(1:29696)) {
    vertex_lh[vert] <- MAD_row[(PFN-1)*59412+vert]
  }
  #assign MAD to right hemi
  vertex_rh <- matrix(0,29716)
  for (vert in c(1:29716)) {
    vertex_rh[vert] <- MAD_row[(PFN-1)*59412+(vert+29696)]
  }
  
  #define map variable using hard parcel cifti file
  MAD_map <- PFN_template
  
  #assign MAD values into hard parcel
  MAD_map$data$cortex_left <- vertex_lh
  MAD_map$data$cortex_right <- vertex_rh
  
  #add back in medial wall
  MAD_map_med <- move_from_mwall(MAD_map)
  
  #define output file for averaged map
  outfile <- paste0(dirout,'PFN_',PFN,'_MAD')
  
  #save out cifti files
  write_cifti(MAD_map_med,outfile)
}


#Sum MAD values across 17 networks
summed_MAD <- matrix(0,59412)
for (vert in c(1:59412)) {
  for (PFN in c(1:17)) {
    summed_MAD[vert] <- summed_MAD[vert] + MAD_row[(PFN-1)*59412+vert]
  }
}
write.csv(summed_MAD,file = paste0(dirout,"Summed_MAD.csv"))

#Avg MAD
avg_MAD <- summed_MAD/17
write.csv(avg_MAD,file = paste0(dirout,"Avg_MAD.csv"))

#assign avg MAD to left hemi
vertex_lh <- matrix(0,29696)
for (vert in c(1:29696)) {
  vertex_lh[vert] <- avg_MAD[vert]
}
#assign avg MAD to right hemi
vertex_rh <- matrix(0,29716)
for (vert in c(1:29716)) {
  vertex_rh[vert] <- avg_MAD[vert+29696]
}

#Mask by low signal regions
vertex_lh[all_zero_lh == 1] <- NA
vertex_rh[all_zero_rh == 1] <- NA

#define map variable using hard parcel cifti file
MAD_avg_map <- PFN_template

#assign MAD values into hard parcel
MAD_avg_map$data$cortex_left <- vertex_lh
MAD_avg_map$data$cortex_right <- vertex_rh

#add back in medial wall
MAD_avg_map_med <- move_from_mwall(MAD_avg_map)

#define output file for averaged map
outfile <- paste0(dirout,'MAD_Avg')

#save out cifti files
write_cifti(MAD_avg_map_med,outfile)



#Max MAD values across 17 networks
max_MAD <- matrix(0,59412)
for (vert in c(1:59412)) {
  max_tracker <- 0
  for (PFN in c(1:17)) {
    if (MAD_row[(PFN-1)*59412+vert] > max_tracker) {
      max_tracker <- MAD_row[(PFN-1)*59412+vert]
    }
  }
  max_MAD[vert] <- max_tracker
}
write.csv(max_MAD,file = paste0(dirout,"Max_MAD.csv"))

#assign summed MAD to left hemi
vertex_lh <- matrix(0,29696)
for (vert in c(1:29696)) {
  vertex_lh[vert] <- max_MAD[vert]
}
#assign summed MAD to right hemi
vertex_rh <- matrix(0,29716)
for (vert in c(1:29716)) {
  vertex_rh[vert] <- max_MAD[vert+29696]
}

#Mask by low signal regions
vertex_lh[all_zero_lh == 1] <- NA
vertex_rh[all_zero_rh == 1] <- NA

#define map variable using hard parcel cifti file
MAD_max_map <- PFN_template

#assign MAD values into hard parcel
MAD_max_map$data$cortex_left <- vertex_lh
MAD_max_map$data$cortex_right <- vertex_rh

#add back in medial wall
MAD_max_map_med <- move_from_mwall(MAD_max_map)

#define output file for averaged map
outfile <- paste0(dirout,'MAD_Max')

#save out cifti files
write_cifti(MAD_max_map_med,outfile)


#Binary map of zero vs non-zero variability
MAD_binary <- matrix(0,59412)
for (vert in c(1:59412)) {
  if (max_MAD[vert] == 0) {
    MAD_binary[vert] <- 1
  }
}

#assign binary MAD to left hemi
vertex_lh <- matrix(0,29696)
for (vert in c(1:29696)) {
  vertex_lh[vert] <- MAD_binary[vert]
}
#assign binary MAD to right hemi
vertex_rh <- matrix(0,29716)
for (vert in c(1:29716)) {
  vertex_rh[vert] <- MAD_binary[vert+29696]
}

#define map variable using hard parcel cifti file
MAD_binary_map <- PFN_template

#assign binary values into hard parcel
MAD_binary_map$data$cortex_left <- vertex_lh
MAD_binary_map$data$cortex_right <- vertex_rh

#add back in medial wall
MAD_binary_map <- move_from_mwall(MAD_binary_map)

#define output file for averaged map
outfile <- paste0(dirout,'MAD_Binary')

#save out cifti files
write_cifti(MAD_binary_map,outfile)