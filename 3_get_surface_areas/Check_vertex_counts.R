library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)
library(dplyr)
library(rhdf5)
library(hdf5r)

#read in ciftis to check vertex counts (for purpose of -cifti-math and medial wall discrepancy)
vertex_area_cifti <- read_cifti('sub-997818717_vertex_area.dscalar.nii')
PFN_1_soft_parcel_cifti <-read_cifti('./Inputs/PFN1_soft_parcel.dscalar.nii')
 
#read in group atlas to look into if all loadings across networks sum to 1 for each vertex
readin <- readMat('../PNC_data/PNC_group_atlas/PNC_group_atlas.mat')
PNC_group_atlas <- readin$initV
vertex_sums <- rowSums(PNC_group_atlas)
vertex_sums <- data.frame(vertex_sums)

vertex_sums_nonzero <- vertex_sums %>% filter(vertex_sums != "0")
average_sum <- mean(vertex_sums_nonzero$vertex_sums)

ggplot(vertex_sums_nonzero, aes(x = vertex_sums)) +
  geom_histogram(binwidth = 0.025,    # Width of the bins
                 fill = "gray", # Fill color of the bars
                 color = "black") + # Border color of the bars
  labs(title = "Distribution of Loading Value Sums in PNC Group Atlas",
       x = "Loading Value Sum Across Networks",
       y = "Number of Vertices") +
  theme(plot.title = element_text(hjust = 0.5))

#including zeroes
ggplot(vertex_sums, aes(x = vertex_sums)) +
  geom_histogram(binwidth = 0.025,    # Width of the bins
                 fill = "gray", # Fill color of the bars
                 color = "black") + # Border color of the bars
  labs(title = "Distribution of Loading Value Sums in PNC Group Atlas",
       x = "Loading Value Sum Across Networks",
       y = "Number of Vertices") +
  theme(plot.title = element_text(hjust = 0.5))

#read in example pNet PFN to look into if all loadings across networks sum to 1 for each vertex
readin <- readMat('C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/Longitudinal_PFNs/19_PFNs_INV10-_101525/Personalized_FN/sub-NDARINV10EP1VM2_ses-baselineYear1Arm1_desc-merged_timeseries.dtseries.ni/1/FN.mat')
sub_NDARINV10EP1VM2_baseline_pNet <- readin$FN
vertex_sums_indiv_pnet <- rowSums(sub_NDARINV10EP1VM2_baseline_pNet)
vertex_sums_indiv_pnet <- data.frame(vertex_sums_indiv_pnet)

vertex_sums_indiv_pnet_nonzero <- vertex_sums_indiv_pnet %>% filter(vertex_sums_indiv_pnet != "0")
average_sum <- mean(vertex_sums_indiv_pnet_nonzero$vertex_sums_indiv_pnet)

ggplot(vertex_sums_indiv_pnet_nonzero, aes(x = vertex_sums_indiv_pnet)) +
  geom_histogram(binwidth = 0.025,    # Width of the bins
                 fill = "gray", # Fill color of the bars
                 color = "black") + # Border color of the bars
  labs(title = "Distribution of Loading Value Sums in ABCD Participant",
       x = "Loading Value Sum Across Networks",
       y = "Number of Vertices") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim =c(0,2))

#including zeroes
ggplot(vertex_sums_indiv_pnet, aes(x = vertex_sums_indiv_pnet)) +
  geom_histogram(binwidth = 0.025,    # Width of the bins
                 fill = "gray", # Fill color of the bars
                 color = "black") + # Border color of the bars
  labs(title = "Distribution of Loading Value Sums in ABCD Participant",
       x = "Loading Value Sum Across Networks",
       y = "Number of Vertices") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim =c(0,2))


#Map vertex sums of group atlas to cortex:

PFNs_hardparcel<-read_cifti('hardparcel_group.dscalar.nii')

#assign summed loadings to left hemi
vertex_sum_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_sum_lh[i,1] <- vertex_sums$vertex_sums[i]
}

#assign summed loadings to right hemi
vertex_sum_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_sum_rh[i,1] <- vertex_sums$vertex_sums[i+29696]
}

#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel

#assign summed loadings into hard parcel
Weight_map$data$cortex_left <- vertex_sum_lh
Weight_map$data$cortex_right <- vertex_sum_rh

#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)

#define output file for averaged map
outfile <- paste0('vertex_loading_sums')

#save out cifti files
write_cifti(Weight_map_med,outfile)



#Map vertex sums of group atlas to cortex:

PFNs_hardparcel<-read_cifti('hardparcel_group.dscalar.nii')

#assign summed loadings to left hemi
vertex_sum_lh <- matrix(0,29696,1)
for (i in c(1:29696)) {
  vertex_sum_lh[i,1] <- vertex_sums_indiv_pnet$vertex_sums_indiv_pnet[i]
}

#assign summed loadings to right hemi
vertex_sum_rh <- matrix(0,29716,1)
for (i in c(1:29716)) {
  vertex_sum_rh[i,1] <- vertex_sums_indiv_pnet$vertex_sums_indiv_pnet[i+29696]
}

#define map variable using hard parcel cifti file
Weight_map <- PFNs_hardparcel

#assign summed loadings into hard parcel
Weight_map$data$cortex_left <- vertex_sum_lh
Weight_map$data$cortex_right <- vertex_sum_rh

#add back in medial wall
Weight_map_med <- move_from_mwall(Weight_map)

#define output file for averaged map
outfile <- paste0('vertex_loading_sums_sub_NDARINV10EP1VM2_baseline_pNet')

#save out cifti files
write_cifti(Weight_map_med,outfile)

