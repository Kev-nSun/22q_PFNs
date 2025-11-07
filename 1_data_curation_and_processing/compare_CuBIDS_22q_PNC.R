library(readr)
library(dplyr)
library(stringr)

DS_22q_CuBIDS <- read_tsv("../v2_summary_edited.tsv")
DS_study_col <- data.frame(Study = rep("22q",nrow(DS_22q_CuBIDS)))
DS_22q_CuBIDS <- cbind(DS_study_col,DS_22q_CuBIDS)

PNC_CuBIDS <- read_tsv("../study-PNC_desc-CuBIDS_summary.tsv")
PNC_study_col <- data.frame(Study = rep("PNC",nrow(PNC_CuBIDS)))
PNC_CuBIDS <- cbind(PNC_study_col,PNC_CuBIDS)

shared_cols <- intersect(names(DS_22q_CuBIDS),names(PNC_CuBIDS))

merged_CuBIDS <- bind_rows(DS_22q_CuBIDS, PNC_CuBIDS)
merged_CuBIDS <- merged_CuBIDS[, shared_cols]
merged_CuBIDS <- merged_CuBIDS[, -c(2,3,4)]

 
T1w <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "T1w"))
write.csv(T1w,"../T1w_22q_vs_PNC_CuBIDS.csv")

idemo <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "idemo"))
write.csv(idemo,"../task-idemo_22q_vs_PNC_CuBIDS.csv")

rest <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "task-rest"))
write.csv(rest,"../task-rest_22q_vs_PNC_CuBIDS.csv")
  
mag1 <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "magnitude1"))
write.csv(mag1,"../magnitude1_22q_vs_PNC_CuBIDS.csv")
  
mag2 <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "magnitude2"))
write.csv(mag2,"../magnitude2_22q_vs_PNC_CuBIDS.csv")
  
phase1 <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "phase1"))
write.csv(phase1,"../phase1_22q_vs_PNC_CuBIDS.csv")
  
phase2 <- merged_CuBIDS %>%
  filter(str_detect(KeyParamGroup, "phase2"))
write.csv(phase2,"../phase2_22q_vs_PNC_CuBIDS.csv")
