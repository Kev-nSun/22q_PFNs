library(readr)
library(dplyr)
library(stringr)

validation<-read_tsv("../v0_validation.tsv")
code_count<-as.data.frame(table(validation$code))
write.csv(code_count,"../v0_validation_codes.csv")

#code breakdown by file
code_count_mod <- validation %>%
  group_by(code) %>%
  summarise(
    fmap_mag1 = sum(str_detect(location, "magnitude1")),
    fmap_mag2 = sum(str_detect(location, "magnitude2")),
    fmap_phase1 = sum(str_detect(location, "phase1")),
    fmap_phase2 = sum(str_detect(location, "phase2")),
    fmap_phasediff = sum(str_detect(location, "phasediff")),
    func_idemo = sum(str_detect(location, "task-idemo")),
    func_jolo = sum(str_detect(location, "task-jolo")),
    func_rest = sum(str_detect(location, "task-rest")),
    anat_T1w = sum(str_detect(location, "T1w")),
    .groups = "drop"
  )
write.csv(code_count_mod,"../v0_codes_by_mod.csv")

#SIDECAR_KEY_RECOMMENDED breakdown
validation.sidecar_rec<-validation[validation$code=="SIDECAR_KEY_RECOMMENDED",]
sidecar_keys_count<-as.data.frame(table(validation.sidecar_rec$subCode))
write.csv(sidecar_keys_count,"../v0_SIDECAR_KEY_RECOMMENDED_subCodes.csv")

#SIDECAR_KEY_RECOMMENDED breakdown by file
sidecar_keys_count_mod <- validation.sidecar_rec %>%
  group_by(subCode) %>%
  summarise(
    fmap_mag1 = sum(str_detect(location, "magnitude1")),
    fmap_mag2 = sum(str_detect(location, "magnitude2")),
    fmap_phase1 = sum(str_detect(location, "phase1")),
    fmap_phase2 = sum(str_detect(location, "phase2")),
    fmap_phasediff = sum(str_detect(location, "phasediff")),
    func_idemo = sum(str_detect(location, "task-idemo")),
    func_jolo = sum(str_detect(location, "task-jolo")),
    func_rest = sum(str_detect(location, "task-rest")),
    anat_T1w = sum(str_detect(location, "T1w")),
    .groups = "drop"
  )
write.csv(sidecar_keys_count_mod,"../v0_SIDECAR_KEY_RECOMMENDED_subCodes_by_mod.csv")

#SIDECAR_FIELD_OVERRIDE breakdown
validation.sidecar_override<-validation[validation$code=="SIDECAR_FIELD_OVERRIDE",]
sidecar_override_count<-as.data.frame(table(validation.sidecar_override$subCode))
write.csv(sidecar_override_count,"../v0_SIDECAR_FIELD_OVERRIDE_subCodes.csv")

#TSV_ADDITIONAL_COLUMNS_UNDEFINED breakdown
validation.tsv_undefined<-validation[validation$code=="TSV_ADDITIONAL_COLUMNS_UNDEFINED",]
tsv_undefined_count<-as.data.frame(table(validation.tsv_undefined$subCode))

