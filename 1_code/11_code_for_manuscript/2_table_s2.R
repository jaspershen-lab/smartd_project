####for RF
setwd(r4projects::get_project_wd())

table2 <- readxl::read_xlsx("4_manuscript/Supplementary_data/Table S1.xlsx")

setwd(
  "2_data/data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/"
)
###normal data
library(randomForest)

##load data
load("../../../sample_data_dis")
load("../../../sample_data_val")
load("../../../sample_data_dis_x")
load("../../../sample_data_val_x")
load("../../../metabolite_tags")

##############################################################################
####random forest
#############################################################################
library(randomForest)
##use boruta method
library(Boruta)
library(tidyverse)

sample_data_dis_y <-
  sample_data_dis %>%
  dplyr::select(GA) %>%
  as.matrix()

sample_data_val_y <-
  sample_data_val %>%
  dplyr::select(GA) %>%
  as.matrix()

marker_rf <-
  readr::read_csv("marker_rf_final.csv")

info <-
  readxl::read_xlsx(
    "../../../../../../../2_data/patient_information/SmartD_ClinicalVariables_PartiallySummarized.xlsx"
  )

info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

info <-
  info %>%
  filter(!is.na(`C/S`))

info <-
  info %>%
  filter(`C/S` == "N")

idx_dis <-
  which(sample_data_dis$Patient_ID %in% info$ID)


idx_val <-
  which(sample_data_val$Patient_ID %in% info$ID)


sample_data_dis <-
  sample_data_dis[idx_dis, ]

sample_data_dis_x <-
  sample_data_dis_x[idx_dis, ]

sample_data_val <-
  sample_data_val[idx_val, ]

sample_data_val_x <-
  sample_data_val_x[idx_val, ]

unique(sample_data_dis$Patient_ID)

unique(sample_data_val$Patient_ID)

table2 <-
  table2 %>%
  dplyr::filter(subject_id %in% c(
    unique(sample_data_dis$Patient_ID),
    unique(sample_data_val$Patient_ID)
  ))

setwd(r4projects::get_project_wd())
setwd("4_manuscript/Supplementary_data/")
openxlsx::write.xlsx(table2, "Table S2.xlsx")