library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

###read data
load("3_data_analysis/1_data_preparation/0_demographic_data/demographic_data.rda")

dir.create("4_manuscript/Supplementary_data/", recursive = TRUE)

setwd("4_manuscript/Supplementary_data/")

library(ggExtra)

## we have
colnames(demographic_data)

data <-
  demographic_data %>%
  dplyr::select(
    subject_id,
    mother_age,
    mother_induction,
    child_sex,
    child_weight,
    mother_height,
    mother_weight,
    mother_bmi,
    mother_ethnicity,
    mother_ethnicity,
    mother_delivery_weeks
  )

openxlsx::write.xlsx(data, "Table S1.xlsx")
