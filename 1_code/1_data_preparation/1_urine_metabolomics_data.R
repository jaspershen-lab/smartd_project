library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

###load data
# load(
#   "3_data_analysis/data_analysis20200108/urine_metabolome/sample_information/sample_information"
# )
#
# sample_info <-
#   sample_information %>%
#   dplyr::rename(
#     sample_id = Sample_ID,
#     subject_id = Patient_ID,
#     collection_date = Date.Acquired,
#     gestational_age = GA,
#     visit_time_point = Visit
#   )

####metabolomics data
load(
  "3_data_analysis/data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data"
)
load(
  "3_data_analysis/data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info"
)
load(
  "3_data_analysis/data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info"
)

library(tidymass)
dim(sample_info)
dim(variable_info)
dim(expression_data)

dir.create("3_data_analysis/1_urine_metabolomics_data", recursive = TRUE)
setwd("3_data_analysis/1_urine_metabolomics_data")

variable_info <-
  variable_info %>%
  dplyr::rename(variable_id = name)

urine_metabolomics_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

###peaks

urine_metabolomics_data

dir.create("peaks")

save(urine_metabolomics_data, file = "peaks/urine_metabolomics_data.rda")

###metabolites
dir.create("metabolites")

urine_metabolomics_data <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(Compound.name))

save(urine_metabolomics_data, file = "metabolites/urine_metabolomics_data.rda")
