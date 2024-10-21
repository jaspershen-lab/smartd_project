library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidymass)
##load data
load("3_data_analysis/1_data_preparation/0_demographic_data/demographic_data.rda")

load("3_data_analysis/1_data_preparation/0_demographic_data/sample_information.rda")

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

dim(sample_info)
dim(variable_info)
dim(expression_data)

dir.create("3_data_analysis/1_data_preparation/1_urine_metabolomics_data", recursive = TRUE)
setwd("3_data_analysis/1_data_preparation/1_urine_metabolomics_data")

variable_info <-
  variable_info %>%
  dplyr::rename(variable_id = name)

sample_info <-
  sample_info[, c("sample_id", "injection.order", "class", "batch", "group")] %>%
  left_join(sample_information, by = "sample_id") %>%
  dplyr::select(
    c(
      sample_id,
      subject_id,
      injection.order,
      class,
      batch,
      group,
      Visit,
      Date.Acquired,
      GA,
      enrollment_date,
      mother_dob:mother_delivery_weeks,
      dplyr::everything()
    )
  )

sample_info$sample_id == colnames(expression_data)

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

dim(urine_metabolomics_data)

save(urine_metabolomics_data, file = "metabolites/urine_metabolomics_data.rda")
