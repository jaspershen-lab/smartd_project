# ##avoid source
# no_function()
# 
# setwd(r4projects::get_project_wd())
# rm(list = ls())
# library(tidyverse)
# source("1_code/100_tools.R")
# 
# ###load data
# load(
#   "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/metabolites/urine_metabolomics_data.rda"
# )
# 
# ##first set the work directory to project folder
# dir.create("3_data_analysis/4_prediction/1_data_preparation",
#            recursive = TRUE)
# setwd("3_data_analysis/4_prediction/1_data_preparation")
# 
# ###remove qc samples and postpartum samples
# urine_metabolomics_data@sample_info$GA
# urine_metabolomics_data2 <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!is.na(GA))
# 
# expression_data <-
#   extract_expression_data(urine_metabolomics_data2)
# sample_info <-
#   extract_sample_info(urine_metabolomics_data2)
# variable_info <-
#   extract_variable_info(urine_metabolomics_data2)
# 
# ##P samples have no GA information
# ####log 10 and scale
# subject_data <-
#   log(expression_data, 10) %>%
#   as.data.frame()
# 
# subject_data <-
#   apply(subject_data, 1, function(x) {
#     (x - mean(x)) / sd(x)
#   })
# 
# 
# #####################################################################################################
# #####get discovery and validation data###############################################################
# #####################################################################################################
# ###randomly divided into discovery and valdiation datasets
# # set.seed(seed = 1)
# sample_info$subject_id %>% unique() %>%
#   length()
# 
# ##we have 36 people in total
# index_dis <- which(sample_info$batch == 1)
# index_val <- which(sample_info$batch == 2)
