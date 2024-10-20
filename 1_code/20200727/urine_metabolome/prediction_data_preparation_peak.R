#to avoid source
no_exist_function()

library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)

rm(list = ls())
sxtTools::setwd_project()
source("R/20200727/tools.R")

##load data peaks
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info")

sxtTools::setwd_project()
##first set the work directory to project folder
setwd("data_analysis20200108/urine_metabolome/prediction/feature/")

subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id[sample_info$class == "Subject"]))

qc_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id[sample_info$class == "QC"]))

####log 10 and scale
subject_data <- 
  log(subject_data, 10) %>% 
  as.data.frame()

subject_data <-
  apply(subject_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

colnames(subject_data) <- variable_info$name

subject_data <- 
  subject_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample_id")

qc_data <- 
  log(qc_data, 10) %>% 
  as.data.frame()

qc_data <-
  apply(qc_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

colnames(qc_data) <- variable_info$name

qc_data <- 
  qc_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample_id")


###add the patient information to the sample_data
subject_data <- 
  inner_join(x = sample_info[,c(1,2,3,4,5,6,7,8,9)], 
             y = subject_data, 
             by = "sample_id")

save(qc_data, file = "qc_data")


###remnove samples with GA == 0
###sample_data_old is the data that contains PP samples

subject_data$GA[is.na(subject_data$GA)] <- 0

subject_data_old <- 
  subject_data

subject_data <- 
  subject_data %>% 
  dplyr::filter(GA != 0)


#####################################################################################################
#####get discovery and validation data###############################################################
#####################################################################################################
###randomly divided into discovery and valdiation datasets
# set.seed(seed = 1)
subject_data$subject_id %>% unique()

##we have 36 peple in total

index_dis <- which(subject_data$sample_id %in% sample_info$sample_id[sample_info$batch == 1])
index_val <- which(subject_data$sample_id %in% sample_info$sample_id[sample_info$batch == 2])

save(index_dis, file = "index_dis")
save(index_val, file = "index_val")

sample_data_dis <- 
  subject_data[index_dis,]

sample_data_val <- 
  subject_data[index_val,]

sample_data_dis_x <- 
  sample_data_dis %>% 
  dplyr::select(-(sample_id:Date.Acquired)) %>% 
  as.matrix()

sum(is.na(sample_data_dis_x))
sample_data_dis_x[,1]


sample_data_val_x <- 
  sample_data_val %>% 
  dplyr::select(-(sample_id:Date.Acquired)) %>% 
  as.matrix()

sum(is.na(sample_data_dis_x))
sample_data_val_x[,1]

save(sample_data_dis, file = "sample_data_dis")
save(sample_data_val, file = "sample_data_val")

save(sample_data_dis_x, file = "sample_data_dis_x")
save(sample_data_val_x, file = "sample_data_val_x")




