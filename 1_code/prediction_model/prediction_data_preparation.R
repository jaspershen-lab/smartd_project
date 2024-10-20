library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
##first set the work directory to project folder
setwd(r4projects::get_project_wd())
setwd("data_analysis20191015/prediction/identification_table/")

#####################################################################################################
#####data preparation###############################################################################
#####################################################################################################

###This is batch data set
smartd_rplc_batch1 <- 
  readr::read_csv("smartd_rplc_batch1.csv")

smartd_rplc_batch2 <- 
  readr::read_csv("ms1_data_rplc_batch2.csv")

###combine batch 1 and batch 2 data
##read the batch 1 and batch 2 match information
rplc_pos_batch1_batch2_matched_info <- 
  readr::read_csv("rplc_pos_batch1_batch2_matched_info.csv")

rplc_neg_batch1_batch2_matched_info <- 
  readr::read_csv("rplc_neg_batch1_batch2_matched_info.csv")

rplc_matched_info <- rbind(rplc_pos_batch1_batch2_matched_info[,-1],
                           rplc_neg_batch1_batch2_matched_info[,-1])


index1 <- 
  rplc_matched_info$batch1 %in%  smartd_rplc_batch1$name %>% 
  which()

index2 <- 
  rplc_matched_info$batch2 %in%  smartd_rplc_batch2$name %>% 
  which()

index <- intersect(index1, index2)

rplc_matched_info <- 
  rplc_matched_info[index,]


smartd_rplc_batch1 <- 
  smartd_rplc_batch1[match(rplc_matched_info$batch1, smartd_rplc_batch1$name),]

smartd_rplc_batch2 <- 
  smartd_rplc_batch2[match(rplc_matched_info$batch2, smartd_rplc_batch2$name),]

##remove QC and blank form batch 2 dataset
smartd_rplc_batch2 <-
  smartd_rplc_batch2 %>% 
  select(-contains("blk")) %>% 
  select(-contains("QC"))

sample_removed <-
  which(is.na(smartd_rplc_batch2), arr.ind = TRUE)[,2] %>% 
  unique() %>% 
  `[`(colnames(smartd_rplc_batch2), .)

##These samples may have positive or negative samples but don't have, so we should remove them
smartd_rplc_batch2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(-sample_removed)

##remove the tags of metabolites
smartd_rplc_batch2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(-(name:rt))

###rename P samples in batch2
colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] <-
  colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] %>% 
  stringr::str_replace("P", "") %>% 
  as.numeric() %>% 
  `+`(177) %>% 
  paste("X", ., sep = "")

##remove the tags of metabolites
sample_data1 <- 
  smartd_rplc_batch1 %>% 
  dplyr::select(., -(name:Database))

sample_data2 <- smartd_rplc_batch2

##remove NA samples in sample_data1 and sample_data2
##what samples have NAs
which(is.na(sample_data1), arr.ind = TRUE)[,2] %>% 
  unique()

colnames(sample_data1)[c(146, 147, 148)]

###"SFU65" "SFU74" "SFU91" may have no positive data, so remove them from the dataset
sample_data1 <- 
  sample_data1 %>% 
  dplyr::select(., -c(SFU65, SFU74, SFU91))

###metabolite tags
metabolite_tags <- 
  smartd_rplc_batch1 %>% 
  dplyr::select(., name:Database)

###batch 1 batch 2 data integration
mean1 <- apply(sample_data1, 1, mean)
mean2 <- apply(sample_data2, 1, mean)
ref <- mean1 / mean2

sample_data2 <- 
  sample_data2 * ref
####combine batch 1 and batch 2 data
sample_data <- 
  cbind(sample_data1, sample_data2)

###patient and sampel information
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

match(colnames(sample_data), sample_info_191021$Sample_ID)
##P samples have no GA information
####log 10 and scale
sample_data <- 
  log(sample_data, 10) %>% 
  as.data.frame()

sample_data <-
  apply(sample_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

sample_data <- 
  sample_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

colnames(sample_data)[-1] <- 
  smartd_rplc_batch1$name

###add the patient information to the sample_data
sample_data <- 
  inner_join(x = sample_info_191021, 
             y = sample_data, 
             by = "Sample_ID")


###remnove samples with GA == 0
###sample_data_old is the data that contains PP samples
sample_data_old <- 
  sample_data

sample_data <- 
  sample_data %>% 
  dplyr::filter(GA != 0)


#####################################################################################################
#####get discovery and validation data###############################################################
#####################################################################################################
###randomly divided into discovery and valdiation datasets
# set.seed(seed = 1)
sample_data$Patient_ID %>% unique()
##we have 36 peple in total
set.seed(3)
dis_patient_id <- 
  sample_data$Patient_ID %>% 
  unique() %>% 
  sample(size = 18, replace = FALSE)

index_dis <-
  (sample_data$Patient_ID %in% dis_patient_id) %>%
  which()

index_val <- 
  setdiff(1:nrow(sample_data), index_dis)

save(index_dis, file = "index_dis")
save(index_val, file = "index_val")

sample_data_dis <- 
  sample_data[index_dis,]

sample_data_val <- 
  sample_data[index_val,]

sample_data_dis_x <- 
  sample_data_dis %>% 
  dplyr::select(-(Patient_ID:`Birth Control at Discharge?`)) %>% 
  as.matrix()
sum(is.na(sample_data_dis_x))
sample_data_dis_x[,1]


sample_data_val_x <- 
  sample_data_val %>% 
  dplyr::select(-(Patient_ID:`Birth Control at Discharge?`)) %>% 
  as.matrix()
sum(is.na(sample_data_dis_x))
sample_data_val_x[,1]


save(sample_data_dis, file = "sample_data_dis")
save(sample_data_val, file = "sample_data_val")

save(sample_data_dis_x, file = "sample_data_dis_x")
save(sample_data_val_x, file = "sample_data_val_x")

save(metabolite_tags, file = "metabolite_tags")




##information for discovery and validation patients

###BMI
patient_id_dis <- 
  sample_data_dis$Patient_ID %>% 
  unique()

patient_id_val <- 
  sample_data_val$Patient_ID %>% 
  unique()

sample_data_dis$GA
sample_data_val$GA

data.frame(ga = c(sample_data_dis$GA, sample_data_val$GA), 
           class = c(rep("dis", nrow(sample_data_dis)), rep("val", nrow(sample_data_val))),
           stringsAsFactors = FALSE) %>% 
  ggplot(aes(x = class, y = ga, colour = class)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()
  
