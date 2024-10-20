library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
rm(list = ls())
sxtTools::setwd_project()
setwd("data_analysis20200108/data_preparation_for_analysis/")
load("peak_table")

sxtTools::setwd_project()
load("data_analysis20200108/data_cleaning/RPLC/POS/rplc_pos_6")
load("data_analysis20200108/data_cleaning/RPLC/NEG/rplc_neg_6")

##first set the work directory to project folder
sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/features/")

sample_info_pos <- 
  rplc_pos_6@sample.info

sample_info_neg <- 
  rplc_neg_6@sample.info


subject_data <- 
  peak_table %>% 
  select(-contains("QC")) %>% 
  select(-contains("name")) %>% 
  select(-contains("mz")) %>% 
  select(-contains("rt"))


qc_data <- 
  peak_table %>% 
  select(contains("QC"))

peak_tags <- 
  peak_table[,c("name", "mz", "rt")]
  
rownames(subject_data) <- 
  rownames(qc_data) <-
  peak_tags$name

sample_info <- 
  sample_info_pos %>% 
  filter(sample.name %in% colnames(subject_data))


sum(sort(sample_info$sample.name) == sort(colnames(subject_data)))

##patient and sampel information
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

match(colnames(subject_data), sample_info_191021$Sample_ID)
##P samples have no GA information
####log 10 and scale
subject_data <- 
  log(subject_data, 10) %>% 
  as.data.frame()

subject_data <-
  apply(subject_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

subject_data <- 
  subject_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")


qc_data <- 
  log(qc_data, 10) %>% 
  as.data.frame()

qc_data <-
  apply(qc_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

qc_data <- 
  qc_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

###add the patient information to the sample_data
subject_data <- 
  inner_join(x = sample_info_191021, 
             y = subject_data, 
             by = "Sample_ID")

save(qc_data, file = "qc_data")


###remnove samples with GA == 0
###sample_data_old is the data that contains PP samples
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
subject_data$Patient_ID %>% unique()
##we have 36 peple in total

index_dis <- which(subject_data$Sample_ID %in% sample_info$sample.name[sample_info$batch == 1])
index_val <- which(subject_data$Sample_ID %in% sample_info$sample.name[sample_info$batch == 2])

save(index_dis, file = "index_dis")
save(index_val, file = "index_val")

sample_data_dis <- 
  subject_data[index_dis,]

sample_data_val <- 
  subject_data[index_val,]

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


  

#####for metabolite identification result
##load data
library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
rm(list=ls())
sxtTools::setwd_project()
setwd("data_analysis20200108/data_preparation_for_analysis/")
load("metabolite_table")
load("metabolite_tags")
sxtTools::setwd_project()
load("data_analysis20200108/data_cleaning/RPLC/POS/rplc_pos_6")
load("data_analysis20200108/data_cleaning/RPLC/NEG/rplc_neg_6")

##first set the work directory to project folder
sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/")

sample_info_pos <- 
  rplc_pos_6@sample.info

sample_info_neg <- 
  rplc_neg_6@sample.info


subject_data <- 
  metabolite_table %>% 
  select(-contains("QC")) %>% 
  select(-contains("name")) %>% 
  select(-contains("mz")) %>% 
  select(-contains("rt"))


qc_data <- 
  metabolite_table %>% 
  select(contains("QC"))

rownames(subject_data) <- 
  rownames(qc_data) <-
  metabolite_tags$name

sample_info <- 
  sample_info_pos %>% 
  filter(sample.name %in% colnames(subject_data))


sum(sort(sample_info$sample.name) == sort(colnames(subject_data)))

##patient and sampel information
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

match(colnames(subject_data), sample_info_191021$Sample_ID)
##P samples have no GA information
####log 10 and scale
subject_data <- 
  log(subject_data, 10) %>% 
  as.data.frame()

subject_data <-
  apply(subject_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

subject_data <- 
  subject_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

sum(is.na(subject_data))

qc_data <- 
  log(qc_data, 10) %>% 
  as.data.frame()

qc_data <-
  apply(qc_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

qc_data <- 
  qc_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

###add the patient information to the sample_data
subject_data <- 
  inner_join(x = sample_info_191021, 
             y = subject_data, 
             by = "Sample_ID")

save(qc_data, file = "qc_data")

###remnove samples with GA == 0
###sample_data_old is the data that contains PP samples
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
subject_data$Patient_ID %>% unique()
##we have 36 peple in total

index_dis <- which(subject_data$Sample_ID %in% sample_info$sample.name[sample_info$batch == 1])
index_val <- which(subject_data$Sample_ID %in% sample_info$sample.name[sample_info$batch == 2])

save(index_dis, file = "index_dis")
save(index_val, file = "index_val")

sample_data_dis <- 
  subject_data[index_dis,]

sample_data_val <- 
  subject_data[index_val,]

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




#####for cytokine data
##load data
library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
rm(list=ls())
sxtTools::setwd_project()
setwd("data_analysis20200108/data_preparation_for_analysis/cytokine/")
load("cytokine_pheno")
load("cytokine_table")
load("cytokine_tags")

##first set the work directory to project folder
sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/cytokine/")

##patient and sampel information
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

match(colnames(cytokine_table), sample_info_191021$Sample_ID)
##P samples have no GA information
####log 10 and scale
cytokine_table <- 
  log(cytokine_table, 10) %>% 
  as.data.frame()

cytokine_table <-
  apply(cytokine_table, 1, function(x){
    (x - mean(x))/sd(x)
  })

cytokine_table <- 
  cytokine_table %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

sum(is.na(cytokine_table))

colnames(cytokine_table)[-1] <- cytokine_tags$name

###add the patient information to the sample_data
cytokine_table <- 
  inner_join(x = sample_info_191021, 
             y = cytokine_table, 
             by = "Sample_ID")

###remnove samples with GA == 0
###sample_data_old is the data that contains PP samples
cytokine_table_old <- 
  cytokine_table

cytokine_table <- 
  cytokine_table %>% 
  dplyr::filter(GA != 0)

cytokine_pheno <- 
  cytokine_pheno %>% 
  filter(sample.id %in% cytokine_table$Sample_ID) %>% 
  arrange(sample.id)

cytokine_table <- 
  cytokine_table %>% 
  arrange(Sample_ID)

cytokine_pheno <-
  cytokine_table %>% 
  select(-one_of(cytokine_tags$name)) %>% 
    cbind(cytokine_pheno) %>% 
  tibble::as_tibble()
  
cytokine_table <- 
  cytokine_table %>% 
  select(one_of(cytokine_tags$name))


save(cytokine_pheno, file = "cytokine_pheno")
save(cytokine_table, file = "cytokine_table")
save(cytokine_tags, file = "cytokine_tags")

#####################################################################################################
#####get discovery and validation data###############################################################
#####################################################################################################
###randomly divided into discovery and valdiation datasets
# set.seed(seed = 1)
cytokine_pheno$Patient_ID %>% unique()
##we have 36 peple in total
#use the same dis and val wiht metabolomics dataset

load("../metabolites/sample_data_dis")
load("../metabolites/sample_data_val")

dis_id <- sample_data_dis$Sample_ID
val_id <- sample_data_val$Sample_ID

idx_dis <- which(cytokine_pheno$Sample_ID %in% dis_id)
idx_val <- which(cytokine_pheno$Sample_ID %in% val_id)


cytokine_pheno_dis <- 
  cytokine_pheno[idx_dis,]
  
cytokine_table_dis <-
  cytokine_table[idx_dis,]


cytokine_pheno_val <- 
  cytokine_pheno[idx_val,]

cytokine_table_val <-
  cytokine_table[idx_val,]


save(cytokine_pheno_dis, file = "cytokine_pheno_dis")
save(cytokine_table_dis, file = "cytokine_table_dis")

save(cytokine_pheno_val, file = "cytokine_pheno_val")
save(cytokine_table_val, file = "cytokine_table_val")

save(cytokine_tags, file = "cytokine_tags")



#####for metabolite and cytokine
##load data
library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
rm(list = ls())
sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolite_cytokine/")

load("../metabolites/sample_data_dis")
load("../metabolites/sample_data_val")
load("../metabolites/metabolite_tags")

load("../cytokine/cytokine_table_dis")
load("../cytokine/cytokine_table_val")

load("../cytokine/cytokine_pheno_dis")
load("../cytokine/cytokine_pheno_val")
load("../cytokine/cytokine_tags")

metabolite_table_dis <- 
sample_data_dis

metabolite_table_val <- 
  sample_data_val

colnames(metabolite_table_dis)


metabolite_pheno_dis <- 
  metabolite_table_dis %>% 
  select(c(Patient_ID:`Birth Control at Discharge?`))


metabolite_table_dis <- 
  metabolite_table_dis %>% 
  select(-c(Patient_ID:`Birth Control at Discharge?`))


metabolite_pheno_val <- 
  metabolite_table_val %>% 
  select(c(Patient_ID:`Birth Control at Discharge?`))


metabolite_table_val <- 
  metabolite_table_val %>% 
  select(-c(Patient_ID:`Birth Control at Discharge?`))


intersect_name_dis <- 
  intersect(metabolite_pheno_dis$Sample_ID, cytokine_pheno_dis$Sample_ID)


intersect_name_val <- 
  intersect(metabolite_pheno_val$Sample_ID, cytokine_pheno_val$Sample_ID)


idx_dis_met <- match(intersect_name_dis, metabolite_pheno_dis$Sample_ID)
idx_dis_cyto <- match(intersect_name_dis, cytokine_pheno_dis$Sample_ID)

idx_val_met <- match(intersect_name_val, metabolite_pheno_val$Sample_ID)
idx_val_cyto <- match(intersect_name_val, cytokine_pheno_val$Sample_ID)

met_cyto_table_dis <-
  cbind(
    metabolite_table_dis[idx_dis_met,],
    cytokine_table_dis[idx_dis_cyto,]
  )
  
met_cyto_pheno_dis <- 
  cbind(
    metabolite_pheno_dis[idx_dis_met,],
    cytokine_pheno_dis[idx_dis_cyto,]  
  )[,union(colnames(metabolite_pheno_dis), colnames(cytokine_pheno_dis))]


met_cyto_table_val <-
  cbind(
    metabolite_table_val[idx_val_met,],
    cytokine_table_val[idx_val_cyto,]
  )

met_cyto_pheno_val <- 
  cbind(
    metabolite_pheno_val[idx_val_met,],
    cytokine_pheno_val[idx_val_cyto,]  
  )[,union(colnames(metabolite_pheno_val), colnames(cytokine_pheno_val))]

##tags
met_cyto_tags <- 
  metabolite_tags %>% 
  dplyr::full_join(cytokine_tags, by = "name")

  
save(met_cyto_table_dis, file = "met_cyto_table_dis")
save(met_cyto_pheno_dis, file = "met_cyto_pheno_dis")

save(met_cyto_table_val, file = "met_cyto_table_val")
save(met_cyto_pheno_val, file = "met_cyto_pheno_val")

save(met_cyto_tags, file = "met_cyto_tags")
