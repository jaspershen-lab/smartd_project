###RPLC data
###align batch 1 and batch 2 data
##positive mode
setwd(r4projects::get_project_wd())
setwd("data_analysis20191125/batch1_2_alignment/RPLC/POS/")
rm(list = ls())
load("smartd_rplc_pos_batch1_5")
load("smartd_rplc_pos_batch2_5")

library(tidyverse)

###batch 1 and batch 2 alignment
batch1_pos_info <- smartd_rplc_pos_batch1_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

batch2_pos_info <- smartd_rplc_pos_batch2_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_pos_info %>% 
  dplyr::select(., mz:rt)

data2 <- 
  batch2_pos_info %>% 
  dplyr::select(., mz:rt)


dim(data1)
dim(data2)


match_result <-
  masstools::sxtMTmatch(data1 = as.matrix(data1), 
                     data2 = as.matrix(data2), 
                     mz.tol = 25, 
                     rt.tol = 120, 
                     rt.error.type = "abs")

match_result <- as.data.frame(match_result)

###remove the 1:multiple
match_result <- 
  match_result %>% 
  as_tibble() %>% 
  group_by(Index1) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup() %>% 
  group_by(Index2) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup()

unique(match_result[,1])
unique(match_result[,2])

name1 <- batch1_pos_info$name[match_result$Index1]
name2 <- batch2_pos_info$name[match_result$Index2]

rplc_pos_batch1_batch2_matched_info <- 
  data.frame(batch1 = name1, batch2 = name2, 
             stringsAsFactors = FALSE)


write.csv(rplc_pos_batch1_batch2_matched_info, 
          'rplc_pos_batch1_batch2_matched_info.csv')

##negative mode
setwd(r4projects::get_project_wd())
setwd("data_analysis20191125/batch1_2_alignment/RPLC/NEG/")
load("smartd_rplc_neg_batch1_5")
load("smartd_rplc_neg_batch2_5")

library(tidyverse)

batch1_neg_info <- smartd_rplc_neg_batch1_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

batch2_neg_info <- smartd_rplc_neg_batch2_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_neg_info %>% 
  dplyr::select(., mz:rt)

data2 <- 
  batch2_neg_info %>% 
  dplyr::select(., mz:rt)

dim(data1)
dim(data2)

match_result <-
  masstools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25, rt.tol = 120, 
                       rt.error.type = "abs")

unique(match_result[,1])

dim(match_result)

match_result <- as.data.frame(match_result)

range(data1$rt)
range(data2$rt)

###remove the 1:multiple
match_result <- 
  match_result %>% 
  as_tibble() %>% 
  group_by(Index1) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup() %>% 
  group_by(Index2) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup()

name1 <- batch1_neg_info$name[match_result$Index1]
name2 <- batch2_neg_info$name[match_result$Index2]

rplc_neg_batch1_batch2_matched_info <- 
  data.frame(batch1 = name1, batch2 = name2, 
             stringsAsFactors = FALSE)


write.csv(rplc_neg_batch1_batch2_matched_info, 
          'rplc_neg_batch1_batch2_matched_info.csv')



## aligned tables for batch 1 and batch 2
setwd(r4projects::get_project_wd())
setwd("data_analysis20191125/batch1_2_alignment/RPLC/POS_NEG/")
load("../POS/smartd_rplc_pos_batch1_5")
load("../POS/smartd_rplc_pos_batch2_5")

load("../NEG/smartd_rplc_neg_batch1_5")
load("../NEG/smartd_rplc_neg_batch2_5")


match_info_pos <- readr::read_csv("../POS/rplc_pos_batch1_batch2_matched_info.csv")
match_info_neg <- readr::read_csv("../NEG/rplc_neg_batch1_batch2_matched_info.csv")

match_info_pos$batch1 <- paste(match_info_pos$batch1, "POS", sep = "_")
match_info_neg$batch1 <- paste(match_info_neg$batch1, "NEG", sep = "_")

match_info_pos$batch2 <- paste(match_info_pos$batch2, "POS", sep = "_")
match_info_neg$batch2 <- paste(match_info_neg$batch2, "NEG", sep = "_")

match_info <- rbind(match_info_pos, match_info_neg)

peak_table_batch1_pos <- smartd_rplc_pos_batch1_5@ms1.data[[1]]
peak_table_batch1_neg <- smartd_rplc_neg_batch1_5@ms1.data[[1]]

peak_table_batch2_pos <- smartd_rplc_pos_batch2_5@ms1.data[[1]]
peak_table_batch2_neg <- smartd_rplc_neg_batch2_5@ms1.data[[1]]

dim(peak_table_batch1_pos)
dim(peak_table_batch1_neg)

dim(peak_table_batch2_pos)
dim(peak_table_batch2_neg)

peak_table_batch1_pos$name <- 
  paste(peak_table_batch1_pos$name, "POS", sep = "_")
peak_table_batch1_neg$name <- 
  paste(peak_table_batch1_neg$name, "NEG", sep = "_")


peak_table_batch2_pos$name <- 
  paste(peak_table_batch2_pos$name, "POS", sep = "_")
peak_table_batch2_neg$name <- 
  paste(peak_table_batch2_neg$name, "NEG", sep = "_")

setdiff(colnames(peak_table_batch1_pos), colnames(peak_table_batch1_neg))
setdiff(colnames(peak_table_batch1_neg), colnames(peak_table_batch1_pos))

intersect_sample <- 
  intersect(colnames(peak_table_batch1_pos), colnames(peak_table_batch1_neg))


intersect_sample <-
  intersect_sample[-grep("BLK", intersect_sample)]


peak_table_batch1_pos <- 
  peak_table_batch1_pos[,intersect_sample]

peak_table_batch1_neg <- 
  peak_table_batch1_neg[,intersect_sample]


# peak_table_batch1 <- dplyr::full_join(peak_table_batch1_pos, 
#                                           peak_table_batch1_neg, 
#                                           by = intersect(colnames(peak_table_batch1_pos), colnames(peak_table_batch1_neg)))


peak_table_batch1 <- rbind(peak_table_batch1_pos, 
                           peak_table_batch1_neg)


sum(is.na(peak_table_batch1[,intersect_sample]))
sum(is.na(peak_table_batch1))


setdiff(colnames(peak_table_batch2_pos), colnames(peak_table_batch2_neg))
setdiff(colnames(peak_table_batch2_neg), colnames(peak_table_batch2_pos))

intersect_sample <- 
  intersect(colnames(peak_table_batch2_pos), colnames(peak_table_batch2_neg))

intersect_sample <-
  intersect_sample[-grep("blk|DL", intersect_sample)]



# peak_table_batch2 <- dplyr::full_join(peak_table_batch2_pos, 
#                                           peak_table_batch2_neg, 
#                                           by = intersect(colnames(peak_table_batch2_pos), colnames(peak_table_batch2_neg)))


peak_table_batch2_pos <- 
  peak_table_batch2_pos[,intersect_sample]

peak_table_batch2_neg <- 
  peak_table_batch2_neg[,intersect_sample]


peak_table_batch2 <- rbind(peak_table_batch2_pos, 
                           peak_table_batch2_neg)

sum(is.na(peak_table_batch2[,intersect_sample]))
sum(is.na(peak_table_batch2))

####only remain the peaks whcih are matched in batch 1 and batch 2
peak_table_batch1 <-
  match(match_info$batch1, peak_table_batch1$name) %>% 
    `[`(peak_table_batch1,.,)

peak_table_batch2 <-
  match(match_info$batch2, peak_table_batch2$name) %>% 
  `[`(peak_table_batch2,.,)
    
write.csv(peak_table_batch1, "peak_table_batch1.csv", row.names = FALSE)
write.csv(peak_table_batch2, "peak_table_batch2.csv", row.names = FALSE)



peak_table_batch2 <- 
  peak_table_batch2 %>% 
  rename(name2 = name, mz2 = mz, rt2 = rt)

peak_table <- 
  cbind(peak_table_batch1, peak_table_batch2) %>% 
  select(name:rt, name2:rt2, everything())


write.csv(peak_table, "peak_table.csv", row.names = FALSE)



###get the table which only have annotation
setwd(r4projects::get_project_wd())
peak_table <- readr::read_csv("data_analysis20191125/batch1_2_alignment/RPLC/POS_NEG/peak_table.csv")
setwd("data_analysis20191125/annotation/RPLC/")
identification.pos <- readr::read_csv("POS/identification.table.new.csv")
identification.neg <- readr::read_csv("NEG/identification.table.new.csv")

identification.pos$name <- 
  paste(identification.pos$name, "POS", sep = "_")

identification.neg$name <- 
  paste(identification.neg$name, "NEG", sep = "_")


identification <- rbind(identification.pos, identification.neg)


metabolite_table <- 
  peak_table %>% 
  dplyr::inner_join(identification, by = c("name", "mz", "rt")) %>% 
  select(name:rt2,MS2.spectrum.name:Database,everything())

write.csv(metabolite_table, "metabolite_table.csv", row.names = FALSE)
  
####output MS2 spectra of metabolite
##positive mode
load("POS/NCE25/result.pRPLC.nce25")
load("POS/NCE50/result.pRPLC.nce50")

load("NEG/NCE25/result.nRPLC.nce25")
load("NEG/NCE50/result.nRPLC.nce50")

peak_name <- 
metabolite_table %>% 
  filter(Level != 3) %>% 
  select(name, Database)

peak_name_pos <- 
  peak_name %>% 
  filter(stringr::str_detect(name, "POS"))

peak_name_neg <- 
  peak_name %>% 
  filter(stringr::str_detect(name, "NEG"))

load("POS/NCE25/hmdbDatabase0.0.1")
load("POS/NCE25/massbankDatabase0.0.1")
load("POS/NCE25/metlinDatabase0.0.1")
load("POS/NCE25/monaDatabase0.0.1")
load("POS/NCE25/msDatabase_rplc0.0.1")
load("POS/NCE25/nistDatabase0.0.1")
load("POS/NCE25/orbitrapDatabase0.0.1")


library(metID)

name <- 
peak_name_pos %>% 
  filter(stringr::str_detect(name, "POS"), 
         stringr::str_detect(Database, "orbitrapDatabase0.0.1")) %>% 
  pull(name) %>% 
  stringr::str_replace("_POS", "")

name

# metID::ms2plot(object = result.pRPLC.nce25[["orbitrapDatabase0.0.1"]],
#                database = orbitrapDatabase0.0.1, 
#                which.peak = name, 
#                path = "POS", one.folder = FALSE)


metID::ms2plot(object = result.pRPLC.nce50[["orbitrapDatabase0.0.1"]],
               database = orbitrapDatabase0.0.1, 
               which.peak = name, 
               path = "POS", one.folder = FALSE)


#------------------------------------------------------------------------------
###HILIC data
###align batch 1 and batch 2 data
##positive mode
setwd(r4projects::get_project_wd())
setwd("smartD_batch1_2/HILIC/POS/")
load("smartd_hilic_pos_batch1_5")
load("smartd_hilic_pos_batch2_5")

library(tidyverse)

batch1_pos_info <- smartd_hilic_pos_batch1_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

batch2_pos_info <- smartd_hilic_pos_batch2_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_pos_info %>% 
  dplyr::select(., mz:rt)

data2 <- 
  batch2_pos_info %>% 
  dplyr::select(., mz:rt)

match_result <-
  masstools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25, rt.tol = 120, 
                       rt.error.type = "abs")

###remove the 1:multiple
match_result <- 
  match_result %>% 
  as_tibble() %>% 
  group_by(Index1) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup()


unique(match_result[,1])

name1 <- batch1_pos_info$name[match_result$Index1]
name2 <- batch2_pos_info$name[match_result$Index2]

hilic_pos_batch1_batch2_matched_info <- 
  data.frame(batch1 = name1, batch2 = name2, 
             stringsAsFactors = FALSE)


write.csv(hilic_pos_batch1_batch2_matched_info, 
          'hilic_pos_batch1_batch2_matched_info.csv')

##negative mode
setwd("smartD_batch1_2/HILIC/NEG/")
load("smartd_hilic_neg_batch1_5")
load("smartd_hilic_neg_batch2_5")

library(tidyverse)

batch1_neg_info <- smartd_hilic_neg_batch1_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

batch2_neg_info <- smartd_hilic_neg_batch2_5@ms1.data[[1]] %>% 
  dplyr::select(., name, mz, rt, contains("QC"))

data1 <- 
  batch1_neg_info %>% 
  dplyr::select(., mz:rt)

data2 <- 
  batch2_neg_info %>% 
  dplyr::select(., mz:rt)

match_result <-
  masstools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25, rt.tol = 120, 
                       rt.error.type = "abs")

unique(match_result[,1])

###remove the 1:multiple
match_result <- 
  match_result %>% 
  as_tibble() %>% 
  group_by(Index1) %>% 
  filter(`rt error` == min(`rt error`)) %>% 
  ungroup()

name1 <- batch1_neg_info$name[match_result$Index1]
name2 <- batch2_neg_info$name[match_result$Index2]

hilic_neg_batch1_batch2_matched_info <- 
  data.frame(batch1 = name1, batch2 = name2, 
             stringsAsFactors = FALSE)


write.csv(hilic_neg_batch1_batch2_matched_info, 
          'hilic_neg_batch1_batch2_matched_info.csv')


####how many metabolites can be found in batch2
setwd("smartD_batch1_2/HILIC/POS_NEG/")
smartd_hilic_batch1 <- readr::read_csv("smartd_hilic_batch1.csv")

name_pos <- paste(hilic_pos_batch1_batch2_matched_info$batch1, "POS", sep = "_")
name_neg <- paste(hilic_neg_batch1_batch2_matched_info$batch1, "NEG", sep = "_")

sum(smartd_hilic_batch1$name %in% c(name_pos, name_neg))



hilic_pos_batch1_batch2_matched_info$batch1 <- 
  paste(hilic_pos_batch1_batch2_matched_info$batch1, "POS", sep = "_")

hilic_pos_batch1_batch2_matched_info$batch2 <- 
  paste(hilic_pos_batch1_batch2_matched_info$batch2, "POS", sep = "_")


hilic_neg_batch1_batch2_matched_info$batch1 <- 
  paste(hilic_neg_batch1_batch2_matched_info$batch1, "NEG", sep = "_")

hilic_neg_batch1_batch2_matched_info$batch2 <- 
  paste(hilic_neg_batch1_batch2_matched_info$batch2, "NEG", sep = "_")

write.csv(hilic_pos_batch1_batch2_matched_info, 
          'hilic_pos_batch1_batch2_matched_info.csv')

write.csv(hilic_neg_batch1_batch2_matched_info, 
          'hilic_neg_batch1_batch2_matched_info.csv')


#######get the identification table for batch 2
setwd("smartD_batch1_2/HILIC/POS_NEG/")
load("../POS/smartd_hilic_pos_batch2_5")
load("../NEG/smartd_hilic_neg_batch2_5")

ms1_data_hilic_batch2_pos <-
  smartd_hilic_pos_batch2_5@ms1.data[[1]]

ms1_data_hilic_batch2_neg <-
  smartd_hilic_neg_batch2_5@ms1.data[[1]]

ms1_data_hilic_batch2_pos$name <- 
  paste(ms1_data_hilic_batch2_pos$name, "POS", sep = "_")

ms1_data_hilic_batch2_neg$name <- 
  paste(ms1_data_hilic_batch2_neg$name, "NEG", sep = "_")


ms1_data_hilic_batch2 <- 
  dplyr::full_join(x = ms1_data_hilic_batch2_pos, 
                   y = ms1_data_hilic_batch2_neg, 
                   by = intersect(colnames(ms1_data_hilic_batch2_pos), colnames(ms1_data_hilic_batch2_neg)))


write.csv(ms1_data_hilic_batch2, "ms1_data_hilic_batch2.csv", 
          row.names = FALSE)


hilic_pos_batch1_batch2_matched_info
hilic_neg_batch1_batch2_matched_info

dim(hilic_pos_batch1_batch2_matched_info)
dim(hilic_neg_batch1_batch2_matched_info)

ms1_data_hilic_batch2 <- 
  ms1_data_hilic_batch2 %>% 
  filter(name %in% c(hilic_pos_batch1_batch2_matched_info$batch2, hilic_neg_batch1_batch2_matched_info$batch2))

write.csv(ms1_data_hilic_batch2, "ms1_data_hilic_batch2.csv", 
          row.names = FALSE)





















