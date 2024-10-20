sxtTools::setwd_project()
setwd("troubleshooting/")

setwd("POS/XCMS_processing/")

library(metflow2)
library(xcms)
library(MSnbase)
metflow2::processData(polarity = 'positive',
                      peakwidth = c(5, 30), 
                      threads = 10, 
                      min.fraction = 0.8, 
                      fill.peaks = FALSE,
                      output.tic = TRUE, 
                      output.bpc = TRUE,
                      output.rt.correction.plot = TRUE)


setwd("../../NEG/XCMS_processing/")

library(metflow2)
library(xcms)
library(MSnbase)
metflow2::processData(polarity = 'negative',
                      peakwidth = c(5, 30), 
                      threads = 10, 
                      min.fraction = 0.8, 
                      fill.peaks = FALSE,
                      output.tic = TRUE, 
                      output.bpc = TRUE,
                      output.rt.correction.plot = TRUE)

####the internal standards must be in the raw data
library(metflow2)
sxtTools::setwd_project()
##positive
setwd("troubleshooting/POS/insternal_standards/")

peak_data_pos <- metflow2::extractPeaks(threads = 3, peak.table = "is.xlsx")

metflow2::showPeak(object = peak_data_pos, peak.index = 1, alpha = 0)
metflow2::showPeak(object = peak_data_pos, peak.index = 2, alpha = 0)
metflow2::showPeak(object = peak_data_pos, peak.index = 3, alpha = 0)
metflow2::showPeak(object = peak_data_pos, peak.index = 4, alpha = 0)
metflow2::showPeak(object = peak_data_pos, peak.index = 5, alpha = 0)
metflow2::showPeak(object = peak_data_pos, peak.index = 6, alpha = 0)
metflow2::showPeak(object = peak_data_pos, peak.index = 7, alpha = 0)


setwd("..")
peak_table_pos <- readr::read_csv("Peak.table.csv")

peak_table_pos %>% 
  select()

is_table_pos <- readxl::read_xlsx("insternal_standards//is.xlsx")

data1 <- data.frame(mz = is_table_pos$mz, rt = 1, 
                    stringsAsFactors = FALSE)

data2 <- peak_table_pos[,c("mzmed", "rtmed")]

data2 <- 
  data2 %>% 
  rename(mz = mzmed, rt = rtmed)

match_result_pos <- sxtTools::sxtMTmatch(
                     data1 = as.matrix(data1), 
                     data2 = as.matrix(data2), 
                     mz.tol = 25, 
                     rt.tol = 100000000000,
                     rt.error.type = "abs")


match_result_pos





##negative
setwd("../NEG/")

peak_data_neg <- metflow2::extractPeaks(threads = 3, peak.table = "is.xlsx")

metflow2::showPeak(object = peak_data_neg, peak.index = 1)#no peaks
metflow2::showPeak(object = peak_data_neg, peak.index = 2)#no peaks
metflow2::showPeak(object = peak_data_neg, peak.index = 3, alpha = 0)
metflow2::showPeak(object = peak_data_neg, peak.index = 4, alpha = 0)#no peaks
metflow2::showPeak(object = peak_data_neg, peak.index = 5, alpha = 0)#no peak
metflow2::showPeak(object = peak_data_neg, peak.index = 6, alpha = 0)#no peak
metflow2::showPeak(object = peak_data_neg, peak.index = 7, alpha = 0)#no peak











load("POS/smartd_rplc_pos_batch1_5")
load("POS/smartd_rplc_pos_batch2_5")

load("NEG/smartd_rplc_neg_batch1_5")
load("NEG/smartd_rplc_neg_batch2_5")

library(tidyverse)
##positive mode
peak_table1_pos <- 
  smartd_rplc_pos_batch1_5@ms1.data[[1]] %>% 
  select(name,contains("QC"))

peak_table2_pos <- 
  smartd_rplc_pos_batch2_5@ms1.data[[1]] %>% 
  select(name,contains("QC"))

peak_table1_neg <- 
  smartd_rplc_neg_batch1_5@ms1.data[[1]] %>% 
  select(name,contains("QC"))

peak_table2_neg <- 
  smartd_rplc_neg_batch2_5@ms1.data[[1]] %>% 
  select(name,contains("QC"))


pos1 <- 
  peak_table1_pos %>% 
  mutate(mean.int = apply(peak_table1_pos[,-1], 1, mean)) %>% 
  arrange(desc(mean.int))

pos2 <- 
  peak_table2_pos %>% 
  mutate(mean.int = apply(peak_table2_pos[,-1], 1, mean)) %>% 
  arrange(desc(mean.int))

neg1 <- 
  peak_table1_neg %>% 
  mutate(mean.int = apply(peak_table1_neg[,-1], 1, mean)) %>% 
  arrange(desc(mean.int))

neg2 <- 
  peak_table2_neg %>% 
  mutate(mean.int = apply(peak_table2_neg[,-1], 1, mean)) %>% 
  arrange(desc(mean.int))

pos1_name <- pos1$name[1:100]
pos2_name <- pos2$name[1:100]

neg1_name <- neg1$name[1:100]
neg2_name <- neg2$name[1:100]



pos1 <- 
smartd_rplc_pos_batch1_5@ms1.data[[1]] %>% 
  select(name:rt) %>% 
  filter(name %in% pos1_name)

pos2 <- 
  smartd_rplc_pos_batch2_5@ms1.data[[1]] %>% 
  select(name:rt) %>% 
  filter(name %in% pos2_name)

neg1 <- 
  smartd_rplc_neg_batch1_5@ms1.data[[1]] %>% 
  select(name:rt) %>% 
  filter(name %in% neg1_name)

neg2 <- 
  smartd_rplc_neg_batch2_5@ms1.data[[1]] %>% 
  select(name:rt) %>% 
  filter(name %in% neg2_name)


result <- 
sxtTools::sxtMTmatch(data1 = as.matrix(pos1[,-1]), 
                     data2 = as.matrix(pos2[,-1]), 
                     mz.tol = 25, rt.tol = 1000000000000000000, 
                     rt.error.type = "abs")


result <- 
  sxtTools::sxtMTmatch(data1 = as.matrix(pos1[,-1]), 
                       data2 = as.matrix(pos2[,-1]), 
                       mz.tol = 25, rt.tol = 120, 
                       rt.error.type = "abs")



















