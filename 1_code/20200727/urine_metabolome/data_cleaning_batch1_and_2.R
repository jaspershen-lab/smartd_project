#to avoind source
no_exist_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list=ls())
source("R/20200727/tools.R")

####Positive mode
setwd("data_analysis20200108/urine_metabolome/data_cleaning/RPLC/POS/")

sample_info <- readr::read_csv("sample_info.csv")
peak_table <- readr::read_csv("Peak_table_for_cleaning.csv")

colnames(peak_table)
sample_info$sample.name

##remove duplicated samples
name <- grep("_[0-9]{6,15}", colnames(peak_table), value = TRUE)


remove_name <- NULL

for(i in seq_along(name)){
  temp_date <- stringr::str_extract(string = name[i], 
                                    pattern = "_[0-9]{6,15}") %>% 
    stringr::str_replace("_", "")
  
  temp_name <- 
    name[i] %>% 
    stringr::str_replace("_[0-9]{6,15}", "")
  
  temp_peak_table <- 
    peak_table[,c(which(paste(colnames(peak_table), temp_date, sep = "_") == name[i]),
                  which(colnames(peak_table) == name[i]))]
  
  
  remove_name <- c(
  remove_name, 
  apply(temp_peak_table, 2, function(x) sum(is.na(x))) %>% 
    which.min() %>% 
    names() %>% 
    `!=`(colnames(temp_peak_table)) %>% 
    which() %>% 
    `[`(colnames(temp_peak_table), .)
  )
}


peak_table <- 
  peak_table %>% 
  select(-remove_name)


##remove wrong names
colnames(peak_table) <- 
colnames(peak_table) %>% 
  sapply(function(x){
    if(stringr::str_detect(x, '_[0-9]{10,20}')){
      x <- stringr::str_replace(x, '_[0-9]{10,20}', "")
      x
    }else{
     return(x) 
    }
  }) %>% 
  unname()


setdiff(colnames(peak_table), sample_info$sample.name)
setdiff(sample_info$sample.name, colnames(peak_table))

peak_table <- 
  peak_table %>% 
  select(-c(QC2.4, QC2.5))

readr::write_csv(peak_table, "peak_table.csv")



####data cleaning
library(metflow2)
library(tidyverse)

(
  rplc_pos_1 <-
    creatMetflowObject(
      ms1.data = "peak_table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(rplc_pos_1, file = "rplc_pos_1")


###remove peaks
# (
#   rplc_pos_2 <-
#     filterPeak(
#       object = rplc_pos_1,
#       min.fraction.qc = 0.8,##QC at least more than 80%
#       min.fraction = 0.2,##Subject no restraction
#       min.subject.qc.ratio = 0,##don't use blank to remove any peaks
#       dl.qc.r2.cutoff = 0
#     )
# )
# 
# 
###remove some samples (outliers)
# (
#   rplc_pos_3 <-
#     filterSample(object = rplc_pos_2,
#                  min.fraction.peak = 0.5)
# # )
# 
# (
#   plot <- 
#     getMVplot4sample(object = rplc_pos_3)
# )
# 
# 
# 
# 
# ###MV imputateion
# (
#   rplc_pos_4 <- 
#     imputeMV(object = rplc_pos_3, method = 'knn')
# )
# 
# save(rplc_pos_4, file = "rplc_pos_4")
# qc_data <- 
#   getData(object = rplc_pos_batch1_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.rplc.pos4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.position = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
# qc_rsd <- calRSD(object = rplc_pos_4, slot = "QC")
# 
# (plot1 <- 
#     ggplot(data = qc_rsd, aes(index, rsd)) +
#     geom_hline(yintercept = 30, linetype = 2, colour = "red") +
#     geom_point() +
#     scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
#                        breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
#     labs(x = "Peak index") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 15),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 15),
#           legend.text = element_text(size = 12),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15))
# )
# 
# 
# 
# (plot2 <- 
#     ggplot(data = qc_rsd, aes(x = rsd)) +
#     geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
#     scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
#     labs(y = "Peak number") +
#     # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
#     annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
#              alpha = 0.1, fill = "#ED00007F") +
#     theme_bw()+
#     theme(axis.title = element_text(size = 15),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 15),
#           legend.text = element_text(size = 12),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15))
# )
# 
# 
# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization
# 
# (
#   rplc_pos_5 <- 
#     normalizeData(object = rplc_pos_4, 
#                   method = "mean")
# )
# 
# save(rplc_pos_5, file = "rplc_pos_5")
# ###batch intergration
# 
# (
#   rplc_pos_6 <- 
#     metflow2:::integrateData(
#       object = rplc_pos_5, 
#       method = "subject.mean")
# )
# 
# qc_rsd <- calRSD(object = rplc_pos_6, 
#                  slot = "QC")
# 
# plot(qc_rsd$rsd)
# save(rplc_pos_6, file = "rplc_pos_6")
load("rplc_pos_6")

##------------------------------------------------------------------------------
####negative mode
# setwd(r4projects::get_project_wd())
# setwd("data_analysis20200108/data_cleaning/RPLC/NEG/")
# library(tidyverse)
# sample_info <- readr::read_csv("sample_info.csv")
# peak_table <- readr::read_csv("Peak_table_for_cleaning.csv")
# 
# colnames(peak_table)
# sample_info$sample.name
# 
# ##remove duplicated samples
# name <- grep("_[0-9]{4,15}", colnames(peak_table), value = TRUE)
# 
# test <- 
#   peak_table %>% 
#   select(dplyr::starts_with(stringr::str_split(name[4], "_")[[1]][1]))
# 
# 
# # plot(test$QCU28, test$QCU28_171203020739)
# 
# remove_name <- NULL
# for(i in seq_along(name)){
#   temp_name <- name[i] %>% 
#     stringr::str_split("_") %>% 
#     `[[`(1) %>% 
#     `[`(1)
#   
#   temp_peak_table <- 
#     peak_table %>% 
#     select(starts_with(temp_name))
#   
#   
#   remove_name <- c(
#     remove_name, 
#     apply(temp_peak_table, 2, function(x) sum(is.na(x))) %>% 
#       which.min() %>% 
#       names() %>% 
#       `!=`(colnames(temp_peak_table)) %>% 
#       which() %>% 
#       `[`(colnames(temp_peak_table), .)
#   )
# }
# 
# 
# peak_table <- 
#   peak_table %>% 
#   select(-remove_name)
# 
# 
# ##remove wrong names
# colnames(peak_table) <- 
#   colnames(peak_table) %>% 
#   sapply(function(x){
#     if(stringr::str_detect(x, '_[0-9]{10,20}')){
#       x <- stringr::str_replace(x, '_[0-9]{10,20}', "")
#       x
#     }else{
#       return(x) 
#     }
#   }) %>% 
#   unname()
# 
# name1 <- 
#   setdiff(colnames(peak_table)[-c(1:3)], sample_info$sample.name)
# 
# peak_table <-
#   peak_table %>% 
#   select(-name1)
# 
# name2 <- 
#   setdiff(sample_info$sample.name, colnames(peak_table))
# 
# if(length(name2) != 0){
#   sample_info <- 
#     sample_info %>% 
#     filter(!sample.name %in% name2)  
# }
# 
# 
# 
# readr::write_csv(peak_table, "peak_table.csv")
# readr::write_csv(sample_info, "sample_info.csv")
# 
# ####data cleaning
# library(metflow2)
# library(tidyverse)
# 
# (
#   rplc_neg_1 <-
#     creatMetflowObject(
#       ms1.data = "peak_table.csv",
#       sample.information = "sample_info.csv",
#       path = "."
#     )
# )
# 
# save(rplc_neg_1, file = "rplc_neg_1")
# 
# 
# ###remove peaks
# (
#   rplc_neg_2 <-
#     filterPeak(
#       object = rplc_neg_1,
#       min.fraction.qc = 0.8,##QC at least more than 80%
#       min.fraction = 0.2,##Subject no restraction
#       min.subject.qc.ratio = 0,##don't use blank to remove any peaks
#       dl.qc.r2.cutoff = 0
#     )
# )
# 
# 
# ###remove some samples (outliers)
# (
#   rplc_neg_3 <-
#     filterSample(object = rplc_neg_2,
#                  min.fraction.peak = 0.5)
# )
# 
# (
#   plot <- 
#     getMVplot4sample(object = rplc_neg_3)
# )
# 
# 
# 
# 
# ###MV imputateion
# (
#   rplc_neg_4 <- 
#     imputeMV(object = rplc_neg_3, method = 'knn')
# )
# 
# save(rplc_neg_4, file = "rplc_neg_4")
# 
# qc_data <- 
#   getData(object = rplc_neg_batch1_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.rplc.neg4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.negition = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
# qc_rsd <- calRSD(object = rplc_neg_4, slot = "QC")
# 
# (plot1 <- 
#     ggplot(data = qc_rsd, aes(index, rsd)) +
#     geom_hline(yintercept = 30, linetype = 2, colour = "red") +
#     geom_point() +
#     scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
#                        breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
#     labs(x = "Peak index") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 15),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 15),
#           legend.text = element_text(size = 12),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15))
# )
# 
# 
# 
# (plot2 <- 
#     ggplot(data = qc_rsd, aes(x = rsd)) +
#     geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
#     scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
#     labs(y = "Peak number") +
#     # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
#     annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
#              alpha = 0.1, fill = "#ED00007F") +
#     theme_bw()+
#     theme(axis.title = element_text(size = 15),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 15),
#           legend.text = element_text(size = 12),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15))
# )
# 
# 
# # g1 = ggplot2::ggplotGrob(plot1)
# # (plot3 <- plot2 +
# #     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# # )
# 
# # export::graph2ppt(plot3, "test", width = 6, heigh = 4)
# ####data normalization
# 
# (
#   rplc_neg_5 <- 
#     normalizeData(object = rplc_neg_4, 
#                   method = "mean")
# )
# 
# save(rplc_neg_5, file = "rplc_neg_5")
# ###batch intergration
# 
# (
#   rplc_neg_6 <- 
#     metflow2:::integrateData(
#       object = rplc_neg_5, 
#       method = "subject.mean")
# )
# 
# qc_rsd <- calRSD(object = rplc_neg_6, 
#                  slot = "QC")
# 
# plot(qc_rsd$rsd)
# save(rplc_neg_6, file = "rplc_neg_6")
load("rplc_pos_6")




####HILIC
setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())
####Positive mode
setwd("data_analysis20200108/data_cleaning/HILIC/POS/")
sample_info <- readr::read_csv("sample_info.csv")
peak_table <- readr::read_csv("Peak_table_for_cleaning.csv")

colnames(peak_table)
sample_info$sample.name

##remove duplicated samples
name <- grep("_[0-9]{6,15}", colnames(peak_table), value = TRUE)


remove_name <- NULL

for(i in seq_along(name)){
  temp_date <- stringr::str_extract(string = name[i], 
                                    pattern = "_[0-9]{6,15}") %>% 
    stringr::str_replace("_", "")
  
  temp_name <- 
    name[i] %>% 
    stringr::str_replace("_[0-9]{6,15}", "")
  
  temp_peak_table <- 
    peak_table[,c(which(paste(colnames(peak_table), temp_date, sep = "_") == name[i]),
                  which(colnames(peak_table) == name[i]))]
  
  if(ncol(temp_peak_table) == 1){next()}
  
  remove_name <- c(
    remove_name, 
    apply(temp_peak_table, 2, function(x) sum(is.na(x))) %>% 
      which.min() %>% 
      names() %>% 
      `!=`(colnames(temp_peak_table)) %>% 
      which() %>% 
      `[`(colnames(temp_peak_table), .)
  )
}


peak_table <- 
  peak_table %>% 
  select(-remove_name)


##remove wrong names
colnames(peak_table) <- 
  colnames(peak_table) %>% 
  sapply(function(x){ 
    if(stringr::str_detect(x, '_[0-9]{10,20}')){
      x <- stringr::str_replace(x, '_[0-9]{10,20}', "")
      x
    }else{
      return(x) 
    }
  }) %>% 
  unname()

####change some sample IDs
urine_id <- 
  readr::read_csv("E:/project/smartD/patient information/SmartD_all346urine.csv")

urine_id$sample_id_HILIC %>% stringr::str_sort(numeric = TRUE) %>% 
  grep("SFU_B", ., value = TRUE)


##change the HILIC ID to RPLC ID
match(colnames(peak_table), urine_id$sample_id_HILIC)



name <- 
match(colnames(peak_table), urine_id$sample_id_HILIC) %>% 
  `[`(urine_id$sample_id_RPLC,.)

name <- 
  data.frame(colnames(peak_table), name, stringsAsFactors = FALSE)


name <- 
  apply(name, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2])){
      return(x[1])
    }else{
      return(x[2])
    }
  })


name[grep("SFU_B", name)] <- 
  name[grep("SFU_B", name)] %>% 
  gsub(pattern = "SFU_B", replacement = "", x = .) %>% 
  paste("X", ., sep = "")


##X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on
temp_idx <- match(paste("X", 178:198, sep = ""), name)
name[temp_idx] <-
  name[temp_idx] %>% 
  stringr:::str_replace("X", "") %>% 
    as.numeric() %>% 
    `-`(177) %>% 
    paste("P", ., sep = "")

cbind(colnames(peak_table), name)  


colnames(peak_table) <- name


name1 <- 
  setdiff(colnames(peak_table), sample_info$sample.name)
name2 <- setdiff(sample_info$sample.name, colnames(peak_table))

name1 <- 
  name1[-grep("QC|Blank|BLK", name1)] 

name1 <- name1[-c(1:3)]

name2 <- 
  name2[-grep("QC|Blank|blk|BLK", name2)]


sort(name1)
sort(name2)


sample_info <- 
  sample_info %>% 
  filter(group == "Subject") %>% 
  filter(!sample.name %in% c("X107", "X108", "X109", "X110", "X111"))


blank <- grep("Blank|BLK",colnames(peak_table), value = TRUE)
qc <- grep("QC",colnames(peak_table), value = TRUE)


blank <- data.frame(sample.name = blank, injection.order = 0, 
                    class = "Blank", batch = 1, group = "Blank", GA = 0,
                    stringsAsFactors = FALSE)

qc <- data.frame(sample.name = qc, injection.order = 0, 
                    class = "QC", batch = 1, group = "QC", GA = 0,
                    stringsAsFactors = FALSE)


sample_info <- rbind(blank, qc, sample_info)

sample_info <- tibble::as_tibble(sample_info)

####give the right batch inforamtion
batch_info_pos <- readr::read_csv("batch_info_pos.csv")

#rename batch_info_pos ID
batch_info_pos$name <-
  batch_info_pos$name %>% 
  stringr::str_replace("_[0-9]{8,20}", "") 
  
batch_info_pos <- 
  batch_info_pos %>% 
  dplyr::distinct()

##change name
name <- 
  match(batch_info_pos$name, urine_id$sample_id_HILIC) %>% 
  `[`(urine_id$sample_id_RPLC, .)

name <- 
data.frame(name1 = batch_info_pos$name, 
           name2 = name, stringsAsFactors = FALSE) %>% 
  apply(1, function(x){
    x <- as.character(x)
    if(is.na(x[2])){
      return(x[1])
    }else{
      return(x[2])
    }
  })


name[grep("SFU_B", name)] <- 
  name[grep("SFU_B", name)] %>% 
  gsub(pattern = "SFU_B", replacement = "", x = .) %>% 
  paste("X", ., sep = "")

##X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on

temp_idx <- match(paste("X", 178:198, sep = ""), name)
name[temp_idx] <-
  name[temp_idx] %>% 
  stringr:::str_replace("X", "") %>% 
  as.numeric() %>% 
  `-`(177) %>% 
  paste("P", ., sep = "")

cbind(batch_info_pos$name, name)  


batch_info_pos$name <- name

batch <- 
match(sample_info$sample.name, batch_info_pos$name) %>% 
  `[`(batch_info_pos$batch, .) 

sample_info$batch <- batch


readr::write_csv(sample_info, "sample_info.csv")

readr::write_csv(peak_table, "peak_table.csv")


####data cleaning
library(metflow2)
library(tidyverse)

(
  hilic_pos_1 <-
    creatMetflowObject(
      ms1.data = "peak_table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(hilic_pos_1, file = "hilic_pos_1")


###remove peaks
(
  hilic_pos_2 <-
    filterPeak(
      object = hilic_pos_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.2,##Subject no restraction
      min.subject.qc.ratio = 0,##don't use blank to remove any peaks
      dl.qc.r2.cutoff = 0
    )
)


###remove some samples (outliers)
(
  hilic_pos_3 <-
    filterSample(object = hilic_pos_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = hilic_pos_3)
)




###MV imputateion
(
  hilic_pos_4 <- 
    imputeMV(object = hilic_pos_3, method = 'knn')
)

save(hilic_pos_4, file = "hilic_pos_4")

# qc_data <- 
#   getData(object = hilic_pos_batch1_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.hilic.pos4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.position = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- calRSD(object = hilic_pos_4, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)



(plot2 <- 
    ggplot(data = qc_rsd, aes(x = rsd)) +
    geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
    scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
    labs(y = "Peak number") +
    # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
    annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
             alpha = 0.1, fill = "#ED00007F") +
    theme_bw()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization

(
  hilic_pos_5 <- 
    normalizeData(object = hilic_pos_4, 
                  method = "mean")
)

save(hilic_pos_5, file = "hilic_pos_5")
###batch intergration

(
  hilic_pos_6 <- 
    metflow2:::integrateData(
      object = hilic_pos_5, 
      method = "subject.mean")
)

qc_rsd <- calRSD(object = hilic_pos_6, 
                 slot = "QC")

plot(qc_rsd$rsd)
save(hilic_pos_6, file = "hilic_pos_6")





####negative mode
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_cleaning/HILIC/neg/")
rm(list=ls())
sample_info <- readr::read_csv("sample_info.csv")
peak_table <- readr::read_csv("Peak_table_for_cleaning.csv")

colnames(peak_table)
sample_info$sample.name

##remove duplicated samples
name <- grep("_[0-9]{6,15}", colnames(peak_table), value = TRUE)


remove_name <- NULL

for(i in seq_along(name)){
  temp_date <- stringr::str_extract(string = name[i], 
                                    pattern = "_[0-9]{6,15}") %>% 
    stringr::str_replace("_", "")
  
  temp_name <- 
    name[i] %>% 
    stringr::str_replace("_[0-9]{6,15}", "")
  
  temp_peak_table <- 
    peak_table[,c(which(paste(colnames(peak_table), temp_date, sep = "_") == name[i]),
                  which(colnames(peak_table) == name[i]))]
  
  if(ncol(temp_peak_table) == 1){next()}
  
  remove_name <- c(
    remove_name, 
    apply(temp_peak_table, 2, function(x) sum(is.na(x))) %>% 
      which.min() %>% 
      names() %>% 
      `!=`(colnames(temp_peak_table)) %>% 
      which() %>% 
      `[`(colnames(temp_peak_table), .)
  )
}


peak_table <- 
  peak_table %>% 
  select(-remove_name)


##remove wrong names
colnames(peak_table) <- 
  colnames(peak_table) %>% 
  sapply(function(x){ 
    if(stringr::str_detect(x, '_[0-9]{10,20}')){
      x <- stringr::str_replace(x, '_[0-9]{10,20}', "")
      x
    }else{
      return(x) 
    }
  }) %>% 
  unname()

####change some sample IDs
urine_id <- 
  readr::read_csv("E:/project/smartD/patient information/SmartD_all346urine.csv")

urine_id$sample_id_HILIC %>% stringr::str_sort(numeric = TRUE) %>% 
  grep("SFU_B", ., value = TRUE)


##change the HILIC ID to RPLC ID
match(colnames(peak_table), urine_id$sample_id_HILIC)



name <- 
  match(colnames(peak_table), urine_id$sample_id_HILIC) %>% 
  `[`(urine_id$sample_id_RPLC,.)

name <- 
  data.frame(colnames(peak_table), name, stringsAsFactors = FALSE)


name <- 
  apply(name, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2])){
      return(x[1])
    }else{
      return(x[2])
    }
  })


name[grep("SFU_B", name)] <- 
  name[grep("SFU_B", name)] %>% 
  gsub(pattern = "SFU_B", replacement = "", x = .) %>% 
  paste("X", ., sep = "")


##X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on
temp_idx <- match(paste("X", 178:198, sep = ""), name)
name[temp_idx] <-
  name[temp_idx] %>% 
  stringr:::str_replace("X", "") %>% 
  as.numeric() %>% 
  `-`(177) %>% 
  paste("P", ., sep = "")

cbind(colnames(peak_table), name)  


colnames(peak_table) <- name


name1 <- 
  setdiff(colnames(peak_table), sample_info$sample.name)
name2 <- setdiff(sample_info$sample.name, colnames(peak_table))

name1 <- 
  name1[-grep("QC|Blank|BLK", name1)] 

name1 <- name1[-c(1:3)]

name2 <- 
  name2[-grep("QC|Blank|blk|BLK", name2)]


sort(name1)
sort(name2)


sample_info <- 
  sample_info %>% 
  filter(group == "Subject") %>% 
  filter(!sample.name %in% c("X128"))


blank <- grep("Blank|BLK",colnames(peak_table), value = TRUE)
qc <- grep("QC",colnames(peak_table), value = TRUE)


blank <- data.frame(sample.name = blank, injection.order = 0, 
                    class = "Blank", batch = 1, group = "Blank", GA = 0,
                    stringsAsFactors = FALSE)

qc <- data.frame(sample.name = qc, injection.order = 0, 
                 class = "QC", batch = 1, group = "QC", GA = 0,
                 stringsAsFactors = FALSE)


sample_info <- rbind(blank, qc, sample_info)

sample_info <- tibble::as_tibble(sample_info)

####give the right batch inforamtion
batch_info_neg <- readr::read_csv("batch_info_neg.csv")

#rename batch_info_neg ID
batch_info_neg$name <-
  batch_info_neg$name %>% 
  stringr::str_replace("_[0-9]{8,20}", "") 

batch_info_neg <- 
  batch_info_neg %>% 
  dplyr::distinct()

##change name
name <- 
  match(batch_info_neg$name, urine_id$sample_id_HILIC) %>% 
  `[`(urine_id$sample_id_RPLC, .)

name <- 
  data.frame(name1 = batch_info_neg$name, 
             name2 = name, stringsAsFactors = FALSE) %>% 
  apply(1, function(x){
    x <- as.character(x)
    if(is.na(x[2])){
      return(x[1])
    }else{
      return(x[2])
    }
  })


name[grep("SFU_B", name)] <- 
  name[grep("SFU_B", name)] %>% 
  gsub(pattern = "SFU_B", replacement = "", x = .) %>% 
  paste("X", ., sep = "")

##X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on

temp_idx <- match(paste("X", 178:198, sep = ""), name)
name[temp_idx] <-
  name[temp_idx] %>% 
  stringr:::str_replace("X", "") %>% 
  as.numeric() %>% 
  `-`(177) %>% 
  paste("P", ., sep = "")

cbind(batch_info_neg$name, name)  


batch_info_neg$name <- name

batch <- 
  match(sample_info$sample.name, batch_info_neg$name) %>% 
  `[`(batch_info_neg$batch, .) 

sample_info$batch <- batch


readr::write_csv(sample_info, "sample_info.csv")

readr::write_csv(peak_table, "peak_table.csv")


####data cleaning
library(metflow2)
library(tidyverse)

(
  hilic_neg_1 <-
    creatMetflowObject(
      ms1.data = "peak_table.csv",
      sample.information = "sample_info.csv",
      path = "."
    )
)

save(hilic_neg_1, file = "hilic_neg_1")


###remove peaks
(
  hilic_neg_2 <-
    filterPeak(
      object = hilic_neg_1,
      min.fraction.qc = 0.8,##QC at least more than 80%
      min.fraction = 0.2,##Subject no restraction
      min.subject.qc.ratio = 0,##don't use blank to remove any peaks
      dl.qc.r2.cutoff = 0
    )
)


###remove some samples (outliers)
(
  hilic_neg_3 <-
    filterSample(object = hilic_neg_2,
                 min.fraction.peak = 0.5)
)

(
  plot <- 
    getMVplot4sample(object = hilic_neg_3)
)




###MV imputateion
(
  hilic_neg_4 <- 
    imputeMV(object = hilic_neg_3, method = 'knn')
)

save(hilic_neg_4, file = "hilic_neg_4")

# qc_data <- 
#   getData(object = hilic_neg_batch1_4, slot = "QC")
# 
# qc_data <- t(apply(qc_data, 1, function(x){
#   (x - mean(x))/sd(x)
# }))
# 
# qc_data <- tibble::as_tibble(qc_data)
# 
# qc_data <- 
#   gather(data = qc_data, colnames(qc_data), 
#          key = "QC", value = "Intensity")
# 
# plot <- 
#   ggplot(qc_data, aes(x = QC, y = Intensity)) +
#   geom_boxplot(colour = "#00468B7F", outlier.size = 0.5) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(y = "Intensity (z-score)") +
#   coord_flip() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15))
# 
# plot
# export::graph2ppt(plot, "test", width = 4, heigh = 8)


###PCA analysis
# plot <- 
#   pcaAnalysis(object = smart.d.hilic.neg4, scale.method = "auto")
# (plot <- plot +
#     theme(legend.negition = c(0, 1),legend.justification = c(0, 1),
#           legend.background = element_blank())
# )

# export::graph2ppt(plot, "test", width = 6, heigh = 4)


###RSD distribution
qc_rsd <- calRSD(object = hilic_neg_4, slot = "QC")

(plot1 <- 
    ggplot(data = qc_rsd, aes(index, rsd)) +
    geom_hline(yintercept = 30, linetype = 2, colour = "red") +
    geom_point() +
    scale_y_continuous(name = "Relative standard deviation (RSD, %)", limits = c(0,100),
                       breaks = c(0, 25, 30, 50, 75, 100), labels = c(c(0, 25, 30, 50, 75, 100)))+
    labs(x = "Peak index") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)



(plot2 <- 
    ggplot(data = qc_rsd, aes(x = rsd)) +
    geom_histogram(binwidth = 10, colour = "white", fill = "#0099B47F", origin = 0) +
    scale_x_continuous(name = "Relative standard deviation (RSD, %)") +
    labs(y = "Peak number") +
    # geom_vline(xintercept = 30, linetype = 2, colour = "red") +
    annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 380,
             alpha = 0.1, fill = "#ED00007F") +
    theme_bw()+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
)


# g1 = ggplot2::ggplotGrob(plot1)
# (plot3 <- plot2 +
#     annotation_custom(grob = g1, xmin = 100, xmax = 500, ymin = 50, ymax = 350)
# )

# export::graph2ppt(plot3, "test", width = 6, heigh = 4)
####data normalization

(
  hilic_neg_5 <- 
    normalizeData(object = hilic_neg_4, 
                  method = "mean")
)

save(hilic_neg_5, file = "hilic_neg_5")
###batch intergration

(
  hilic_neg_6 <- 
    metflow2:::integrateData(
      object = hilic_neg_5, 
      method = "subject.mean")
)

qc_rsd <- calRSD(object = hilic_neg_6, 
                 slot = "QC")

plot(qc_rsd$rsd)
save(hilic_neg_6, file = "hilic_neg_6")
