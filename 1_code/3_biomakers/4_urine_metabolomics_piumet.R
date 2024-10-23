##avoid source
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
source("1_code/100_tools.R")

##load data
load(
  "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/peaks/urine_metabolomics_data.rda"
)

dir.create("3_data_analysis/3_biomakers/4_urine_metabolomics_piumet",
           recursive = TRUE)
setwd("3_data_analysis/3_biomakers/4_urine_metabolomics_piumet")

##remove qc and blank samples
urine_metabolomics_data2 <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject")

expression_data <-
  extract_expression_data(urine_metabolomics_data2)

sample_info <-
  extract_sample_info(urine_metabolomics_data2)

variable_info <-
  extract_variable_info(urine_metabolomics_data2)

dim(expression_data)

###find marker which are change according to pregnancy

#log transformation
subject_data <-
  log(expression_data + 1, 10)

sample_info$GA[is.na(sample_info$GA)] <- 45

###remove lm_fdr > 0.05
text_colour <-
  c(colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13), "red")

colnames(subject_data) == sample_info$sample_id

###combine different samples in one ga range together
subject_data <-
  apply(subject_data, 1, function(x) {
    (x) / sd(x)
  })

load("../2_urine_metabolomics_makers_compared2baseline/subject_data_mean")
load("../2_urine_metabolomics_makers_compared2baseline/subject_data_sd")
load("../2_urine_metabolomics_makers_compared2baseline/subject_data_sem")
load("../2_urine_metabolomics_makers_compared2baseline/fc_p_value")

###PIUMet analysis
####output the table for PIUMet analysis
# purrr::walk2(
#   .x = fc_p_value,
#   .y = names(fc_p_value),
#   .f = function(x, y) {
#     idx <- which(x$p_value < 0.05)
#     if(length(idx) == 0){
#       return()
#     }
#     x <-
#       cbind(variable_info, x)[idx, ] %>%
#       as.data.frame() %>%
#       dplyr::mutate(polarity =
#                       case_when(
#                         stringr::str_detect(variable_id, "POS") ~ "positive",
#                         stringr::str_detect(variable_id, "NEG") ~ "negative"
#                       )) %>%
#       dplyr::select(mz, polarity, p_value) %>%
#       dplyr::mutate(p_value = -log(p_value, 10))
#     write.csv(x, file.path(file.path(y, "marker.csv")))
#
#     write.table(
#       x,
#       file.path(file.path(y, "marker.txt")),
#       sep = "\t",
#       quote = FALSE,
#       row.names = FALSE,
#       col.names = FALSE
#     )
#
#   }
# )

###20_22
# readPIUMet(path = "(20,22]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(20,22]`)
#
###22_24
# readPIUMet(path = "(22,24]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(22,24]`)
#
# ###28_30
# readPIUMet(path = "(24,26]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(24,26]`)
#
# ###26_28
# readPIUMet(path = "(26,28]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(26,28]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
# ###28_30
# readPIUMet(path = "(28,30]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(28,30]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
#
# ###30_32
# readPIUMet(path = "(30,32]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(30,32]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
# ###32_34
# readPIUMet(path = "(32,34]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(32,34]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
# ###34_36
# readPIUMet(path = "(34,36]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(34,36]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
#
# ###36_38
# readPIUMet(path = "(36,38]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(36,38]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
# ###38_42
# readPIUMet(path = "(38,42]/piumet_output",
#            variable_info = variable_info,
#            fc_p_table = fc_p_value$`(38,42]`,
#            text = FALSE,
#            layout = "kk",
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
#
# # ###PP
# readPIUMet(
#   path = "PP/piumet_output",
#   variable_info = variable_info,
#   fc_p_table = fc_p_value$`PP`,
#   text = FALSE,
#   layout = "kk",
#   size_range = c(1, 3),
#   width_range = c(0.2, 0.5)
# )


##pathway enrichment for each data point
# #22_24
# load("(22,24]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
# lapply(hmdb_id, function(x){
#   metflow2::transID(query = x,
#                     from = "Human Metabolome Database",
#                     to = "KEGG", top = 1)
# }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
# apply(kegg_id, 1, function(x){
#   x <- as.character(x)
#   if(is.na(x[2]) & is.na(x[3])){
#     return(NA)
#   }
#
#   if(!is.na(x[2])){
#     return(x[2])
#   }
#
#   if(!is.na(x[3])){
#     return(x[3])
#   }
#
#
# })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
#
# node_data <-
# node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(22,24]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# load("hsa_pathway")
#
# #remove null
# remove_idx <-
# lapply(hsa_pathway, is.null) %>%
#   unlist() %>%
#   which()
#
# if(length(remove_idx) > 0){
#   hsa_pathway <- hsa_pathway[-remove_idx]
# }
#
# save(hsa_pathway, file = "hsa_pathway")
#
# path_result <-
# enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result <-
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
#
# save(path_result, file = "(22,24]/path_result")
#
#
# ##28_30
# load("(24,26]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(24,26]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
#   path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
#
# save(path_result, file = "(24,26]/path_result")
#
#
#
# ##26_28
# load("(26,28]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(26,28]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(26,28]/path_result")
#
#
#
# ##28_30
# load("(28,30]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# # load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(28,30]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(28,30]/path_result")
#
#
#
#
#
#
# ##30_32
# load("(30,32]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
# kegg_id
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(30,32]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# kegg_id
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(30,32]/path_result")
#
#
# ##32_34
# load("(32,34]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
# kegg_id
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(32,34]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# kegg_id
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(32,34]/path_result")
#
#
#
#
#
# ##34_36
# load("(34,36]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
# kegg_id
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(34,36]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# kegg_id
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(34,36]/path_result")
#
#
#
# ##36_38
# load("(36,38]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
# kegg_id
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(36,38]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# kegg_id
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(36,38]/path_result")
#
#
#
#
#
# ##38_42
# load("(38,42]/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
# kegg_id
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "(38,42]/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# kegg_id
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "(38,42]/path_result")
#
#
#
#
#
#
#
# ##PP
# load("PP/piumet_output/Result/node_data")
# hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])
#
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human Metabolome Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
#
# kegg_id
#
# load("hmdbMS1Database0.0.1")
#
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
#
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
#
#
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
#
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
#
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
#
#     if(!is.na(x[2])){
#       return(x[2])
#     }
#
#     if(!is.na(x[3])){
#       return(x[3])
#     }
#
#
#   })
#
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
#
# kegg_id
#
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
#
# annotation_result <- node_data
#
# save(annotation_result, file = "PP/annotation_result")
#
# kegg_id <-
#   kegg_id$KEGG_ID
#
# kegg_id <-
#   unique(kegg_id[!is.na(kegg_id)])
#
# kegg_id
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = kegg_id, database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)
#
# save(path_result, file = "PP/path_result")
#



##-----------------------------------------------------------------------------
####heatmap for pathway in different data points
load("(22,24]/path_result")
load("(22,24]/annotation_result")
annotation_result_22_24 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(22,24]/piumet_output/Result/annotation_result")
annotation_result_22_24 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_22_24, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))

path_result_22_24 <- path_result

#####
load("(24,26]/path_result")
load("(24,26]/annotation_result")
annotation_result_24_26 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(24,26]/piumet_output/Result/annotation_result")
annotation_result_24_26 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_24_26, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))

path_result_24_26 <- path_result

##----
load("(26,28]/path_result")
load("(26,28]/annotation_result")
annotation_result_26_28 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(26,28]/piumet_output/Result/annotation_result")
annotation_result_26_28 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_26_28, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))

path_result_26_28 <- path_result

load("(28,30]/path_result")
load("(28,30]/annotation_result")
annotation_result_28_30 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(28,30]/piumet_output/Result/annotation_result")
annotation_result_28_30 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_28_30, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_28_30 <- path_result

load("(30,32]/path_result")
load("(30,32]/annotation_result")
annotation_result_30_32 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(30,32]/piumet_output/Result/annotation_result")
annotation_result_30_32 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_30_32, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_30_32 <- path_result

load("(32,34]/path_result")
load("(32,34]/annotation_result")
annotation_result_32_34 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(32,34]/piumet_output/Result/annotation_result")
annotation_result_32_34 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_32_34, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_32_34 <- path_result

load("(34,36]/path_result")
load("(34,36]/annotation_result")
annotation_result_34_36 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(34,36]/piumet_output/Result/annotation_result")
annotation_result_34_36 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_34_36, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_34_36 <- path_result

load("(36,38]/path_result")
load("(36,38]/annotation_result")
annotation_result_36_38 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(36,38]/piumet_output/Result/annotation_result")
annotation_result_36_38 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_36_38, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_36_38 <- path_result

load("(38,42]/path_result")
load("(38,42]/annotation_result")
annotation_result_38_42 <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("(38,42]/piumet_output/Result/annotation_result")
annotation_result_38_42 <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_38_42, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_38_42 <- path_result

load("PP/path_result")
load("PP/annotation_result")
annotation_result_pp <- annotation_result %>%
  dplyr::select(KEGG_ID, node) %>%
  dplyr::filter(!is.na(KEGG_ID))

load("PP/piumet_output/Result/annotation_result")
annotation_result_pp <-
  annotation_result %>%
  dplyr::select(variable_id, Metabolite.Name) %>%
  dplyr::left_join(annotation_result_pp, by = c("Metabolite.Name" = "node")) %>%
  dplyr::filter(!is.na(KEGG_ID))
path_result_pp <- path_result

###pahtway
path_result_22_24 <-
  path_result_22_24 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_24_26 <-
  path_result_24_26 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_26_28 <-
  path_result_26_28 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_28_30 <-
  path_result_28_30 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_30_32 <-
  path_result_30_32 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_32_34 <-
  path_result_32_34 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_34_36 <-
  path_result_34_36 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_36_38 <-
  path_result_36_38 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_38_42 <-
  path_result_38_42 %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_result_pp <-
  path_result_pp %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

# all_path_name <- unique(c(path_result_22_24$Pathway.name,
#                           path_result_28_30$Pathway.name,
#                           path_result_26_28$Pathway.name,
#                           path_result_28_30$Pathway.name,
#                           path_result_30_32$Pathway.name,
#                           path_result_32_34$Pathway.name,
#                           path_result_34_36$Pathway.name,
#                           path_result_36_38$Pathway.name,
#                           path_result_38_42$Pathway.name,
#                           path_result_pp$Pathway.name
#                           ))
#
# all_path_id <- unique(c(path_result_22_24$Pathway.ID,
#                           path_result_28_30$Pathway.ID,
#                           path_result_26_28$Pathway.ID,
#                           path_result_28_30$Pathway.ID,
#                           path_result_30_32$Pathway.ID,
#                           path_result_32_34$Pathway.ID,
#                           path_result_34_36$Pathway.ID,
#                           path_result_36_38$Pathway.ID,
#                           path_result_38_42$Pathway.ID,
#                           path_result_pp$Pathway.ID
# ))
#
# path_result_22_24
#
# save(all_path_name, file = "all_path_name")
# save(all_path_id, file = "all_path_id")

load("all_path_name")
load("all_path_id")

###
###annotation result
# library(plyr)
#
# annotation_result_raw <-
#   rbind(
#     annotation_result_22_24,
#     # annotation_result_24_26,
#     annotation_result_26_28,
#     annotation_result_28_30,
#     annotation_result_30_32,
#     annotation_result_32_34,
#     annotation_result_34_36,
#     annotation_result_36_38,
#     annotation_result_38_42,
#     annotation_result_pp
#   ) %>%
#   dplyr::arrange(KEGG_ID) %>%
#   dplyr::mutate(KEGG_ID = stringr::str_replace(KEGG_ID, "\\\t ", ""))
#
# ##remove wrong annoation
# library(plyr)
# annotation_result_raw <-
# annotation_result_raw %>%
# plyr::dlply(.variables = .(KEGG_ID))  %>%
#   purrr::map(.f = function(x){
#     if(nrow(x) == 1){
#       return(x)
#     }
#
#     temp_x <-
#     x %>%
#       dplyr::group_by(KEGG_ID, name) %>%
#       dplyr::summarise(n = n()) %>%
#       dplyr::ungroup()
#
#     if(nrow(temp_x) == 1){
#       return(x)
#     }
#
#     if(any(temp_x$n >= 3) & any(temp_x$n < 2)){
#       temp_x <-
#         temp_x %>%
#         dplyr::filter(n >= 3)
#       x <-
#         x %>%
#         dplyr::filter(name %in% temp_x$name)
#     }
#     x
#   }) %>%
#   do.call(rbind, .)
#
# annotation_result_raw %>%
#   dplyr::group_by(name, KEGG_ID) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::arrange(desc(n))
#
# annotation_result <-
#   annotation_result_raw %>%
#   dplyr::distinct() %>%
#   plyr::dlply(.variables = .(KEGG_ID)) %>%
#   lapply(function(x){
#     x$name <- paste(x$name, collapse = ";")
#     x <-
#     x %>%
#       dplyr::select(KEGG_ID, name,Metabolite.Name) %>%
#       dplyr::distinct()
#     x
#   }) %>%
#   do.call(rbind, .)
#
# save(annotation_result, file = "annotation_result")

load("annotation_result")

dim(subject_data_mean)

plot(density(as.numeric(subject_data_mean[6, ])))

###from metabolite to pathway
temp_data <-
  purrr::map(
    as.data.frame(t(annotation_result)),
    .f = function(x) {
      peak_name <- stringr::str_split(x[2], ";")[[1]]
      subject_data_mean[peak_name, , drop = FALSE] %>%
        apply(2, mean)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

load("hsa_pathway")

###heatmap for each metabolite and pathway
data_for_each_path <-
  purrr::map2(
    .x = all_path_id,
    .y = all_path_name,
    .f = function(x, y) {
      temp_idx <- grep(x, names(hsa_pathway))
      temp_idx <- which(rownames(temp_data) %in% hsa_pathway[[temp_idx]])
      temp_data[temp_idx, , drop = FALSE] %>%
        apply(1, function(x) {
          (x - mean(x)) / sd(x)
        }) %>%
        t() %>%
        as.data.frame() %>%
        data.frame(
          path_name = y,
          .,
          stringsAsFactors = FALSE ,
          check.names = FALSE
        )
    }
  )

names(data_for_each_path) <- all_path_name

library(ComplexHeatmap)
##top 6 pathway are consistant
all_path_name[1:6]
Heatmap(
  as.matrix(data_for_each_path[[5]][, -1]),
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

anno_path1 <-
  annotation_result %>%
  dplyr::filter(KEGG_ID %in% rownames(data_for_each_path[[1]]))

# idx <- 16
# variable_info %>%
#   dplyr::filter(name %in% stringr::str_split(anno_path1$name[[idx]], ";")[[1]]) %>%
#   dplyr::select(Compound.name, KEGG.ID)
#
# anno_path1[idx,]


###times for each pathway
enrichment <-
  rbind(
    path_result_22_24[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(22,24]"),
    
    # path_result_24_26[,c("Pathway.name", "Overlap","p.value.fdr")] %>%
    #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
    #   data.frame(., ga = "(24,26]")
    
    path_result_26_28[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(26,28]"),
    
    path_result_28_30[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(28,30]"),
    
    path_result_30_32[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(30,32]"),
    
    path_result_32_34[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(32,34]"),
    
    path_result_34_36[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(34,36]"),
    
    path_result_36_38[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(36,38]"),
    
    path_result_38_42[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(38,42]"),
    
    path_result_pp[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "PP")
  ) %>%
  as.data.frame() %>%
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", ""))


enrichment %>%
  dplyr::group_by(Pathway.name) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

col_fun = colorRamp2(c(-3, 0, 3), c("#4292C6", "white", "red"))

data_for_each_path2 <-
  data_for_each_path[1:6] %>%
  do.call(rbind, .)

data_for_each_path2$path_name <-
  data_for_each_path2$path_name %>%
  stringr::str_replace(" - Homo sapiens \\(human\\)", "")

data_for_each_path2 <-
  data_for_each_path2 %>%
  tibble::rownames_to_column(var = "kegg_id")

plot <-
  Heatmap(
    matrix = t(as.matrix(data_for_each_path2[, -c(1, 2)])),
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    top_annotation = columnAnnotation(
      pathway = factor(
        data_for_each_path2[, 2],
        levels = c(
          "Steroid hormone biosynthesis",
          "Ovarian steroidogenesis",
          "Cortisol synthesis and secretion",
          "Prolactin signaling pathway",
          "Aldosterone synthesis and secretion",
          "Phenylalanine, tyrosine and tryptophan biosynthesis"
        )
      ),
      col = list(
        pathway = c(
          "Steroid hormone biosynthesis" = ggsci::pal_lancet()(10)[1],
          "Ovarian steroidogenesis" = ggsci::pal_lancet()(10)[2],
          "Cortisol synthesis and secretion" = ggsci::pal_lancet()(10)[3],
          "Prolactin signaling pathway" = ggsci::pal_lancet()(10)[4],
          "Aldosterone synthesis and secretion" = ggsci::pal_lancet()(10)[5],
          "Phenylalanine, tyrosine and tryptophan biosynthesis" = ggsci::pal_lancet()(10)[6]
        )
      )
    ),
    # row_title_rot = 0
    border = TRUE,
    gap = TRUE,
    col = col_fun,
    rect_gp = gpar(col = "white"),
    row_names_rot = 45,
    row_names_gp = gpar(col = text_colour)
  )

plot <- as.ggplot(plot)
plot
# ggsave(plot, filename = "pathway_heatmap.pdf", width = 20, height = 7)


###for each pathway
names(data_for_each_path)

# data_for_each_path[[4]] %>%
#   dplyr::select(-path_name) %>%
#   apply(1, function(x){
#     x <- x - x[1]
#   }) %>%
#   t() %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "metabolite") %>%
#   tidyr::pivot_longer(-metabolite, names_to = "ga", values_to = "value") %>%
#   ggplot(aes(ga, value)) +
#   geom_line(aes(group = metabolite, color = metabolite)) +
#   theme_bw() +
#   ggsci::scale_color_aaas()


path_data <-
  purrr::map(
    all_path_id,
    .f = function(x) {
      temp_idx <- grep(x, names(hsa_pathway))
      temp_idx <- which(rownames(temp_data) %in% hsa_pathway[[temp_idx]])
      temp_data[temp_idx, , drop = FALSE] %>%
        apply(2, mean)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(path_data) <- all_path_name

library(ComplexHeatmap)

path_data <-
  apply(path_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t()

rownames(path_data) <- all_path_name
# stringr::str_replace(" - Homo sapiens \\(human\\)", "")

##calculate p_value
temp_data2 <-
  purrr::map(
    as.data.frame(t(annotation_result)),
    .f = function(x) {
      peak_name <- stringr::str_split(x[2], ";")[[1]]
      # peak_name <- stringr::str_split(x[2], ";")[[1]]
      subject_data[, peak_name , drop = FALSE] %>%
        apply(1, mean)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

path_data2 <-
  purrr::map(
    all_path_id,
    .f = function(x) {
      temp_idx <- grep(x, names(hsa_pathway))
      temp_idx <- which(rownames(temp_data2) %in% hsa_pathway[[temp_idx]])
      temp_data2[temp_idx, , drop = FALSE] %>%
        apply(2, mean)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(path_data) <- all_path_name

ga <- sample_info$GA[match(colnames(path_data2), sample_info$sample_id)]
ga_range <- sample_info$ga_range[match(colnames(path_data2), sample_info$sample_id)]
ga_level <- unique(ga_range) %>% sort()

fdr <-
  apply(path_data2, 1, function(x) {
    x <- as.numeric(x)
    p <- lapply(ga_level[-1], function(y) {
      t.test(x[which(ga_range == ga_level[1])], x[which(ga_range == y)])$p.value
    }) %>% unlist()
    
    fdr <- c(1.3, -log(p.adjust(p, method = "fdr"), 10))
    fdr
  }) %>%
  t() %>%
  as.data.frame()

colnames(fdr) <- ga_level
rownames(fdr) <- all_path_name

fdr <-
  fdr %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "path_id") %>%
  tidyr::pivot_longer(cols = -path_id,
                      names_to = "ga",
                      values_to = "fdr")

path_data <-
  path_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "path_id") %>%
  tidyr::pivot_longer(cols = -path_id,
                      names_to = "ga",
                      values_to = "value")

path_data <-
  path_data %>%
  dplyr::left_join(fdr, by = c("path_id", "ga"))


##add pathway enrichment P value to them
enrichment <-
  rbind(
    path_result_22_24[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(22,24]"),
    
    # path_result_24_26[,c("Pathway.name", "Overlap","p.value.fdr")] %>%
    #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
    #   data.frame(., ga = "(24,26]")
    
    path_result_26_28[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(26,28]"),
    
    path_result_28_30[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(28,30]"),
    
    path_result_30_32[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(30,32]"),
    
    path_result_32_34[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(32,34]"),
    
    path_result_34_36[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(34,36]"),
    
    path_result_36_38[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(36,38]"),
    
    path_result_38_42[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "(38,42]"),
    
    path_result_pp[, c("Pathway.name", "Overlap", "p.value.fdr")] %>%
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>%
      data.frame(., ga = "PP")
  ) %>%
  as.data.frame() %>%
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", ""))

path_data$path_id <-
  path_data$path_id %>%
  stringr::str_replace(" - Homo sapiens \\(human\\)", "")

path_data <-
  path_data %>%
  dplyr::left_join(enrichment, by = c("path_id" = "Pathway.name", "ga" = "ga"))

path_data <-
  path_data %>%
  dplyr::mutate(enriched = case_when(is.na(p.value.fdr) ~ "No", TRUE ~ "Yes"))

###
temp_data <-
  enrichment %>%
  group_by(Pathway.name) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

plot <-
  temp_data %>%
  mutate(Pathway.name = factor(Pathway.name, levels = rev(Pathway.name))) %>%
  ggplot(aes(x = n, y = Pathway.name)) +
  geom_bar(stat = "identity", fill = ggsci::pal_aaas(alpha = 0.7)(10)[5]) +
  theme_classic() +
  labs(x = "#Number", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

plot

# ggsave(
#   plot,
#   file = file.path("pathway_number.pdf"),
#   width = 5,
#   height = 7,
#   bg = "transparent"
# )

text_colour <- c(colorRampPalette(colors = c(
  alpha("#155F83FF", 1),
  alpha("#155F83FF", 0.4),
  alpha("#FFA319FF", 0.4),
  alpha("#FFA319FF", 1)
))(13), "red")

names(text_colour) <-
  c(
    "(10,16]",
    "(16,18]",
    "(18,20]",
    "(20,22]",
    "(22,24]",
    "(24,26]",
    "(26,28]",
    "(28,30]",
    "(30,32]",
    "(32,34]",
    "(34,36]",
    "(36,38]",
    "(38,42]",
    "PP"
  )

range(path_data$value)
path_data$value[path_data$value < -1.95] <- -1.95

plot <-
  path_data %>%
  dplyr::filter(ga != "(10,16]") %>%
  dplyr::filter(ga != "(16,18]") %>%
  dplyr::mutate(path_id = factor(path_id, levels = rev(temp_data$Pathway.name))) %>%
  ggplot(aes(ga, path_id)) +
  geom_point(aes(fill = value, size = fdr),
             shape = 21,
             color = "black") +
  geom_point(aes(shape = enriched), size = 3) +
  scale_shape_manual(values = c(Yes = 3, No = NA)) +
  geom_line(aes(ga, path_id)) +
  geom_hline(yintercept = c(2:length(unique(path_data$path_id))) - 0.5, color = "grey") +
  geom_vline(xintercept = c(2:length(unique(path_data$ga))) - 0.5, color = "grey") +
  # scale_size_continuous(range = c()) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_fill_gradient2(
    low = "#3B4992FF",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(8, 15)) +
  scale_y_discrete(position = "right") +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = text_colour[-c(1, 2)]
    ),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(
      angle = -45,
      hjust = -1,
      vjust = 0
    ),
    legend.position = "right"
  )

plot

# ggsave(
#   plot,
#   file = file.path("pathway.pdf"),
#   width = 20,
#   height = 7,
#   bg = "transparent"
# )

#####pathway way for the decrease peaks
# annotation_result_raw
# down_marker <-
#   lapply(marker_each_point, function(x){
#     if(is.null(x)){
#       return(NULL)
#     }else{
#       x %>%
#         dplyr::filter(class == "decrease") %>%
#         pull(gene_id)
#     }
#   }) %>%
#   unlist() %>%
#   unname() %>%
#   unique()
#
#
# up_marker <-
#   lapply(marker_each_point, function(x){
#     if(is.null(x)){
#       return(NULL)
#     }else{
#       x %>%
#         dplyr::filter(class == "increase") %>%
#         pull(gene_id)
#     }
#   }) %>%
#   unlist() %>%
#   unname() %>%
#   unique()
#
# library(VennDiagram)
#
# plot <-
#   venn.diagram(x = list(up = up_marker, down = down_marker), filename = NULL)
# grid.draw(plot)
#
#
# idx <- match(setdiff(down_marker, up_marker), annotation_result_raw$name)
# idx <- idx[!is.na(idx)]
#
# annotation_result_raw$KEGG_ID[idx]
#
# load("hsa_pathway")
#
# path_result <-
#   enrichPathway(id = annotation_result_raw$KEGG_ID[idx],
#                 database = hsa_pathway)
#
# path_result %>%
#   dplyr::filter(Overlap > 3 & p.value.fdr < 0.05)

####
library(ggridges)

###from metabolite to pathway
temp_data <-
  purrr::map(
    as.data.frame(t(annotation_result)),
    .f = function(x) {
      peak_name <- stringr::str_split(x[2], ";")[[1]]
      subject_data[, peak_name , drop = FALSE] %>%
        apply(1, mean)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

load("hsa_pathway")

data_for_each_path <-
  purrr::map2(
    .x = all_path_id[1:6],
    .y = all_path_name[1:6],
    .f = function(x, y) {
      temp_idx <- grep(x, names(hsa_pathway))
      temp_idx <- which(rownames(temp_data) %in% hsa_pathway[[temp_idx]])
      temp_data[temp_idx, , drop = FALSE] %>%
        apply(1, function(x) {
          (x - mean(x)) / sd(x)
        }) %>%
        t() %>%
        as.data.frame() %>%
        data.frame(
          path_name = y,
          .,
          stringsAsFactors = FALSE ,
          check.names = FALSE
        )
    }
  )

names(data_for_each_path) <- all_path_name[1:6]

###for pathway 1
temp_data <-
  # data_for_each_path[[2]] %>%
  data_for_each_path %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = -path_name,
                      names_to = "sample_id",
                      values_to = "value") %>%
  dplyr::left_join(sample_info[, c("sample_id", "ga_range")], by = "sample_id")

library(ggplot2)
library(ggridges)

plot <-
  temp_data %>%
  dplyr::mutate(ga_range = factor(ga_range, levels = rev(sort(
    unique(temp_data$ga_range)
  )))) %>%
  dplyr::mutate(path_name = stringr::str_replace(path_name, " - Homo sapiens \\(human\\)", "")) %>%
  dplyr::mutate(path_name = factor(
    path_name,
    levels =
      c(
        "Steroid hormone biosynthesis",
        "Ovarian steroidogenesis",
        "Cortisol synthesis and secretion",
        "Prolactin signaling pathway",
        "Aldosterone synthesis and secretion",
        "Phenylalanine, tyrosine and tryptophan biosynthesis"
      )
  )) %>%
  dplyr::filter(
    ga_range %in% c(
      "(10,16]",
      "(18,20]",
      "(22,24]",
      "(26,28]",
      "(30,32]",
      "(34,36]",
      "(38,42]",
      "PP"
    )
  ) %>%
  ggplot(aes(x = value, y = ga_range)) +
  # geom_vline(xintercept = 0) +
  labs(x = "", y = "") +
  geom_density_ridges_gradient(
    aes(fill = ga_range),
    scale = 3,
    size = 0.1,
    show.legend = FALSE
  ) +
  theme(legend.position = "none") +
  theme_bw() +
  scale_fill_manual(values = text_colour) +
  facet_wrap(facets = vars(path_name),
             nrow = 3,
             scales = "free_x")

plot

ggsave(
  plot,
  file = file.path("pathway_ridges.pdf"),
  width = 6,
  height = 10,
  bg = "transparent"
)

######correlation network to find the module
dim(annotation_result)

##remove duplicated KEGG OD
annotation_result1 <-
  annotation_result %>%
  dplyr::distinct(KEGG_ID, .keep_all = TRUE)

##from peak to metabolite
temp_subject_data <-
  t(subject_data) %>%
  as.data.frame()

temp_data <-
  annotation_result1$name %>%
  purrr::map(
    .f = function(x) {
      x <- stringr::str_split(x, ";")[[1]]
      temp_subject_data[x, , drop = FALSE] %>%
        apply(2, mean)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(temp_data) <- annotation_result1$KEGG_ID

# ##calculate correlation
# library(corrr)
# cor_value <- corrr::correlate(x = t(temp_data), method = "spearman")
# cor_value <-
#   cor_value %>%
#   shave() %>%
#   stretch() %>%
#   dplyr::filter(!is.na(r))
#
# p <-
#   as.data.frame(t(cor_value)) %>%
#   purrr::map(.f = function(x){
#     cor.test(as.numeric(temp_data[x[1],]),
#              as.numeric(temp_data[x[2],]), method = "spearman"
#     )$p.value
#   }) %>%
#   unlist()
#
# fdr <- p.adjust(p, method = "fdr")
# fdr[fdr == 0] <- min(fdr[fdr != 0])
#
# cor_value <-
#   data.frame(cor_value, p, fdr, stringsAsFactors = FALSE) %>%
#   dplyr::filter(fdr < 0.05)
#
# save(cor_value, file = "cor_value")
load("cor_value")

cor_value <- cor_value %>%
  dplyr::filter(abs(r) > 0.5)

edge_data <-
  cor_value %>%
  dplyr::mutate(
    from = x,
    to = y,
    fdr = -log(fdr, 10),
    cor = r,
    abs.cor = abs(r)
  ) %>%
  dplyr::select(from, to, cor, abs.cor, fdr)

node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)), stringsAsFactors = FALSE)

node_data <-
  node_data %>%
  dplyr::left_join(annotation_result1, by = c("node" = "KEGG_ID")) %>%
  dplyr::rename(peak_name = name)

library(igraph)
library(tidygraph)

graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


graph <- tidygraph::as.igraph(x = graph)

# subnetworks <-
#   igraph::cluster_edge_betweenness(graph = graph,
#                                    weights = abs(edge_attr(graph,
#                                                            "cor")))
# save(subnetworks, file = "subnetworks")
load("subnetworks")

table(membership(subnetworks))

plot <-
  ggplot(
    data.frame(
      index = 1:length(subnetworks$modularity),
      modu = subnetworks$modularity,
      stringsAsFactors = FALSE
    ),
    aes(index, modu)
  ) +
  geom_vline(
    xintercept = which.max(subnetworks$modularity),
    linetype = 2,
    colour = "#800000B2"
  ) +
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line(colour = "black") +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

plot <-
  plot +
  ggplot2::annotate(
    geom = "point",
    x = which.max(subnetworks$modularity),
    y = max(subnetworks$modularity),
    size = 3,
    colour = "#FFA319FF"
  ) +
  annotate(
    geom = "text",
    x = which.max(subnetworks$modularity),
    y = max(subnetworks$modularity),
    label = paste(
      "(",
      which.max(subnetworks$modularity),
      ",",
      max(subnetworks$modularity) %>% round(3),
      ")"
    ),
    size = 5,
    colour = "#FFA319FF"
  )

plot

# ggsave(plot, filename = "modularity.pdf", width = 7, height = 7)

table(membership(communities = subnetworks))
which(table(membership(communities = subnetworks)) >= 5)
##cluster 1,2,3,8,10,15
node <- vertex_attr(graph = graph, name = "node")

membership <-
  data.frame(node,
             membership = as.numeric(membership(communities = subnetworks)),
             stringsAsFactors = FALSE) %>%
  # dplyr::arrange(membership) %>%
  # dplyr::filter(membership %in% c(1,2,3,4,5,7,8)) %>%
  dplyr::mutate(membership = paste("Cluster", membership, sep = "")) %>%
  dplyr::mutate(membership = case_when(membership %in% c(paste(
    "Cluster", c(1, 2, 3, 8, 10, 15), sep = ""
  )) ~ membership, TRUE ~ "Other"))

graph <-
  igraph::set_vertex_attr(graph = graph,
                          name = "membership",
                          value = membership$membership)

# save(graph, file = "graph")
load("graph")
library(ggraph)

cluster_color <-
  c(ggsci::pal_futurama()(11)[1:6], "grey")

names(cluster_color) <-
  unique(igraph::vertex_attr(graph = graph, name = "membership")) %>%
  stringr::str_sort(numeric = TRUE)

plot <-
  ggraph(graph, layout = 'fr') +
  geom_edge_link(aes(edge_width = fdr, color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(
    aes(size = Degree, fill = membership),
    alpha = 1,
    shape = 21,
    show.legend = TRUE
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = "#3B4992FF",
    mid = "white",
    high = "#EE0000FF",
    midpoint = 0
  ) +
  ggraph::scale_edge_width(range = c(0.3, 1)) +
  # geom_node_text(aes(label = Metabolite.Name), repel = TRUE) +
  scale_size_continuous(range = c(1.5, 5)) +
  # ggsci::scale_fill_futurama(alpha = 0.7) +
  scale_fill_manual(values = cluster_color) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
# coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-1.4,1.4))

plot

# extrafont::font_import()
extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = "cor_network_all.pdf",
#   width = 8,
#   height = 7,
#   bg = "transparent"
# )

###from metabolite to module
membership1 <-
  membership %>%
  dplyr::filter(membership != "Other")

temp_data <-
  lapply(unique(membership1$membership), function(x) {
    x <- membership1$node[membership1$membership == x]
    peak_name <- annotation_result1$name[match(x, annotation_result1$KEGG_ID)]
    peak_name <- stringr::str_split(peak_name, ";") %>%
      unlist() %>%
      unique()
    subject_data_mean[peak_name, ] %>%
      apply(2, mean)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(temp_data) <- unique(membership1$membership)

temp_data <-
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c(ggsci::pal_aaas()(10)[5], "white", ggsci::pal_aaas()(10)[6]))

range(temp_data)

temp_data[temp_data < -2.16] <- -2.16

plot <-
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    border = TRUE,
    col = col_fun
  )

plot

library(ggplotify)

plot <- as.ggplot(plot)

# ggsave(
#   plot,
#   file = "module_heatmap.pdf",
#   width = 7,
#   height = 7,
#   bg = "transparent"
# )

####
##module quantative
ga <- sample_info$GA

ga_range <-
  sample_info$ga_range

table(ga_range)

ga_level <- sort(unique(ga_range))

temp_data <-
  lapply(unique(membership1$membership), function(x) {
    x <- membership1$node[membership1$membership == x]
    peak_name <- annotation_result1$name[match(x, annotation_result1$KEGG_ID)]
    
    temp <-
      lapply(peak_name, function(y) {
        y <- stringr::str_split(y, ";")[[1]]
        apply(subject_data_mean[y, , drop = FALSE], 2, mean)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    rownames(temp) <- x
    temp
  })

names(temp_data) <- unique(membership1$membership)

temp_data <-
  purrr::map2(
    .x = temp_data,
    .y = names(temp_data),
    .f = function(x, y) {
      x <-
        x %>%
        tibble::rownames_to_column(var = "KEGG_ID")
      x <- data.frame(
        x,
        cluster = y,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        check.rows = FALSE
      ) %>%
        dplyr::select(KEGG_ID, cluster, everything())
      x
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

cluster_color <-
  ggsci::pal_futurama()(7)[1:6]

names(cluster_color) <- unique(temp_data$cluster)

temp_data2 <-
  temp_data[, -c(1:2)] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t()

plot <-
  Heatmap(
    matrix = as.matrix(temp_data2),
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    left_annotation = rowAnnotation(
      cluster = factor(temp_data[, 2], levels = unique(temp_data$cluster)),
      col = list(cluster = cluster_color)
    ),
    # row_title_rot = 0
    border = TRUE,
    gap = TRUE,
    col = col_fun,
    rect_gp = gpar(col = "white"),
    column_names_rot = 45,
    column_names_gp = gpar(col = text_colour)
  )

plot

plot <- as.ggplot(plot)

# ggsave(plot,
#        filename = "heatmap_for_each_cluster_metabolite_level.pdf",
#        width = 10, height = 7)

temp_data <-
  temp_data %>%
  tidyr::pivot_longer(
    cols = -c(cluster, KEGG_ID),
    names_to = "ga_range",
    values_to = "value"
  )

plot <-
  temp_data %>%
  dplyr::filter(cluster %in% c("Cluster2", "Cluster3")) %>%
  group_by(ga_range, cluster) %>%
  dplyr::mutate(mean = mean(value), sem = sd(value) / sqrt(nrow(temp_data))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ga_range = factor(ga_range, levels = ga_level)) %>%
  ggplot(aes(ga_range, mean)) +
  geom_line(aes(group = cluster, color = cluster), show.legend = FALSE) +
  geom_errorbar(
    mapping = aes(
      ymin = mean - sem,
      ymax = mean + sem,
      color = cluster
    ),
    width = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_color) +
  labs(x = "", y = "PIUMet module intensity") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text.x = element_text(
      size = 12,
      color = text_colour,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 12)
  )  +
  facet_wrap(~ cluster, ncol = 1, scales = "free_y")

plot

# ggsave(plot,
#        filename = "module_in_different_point.pdf",
#        width = 7,
#        height = 7)

#####grah cluster 2 and 3
membership <- igraph::vertex_attr(graph = graph)$membership
node <- igraph::vertex_attr(graph = graph)$node

igraph::induced_subgraph(graph = graph, v = which(membership %in% c("Cluster2")))

graph2 <- igraph::induced_subgraph(graph = graph, v = which(membership %in% c("Cluster2")))

node2 <-  igraph::vertex_attr(graph = graph2, name = "node")

graph3 <- igraph::induced_subgraph(graph = graph, v = which(membership %in% c("Cluster3")))

node3 <-  igraph::vertex_attr(graph = graph3, name = "node")

##class of node
load("hmdbMS1Database0.0.1")
hmdb_data <- hmdbMS1Database0.0.1@spectra.info

node_class2 <- hmdb_data[match(node2, hmdb_data$KEGG.ID), "Super.class"]

node_pathway2 <-
  lapply(node2, function(x) {
    lapply(hsa_pathway, function(y) {
      x %in% y
    }) %>%
      unlist() %>%
      which() %>%
      names()
  })

graph2 <-
  igraph::set_vertex_attr(graph = graph2,
                          name = "class",
                          value = node_class2)

node_class3 <- hmdb_data[match(node3, hmdb_data$KEGG.ID), "Super.class"]

node_pathway3 <-
  lapply(node3, function(x) {
    lapply(hsa_pathway, function(y) {
      x %in% y
    }) %>%
      unlist() %>%
      which() %>%
      names()
  })

graph3 <-
  igraph::set_vertex_attr(graph = graph3,
                          name = "class",
                          value = node_class3)

plot2 <-
  ggraph(graph2, layout = 'kk') +
  geom_edge_link(aes(edge_width = fdr, color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(
    aes(size = Degree, fill = class),
    alpha = 1,
    shape = 21,
    show.legend = TRUE
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = "#3B4992FF",
    mid = "white",
    high = "#EE0000FF",
    midpoint = 0
  ) +
  
  ggraph::scale_edge_width(range = c(0.1, 1)) +
  geom_node_text(aes(label = Metabolite.Name),
                 repel = TRUE,
                 size = 3) +
  scale_size_continuous(range = c(2, 10)) +
  ggsci::scale_fill_aaas() +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot2

# extrafont::font_import()
extrafont::loadfonts()

# ggsave(
#   plot2,
#   filename = "cluster_plot2.pdf",
#   width = 8,
#   height = 7,
#   bg = "transparent"
# )

plot3 <-
  ggraph(graph3, layout = 'kk') +
  geom_edge_link(aes(edge_width = fdr, color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(
    aes(size = Degree, fill = class),
    alpha = 1,
    shape = 21,
    show.legend = TRUE
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = "#3B4992FF",
    mid = "white",
    high = "#EE0000FF",
    midpoint = 0
  ) +
  
  ggraph::scale_edge_width(range = c(0.1, 1)) +
  geom_node_text(aes(label = Metabolite.Name),
                 repel = TRUE,
                 size = 3) +
  scale_size_continuous(range = c(2, 10)) +
  ggsci::scale_fill_aaas() +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot3

# ggsave(
#   plot3,
#   filename = "cluster_plot3.pdf",
#   width = 8,
#   height = 7,
#   bg = "transparent"
# )

node_pathway2 <-
  node_pathway2 %>%
  unlist() %>%
  table() %>%
  data.frame()

colnames(node_pathway2) <-
  c("pathway", "frequency")

node_pathway2$pathway <-
  node_pathway2$pathway %>%
  stringr::str_split(pattern = ";") %>%
  purrr::map(
    .f = function(x) {
      x[1]
    }
  ) %>%
  unlist() %>%
  stringr::str_replace(" - Homo sapiens \\(human\\)", "")

node_pathway2 <-
  node_pathway2 %>%
  dplyr::arrange(desc(frequency))

pathway2 <-
  node_pathway2 %>%
  dplyr::filter(frequency >= 3) %>%
  dplyr::mutate(pathway = factor(pathway, levels = rev(pathway))) %>%
  ggplot(aes(frequency, pathway)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )


pathway2

# ggsave(pathway2,
#        filename = "pathway2.pdf",
#        width = 7,
#        height = 7)

node_pathway3 <-
  node_pathway3 %>%
  unlist() %>%
  table() %>%
  data.frame()

colnames(node_pathway3) <-
  c("pathway", "frequency")

node_pathway3$pathway <-
  node_pathway3$pathway %>%
  stringr::str_split(pattern = ";") %>%
  purrr::map(
    .f = function(x) {
      x[1]
    }
  ) %>%
  unlist() %>%
  stringr::str_replace(" - Homo sapiens \\(human\\)", "")

node_pathway3 <-
  node_pathway3 %>%
  dplyr::arrange(desc(frequency))

pathway3 <-
  node_pathway3 %>%
  # dplyr::filter(frequency >=3) %>%
  dplyr::mutate(pathway = factor(pathway, levels = rev(pathway))) %>%
  ggplot(aes(frequency, pathway)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )


pathway3

# ggsave(pathway3,
#        filename = "pathway3.pdf",
#        width = 7,
#        height = 7)
