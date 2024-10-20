##avoid source 
no_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)
source("R/20200727/tools.R")

##load data
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info")

setwd("data_analysis20200108/urine_metabolome/DEG_analysis/")

# load("data_preparation_for_analysis/metabolite_table")
# load("data_preparation_for_analysis/metabolite_tags")
# load("data_preparation_for_analysis/peak_table")

# info <-
#   readxl::read_xlsx("/Users/shenxt/projects/smartD/patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
# info <-
#   info %>%
#   mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
#   mutate(ID = paste("SF", ID, sep = ""))
# 
# sample_info <-
#   readr::read_csv("/Users/shenxt/projects/smartD/patient information/sample_info_191021.csv")

# sxtTools::setwd_project()
# marker1 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/GA_prediction/marker_rf_final.csv")
# marker2 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/marker_rf_final.csv")

dim(expression_data)

rownames(expression_data)
variable_info$name

rownames(expression_data) == variable_info$name
colnames(expression_data) == sample_info$sample_id

###find marker which are change according to pregnancy
##remove QC and blank samples
sample_info <-
  sample_info %>%
  dplyr::filter(!stringr::str_detect(sample_id, "QC")) 

expression_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

expression_data[5,] %>% 
  as.numeric() %>% 
  density() %>%
  plot()

subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

#log transformation
subject_data <-
  log(subject_data + 1, 10)

subject_data[5,] %>% 
  as.numeric() %>% 
  density() %>%
  plot()

dim(subject_data)
dim(variable_info)

rownames(subject_data) <- variable_info$name

sample_info$GA[is.na(sample_info$GA)] <- 45

# ###remove samples which are after birth
# remove_idx <- which(sample_info$GA == 45)
# # remove_name <- sample_info$sample_id[remove_idx]
# 
# #SAM analysi
# # sam_test <-
# # samr::SAM(x = subject_data[,-remove_idx],
# #           y = (sample_info$GA)[-remove_idx],
# #           resp.type = "Quantitative",
# #           geneid = rownames(subject_data),
# #           genenames = rownames(subject_data),
# #           return.x = FALSE,
# #           fdr.output = 0.05,
# #           regression.method = "ranks",
# #           random.seed = "123",
# #           nperms = 1000)
# # 
# # save(sam_test, file = "sam_test")
# 
# load("sam_test")
# 
# gene_up <- sam_test$siggenes.table$genes.up
# gene_down <- sam_test$siggenes.table$genes.lo
# 
# peak_marker <- c(gene_up[,1], gene_down[,1])
# 
# plot(sample_info$GA, subject_data[peak_marker[1],])
# 
# plot(sample_info$GA, expression_data[peak_marker[1],])
# 
# plot(sample_info$GA, expression_data[gene_down[4,1],])
# 
# ###linear regression
# subject_data <- 
#   subject_data[peak_marker,]
# 
# rownames(subject_data) <- peak_marker
# 
# colnames(subject_data) == sample_info$sample_id
# 
# ga <- sample_info$GA
# batch <- sample_info$batch
# bmi <- sample_info$bmi
# age <- sample_info$Age
# parity <- sample_info$parity
# ethnicity <- sample_info$ethnicity
# 
# remove_idx <- which(ga == 45)
# 
# lm_p_value <- 
# apply(subject_data[,-remove_idx], 1, function(x){
#   x <- as.numeric(x)
#   lm_reg <- lm(formula = ga[-remove_idx] ~ x + batch[-remove_idx] + 
#                  bmi[-remove_idx] + age[-remove_idx] + parity[-remove_idx] +
#                  ethnicity[-remove_idx])
#   temp <- summary(lm_reg)$coefficients  %>% as.data.frame()
#   as.numeric(temp["x",4])
# })
# 
# lm_fdr <- p.adjust(lm_p_value, method = "fdr")
# 
# ###correlation
# 
# # cor_value <-
# #   purrr::map(.x = as.data.frame(t(subject_data[,-remove_idx])), .f = function(x){
# #     temp1 <- cor.test(x, sample_info$GA[-remove_idx], method = "spearman")
# #     temp2 <- cor.test(sxt_rank(x), sample_info$GA[-remove_idx], method = "spearman")
# #     c(temp1$estimate, temp1$p.value, temp2$estimate, temp2$p.value)
# #   })
# # 
# # cor_value <-
# #   do.call(rbind, cor_value)
# # 
# # colnames(cor_value) <- c("correlation1", "p1", "correlation2", "p2")
# # 
# # cor_value <- data.frame(name = rownames(subject_data),
# #                         cor_value, stringsAsFactors = FALSE)
# # 
# # 
# # cor_value$p1_adj <- p.adjust(p = cor_value$p1, method = "fdr")
# # 
# # cor_value$p2_adj <- p.adjust(p = cor_value$p2, method = "fdr")
# # 
# # plot(cor_value$correlation1)
# # 
# # plot(cor_value$correlation2)
# # 
# # save(cor_value, file = "cor_value")
# 
# load("cor_value")
# 
# peak_marker <-
#   rbind(gene_up, gene_down) %>%
#   as.data.frame() %>%
#   dplyr::mutate(score = as.numeric(`Score(d)`)) %>%
#   dplyr::left_join(cor_value, by = c("Gene ID" = "name")) %>%
#   # dplyr::arrange(score) %>%
#   dplyr::mutate(index = 1:(nrow(gene_up) + nrow(gene_down))) %>%
#   dplyr::mutate(class = case_when(
#     score > 0 & correlation1 > 0 ~ "up",
#     score < 0 & correlation1 < 0 ~ "down",
#     TRUE ~ "no"
#   ))
# 
# 
# peak_marker$lm_fdr <- lm_fdr
# 
# save(peak_marker, file = "peak_marker")
# load("peak_marker")
# 
# 
# plot <- peak_marker %>% 
#   ggplot(aes(score, correlation1)) +
#   geom_hline(yintercept = 0, linetype = 2, color = "black") +
#   geom_vline(xintercept = 0, linetype = 2, color = "black") +
#   geom_point(aes(color = class, size = -log(p1_adj, 10)), 
#              show.legend = TRUE, alpha = 0.8) +
#   labs(x = "Score (SAM test)", y = "Correlation (Spearman)") +
#   scale_color_manual(values = c("up" = ggsci::pal_aaas()(10)[2],
#                                 "down" = ggsci::pal_aaas()(10)[1],
#                                 "no" = "#D9D9D9")) +
#   scale_size_continuous(range = c(0.05, 5)) +
#   guides(size = guide_legend(override.aes = list(color = "black"))) +
#   # scale_colour_brewer() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,0), 
#         legend.justification = c(1,0),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 13)) 
# 
# plot
# 
# ggsave(plot, filename = "dem_plot_light.pdf",
#        width = 7, height = 7, bg = "transparent")
# 
# ggsave(plot, filename = "dem_plot_light.png",
#        width = 7, height = 7, bg = "transparent")
# 
# table(peak_marker$class)
# 
# variable_info <- 
#   variable_info[match(rownames(subject_data), variable_info$name),]
# 
# 
# ###remove lm_fdr > 0.05
# remain_idx <- which(peak_marker$lm_fdr < 0.05)
# 
# peak_marker <- peak_marker[remain_idx,]
# variable_info <- variable_info[remain_idx,]
# subject_data <- subject_data[remain_idx,]

##ga range
dim(sample_info)
dim(subject_data)

sample_info$sample_id == colnames(subject_data)

ga_range <- 
as.character(cut_width(x = sample_info$GA, width = 2, center = 3)) %>%
  stringr::str_replace("\\[", '(')

stringr::str_sort(unique(ga_range))

table(ga_range)

## for each person, just combine the samples in the sample ga range
sample_info$ga_range <-
  ga_range

sample_info$ga_range %>% unique() %>% sort()

sample_info$ga_range %>% table()

sample_info$ga_range[sample_info$ga_range == "(10,12]"] <- "(10,16]"
sample_info$ga_range[sample_info$ga_range == "(12,14]"] <- "(10,16]"
sample_info$ga_range[sample_info$ga_range == "(14,16]"] <- "(10,16]"

sample_info$ga_range[sample_info$ga_range == "(38,40]"] <- "(38,42]"
sample_info$ga_range[sample_info$ga_range == "(40,42]"] <- "(38,42]"

sample_info$ga_range[sample_info$ga_range == "(44,46]"] <- "PP"

###combine different samples in one ga range together
library(plyr)

subject_data <- 
  apply(subject_data, 1, function(x){
    (x) / sd(x)
  })

subject_data2 <-
  subject_data %>%
  data.frame(., ga_range = sample_info$ga_range, stringsAsFactors = FALSE) %>%
  mutate(ga_range = factor(
    ga_range,
    levels = sample_info$ga_range %>% unique() %>% stringr::str_sort(numeric = TRUE)
  )) %>%
  plyr::dlply(.variables = .(ga_range))

# subject_data_mean <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, mean)
#   }) %>%
#   do.call(cbind, . )
# 
# subject_data_sd <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, sd)
#   }) %>%
#   do.call(cbind, . )
# 
# subject_data_sem <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, function(y){
#       sd(y)/(nrow(x) - 1)
#     })
#   }) %>%
#   do.call(cbind, . )
# 
# save(subject_data_mean, file = "subject_data_mean")
# save(subject_data_sd, file = "subject_data_sd")
# save(subject_data_sem, file = "subject_data_sem")

load("subject_data_mean")
load("subject_data_sd")
load("subject_data_sem")

subject_data2 <- 
  lapply(subject_data2, function(x){
    x <-
      x %>% 
      dplyr::select(-ga_range)
  })

# plot(subject_data_mean[peak_marker$`Gene ID`[300],])
# temp_data <- 
#   purrr::map2(subject_data2, names(subject_data2), .f = function(x,y){
#     data.frame(value = x[,peak_marker$`Gene ID`[1]], 
#                point = y, stringsAsFactors = FALSE)
#   }) %>% 
#   do.call(rbind, .)
# 
# temp_data %>% 
#   mutate(point = factor(point, levels = temp_data$point %>% unique())) %>% 
#   ggplot(aes(point, value)) +
#   geom_boxplot(fill = "white") +
#   geom_jitter(color = "black") +
#   theme_classic()

## find all the peaks in different time points
# fc_p_value <-
#   pbapply::pblapply(subject_data2[-1], function(x){
#     p_value <- lapply(1:ncol(x), function(idx){
#       t.test(x[,idx], subject_data2[[1]][,idx])$p.value
#     }) %>%
#       unlist() %>%
#       p.adjust(method = "fdr")
# 
#     fc <- lapply(1:ncol(x), function(idx){
#       mean(x[,idx]) /mean(subject_data2[[1]][,idx])
#     }) %>%
#       unlist()
# 
#     fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
# 
#     data.frame(p_value, fc, stringsAsFactors = FALSE)
#   })
# 
# save(fc_p_value, file = "fc_p_value")

load("fc_p_value")

names(fc_p_value)

fc_p_value$`(14,16]`

dir.create("marker_in_different_points")

# for(idx in 1:length(fc_p_value)){
#   cat(idx, " ")
#   plot <-
#     volcano_plot(fc = fc_p_value[[idx]][,2],
#                  p_value = fc_p_value[[idx]][,1],
#                  p.cutoff = 0.05, fc.cutoff = 1,
#                  theme = "light")
#   
#   plot <-
#     plot +
#     labs(title = paste(names(fc_p_value)[idx], "(4,12]", sep = "/"))
#   
#   plot
#   
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_light.pdf", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
#   
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_light.png", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
# 
# }
# 
# 
##find markers for each time points
marker_each_point <- 
  lapply(fc_p_value, function(x){
    idx1 <- which(x$p_value < 0.05 & x$fc > 1)
    idx2 <- which(x$p_value < 0.05 & x$fc < 1)
    
    gene1 <- 
      try(
        data.frame(gene_id = variable_info$name[idx1],
                   x[idx1,],
                   class = "increase",
                   stringsAsFactors = FALSE
        ),silent = TRUE 
      )
    
    if(class(gene1) == "try-error"){
      gene1 <- NULL
    }
    
    gene2 <- 
      try(
        data.frame(gene_id = variable_info$name[idx2],
                   x[idx2,],
                   class = "decrease",
                   stringsAsFactors = FALSE
        ),silent = TRUE
      )
    
    if(class(gene2) == "try-error"){
      gene2 <- NULL
    }
    
    rbind(gene1, gene2)
  })


marker_each_point[[10]]

names(marker_each_point)

#####a sankey 
marker_each_point %>% 
  lapply(nrow) %>% 
  unlist()

all_marker_name <- 
  lapply(marker_each_point, function(x){
    x$gene_id
  }) %>% 
  unlist() %>% 
  unique()


##disease pathway enrichment for each data point
#22_24
load("marker_in_different_points/(22,24]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
lapply(hmdb_id, function(x){
  metflow2::transID(query = x,
                    from = "Human Metabolome Database",
                    to = "KEGG", top = 1)
}) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
apply(kegg_id, 1, function(x){
  x <- as.character(x)
  if(is.na(x[2]) & is.na(x[3])){
    return(NA)
  }

  if(!is.na(x[2])){
    return(x[2])
  }

  if(!is.na(x[3])){
    return(x[3])
  }


})

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

node_data <-
node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(22,24]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

load("hsa_disease_pathway")


path_disease_result <-
enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result <-
path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)


save(path_disease_result, file = "marker_in_different_points/(22,24]/path_disease_result")


##28_30
load("marker_in_different_points/(24,26]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))


node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(24,26]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

  path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)


save(path_disease_result, file = "marker_in_different_points/(24,26]/path_disease_result")



##26_28
load("marker_in_different_points/(26,28]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))


node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(26,28]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(26,28]/path_disease_result")



##28_30
load("marker_in_different_points/(28,30]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

# load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))


node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(28,30]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(28,30]/path_disease_result")



##30_32
load("marker_in_different_points/(30,32]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(30,32]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

kegg_id

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(30,32]/path_disease_result")


##32_34
load("marker_in_different_points/(32,34]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]

kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(32,34]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

kegg_id

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(32,34]/path_disease_result")





##34_36
load("marker_in_different_points/(34,36]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(34,36]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

kegg_id

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(34,36]/path_disease_result")



##36_38
load("marker_in_different_points/(36,38]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(36,38]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

kegg_id

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(36,38]/path_disease_result")





##38_42
load("marker_in_different_points/(38,42]/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/(38,42]/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

kegg_id

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/(38,42]/path_disease_result")







##PP
load("marker_in_different_points/PP/piumet_output/Result/node_data")
hmdb_id <- unique(node_data$HMDB_ID[node_data$HMDB_ID != " "])

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

# save(annotation_result, file = "marker_in_different_points/PP/annotation_result")

kegg_id <-
  kegg_id$KEGG_ID

kegg_id <-
  unique(kegg_id[!is.na(kegg_id)])

kegg_id

load("hsa_disease_pathway")

path_disease_result <-
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)

path_disease_result %>%
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_disease_result, file = "marker_in_different_points/PP/path_disease_result")


##-----------------------------------------------------------------------------
####heatmap for pathway in different data points
load("marker_in_different_points/(22,24]/path_disease_result")
load("marker_in_different_points/(22,24]/annotation_result")
annotation_result_22_24 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(22,24]/piumet_output/Result/annotation_result")
annotation_result_22_24 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_22_24, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))

path_disease_result_22_24 <- path_disease_result


#####
load("marker_in_different_points/(24,26]/path_disease_result")
load("marker_in_different_points/(24,26]/annotation_result")
annotation_result_24_26 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(24,26]/piumet_output/Result/annotation_result")
annotation_result_24_26 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_24_26, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))

path_disease_result_24_26 <- path_disease_result


##----
load("marker_in_different_points/(26,28]/path_disease_result")
load("marker_in_different_points/(26,28]/annotation_result")
annotation_result_26_28 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(26,28]/piumet_output/Result/annotation_result")
annotation_result_26_28 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_26_28, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))

path_disease_result_26_28 <- path_disease_result



load("marker_in_different_points/(28,30]/path_disease_result")
load("marker_in_different_points/(28,30]/annotation_result")
annotation_result_28_30 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(28,30]/piumet_output/Result/annotation_result")
annotation_result_28_30 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_28_30, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_28_30 <- path_disease_result



load("marker_in_different_points/(30,32]/path_disease_result")
load("marker_in_different_points/(30,32]/annotation_result")
annotation_result_30_32 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(30,32]/piumet_output/Result/annotation_result")
annotation_result_30_32 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_30_32, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_30_32 <- path_disease_result




load("marker_in_different_points/(32,34]/path_disease_result")
load("marker_in_different_points/(32,34]/annotation_result")
annotation_result_32_34 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(32,34]/piumet_output/Result/annotation_result")
annotation_result_32_34 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_32_34, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_32_34 <- path_disease_result




load("marker_in_different_points/(34,36]/path_disease_result")
load("marker_in_different_points/(34,36]/annotation_result")
annotation_result_34_36 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(34,36]/piumet_output/Result/annotation_result")
annotation_result_34_36 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_34_36, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_34_36 <- path_disease_result


load("marker_in_different_points/(36,38]/path_disease_result")
load("marker_in_different_points/(36,38]/annotation_result")
annotation_result_36_38 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(36,38]/piumet_output/Result/annotation_result")
annotation_result_36_38 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_36_38, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_36_38 <- path_disease_result

load("marker_in_different_points/(38,42]/path_disease_result")
load("marker_in_different_points/(38,42]/annotation_result")
annotation_result_38_42 <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/(38,42]/piumet_output/Result/annotation_result")
annotation_result_38_42 <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_38_42, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_38_42 <- path_disease_result

load("marker_in_different_points/PP/path_disease_result")
load("marker_in_different_points/PP/annotation_result")
annotation_result_pp <- annotation_result %>% 
  dplyr::select(KEGG_ID, node) %>% 
  dplyr::filter(!is.na(KEGG_ID))

load("marker_in_different_points/PP/piumet_output/Result/annotation_result")
annotation_result_pp <- 
  annotation_result %>% 
  dplyr::select(name, Metabolite.Name) %>% 
  dplyr::left_join(annotation_result_pp, by = c("Metabolite.Name" = "node")) %>% 
  dplyr::filter(!is.na(KEGG_ID))
path_disease_result_pp <- path_disease_result



###pahtway
path_disease_result_22_24 <- 
path_disease_result_22_24 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_24_26 <- 
  path_disease_result_24_26 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_26_28 <- 
  path_disease_result_26_28 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_28_30 <- 
  path_disease_result_28_30 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_30_32 <- 
  path_disease_result_30_32 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_32_34 <- 
  path_disease_result_32_34 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)


path_disease_result_34_36 <- 
  path_disease_result_34_36 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_36_38 <- 
  path_disease_result_36_38 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_38_42 <- 
  path_disease_result_38_42 %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

path_disease_result_pp <- 
  path_disease_result_pp %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

all_path_name <- unique(c(path_disease_result_22_24$Pathway.name,
                          path_disease_result_28_30$Pathway.name,
                          path_disease_result_26_28$Pathway.name,
                          path_disease_result_28_30$Pathway.name,
                          path_disease_result_30_32$Pathway.name,
                          path_disease_result_32_34$Pathway.name,
                          path_disease_result_34_36$Pathway.name,
                          path_disease_result_36_38$Pathway.name,
                          path_disease_result_38_42$Pathway.name,
                          path_disease_result_pp$Pathway.name
                          ))


all_path_id <- unique(c(path_disease_result_22_24$Pathway.ID,
                          path_disease_result_28_30$Pathway.ID,
                          path_disease_result_26_28$Pathway.ID,
                          path_disease_result_28_30$Pathway.ID,
                          path_disease_result_30_32$Pathway.ID,
                          path_disease_result_32_34$Pathway.ID,
                          path_disease_result_34_36$Pathway.ID,
                          path_disease_result_36_38$Pathway.ID,
                          path_disease_result_38_42$Pathway.ID,
                          path_disease_result_pp$Pathway.ID
))

all_path_name <- all_path_name[1]
all_path_id <- all_path_id[1]

path_disease_result_22_24

###
###annotation result
library(plyr)

annotation_result_raw <- 
  rbind(
    annotation_result_22_24,
    # annotation_result_24_26,
    annotation_result_26_28,
    annotation_result_28_30,
    annotation_result_30_32,
    annotation_result_32_34,
    annotation_result_34_36,
    annotation_result_36_38,
    annotation_result_38_42,
    annotation_result_pp
  ) %>% 
  dplyr::arrange(KEGG_ID) %>% 
  dplyr::mutate(KEGG_ID = stringr::str_replace(KEGG_ID, "\\\t ", ""))

##remove wrong annoation
library(plyr)
annotation_result_raw <- 
annotation_result_raw %>% 
plyr::dlply(.variables = .(KEGG_ID))  %>% 
  purrr::map(.f = function(x){
    if(nrow(x) == 1){
      return(x)
    }
    
    temp_x <- 
    x %>% 
      dplyr::group_by(KEGG_ID, name) %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::ungroup()
    
    if(nrow(temp_x) == 1){
      return(x)
    } 
    
    if(any(temp_x$n >= 3) & any(temp_x$n < 2)){
      temp_x <- 
        temp_x %>% 
        dplyr::filter(n >= 3)
      x <- 
        x %>% 
        dplyr::filter(name %in% temp_x$name)
    }
    x
  }) %>% 
  do.call(rbind, .)


annotation_result_raw %>% 
  dplyr::group_by(name, KEGG_ID) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n))

annotation_result <- 
  annotation_result_raw %>% 
  dplyr::distinct() %>% 
  plyr::dlply(.variables = .(KEGG_ID)) %>% 
  lapply(function(x){
    x$name <- paste(x$name, collapse = ";")
    x <- 
    x %>% 
      dplyr::select(KEGG_ID, name,Metabolite.Name) %>% 
      dplyr::distinct()
    x
  }) %>% 
  do.call(rbind, .)


dim(subject_data_mean)

plot(density(as.numeric(subject_data_mean[6,])))

###from metabolite to pathway
temp_data <- 
purrr::map(
  as.data.frame(t(annotation_result)),
  .f = function(x) {
    peak_name <- stringr::str_split(x[2], ";")[[1]]
    subject_data_mean[peak_name, ,drop = FALSE] %>% 
      apply(2, mean)
  }
) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

load("hsa_disease_pathway")

###heatmap for each metabolite and pathway
data_for_each_path <- 
purrr::map2(.x = all_path_id, .y = all_path_name,.f = function(x, y){
  temp_idx <- grep(x, names(hsa_disease_pathway))
  temp_idx <- which(rownames(temp_data) %in% hsa_disease_pathway[[temp_idx]])
  temp_data[temp_idx,,drop = FALSE] %>% 
    apply(1, function(x){
      (x - mean(x))/sd(x)
    }) %>% 
    t() %>% 
    as.data.frame() %>% 
    data.frame(path_name = y, ., stringsAsFactors = FALSE , check.names = FALSE)
}
)

names(data_for_each_path) <- all_path_name

library(ComplexHeatmap)
##top 6 pathway are consistant
all_path_name[1]
Heatmap(as.matrix(data_for_each_path[[1]][,-1]), cluster_columns = FALSE, cluster_rows = FALSE)

densityHeatmap(data = as.matrix(data_for_each_path[[1]][,-1]))

anno_path1 <- 
  annotation_result %>% 
  dplyr::filter(KEGG_ID %in% rownames(data_for_each_path[[1]]))

idx <- 1
variable_info %>% 
  dplyr::filter(name %in% stringr::str_split(anno_path1$name[[idx]], ";")[[1]]) %>% 
  dplyr::select(Compound.name, KEGG.ID)

anno_path1[idx,]


###times for each pathway
enrichment <- 
  rbind(
    # path_disease_result_22_24[,c("Pathway.name", "Overlap", "p.value.fdr")] %>% 
    #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    #   data.frame(., ga = "(22,24]"),
    # 
    # path_disease_result_24_26[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    #   data.frame(., ga = "(24,26]")
    
    path_disease_result_26_28[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(26,28]"),
    
    path_disease_result_28_30[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(28,30]"),
    
    path_disease_result_30_32[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(30,32]"),
    
    path_disease_result_32_34[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(32,34]"),
    
    path_disease_result_34_36[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(34,36]"),
    
    path_disease_result_36_38[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(36,38]"),
    
    path_disease_result_38_42[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(38,42]")
    
    # path_disease_result_pp[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    #   data.frame(., ga = "PP")
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

text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13),
  "red"
)

rownames(data_for_each_path2)

plot <- 
  Heatmap(matrix = t(as.matrix(data_for_each_path2[,-c(1,2)])),
          cluster_row_slices = FALSE, 
          cluster_column_slices = FALSE, 
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    # top_annotation = columnAnnotation(
    #   pathway = factor(data_for_each_path2[,2], 
    #                    levels = c("Steroid hormone biosynthesis", 
    #                               "Ovarian steroidogenesis", 
    #                               "Cortisol synthesis and secretion",
    #                               "Prolactin signaling pathway",
    #                               "Aldosterone synthesis and secretion",
    #                               "Phenylalanine, tyrosine and tryptophan biosynthesis")
    #                    ),
    #   col = list(
    #     pathway = c(
    #       "Steroid hormone biosynthesis" = ggsci::pal_lancet()(10)[1], 
    #       "Ovarian steroidogenesis" = ggsci::pal_lancet()(10)[2], 
    #       "Cortisol synthesis and secretion" = ggsci::pal_lancet()(10)[3],
    #       "Prolactin signaling pathway" = ggsci::pal_lancet()(10)[4],
    #       "Aldosterone synthesis and secretion" = ggsci::pal_lancet()(10)[5],
    #       "Phenylalanine, tyrosine and tryptophan biosynthesis" = ggsci::pal_lancet()(10)[6]
    #     )
    #   )
    # ),
    # row_title_rot = 0
    border = TRUE,
    gap = TRUE, 
    col = col_fun,
    rect_gp = gpar(col= "white"), 
    row_names_rot = -45, 
    row_names_gp = gpar(col = text_colour)
  )


plot <- as.ggplot(plot)
plot
ggsave(plot, filename = "pathway_disease_heatmap.pdf", width = 3, height = 7)


###for each pathway
names(data_for_each_path)

data_for_each_path[[1]] %>% 
  dplyr::select(-path_name) %>% 
  apply(1, function(x){
    x <- x - x[1]
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "metabolite") %>% 
  tidyr::pivot_longer(-metabolite, names_to = "ga", values_to = "value") %>% 
  ggplot(aes(ga, value)) +
  geom_line(aes(group = metabolite, color = metabolite)) +
  theme_bw() +
  ggsci::scale_color_aaas()


path_data <- 
purrr::map(all_path_id, .f = function(x){
  temp_idx <- grep(x, names(hsa_disease_pathway))
  temp_idx <- which(rownames(temp_data) %in% hsa_disease_pathway[[temp_idx]])
  temp_data[temp_idx,,drop = FALSE] %>% 
    apply(2, mean)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(path_data) <- all_path_name

library(ComplexHeatmap)

path_data <- 
  apply(path_data, 1, function(x){
    (x - mean(x))/sd(x)
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
      subject_data[,peak_name ,drop = FALSE] %>% 
        apply(1, mean)
    }
  ) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()



path_data2 <- 
  purrr::map(all_path_id, .f = function(x){
    temp_idx <- grep(x, names(hsa_disease_pathway))
    temp_idx <- which(rownames(temp_data2) %in% hsa_disease_pathway[[temp_idx]])
    temp_data2[temp_idx,,drop = FALSE] %>% 
      apply(2, mean)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(path_data) <- all_path_name

ga <- sample_info$GA[match(colnames(path_data2), sample_info$sample_id)]
ga_range <- sample_info$ga_range[match(colnames(path_data2), sample_info$sample_id)]
ga_level <- unique(ga_range) %>% sort()

fdr <- 
apply(path_data2, 1, function(x){
x <- as.numeric(x)  
p <- lapply(ga_level[-1], function(y){
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
  tidyr::pivot_longer(cols = -path_id, names_to = "ga", values_to = "fdr")

path_data <- 
path_data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "path_id") %>% 
  tidyr::pivot_longer(cols = -path_id, names_to = "ga", values_to = "value")

path_data <- 
path_data %>% 
  dplyr::left_join(fdr, by = c("path_id", "ga"))


##add pathway enrichment P value to them
enrichment <- 
rbind(
  # path_disease_result_22_24[,c("Pathway.name", "Overlap", "p.value.fdr")] %>% 
  #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
  #   data.frame(., ga = "(22,24]"),
  # 
  # path_disease_result_24_26[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
  #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
  #   data.frame(., ga = "(24,26]")
  
  path_disease_result_26_28[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(26,28]"),
  
  path_disease_result_28_30[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(28,30]"),
  
  path_disease_result_30_32[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(30,32]"),
  
  path_disease_result_32_34[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(32,34]"),
  
  path_disease_result_34_36[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(34,36]"),
  
  path_disease_result_36_38[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(36,38]"),
  
  path_disease_result_38_42[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(38,42]")
  
  # path_disease_result_pp[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
  #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
  #   data.frame(., ga = "PP")
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
  dplyr::mutate(enriched = case_when(
    is.na(p.value.fdr) ~ "No",
    TRUE ~ "Yes"
  ))


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
  geom_bar(stat = "identity", 
           fill = ggsci::pal_aaas(alpha = 0.7)(10)[5]) +
  theme_classic() +
  labs(x = "#Number", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0,0))) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10)) 

plot




text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13),
  "red"
)

names(text_colour) <- 
  c("(10,16]", "(16,18]", "(18,20]", "(20,22]", "(22,24]", "(24,26]", "(26,28]", 
    "(28,30]", "(30,32]", "(32,34]", "(34,36]","(36,38]","(38,42]","PP")


range(path_data$value)
path_data$value[path_data$value > 1.23] <- 1.23

plot <- 
path_data %>%
  dplyr::filter(ga != "(10,16]") %>% 
  dplyr::filter(ga != "(16,18]") %>% 
  dplyr::mutate(path_id = factor(path_id, levels = rev(temp_data$Pathway.name))) %>% 
ggplot(aes(ga, path_id)) +
  geom_point(aes(fill = value, size = fdr), shape = 21, color = "black") +
  geom_point(aes(shape = enriched), size = 3) +
  scale_shape_manual(values = c(Yes = 3, No = NA)) +
  geom_line(aes(ga, path_id)) +
  geom_hline(yintercept = c(2:length(unique(path_data$path_id)))-0.5, color = "grey") +
  geom_vline(xintercept = c(2:length(unique(path_data$ga)))-0.5, color = "grey") +
  # scale_size_continuous(range = c()) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_fill_gradient2(low = "#3B4992FF", 
                       mid = "white", 
                       high = "red", 
                       midpoint = 0) +
  scale_size_continuous(range = c(8,15)) +
  scale_y_discrete(position = "right") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   color = text_colour[-c(1,2)]),
        axis.text = element_text(size = 10),
        axis.text.y = element_text(angle = -45, hjust = -1, vjust = 0),
        legend.position = "right"
        )

plot

ggsave(
  plot,
  file = file.path("marker_in_different_points", "pathway.pdf"),
  width = 20,
  height = 7,
  bg = "transparent"
)

####
library(ggridges)

###from metabolite to pathway
temp_data <- 
  purrr::map(
    as.data.frame(t(annotation_result)),
    .f = function(x) {
      peak_name <- stringr::str_split(x[2], ";")[[1]]
      subject_data[,peak_name ,drop = FALSE] %>% 
        apply(1, mean)
    }
  ) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

load("hsa_disease_pathway")

data_for_each_path <- 
  purrr::map2(.x = all_path_id[1], 
              .y = all_path_name[1],.f = function(x, y){
    temp_idx <- grep(x, names(hsa_disease_pathway))
    temp_idx <- which(rownames(temp_data) %in% hsa_disease_pathway[[temp_idx]])
    temp_data[temp_idx,,drop = FALSE] %>% 
      apply(1, function(x){
        (x - mean(x))/sd(x)
      }) %>% 
      t() %>% 
      as.data.frame() %>% 
      data.frame(path_name = y, ., stringsAsFactors = FALSE , check.names = FALSE)
  }
  )

names(data_for_each_path) <- all_path_name[1]

###for pathway 1
temp_data <- 
# data_for_each_path[[2]] %>%
  data_for_each_path %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  tidyr::pivot_longer(cols = -path_name, names_to = "sample_id", values_to = "value") %>% 
  dplyr::left_join(sample_info[,c("sample_id", "ga_range")], by = "sample_id")
library(ggplot2)
library(ggridges)

plot <- 
temp_data %>% 
  dplyr::mutate(ga_range = factor(ga_range, 
                                  levels = rev(sort(unique(temp_data$ga_range))))) %>% 
  dplyr::mutate(path_name = stringr::str_replace(path_name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(path_name = factor(path_name, levels =
                                     c("Cushing syndrome - Homo sapiens (human)"))) %>%
  # dplyr::filter(
  #   ga_range %in% c("(10,16]", "(18,20]", 
  #                   "(22,24]", 
  #                   "(26,28]", 
  #                   "(30,32]", 
  #                   "(34,36]", 
  #                   "(38,42]", 
  #                   "PP")
  # ) %>% 
ggplot(aes(x = value, y = ga_range)) +
  # geom_vline(xintercept = 0) +
  labs(x = "", y = "") +
  geom_density_ridges_gradient(aes(fill = ga_range), 
                               scale = 3, size = 0.1,
                               show.legend = FALSE) +  
  theme(legend.position = "none") +
  theme_bw() +
  scale_fill_manual(values = text_colour) 
  # facet_grid(cols = vars(path_name), rows = NULL, scales = "free")
  # 
  # facet_wrap(facets = vars(path_name)
  #            # nrow = 1, 
  #            # scales = "free_x"
  #            )
plot

ggsave(
  plot,
  file = file.path("marker_in_different_points", "pathway_disease_ridges.pdf"),
  width = 6,
  height = 10,
  bg = "transparent"
)










    