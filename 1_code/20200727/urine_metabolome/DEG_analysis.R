##avoid source 
no_function()

setwd(r4projects::get_project_wd())
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
#   readxl::read_xlsx("/Users/shenxt/projects/smartD/patient_information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
# info <-
#   info %>%
#   mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
#   mutate(ID = paste("SF", ID, sep = ""))
# 
# sample_info <-
#   readr::read_csv("/Users/shenxt/projects/smartD/patient_information/sample_info_191021.csv")

# setwd(r4projects::get_project_wd())
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

text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13),
  "red"
)

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

sample_info %>% 
  dplyr::group_by(ga_range) %>% 
  dplyr::summarise(n = n())

library(plyr)
temp_data <- 
sample_info  %>% 
  plyr::dlply(.(ga_range)) %>% 
  lapply(function(x){
    c(length(unique(x$subject_id)), length(unique(x$sample_id)))
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

library(ggradar)

colnames(temp_data) <- c("subject", "sample")

plot <- 
temp_data %>% 
  tibble::rownames_to_column(var = "ga_range") %>% 
  tidyr::pivot_longer(cols = -ga_range, 
                      names_to = "class", values_to = "number") %>% 
  ggplot(aes(ga_range, number)) +
  geom_line(aes(group = class, 
                color = class),
           size = 1) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point(aes(group = class, fill = class),
             size = 5, shape = 21) +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm(alpha = 0.8) +
  scale_y_continuous(limits = c(5, 30)) +
  theme_bw() +
  labs(x = "", y = "#Subject/Sample") +
  theme(axis.text.x = element_text(size = 12, color = text_colour),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_blank()) +
  coord_polar(start = -0.2) 

# ggsave(plot, filename = "subject_sample_number_in_ga_range.pdf", 
#        width = 8, height = 7)

  


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


##increase
temp_data_increase <- 
lapply(marker_each_point, function(x){
  if(is.null(x)){
   x <- rep(NA, length(marker_each_point))
   names(x) <- names(marker_each_point)
   return(x)
  }
  peak_name1 <- 
    x %>% 
    dplyr::filter(class == "increase") %>% 
    dplyr::pull(gene_id)
  
  if(length(peak_name1) == 0){
    x <- rep(NA, length(marker_each_point))
    names(x) <- names(marker_each_point)
    return(x)
  }
  
  lapply(marker_each_point, function(y){
    if(is.null(y)){
      return(NA)
    }
    
    peak_name2 <- 
      y %>% 
      dplyr::filter(class == "increase") %>% 
      dplyr::pull(gene_id)
    
    if(length(peak_name2) == 0){
      return(NA)
    }
    
    length(intersect(peak_name1, peak_name2))/length(peak_name1)
    
  }) %>% unlist()
  
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

temp_data_increase <- temp_data_increase[-13,-13]





##increase
temp_data_decrease <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    peak_name1 <- 
      x %>% 
      dplyr::filter(class == "decrease") %>% 
      dplyr::pull(gene_id)
    
    if(length(peak_name1) == 0){
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    
    lapply(marker_each_point, function(y){
      if(is.null(y)){
        return(NA)
      }
      
      peak_name2 <- 
        y %>% 
        dplyr::filter(class == "decrease") %>% 
        dplyr::pull(gene_id)
      
      if(length(peak_name2) == 0){
        return(NA)
      }
      
      length(intersect(peak_name1, peak_name2))/length(peak_name1)
      
    }) %>% unlist()
    
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()





temp_data <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    peak_name1_increase <- 
      x %>% 
      dplyr::filter(class == "increase") %>% 
      dplyr::pull(gene_id)
    
    peak_name1_decrease <- 
      x %>% 
      dplyr::filter(class == "decrease") %>% 
      dplyr::pull(gene_id)
    
    if(length(peak_name1_decrease) == 0 & length(peak_name1_increase) == 0){
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    
    lapply(marker_each_point, function(y){
      if(is.null(y)){
        return(NA)
      }
      
      peak_name2_decrease <- 
        y %>% 
        dplyr::filter(class == "decrease") %>% 
        dplyr::pull(gene_id)
      
      peak_name2_increase <- 
        y %>% 
        dplyr::filter(class == "increase") %>% 
        dplyr::pull(gene_id)
      
      if(length(peak_name2_decrease) == 0 & length(peak_name2_increase) == 0){
        return(NA)
      }
      
      (length(intersect(peak_name1_decrease, peak_name2_decrease)) + 
        length(intersect(peak_name1_increase, peak_name2_increase)))/(length(peak_name1_decrease) + length(peak_name1_increase))
      
    }) %>% unlist()
    
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()


temp_data_increase <- temp_data_increase[-13,-13]

temp_data_decrease <- temp_data_decrease[-13,-13]

temp_data <- temp_data[-13,-13]



# temp_data_increase[which(is.na(temp_data_increase), arr.ind = TRUE)] <- 0

library(corrplot)

corrplot(corr = as.matrix(temp_data_increase),
         method = "circle", type = "upper", is.corr = FALSE, tl.col = text_colour,
         diag = FALSE,
         na.label = "NA",
         addCoef.col = "grey",
         col =   colorRampPalette(colors = c(
           ggsci::pal_aaas()(10)[5],
           "white",
           ggsci::pal_aaas()(10)[6]
         ))(100)
         )



corrplot(corr = as.matrix(temp_data_decrease),
         method = "circle", type = "upper", is.corr = FALSE, tl.col = text_colour,
         diag = FALSE,
         na.label = "NA",
         addCoef.col = "grey",
         col =   colorRampPalette(colors = c(
           ggsci::pal_aaas()(10)[6],
           "white",
           ggsci::pal_aaas()(10)[5]
         ))(100)
)


corrplot(corr = as.matrix(temp_data),
         method = "circle", type = "upper", is.corr = FALSE, tl.col = text_colour,
         diag = FALSE,
         na.label = "NA",
         addCoef.col = "grey"
)










# 
# 
# temp_data <- 
#   lapply(1:length(marker_each_point), function(i){
#     x <- marker_each_point[[i]]
#     if(is.null(x)){
#       return(NA)
#     }
#     peak_name1_increase <- 
#       x %>% 
#       dplyr::filter(class == "increase") %>% 
#       dplyr::pull(gene_id)
#     
#     peak_name1_decrease <- 
#       x %>% 
#       dplyr::filter(class == "decrease") %>% 
#       dplyr::pull(gene_id)
#     
#     if(length(peak_name1_decrease) == 0 & length(peak_name1_increase) == 0){
#       return(NA)
#     }
#     
#     peak_name2_decrease <- 
#     lapply(marker_each_point[(i+1):(length(marker_each_point) - 1)], function(y){
#       if(is.null(y)){
#         return(NULL)
#       }
#         y %>% 
#         dplyr::filter(class == "decrease") %>% 
#         dplyr::pull(gene_id)
#     }) %>% unlist() %>% 
#       unname()
#     
#     peak_name2_increase <- 
#       lapply(marker_each_point[(i+1):(length(marker_each_point) - 1)], function(y){
#         if(is.null(y)){
#           return(NULL)
#         }
#         y %>% 
#           dplyr::filter(class == "increase") %>% 
#           dplyr::pull(gene_id)
#       }) %>% unlist() %>% 
#       unname()
#     
#     decrease <- 
#     lapply(peak_name1_decrease, function(z){
#       sum(z == peak_name2_decrease)
#     }) %>% 
#       unlist()
#     
#     increase <- 
#     lapply(peak_name1_increase, function(z){
#       sum(z == peak_name2_increase)
#     }) %>% 
#       unlist()
#     
#     all <- c(decrease, increase)
#     all >= length((i+1):(length(marker_each_point) - 1)) - 1
#   })
# 
# temp_data <- unlist(temp_data)
# temp_data <- temp_data[!is.na(temp_data)]
# sum(unlist(temp_data), na.rm = TRUE)/length(temp_data)





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

library(ggalluvial)

temp_data <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      return(NULL)
    }
    x <- 
      data.frame(gene_id = all_marker_name,
                 stringsAsFactors = FALSE) %>% 
      left_join(x, by = "gene_id") %>% 
      dplyr::select(gene_id, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
    
  })

temp_data <-
  purrr::map2(.x = temp_data, .y = names(temp_data), .f = function(x,y){
    if(is.null(x)){
      return(NULL)
    }
    data.frame(x, point = y, stringsAsFactors = FALSE)
  })

temp_data <- 
  do.call(rbind, temp_data)


temp_data %>% 
  dplyr::group_by(point, class) %>% 
  dplyr::summarise(n = n())


plot1 <- 
  ggplot(temp_data,
         aes(x = point, 
             y = freq,
             stratum = class, 
             alluvium = gene_id,
             fill = class, 
             label = class)) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c(
    "increase" = ggsci::pal_aaas()(10)[6],
    "decrease" = ggsci::pal_aaas()(10)[5],
    "no" = "grey"
  )) +
  ggalluvial::geom_stratum(alpha = 1) +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot1

# ggsave(
#   plot1,
#   file = file.path("marker_in_different_points", "gene_sankey_light.pdf"),
#   width = 14,
#   height = 7,
#   bg = "transparent"
# )

# ggsave(
#   plot1,
#   file = file.path("marker_in_different_points", "gene_sankey_light2.pdf"),
#   width = 20,
#   height = 7,
#   bg = "transparent"
# )


plot2 <- 
  marker_each_point %>% 
  lapply(function(x){
    if(is.null(x)){
      c(0, 0)
    }else{
      x$class %>% table() %>% as.numeric() 
    }
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "point") %>% 
  mutate(decrease = -V1, increase = V2) %>% 
  mutate(point = factor(point, levels = point)) %>% 
  ggplot(aes(point)) +
  geom_bar(aes(point, increase), stat = "identity", fill = ggsci::pal_aaas()(10)[6]) +
  geom_bar(aes(point, decrease), stat = "identity", fill = ggsci::pal_aaas()(10)[5]) +
  theme_bw() +
  labs(x = "", y = "Peak number") +
  theme(legend.position = "top", 
        # panel.border = element_blank(), 
        # panel.grid = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.text.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot2

# ggsave(
#   plot2,
#   file = file.path("marker_in_different_points", "barplot.pdf"),
#   width = 14,
#   height = 7,
#   bg = "transparent"
# )


##heatmap
all_marker_name
subject_data_mean
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c(ggsci::pal_aaas()(10)[5], "white", ggsci::pal_aaas()(10)[6]))
temp_data <- 
  subject_data_mean[all_marker_name,]

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

range(temp_data)

plot3 <- 
Heatmap(temp_data, 
        cluster_columns = FALSE, 
        cluster_rows = TRUE, 
        show_row_names = FALSE, 
        show_column_names = FALSE,
        border = TRUE, 
        col = col_fun,
        row_km = 2)

plot3

library(ggplotify)

plot3 <- as.ggplot(plot3)

ggsave(
  plot3,
  file = file.path("marker_in_different_points", "heatmap.pdf"),
  width = 14,
  height = 7,
  bg = "transparent"
)



###upset plot show overlap of markers
marker_each_point_up <-
  purrr::map(
    .x = marker_each_point,
    # .y = names(marker_each_point),
    .f = function(x) {
      if(is.null(x)){
        return(NULL)
      }
      x <- 
        x %>% 
        dplyr::filter(class == "increase")
      
      if(nrow(x) == 0){
        return(NULL)
      }
      x$gene_id
      # data.frame(name = x$gene_id, time = y, stringsAsFactors = FALSE)
    }
  ) 
  # do.call(rbind, .)

remove_idx <- lapply(marker_each_point_up, is.null) %>% 
  unlist() %>% 
  which()

if(length(remove_idx) > 0){
  marker_each_point_up <- marker_each_point_up[-remove_idx]
}


marker_each_point_down <-
  purrr::map(
    .x = marker_each_point,
    # .y = names(marker_each_point),
    .f = function(x) {
      if(is.null(x)){
        return(NULL)
      }
      x <- 
        x %>% 
        dplyr::filter(class == "decrease")
      
      if(nrow(x) == 0){
        return(NULL)
      }
      x$gene_id
      # data.frame(name = x$gene_id, time = y, stringsAsFactors = FALSE)
    }
  ) 
# do.call(rbind, .)

remove_idx <- lapply(marker_each_point_down, is.null) %>% 
  unlist() %>% 
  which()

if(length(remove_idx) > 0){
  marker_each_point_down <- marker_each_point_down[-remove_idx]
}

library(UpSetR)
library(ComplexHeatmap)

text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13),
  "red"
)

m = make_comb_mat(marker_each_point_up)

plot <- 
UpSet(m, set_order = names(marker_each_point_up), 
      # comb_order = order(comb_size(m)),
      bg_col = c(ggsci::pal_aaas(alpha = 0.5)(10)[6], "white"),
      right_annotation = upset_right_annotation(m, 
                                                # ylim = c(0, 30),
                                                gp = gpar(fill = text_colour[4:14]),
                                                annotation_name_side = "bottom",
                                                axis_param = list(side = "bottom"))
      )

plot <- ggplotify::as.ggplot(plot = plot)

ggsave(plot, filename = "up_marker_upset_plot.pdf", width = 14, height = 10)



m = make_comb_mat(marker_each_point_down)

plot <- 
  UpSet(m, set_order = names(marker_each_point_down), 
        # comb_order = order(comb_size(m)),
        bg_col = c(ggsci::pal_aaas(alpha = 0.5)(10)[5], "white"),
        right_annotation = upset_right_annotation(m, 
                                                  # ylim = c(0, 30),
                                                  gp = gpar(fill = text_colour[3:14]),
                                                  annotation_name_side = "bottom",
                                                  axis_param = list(side = "bottom"))
  )

plot <- ggplotify::as.ggplot(plot = plot)
plot

ggsave(plot, filename = "down_marker_upset_plot.pdf", width = 14, height = 10)


# temp_data <- 
#   expression_data
# 
# colnames(temp_data) == sample_info$sample_id
# 
# temp_data <- temp_data[all_marker_name,]
# 
# temp_data <- 
#   log(temp_data + 1, 10) %>% 
#   apply(1, function(x){
#     (x - mean(x))/sd(x)
#   }) %>% 
#   t() %>% 
#   as.data.frame()
# 
# range(temp_data)
# 
# temp_data[temp_data > 8.13] <- 8.13
# 
# Heatmap(temp_data, 
#         cluster_columns = FALSE, 
#         cluster_rows = TRUE, 
#         show_row_names = FALSE, 
#         show_column_names = FALSE,
#         border = FALSE, 
#         col = col_fun)


####PIUMet analysis
# purrr::walk2(
#   .x = fc_p_value,
#   .y = names(fc_p_value),
#   .f = function(x, y) {
#     idx <- which(x$p_value < 0.05)
#     x <- 
#     cbind(variable_info, x)[idx,] %>% 
#       as.data.frame() %>% 
#       dplyr::mutate(polarity = 
#                       case_when(
#                         stringr::str_detect(name, "POS") ~ "positive",
#                         stringr::str_detect(name, "NEG") ~ "negative"
#                       )) %>% 
#       dplyr::select(mz, polarity, p_value) %>% 
#       dplyr::mutate(p_value = -log(p_value, 10))
#     dir.create(file.path("marker_in_different_points", y))
#     write.csv(x, file.path(file.path("marker_in_different_points", 
#                                      y, "marker.csv")))
#     
#     write.table(x, file.path(file.path("marker_in_different_points", 
#                                      y, "marker.txt")), sep = "\t",
#                 quote = FALSE, row.names = FALSE, col.names = FALSE)
#     
#   }
# )

###20_22
# readPIUMet(path = "marker_in_different_points/(20,22]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(20,22]`)
# 
# ###22_24
# readPIUMet(path = "marker_in_different_points/(22,24]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(22,24]`)
# 
# ###28_30
# readPIUMet(path = "marker_in_different_points/(24,26]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(24,26]`)
# 
# ###26_28
# readPIUMet(path = "marker_in_different_points/(26,28]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(26,28]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# ###28_30
# readPIUMet(path = "marker_in_different_points/(28,30]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(28,30]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# 
# ###30_32
# readPIUMet(path = "marker_in_different_points/(30,32]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(30,32]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# ###32_34
# readPIUMet(path = "marker_in_different_points/(32,34]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(32,34]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# ###34_36
# readPIUMet(path = "marker_in_different_points/(34,36]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(34,36]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# 
# ###36_38
# readPIUMet(path = "marker_in_different_points/(36,38]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(36,38]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# ###38_42
# readPIUMet(path = "marker_in_different_points/(38,42]/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`(38,42]`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))
# 
# ###PP
# readPIUMet(path = "marker_in_different_points/PP/piumet_output", 
#            variable_info = variable_info, 
#            fc_p_table = fc_p_value$`PP`, 
#            text = FALSE, 
#            layout = "kk", 
#            size_range = c(1,3),
#            width_range = c(0.2,0.5))


##pathway enrichment for each data point
# #22_24
# load("marker_in_different_points/(22,24]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(22,24]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(22,24]/path_result")
# 
# 
# ##28_30
# load("marker_in_different_points/(24,26]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(24,26]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(24,26]/path_result")
# 
# 
# 
# ##26_28
# load("marker_in_different_points/(26,28]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(26,28]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(26,28]/path_result")
# 
# 
# 
# ##28_30
# load("marker_in_different_points/(28,30]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(28,30]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(28,30]/path_result")
# 
# 
# 
# 
# 
# 
# ##30_32
# load("marker_in_different_points/(30,32]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(30,32]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(30,32]/path_result")
# 
# 
# ##32_34
# load("marker_in_different_points/(32,34]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(32,34]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(32,34]/path_result")
# 
# 
# 
# 
# 
# ##34_36
# load("marker_in_different_points/(34,36]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(34,36]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(34,36]/path_result")
# 
# 
# 
# ##36_38
# load("marker_in_different_points/(36,38]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(36,38]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(36,38]/path_result")
# 
# 
# 
# 
# 
# ##38_42
# load("marker_in_different_points/(38,42]/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/(38,42]/annotation_result")
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
# save(path_result, file = "marker_in_different_points/(38,42]/path_result")
# 
# 
# 
# 
# 
# 
# 
# ##PP
# load("marker_in_different_points/PP/piumet_output/Result/node_data")
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
# save(annotation_result, file = "marker_in_different_points/PP/annotation_result")
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
# save(path_result, file = "marker_in_different_points/PP/path_result")
# 



##-----------------------------------------------------------------------------
####heatmap for pathway in different data points
load("marker_in_different_points/(22,24]/path_result")
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

path_result_22_24 <- path_result

#####
load("marker_in_different_points/(24,26]/path_result")
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

path_result_24_26 <- path_result


##----
load("marker_in_different_points/(26,28]/path_result")
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

path_result_26_28 <- path_result

load("marker_in_different_points/(28,30]/path_result")
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
path_result_28_30 <- path_result

load("marker_in_different_points/(30,32]/path_result")
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
path_result_30_32 <- path_result

load("marker_in_different_points/(32,34]/path_result")
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
path_result_32_34 <- path_result

load("marker_in_different_points/(34,36]/path_result")
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
path_result_34_36 <- path_result

load("marker_in_different_points/(36,38]/path_result")
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
path_result_36_38 <- path_result

load("marker_in_different_points/(38,42]/path_result")
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
path_result_38_42 <- path_result

load("marker_in_different_points/PP/path_result")
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
# save(all_path_name, file = "marker_in_different_points/all_path_name")
# save(all_path_id, file = "marker_in_different_points/all_path_id")

load("marker_in_different_points/all_path_name")
load("marker_in_different_points/all_path_id")

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
# save(annotation_result, file = "marker_in_different_points/annotation_result")

load("marker_in_different_points/annotation_result")

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

load("hsa_pathway")

###heatmap for each metabolite and pathway
data_for_each_path <- 
purrr::map2(.x = all_path_id, .y = all_path_name,.f = function(x, y){
  temp_idx <- grep(x, names(hsa_pathway))
  temp_idx <- which(rownames(temp_data) %in% hsa_pathway[[temp_idx]])
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
all_path_name[1:6]
Heatmap(as.matrix(data_for_each_path[[5]][,-1]), cluster_columns = FALSE, cluster_rows = FALSE)

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
    path_result_22_24[,c("Pathway.name", "Overlap", "p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(22,24]"),
    
    # path_result_24_26[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    #   data.frame(., ga = "(24,26]")
    
    path_result_26_28[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(26,28]"),
    
    path_result_28_30[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(28,30]"),
    
    path_result_30_32[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(30,32]"),
    
    path_result_32_34[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(32,34]"),
    
    path_result_34_36[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(34,36]"),
    
    path_result_36_38[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(36,38]"),
    
    path_result_38_42[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
      dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
      data.frame(., ga = "(38,42]"),
    
    path_result_pp[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
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

text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13),
  "red"
)

plot <- 
  Heatmap(matrix = t(as.matrix(data_for_each_path2[,-c(1,2)])),
          cluster_row_slices = FALSE, 
          cluster_column_slices = FALSE, 
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    top_annotation = columnAnnotation(
      pathway = factor(data_for_each_path2[,2], 
                       levels = c("Steroid hormone biosynthesis", 
                                  "Ovarian steroidogenesis", 
                                  "Cortisol synthesis and secretion",
                                  "Prolactin signaling pathway",
                                  "Aldosterone synthesis and secretion",
                                  "Phenylalanine, tyrosine and tryptophan biosynthesis")
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
    rect_gp = gpar(col= "white"), 
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
purrr::map(all_path_id, .f = function(x){
  temp_idx <- grep(x, names(hsa_pathway))
  temp_idx <- which(rownames(temp_data) %in% hsa_pathway[[temp_idx]])
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
    temp_idx <- grep(x, names(hsa_pathway))
    temp_idx <- which(rownames(temp_data2) %in% hsa_pathway[[temp_idx]])
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
  path_result_22_24[,c("Pathway.name", "Overlap", "p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(22,24]"),
  
  # path_result_24_26[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
  #   dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
  #   data.frame(., ga = "(24,26]")
  
  path_result_26_28[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(26,28]"),
  
  path_result_28_30[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(28,30]"),
  
  path_result_30_32[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(30,32]"),
  
  path_result_32_34[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(32,34]"),
  
  path_result_34_36[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(34,36]"),
  
  path_result_36_38[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(36,38]"),
  
  path_result_38_42[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
    dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05) %>% 
    data.frame(., ga = "(38,42]"),
  
  path_result_pp[,c("Pathway.name", "Overlap","p.value.fdr")] %>% 
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

# ggsave(
#   plot,
#   file = file.path("marker_in_different_points", "pathway_number.pdf"),
#   width = 5,
#   height = 7,
#   bg = "transparent"
# )

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
path_data$value[path_data$value < -1.95] <- -1.95

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

# ggsave(
#   plot,
#   file = file.path("marker_in_different_points", "pathway.pdf"),
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
      subject_data[,peak_name ,drop = FALSE] %>% 
        apply(1, mean)
    }
  ) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

load("hsa_pathway")

data_for_each_path <- 
  purrr::map2(.x = all_path_id[1:6], 
              .y = all_path_name[1:6],.f = function(x, y){
    temp_idx <- grep(x, names(hsa_pathway))
    temp_idx <- which(rownames(temp_data) %in% hsa_pathway[[temp_idx]])
    temp_data[temp_idx,,drop = FALSE] %>% 
      apply(1, function(x){
        (x - mean(x))/sd(x)
      }) %>% 
      t() %>% 
      as.data.frame() %>% 
      data.frame(path_name = y, ., stringsAsFactors = FALSE , check.names = FALSE)
  }
  )

names(data_for_each_path) <- all_path_name[1:6]

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
                                     c("Steroid hormone biosynthesis", 
                                       "Ovarian steroidogenesis", 
                                       "Cortisol synthesis and secretion",
                                       "Prolactin signaling pathway",
                                       "Aldosterone synthesis and secretion",
                                       "Phenylalanine, tyrosine and tryptophan biosynthesis"))) %>% 
  dplyr::filter(
    ga_range %in% c("(10,16]", "(18,20]", 
                    "(22,24]", 
                    "(26,28]", 
                    "(30,32]", 
                    "(34,36]", 
                    "(38,42]", 
                    "PP")
  ) %>% 
ggplot(aes(x = value, y = ga_range)) +
  # geom_vline(xintercept = 0) +
  labs(x = "", y = "") +
  geom_density_ridges_gradient(aes(fill = ga_range), 
                               scale = 3, size = 0.1,
                               show.legend = FALSE) +  
  theme(legend.position = "none") +
  theme_bw() +
  scale_fill_manual(values = text_colour) +
  facet_wrap(facets = vars(path_name), nrow = 3, scales = "free_x")

plot

ggsave(
  plot,
  file = file.path("marker_in_different_points", "pathway_ridges.pdf"),
  width = 6,
  height = 10,
  bg = "transparent"
)



######correlation network to find the module
dim(annotation_result)

##remove duplicated KEGG OD
annotation_result1 <- 
  annotation_result %>% 
  dplyr::distinct(KEGG_ID,.keep_all = TRUE)

##from peak to metabolite
temp_subject_data <- 
  t(subject_data) %>% 
  as.data.frame()

temp_data <- 
  annotation_result1$name %>% 
  purrr::map(
    .f = function(x){
      x <- stringr::str_split(x, ";")[[1]]
      temp_subject_data[x,,drop = FALSE] %>% 
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
# save(cor_value, file = "marker_in_different_points/cor_value")
load("marker_in_different_points/cor_value")

cor_value <- cor_value %>% 
  dplyr::filter(abs(r) > 0.5)

edge_data <- 
  cor_value %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE)

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
    data.frame(index = 1:length(subnetworks$modularity),
               modu = subnetworks$modularity, stringsAsFactors = FALSE),
    aes(index, modu) 
  ) +
  geom_vline(xintercept = which.max(subnetworks$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line(colour = "black") +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

plot <-
  plot + 
  ggplot2::annotate(geom = "point", 
                    x = which.max(subnetworks$modularity),
                    y = max(subnetworks$modularity), 
                    size = 3, 
                    colour = "#FFA319FF") +
  annotate(geom = "text", 
           x = which.max(subnetworks$modularity),
           y = max(subnetworks$modularity), 
           label = paste("(",  which.max(subnetworks$modularity),
                         ",", 
                         max(subnetworks$modularity) %>% round(3),
                         ")"),
           size = 5,
           colour = "#FFA319FF"
  )

plot

# ggsave(plot, filename = "modularity.pdf", width = 7, height = 7)
# ggsave(plot, filename = "modularity.png", width = 7, height = 7)

table(membership(communities = subnetworks))
which(table(membership(communities = subnetworks)) >= 5)
##cluster 1,2,3,8,10,15 
node <- vertex_attr(graph = graph, name = "node")



membership <- 
  data.frame(node, membership = as.numeric(membership(communities = subnetworks)), 
             stringsAsFactors = FALSE) %>% 
  # dplyr::arrange(membership) %>% 
  # dplyr::filter(membership %in% c(1,2,3,4,5,7,8)) %>% 
  dplyr::mutate(membership = paste("Cluster", membership, sep = "")) %>% 
  dplyr::mutate(membership = case_when(
    membership %in% c(paste("Cluster", c(1,2,3,8,10,15), 
                            sep = "")) ~ membership,
    TRUE ~ "Other"
  ))

graph <- 
  igraph::set_vertex_attr(graph = graph, name = "membership", 
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
  ggraph(graph,
         layout = 'fr') +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE
  ) +
  geom_node_point(aes(size = Degree,
                      fill = membership),
                  alpha = 1, 
                  shape = 21,
                  show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(low = "#3B4992FF", mid = "white",
                                     high = "#EE0000FF", midpoint = 0) +
  ggraph::scale_edge_width(range = c(0.3,1)) +
  # geom_node_text(aes(label = Metabolite.Name), repel = TRUE) +
  scale_size_continuous(range = c(1.5, 5)) +
  # ggsci::scale_fill_futurama(alpha = 0.7) +
  scale_fill_manual(values = cluster_color) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))
# coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-1.4,1.4))

plot

# extrafont::font_import()
extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = "marker_in_different_points/cor_network_all.pdf",
#   width = 8,
#   height = 7,
#   bg = "transparent"
# )



###from metabolite to module
membership1 <- 
  membership %>% 
  dplyr::filter(membership != "Other")

temp_data <- 
  lapply(unique(membership1$membership), function(x){
    x <- membership1$node[membership1$membership == x]
    peak_name <- annotation_result1$name[match(x, annotation_result1$KEGG_ID)]
    peak_name <- stringr::str_split(peak_name, ";") %>%
      unlist() %>% 
      unique()
    subject_data_mean[peak_name,] %>% 
      apply(2, mean)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(temp_data) <- unique(membership1$membership)

temp_data <- 
  temp_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c(ggsci::pal_aaas()(10)[5], "white", ggsci::pal_aaas()(10)[6]))

range(temp_data)

temp_data[temp_data < -2.16] <- -2.16

plot <- 
  Heatmap(temp_data, 
          cluster_columns = FALSE, 
          cluster_rows = TRUE, 
          show_row_names = TRUE, 
          show_column_names = TRUE,
          border = TRUE, 
          col = col_fun)

plot

library(ggplotify)

plot <- as.ggplot(plot)


ggsave(
  plot,
  file = "marker_in_different_points//module_heatmap.pdf",
  width = 7,
  height = 7,
  bg = "transparent"
)

####
##module quantative
ga <- sample_info$GA

ga_range <- 
  as.character(cut_width(x = ga, width = 2, center = 3)) %>%
  stringr::str_replace("\\[", '(')

table(ga_range)

ga_range[ga_range == "(10,12]"] <- "(10,16]"
ga_range[ga_range == "(12,14]"] <- "(10,16]"
ga_range[ga_range == "(14,16]"] <- "(10,16]"

ga_range[ga_range == "(38,40]"] <- "(38,42]"
ga_range[ga_range == "(40,42]"] <- "(38,42]"

ga_range[ga_range == "(44,46]"] <- "PP"


table(ga_range)

ga_level <- sort(unique(ga_range))

temp_data <- 
  lapply(unique(membership1$membership), function(x){
    x <- membership1$node[membership1$membership == x]
    peak_name <- annotation_result1$name[match(x, annotation_result1$KEGG_ID)]
  
    temp <-   
    lapply(peak_name, function(y){
      y <- stringr::str_split(y, ";")[[1]]
      apply(subject_data_mean[y,,drop = FALSE], 2, mean)
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
    x <- data.frame(x, cluster = y, 
                    stringsAsFactors = FALSE, 
                    check.names = FALSE, check.rows = FALSE) %>% 
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
  temp_data[,-c(1:2)] %>% 
  apply(1, function(x){
    (x-mean(x))/sd(x)
  }) %>% 
  t()

plot <- 
  Heatmap(matrix = as.matrix(temp_data2),
          cluster_row_slices = FALSE, 
          cluster_column_slices = FALSE, 
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = FALSE,
          show_column_names = TRUE,
          left_annotation = rowAnnotation(
            cluster = factor(temp_data[,2],
                             levels = unique(temp_data$cluster)
            ),
            col = list(cluster = cluster_color
            )
          ),
          # row_title_rot = 0
          border = TRUE,
          gap = TRUE, 
          col = col_fun,
          rect_gp = gpar(col= "white"), 
          column_names_rot = 45, 
          column_names_gp = gpar(col = text_colour)
  )

plot

plot <- as.ggplot(plot)

# ggsave(plot, 
#        filename = "marker_in_different_points/heatmap_for_each_cluster_metabolite_level.pdf", 
#        width = 10, height = 7)

temp_data <- 
  temp_data %>% 
  tidyr::pivot_longer(cols = -c(cluster, KEGG_ID), 
                      names_to = "ga_range", 
                      values_to = "value")

plot <- 
  temp_data %>%
  dplyr::filter(cluster %in% c(
     "Cluster2","Cluster3")) %>% 
  group_by(ga_range,
           cluster) %>%
  dplyr::mutate(mean = mean(value),
                sem = sd(value) / sqrt(nrow(temp_data))) %>%
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
    width = 0, show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_color) +
  labs(x = "", y = "PIUMet module intensity") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 12, 
                                   color = text_colour, 
                                   angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12))  +
  facet_wrap(~cluster, ncol = 1, scales = "free_y")

plot

ggsave(plot, filename = "marker_in_different_points/module_in_different_point.pdf", 
       width = 7, height = 7)


#####grah cluster 2 and 3
membership <- igraph::vertex_attr(graph = graph)$membership
node <- igraph::vertex_attr(graph = graph)$node



igraph::induced_subgraph(graph = graph,
                         v = which(membership %in% c("Cluster2")))


graph2 <- igraph::induced_subgraph(graph = graph,
                           v = which(membership %in% c("Cluster2")))

node2 <-  igraph::vertex_attr(graph = graph2, name = "node")

graph3 <- igraph::induced_subgraph(graph = graph,
                                   v = which(membership %in% c("Cluster3")))

node3 <-  igraph::vertex_attr(graph = graph3, name = "node")

##class of node
load("hmdbMS1Database0.0.1")
hmdb_data <- hmdbMS1Database0.0.1@spectra.info

node_class2 <- hmdb_data[match(node2, hmdb_data$KEGG.ID),"Super.class"]

node_pathway2 <- 
lapply(node2, function(x){
  lapply(hsa_pathway, function(y){
    x %in% y
  }) %>% 
    unlist() %>% 
    which() %>% 
    names()
})

graph2 <- 
igraph::set_vertex_attr(graph = graph2, name = "class", 
                        value = node_class2)


node_class3 <- hmdb_data[match(node3, hmdb_data$KEGG.ID),"Super.class"]

node_pathway3 <- 
  lapply(node3, function(x){
    lapply(hsa_pathway, function(y){
      x %in% y
    }) %>% 
      unlist() %>% 
      which() %>% 
      names()
  })

graph3 <- 
  igraph::set_vertex_attr(graph = graph3, name = "class", 
                          value = node_class3) 


plot2 <-
  ggraph(graph2,
         layout = 'kk') +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(
    aes(size = Degree,
        fill = class),
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

ggsave(
  plot2,
  filename = "marker_in_different_points/cluster_plot2.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)





plot3 <-
  ggraph(graph3,
         layout = 'kk') +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(
    aes(size = Degree,
        fill = class),
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

ggsave(
  plot3,
  filename = "marker_in_different_points/cluster_plot3.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)

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
  purrr::map(.f = function(x){
    x[1]
  }) %>% 
  unlist() %>% 
  stringr::str_replace(" - Homo sapiens \\(human\\)", "")

node_pathway2 <- 
node_pathway2 %>% 
  dplyr::arrange(desc(frequency))

pathway2 <- 
node_pathway2 %>% 
  dplyr::filter(frequency >=3) %>% 
  dplyr::mutate(pathway = factor(pathway, levels = rev(pathway))) %>% 
  ggplot(aes(frequency, pathway)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0,0))) +
  theme_classic() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))


pathway2


ggsave(pathway2, filename = "marker_in_different_points/pathway2.pdf",
       width = 7, height = 7)





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
  purrr::map(.f = function(x){
    x[1]
  }) %>% 
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
  scale_x_continuous(expand = expansion(mult = c(0,0))) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))


pathway3


ggsave(pathway3, filename = "marker_in_different_points/pathway3.pdf",
       width = 7, height = 7)








    
