#to avoind source
no_exist_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list=ls())
source("R/20200727/tools.R")

#get the identification results
#####for metabolite identification result
##load data
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/urine_metabolome/data_cleaning/RPLC/POS/")
load("rplc_pos_6")

setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/urine_metabolome/data_cleaning/RPLC/NEG/")
load("rplc_neg_6")

##get the identification results
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/urine_metabolome/metabolite_annotation/RPLC/POS/")
identification_table_pos <- readr::read_csv("identification.table.new.csv")

setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/urine_metabolome/metabolite_annotation/RPLC/NEG/")
identification_table_neg <- readr::read_csv("identification.table.new.csv")

##first set the work directory to project folder
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/")

identification_table_pos <- 
  identification_table_pos %>% 
  select(name:rt, MS2.spectrum.name:Database)


identification_table_neg <- 
  identification_table_neg %>% 
  select(name:rt, MS2.spectrum.name:Database)

metabolite_table <- 
  rbind(identification_table_pos, identification_table_neg)


ms1_pos <- rplc_pos_6@ms1.data[[1]][,c(1,2,3)]
ms1_neg <- rplc_neg_6@ms1.data[[1]][,c(1,2,3)]

ms1 <- rbind(ms1_pos, ms1_neg)

ms1 <- 
  ms1 %>% 
  left_join(metabolite_table, by = "name")


ms1 <- 
  ms1 %>% 
  mutate(polarity = case_when(
    stringr::str_detect(name, "POS") ~ "POS",
    stringr::str_detect(name, "NEG") ~ "NEG",
    TRUE ~ "NA"
  ))


ms1 <-
  ms1 %>% 
  mutate(Level2 = 
           case_when(
             Level == 1 ~ "Yes",
             Level == 2 ~ "Yes",
             Level == 3 ~ "No",
             is.na(Level) ~ "Yes",
             TRUE ~ "NA"
           ))


ms1$Level[is.na(ms1$Level)] <- 4

ms1$Level[ms1$Level == 1] <- 2


require(webr)
require(moonBook)


plot <-
  PieDonut(
    ms1,
    aes(pies = polarity, donuts = Level),
    # explode = 1,
    # explodePie = FALSE,
    selected = c(1, 3, 4, 6),
    explodeDonut = TRUE,
    r0 = 0.4
    # color = "grey"
  )

ggsave(plot, filename = "annotation_level_distributation.pdf", width = 7, height = 7)


###we only use the level 1 and level 2 identifications
#1-3 means identification confidence level, 
#4 means MS2 spectra match si not good enough, and 5 means this compound
#may be not in human body
metabolite_table <- 
  metabolite_table %>% 
  filter(Level == 1 | Level == 2)

dim(metabolite_table)

####remove the duplicated metabolites according to score
metabolite_table <- 
  metabolite_table %>% 
  group_by(., Compound.name) %>% 
  filter(., SS == max(SS)) %>% 
  ungroup()

dim(metabolite_table)


metabolite_tags <- 
  metabolite_table %>% 
  select(name:Database)

# ####give new information for metabolites
# kegg.id <-
#   pbapply::pblapply(metabolite_tags$Compound.name, function(x){
#     metflow2::transID(query = x, from = "Chemical name", to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# hmdb.id <-
#   pbapply::pblapply(metabolite_tags$Compound.name, function(x){
#     metflow2::transID(query = x, from = "Chemical name",
#                                to = "Human Metabolome Database", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# ID_trans <-
#   data.frame(
#     KEGG = kegg.id$KEGG,
#     HMDB = hmdb.id$`Human Metabolome Database`,
#     stringsAsFactors = FALSE
#   )
# 
# hmdb_kegg <-
#   metabolite_tags %>%
#   select(KEGG.ID, HMDB.ID) %>%
#   data.frame(ID_trans, stringsAsFactors = FALSE)
# 
# hmdb_kegg <-
#   apply(hmdb_kegg, 1, function(x){
#     hmdb <- x[c(2,4)]
#     hmdb <- hmdb[!is.na(hmdb)]
#     kegg <- x[c(1,3)]
#     kegg <- kegg[!is.na(kegg)]
# 
#     hmdb <- ifelse(length(hmdb) == 0, NA, hmdb[1])
#     kegg <- ifelse(length(kegg) == 0, NA, kegg[1])
# 
#     c(hmdb, kegg)
# 
#   }) %>%
#   t() %>%
#   as_tibble()
# 
# 
# colnames(hmdb_kegg) <-
#   c("HMDB.ID", "KEGG.ID")
# 
# hmdb_kegg$HMDB.ID <-
#   hmdb_kegg$HMDB.ID %>%
#   sapply(function(x) {
#     if(is.na(x)) {
#       return(NA)
#     }
#     if(nchar(x) == 9){
#       pre <- stringr::str_extract(x, "[A-Za-z]{1,}")
#       end <- stringr::str_extract(x, "[0-9]{1,}")
#       end <- paste("00", end, sep = "")
#       x <- paste(pre, end, sep = "")
#     }
#     x
#   }
#   ) %>%
#   unname()
# 
# cbind(metabolite_tags$HMDB.ID, hmdb_kegg[,1])
# 
# metabolite_tags[,c("HMDB.ID", 'KEGG.ID')] <-
#   hmdb_kegg
# 
# 
# ###should use classyfire to calculate the class of all the metabolites
# inchikey1 <-
#   metabolite_tags$HMDB.ID %>%
#   pbapply::pblapply(function(x){
#     metflow2::transID(query = x, from = "Human Metabolome Database", to = "InChIKey", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# inchikey2 <-
#   metabolite_tags$KEGG.ID %>%
#   pbapply::pblapply(function(x){
#     metflow2::transID(query = x, from = "KEGG", to = "InChIKey", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# inchikey3 <-
#   metabolite_tags$name %>%
#   pbapply::pblapply(function(x){
#     metflow2::transID(query = x, from = "Chemical name", to = "InChIKey", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# inchikey <- cbind(inchikey1 = inchikey1$InChIKey,
#                   inchikey2 = inchikey2$InChIKey,
#                   inchikey3 = inchikey3$InChIKey) %>%
#   as_tibble() %>%
#   apply(1, function(x){
#     x <- x[!is.na(x)]
#     ifelse(length(x) == 0, NA, x[1])
#   })
# 
# metabolite_class <- vector(mode = "list", length = length(inchikey))
# 
# 
# metabolite_class <-
#   pbapply::pblapply(inchikey, function(x){
#     Sys.sleep(time = 5)
#     result <- metflow2::get_metclass(inchikey = x, sleep = 5)
#   })
# 
# 
# idx <-
#   lapply(metabolite_class, class) %>% unlist() %>%
#   `==`("logical") %>%
#   which()
# 
# inchikey[idx]
# 
# 
# super_class <- lapply(metabolite_class, function(x){
#   if(is.na(x)) return(NA)
#   x@classification_info %>%
#     dplyr::filter(name == "Superclass") %>%
#     pull(value)
# }) %>%
#   unlist()
# 
# 
# class <- lapply(metabolite_class, function(x){
#   if(is.na(x)) return(NA)
#   x@classification_info %>%
#     dplyr::filter(name == "Class") %>%
#     pull(value)
# }) %>%
#   unlist()
# 
# sub_class <- lapply(metabolite_class, function(x){
#   if(is.na(x)) return(NA)
#   x@classification_info %>%
#     dplyr::filter(name == "Subclass") %>%
#     pull(value)
# }) %>%
#   unlist()
# 
# sub_class[which(sub_class == "Not available")] <- NA
# 
# 
# metabolite_tags <- data.frame(metabolite_tags,
#                               super_class,
#                               class,
#                               sub_class,
#                               stringsAsFactors = FALSE)
# 
# ###get the metabolite table
# peak_table_pos <- rplc_pos_6@ms1.data[[1]]
# peak_table_neg <- rplc_neg_6@ms1.data[[1]]
# 
# ##remove Blank and QC_DL
# peak_table_pos <- 
#   peak_table_pos %>% 
#   select(-dplyr::starts_with('BLK')) %>% 
#   select(-dplyr::starts_with('blk'))
# 
# 
# peak_table_neg <- 
#   peak_table_neg %>% 
#   select(-dplyr::starts_with('BLK')) %>% 
#   select(-dplyr::starts_with('blk'))
# 
# 
# setdiff(colnames(peak_table_pos), colnames(peak_table_neg))
# 
# setdiff(colnames(peak_table_neg), colnames(peak_table_pos))
# 
# intersect_name <- 
#   intersect(colnames(peak_table_pos), colnames(peak_table_neg))
# 
# peak_table_pos <- 
#   peak_table_pos %>% 
#   select(one_of(intersect_name))
# 
# peak_table_neg <- 
#   peak_table_neg %>% 
#   select(one_of(intersect_name))
# 
# 
# colnames(peak_table_pos) == colnames(peak_table_neg)
# 
# peak_table <- 
#   rbind(peak_table_pos, peak_table_neg)
# 
# 
# readr::write_csv(peak_table, "peak_table.csv")
# save(peak_table, file = "peak_table")
# load("peak_table")
# 
# metabolite_table <- 
#   peak_table %>% 
#   filter(name %in% metabolite_tags$name) %>% 
#   arrange(name)
# 
# metabolite_tags <- 
#   metabolite_tags %>% 
#   filter(name %in% metabolite_table$name) %>% 
#   arrange(name)
# 
# metabolite_table$name == metabolite_tags$name
# 
# 
# metabolite_tags$Compound.name[grep("CE[0-9]{1,2}", x = metabolite_tags$Compound.name)] <-
# metabolite_tags$Compound.name %>% 
#   grep("CE[0-9]{1,2}", x = .) %>% 
#   `[`(metabolite_tags$Compound.name, .) %>% 
#   stringr::str_split(";") %>% 
#   lapply(function(x) x[[1]]) %>% 
#   unlist()
# 
# 
# metabolite_tags %>% 
#   filter(stringr::str_detect(Compound.name, "MMV"))
# 
# metabolite_tags <- 
#   metabolite_tags %>% 
#   filter(!stringr::str_detect(Compound.name, "MMV"))
# 
# 
# metabolite_tags %>% 
#   # filter(Database == "monaDatabase0.0.1") %>% 
#   pull(Compound.name) %>% 
#   View()
# 
# 
# metabolite_tags$Compound.name[792] <- 
#   "Lubiprostone"
# 
# 
# ##654 697 791 802 808 should be removed
# 
# metabolite_tags <- metabolite_tags[-c(654, 697, 791, 802, 808),]
# 
# 
# ##remove duplciated metabolites
# metabolite_tags <-
# metabolite_tags %>% 
#   group_by(Compound.name) %>% 
#   filter(SS == max(SS)) %>% 
#   ungroup()
#   
# 
# metabolite_table <- 
#   metabolite_table %>% 
#   dplyr::filter(name %in% metabolite_tags$name)
# 
# metabolite_table$name == metabolite_tags$name
# 
# 
# 
# 
# ###remove some isotopes
# ##creteria
# ##rt <5 s and cor > 0.8 and mz error meet the isotope distribution
# cor_value <- 
#   metabolite_table %>% 
#   select(-c(name, mz, rt)) %>% 
#   select(-contains("QC")) %>% 
#   as.data.frame() %>% 
#   t() %>% 
#   cor(method = "spearman")
# 
# colnames(cor_value) <-
#   rownames(cor_value) <-
#   metabolite_table$name
# 
# cor_value[lower.tri(cor_value)] <- 2
# 
# cor_value <- 
#   cor_value %>% 
#   as.data.frame() %>% 
#   rownames_to_column(., var = "name1") %>% 
#   tidyr::pivot_longer(cols = -name1, 
#                       names_to = "name2", values_to = "Correlation") %>% 
#   distinct() %>% 
#   filter(., Correlation != 1 & Correlation != 2) %>% 
#   arrange(., desc(abs(Correlation)))
# 
# 
# cor_value <- 
#   cor_value %>% 
#   filter(Correlation > 0.8)
# 
# rp_pos <- metID::rp.pos
# rp_neg <- metID::rp.neg
# 
# cor_value_pos <- 
#   cor_value %>% 
#   filter(stringr::str_detect(name1, "_POS") & stringr::str_detect(name2, "_POS"))
# 
# cor_value_neg <- 
#   cor_value %>% 
#   filter(stringr::str_detect(name1, "_NEG") & stringr::str_detect(name2, "_NEG"))
# 
# 
# result <- NULL
# 
# for(i in 1:nrow(cor_value_pos)){
#   x <- cor_value_pos[i,]
#   cat(i, " ")
#   x <- as.character(x)
#   
#   mz <- metabolite_tags %>% 
#     filter(name %in% x[1:2]) %>% 
#     pull(mz)
#   
#   rt <- metabolite_tags %>% 
#     filter(name %in% x[1:2]) %>% 
#     pull(rt)
#   
#   mz_diff <- abs(mz[1] - mz[2])
#   rt_diff <- abs(rt[1] - rt[2])
#   cor <- as.numeric(x[3])
#   
#   if(cor > 0.9 & rt_diff < 2 & mz_diff < 5){
#     result <- c(result, list(c(x, mz_diff, rt_diff)))
#   }
# }
# 
# 
# result <- 
#   do.call(rbind, result)
# 
# colnames(result) <-
#   c("name1", "name2", "cor", "mz_diff", "rt_diff")
# 
# result <- tibble::as_tibble(result)
# 
# idx =  5
# 
# 
# metabolite_tags %>% 
#   filter(name %in% c(as.character(result[idx,1]), as.character(result[idx,2]))) %>% 
#   select(c(name:rt)) %>% 
#   as.data.frame()
#   
# 
# 
# int1 <- metabolite_table %>% 
#   filter(name == as.character(result[idx,1])) %>% 
#   select(-c(name:rt)) %>% 
#   as.numeric()
# 
# 
# int2 <- metabolite_table %>% 
#   filter(name == as.character(result[idx,2])) %>% 
#   select(-c(name:rt)) %>% 
#   as.numeric()
# 
# plot(int1, int2)
# abline(0, 1)
# 
# 
# metabolite_table <- 
#   metabolite_table %>% 
#   filter(!name %in% result$name2)
# 
# 
# metabolite_tags <- 
#   metabolite_tags %>% 
#   filter(!name %in% result$name2)
# 
# metabolite_table$name == metabolite_tags$name
# 
# readr::write_csv(metabolite_table, "metabolite_table.csv")
# save(metabolite_table, file = "metabolite_table")

load("metabolite_table")

# readr::write_csv(metabolite_tags, "metabolite_tags.csv")
# save(metabolite_tags, file = "metabolite_tags") 
##sample sampels may have NAs
load("metabolite_tags")








###Cytokine data
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_preparation_for_analysis/cytokine/")
rm(list = ls())
cytokine_data1 <- 
  readxl::read_xlsx("SmartD_urine_cytokines.xlsx", sheet = 1)

cytokine_data2 <- 
  readxl::read_xlsx("SmartD_urine_cytokines.xlsx", sheet = 2)

cytokine_data3 <- 
  readxl::read_xlsx("SmartD_urine_cytokines.xlsx", sheet = 3)

cytokine_data4 <- 
  readxl::read_xlsx("SmartD_urine_cytokines.xlsx", sheet = 4)

cytokine_data5 <- 
  readxl::read_xlsx("SmartD_urine_cytokines.xlsx", sheet = 5)


cytokine_data1 <- cytokine_data1[-1,]
cytokine_data2 <- cytokine_data2[-1,]
cytokine_data3 <- cytokine_data3[-1,]
cytokine_data4 <- cytokine_data4[-1,]
cytokine_data5 <- cytokine_data5[-1,]


cytokine_data2 <- 
  cytokine_data2 %>% 
  select(colnames(cytokine_data1))

cytokine_data3 <- 
  cytokine_data3 %>% 
  select(colnames(cytokine_data1))

cytokine_data4 <- 
  cytokine_data4 %>% 
  select(colnames(cytokine_data1))

cytokine_data5 <- 
  cytokine_data5 %>% 
  select(colnames(cytokine_data1))


cytokine_data <- 
  rbind(cytokine_data1, cytokine_data2, 
        cytokine_data3, cytokine_data4,
        cytokine_data5)

##remove cotrol samples
cytokine_data <-
  cytokine_data %>% 
  filter(stringr::str_detect(Name, "LL-"))


cytokine_data$Name <-
  cytokine_data$Name %>% 
  stringr::str_extract("SFU.{1,8}")

##get the ID of all the informations
urine_id <- 
  readr::read_csv("E:/project/smartD/patient information/SmartD_all346urine.csv")

cyto_id <- 
  match(cytokine_data$Name, 
        urine_id$sample_id_HIMC) %>% 
  `[`(urine_id$sample_id_RPLC,.)

cyto_id[grep("SFU_B", cyto_id)] <- 
  cyto_id[grep("SFU_B", cyto_id)] %>% 
  gsub(pattern = "SFU_B", replacement = "", x = .) %>% 
  paste("X", ., sep = "")

##X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on
temp_idx <- match(paste("X", 178:198, sep = ""), cyto_id)
cyto_id[temp_idx] <-
  cyto_id[temp_idx] %>% 
  stringr:::str_replace("X", "") %>% 
  as.numeric() %>% 
  `-`(177) %>% 
  paste("P", ., sep = "")


cytokine_data$sample.id <- cyto_id

cytokine_pheno <- 
  cytokine_data %>% 
  select(Sample, Well, Name, Type, sample.id)


cytokine_table <- 
  cytokine_data %>% 
  select(-c(Sample, Well, Name, Type, sample.id))


cytokine_tags <- 
  data.frame(name = colnames(cytokine_table), 
             stringsAsFactors = FALSE)

rownames(cytokine_table) <- 
  cytokine_pheno$sample.id


cytokine_table <- 
  t(cytokine_table)


cytokine_table <- 
  tibble::as_tibble(cytokine_table)

cytokine_table <- 
  cytokine_table %>% 
  dplyr::mutate_all(as.numeric)

save(cytokine_pheno, file = "cytokine_pheno")
save(cytokine_tags, file = "cytokine_tags")
save(cytokine_table, file = "cytokine_table")








