###use metabolites
##first set the work directory to project folder
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/prediction/metabolites/")
rm(list = ls())
##load dataa
load("sample_data_dis")
load("sample_data_val")
load("sample_data_dis_x")
load("sample_data_val_x")
load("metabolite_tags")

setwd("RF/preterm/")
library(randomForest)

library(Boruta)
library(tidyverse)


##remove the C/S
sample_data_dis <- 
sample_data_dis %>% 
  mutate(GA2 = case_when(is.na(GA) ~ 45,
                         !is.na(GA) ~ GA)) %>% 
  mutate(class = case_when(is.na(GA) ~ "PP",
                           !is.na(GA) ~ "Normal")) %>% 
  mutate(Begin.Date = as.Date(EDD) - 280) %>% 
  mutate(term.date = as.Date(DD) - Begin.Date) %>% 
  mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>% 
  mutate(birth_g_stage = 40 + (as.Date(DD) - as.Date(EDD))/7) %>% 
  arrange(diff_day) %>% 
  dplyr::filter(`C/S` == "N")
# dplyr::filter(Spont == "Y" | birth_g_stage >= 39)


sample_data_dis %>% 
  dplyr::filter(diff_day >= 7) %>% 
  pull(Patient_ID) %>% 
  unique()


sample_data_val %>% 
  dplyr::filter(diff_day >= 7) %>% 
  pull(Patient_ID) %>% 
  unique()

sample_data_val <- 
  sample_data_val %>% 
  mutate(GA2 = case_when(is.na(GA) ~ 45,
                         !is.na(GA) ~ GA)) %>% 
  mutate(class = case_when(is.na(GA) ~ "PP",
                           !is.na(GA) ~ "Normal")) %>% 
  mutate(Begin.Date = as.Date(EDD) - 280) %>% 
  mutate(term.date = as.Date(DD) - Begin.Date) %>% 
  mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>% 
  mutate(birth_g_stage = 40 + (as.Date(DD) - as.Date(EDD))/7) %>% 
  arrange(diff_day) %>% 
  dplyr::filter(`C/S` == "N")
  # dplyr::filter(Spont == "Y" | birth_g_stage >= 39)

sample_data_val %>% 
  dplyr::filter(diff_day >= 7) %>% 
  pull(Patient_ID) %>% 
  unique()

colnames(sample_data_val) == colnames(sample_data_dis)

sample_data_dis <- 
  sample_data_dis %>% 
  dplyr::mutate(batch = "1")


sample_data_val <- 
  sample_data_val %>% 
  dplyr::mutate(batch = "2")


sample_data <-
  rbind(sample_data_dis, sample_data_val)

phenotype_data <- 
  sample_data %>% 
  dplyr::select(-matches("POS$")) %>% 
  dplyr::select(-matches("NEG$"))


expression_data <- 
  sample_data %>% 
  dplyr::select(dplyr::matches("POS$|NEG$"))

colnames(expression_data) == metabolite_tags$name

variable_data <- 
  metabolite_tags

# class samples to different time points
phenotype_data %>% 
  dplyr::select(Patient_ID, Sample_ID, GA, diff_day) %>% 
  group_by(Patient_ID) %>% 
  dplyr::arrange(GA) %>% 
  ungroup() %>% 
  plyr::dlply(.(Patient_ID))


phenotype_data <- 
phenotype_data %>% 
  dplyr::mutate(Time_range = case_when(
    GA > 10 & GA <= 15 ~ "10_15",
    GA > 15 & GA <= 20 ~ "15_20",
    GA > 20 & GA <= 25 ~ "20_25",
    GA > 25 & GA <= 30 ~ "25_30",
    GA > 30 & GA <= 35 ~ "30_35",
    GA > 35 & GA <= 40 ~ "35_40",
    GA > 35 & GA <= 45 ~ "40_45"
  )) 

test <- 
phenotype_data %>% 
  dplyr::filter(diff_day >= 7) %>% 
  dplyr::select(Patient_ID, Time_range, diff_day)

unique(test$Time_range)

test <- 
  phenotype_data %>% 
  dplyr::filter(diff_day >= 21) %>% 
  dplyr::select(Patient_ID, Time_range, diff_day)

unique(test$Time_range)



##use the time range 10-20 to predict preterm before 7
##get the subject  have 10-20 time point
temp_idx <- grep("10_15|15_20|20_25|25_30_30_35_35_40", phenotype_data$Time_range)

temp_phenotype_data <- 
  phenotype_data[temp_idx, ]

phenotype_data %>% 
  dplyr::filter(diff_day >= 7) %>% 
  dplyr::select(Patient_ID, Sample_ID)

temp_phenotype_data %>% 
  dplyr::filter(diff_day >= 7) %>% 
  dplyr::select(Patient_ID, Sample_ID) %>% 
  pull(Patient_ID) %>% 
  unique()

temp_expression_data <-
  expression_data[temp_idx,]


sample_data_y <- 
  temp_phenotype_data %>% 
  dplyr::select(diff_day) %>% 
  mutate(y = case_when(
    diff_day >=7 ~ 1,
    TRUE ~ 0
  )) %>% 
  dplyr::select(y) %>% 
  as.matrix()


###get discovery and validation dataset, according to patients
# patient_id_preterm <-
#   temp_phenotype_data$Patient_ID[which(sample_data_y[,1] == 1)] %>%
#   unique()
# 
# patient_id_normal <-
#   temp_phenotype_data$Patient_ID[which(sample_data_y[,1] == 0)] %>%
#   unique()
# 
# patient_id_preterm_dis <-
#   sample(patient_id_preterm, round(length(patient_id_preterm)/2))
# 
# patient_id_normal_dis <-
#   sample(patient_id_normal, round(length(patient_id_normal)/2))
# 
# 
# patient_id_preterm_val <-
#   setdiff(patient_id_preterm, patient_id_preterm_dis)
# 
# 
# patient_id_normal_val <-
#   setdiff(patient_id_normal, patient_id_normal_dis)
# 
# 
# index_dis <-
#   which(temp_phenotype_data$Patient_ID %in% c(patient_id_normal_dis,
#                                              patient_id_preterm_dis
#                                              ))
# index_val <-
#   setdiff(1:nrow(sample_data_y), index_dis)


index_dis <-
  which(temp_phenotype_data$batch == 1)

index_val <-
  which(temp_phenotype_data$batch == 2)

sample_data_dis_x <- expression_data[index_dis,]
sample_data_dis_y <- sample_data_y[index_dis,,drop = FALSE]

sample_data_val_x <- expression_data[index_val,]
sample_data_val_y <- sample_data_y[index_val,,drop = FALSE]

  marker_rf <- 
  vector(mode = "list", length = 10)

set.seed(200)
for(i in 1:10) {
  cat(i, " ")
  boruta_test <- 
    Boruta(x = sample_data_dis_x,
           y = sample_data_dis_y, 
           doTrace = 3, 
           holdHistory = TRUE)  
  marker_rf[[i]] <- 
    boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>% 
    names() %>% 
    sort()
}

save(marker_rf, file = "marker_rf")

marker_rf_new <- 
  unlist(marker_rf) %>% 
  data.frame(name = ., stringsAsFactors = FALSE) %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(n = n()) %>% 
  arrange(desc(n))


print(marker_rf_new, n = 100)


marker_rf <- 
  marker_rf_new %>% 
  dplyr::filter(n >= 5) %>% 
  pull(name)

marker_rf

# ###discard the peaks with bad peak shapes
# is_table <- variable_data %>% 
#   filter(name %in% marker_rf) %>% 
#   select(name, mz, rt)
# 
# is_table_pos <- 
#   is_table %>% 
#   filter(stringr::str_detect(name, "POS"))
# 
# is_table_neg <- 
#   is_table %>% 
#   filter(stringr::str_detect(name, "NEG"))
# 
# 
# xlsx::write.xlsx2(as.data.frame(is_table_pos), 
#                  file = "peak_shape/is_table_pos.xlsx",
#                  row.names = FALSE)
# 
# xlsx::write.xlsx2(as.data.frame(is_table_neg), 
#                   file = "peak_shape/is_table_neg.xlsx",
#                   row.names = FALSE)
# 
# peak_data_pos <- 
#   metflow2::extractPeaks(path = "./peak_shape/POS/", ppm = 15, 
#                          threads = 4,
#                          is.table = "is_table_pos.xlsx")
# 
# 
# for(i in 1:nrow(is_table_pos)){
#   cat(i, " ")
#   plot <- metflow2::showPeak(object = peak_data_pos, 
#                              peak.index = i, 
#                              alpha = 0.5, 
#                              interactive = FALSE)
#   
#   ggsave(plot, filename = file.path("./peak_shape/POS", 
#                                     paste(is_table_pos$name[i], ".png", sep = "")), 
#          width = 8, height = 6)
# }
# 
# 
# is_table_pos <- readxl::read_xlsx("./peak_shape/POS/is_table_pos.xlsx")
# 
# marker_rf_pos <-
#   is_table_pos %>% 
#   filter(peak_shape == "G") %>% 
#   pull(name)
# 
# 
# peak_data_neg <- 
#   metflow2::extractPeaks(path = "./peak_shape/NEG/", ppm = 15, 
#                          threads = 4, 
#                          is.table = "is_table_neg.xlsx")
# 
# for(i in 1:nrow(is_table_neg)){
#   cat(i, " ")
#   plot <- metflow2::showPeak(object = peak_data_neg, 
#                              peak.index = i, 
#                              alpha = 0.5, 
#                              interactive = FALSE)
#   
#   ggsave(plot, filename = file.path("./peak_shape/NEG", 
#                                     paste(is_table_neg$name[i], ".png", sep = "")), 
#          width = 8, height = 6)
# }
# 
# 
# is_table_neg <- readxl::read_xlsx("./peak_shape/NEG/is_table_neg.xlsx")
# marker_rf_neg <-
#   is_table_neg %>% 
#   filter(peak_shape == "G") %>% 
#   pull(name)
# 
# marker_rf <- c(marker_rf_pos, marker_rf_neg)

marker_rf <- 
  variable_data %>% 
  dplyr::filter(name %in% marker_rf)

readr::write_csv(marker_rf, "marker_rf.csv")


# ###check identifiction
# load("./annotation_confirmation/massbankDatabase0.0.1")
# load("./annotation_confirmation/metlinDatabase0.0.1")
# load("./annotation_confirmation/monaDatabase0.0.1")
# load("./annotation_confirmation/msDatabase_rplc0.0.1")
# load("./annotation_confirmation/nistDatabase0.0.1")
# load("./annotation_confirmation/orbitrapDatabase0.0.1")
# load("./annotation_confirmation/result.pRPLC.nce50")
# load("./annotation_confirmation/result.pRPLC.nce25")
# load("./annotation_confirmation/result.nRPLC.nce50")
# load("./annotation_confirmation/result.nRPLC.nce25")
# 
# 
# marker_rf <- 
# dplyr::filter(variable_data, name %in% marker_rf)
# 
# ###POS
# database_name <- 
#   marker_rf %>% 
#   filter(stringr::str_detect(name, "POS")) %>% 
#   pull(Database) %>% 
#   unique()
# 
# for(i in 1:length(database_name)){
#   cat(i, " ")
#   temp_database <- database_name[i]
#   idx <- lapply(result.pRPLC.nce25, function(x) x@database) %>% 
#     unlist() %>% 
#     unname() %>% 
#     `==`(temp_database) %>% 
#     which()
#   temp_name <- marker_rf %>% 
#     filter(Database == temp_database) %>% 
#     filter(stringr::str_detect(name, "POS")) %>% 
#     pull(name)
#   
#   if(length(temp_name) == 0){
#     next()
#   }else{
#     metID::ms2plot(object = result.pRPLC.nce25[[idx]], 
#                    which.peak = temp_name,
#                    database = get(temp_database), 
#                    path = file.path("annotation_confirmation", 
#                                     paste("POS",temp_database, "NCE25", sep= "_")), 
#                    show.plot = FALSE) 
#     
#     metID::ms2plot(object = result.pRPLC.nce50[[idx]], 
#                    which.peak = temp_name,
#                    database = get(temp_database), 
#                    path =  file.path("annotation_confirmation",
#                                      paste("POS",temp_database, "NCE50", sep= "_")),
#                    , show.plot = FALSE)
#   }
#   
# }
# 
# 
# 
# ###neg
# database_name <- 
#   marker_rf %>% 
#   filter(stringr::str_detect(name, "NEG")) %>% 
#   pull(Database) %>% 
#   unique()
# 
# for(i in 1:length(database_name)){
#   cat(i, " ")
#   temp_database <- database_name[i]
#   idx <- lapply(result.nRPLC.nce25, function(x) x@database) %>% 
#     unlist() %>% 
#     unname() %>% 
#     `==`(temp_database) %>% 
#     which()
#   temp_name <- marker_rf %>% 
#     filter(Database == temp_database) %>% 
#     filter(stringr::str_detect(name, "NEG")) %>% 
#     pull(name)
#   
#   if(length(temp_name) == 0){
#     next()
#   }else{
#     metID::ms2plot(object = result.nRPLC.nce25[[idx]], 
#                    which.peak = temp_name,
#                    database = get(temp_database), 
#                    path = file.path("annotation_confirmation", 
#                                     paste("NEG",temp_database, "NCE25", sep= "_")), 
#                    show.plot = FALSE) 
#     
#     metID::ms2plot(object = result.nRPLC.nce50[[idx]], 
#                    which.peak = temp_name,
#                    database = get(temp_database), 
#                    path =  file.path("annotation_confirmation",
#                                      paste("NEG",temp_database, "NCE50", sep= "_")),
#                    , show.plot = FALSE)
#   }
#   
# }


marker_rf <- readr::read_csv("marker_rf.csv")

marker_rf <-
  marker_rf %>% 
  # filter(MS2_match == "G") %>% 
  pull(name)


no_infi <- function(x) all(!is.infinite(x))

temp_data <- 
  boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  dplyr::filter(name %in% marker_rf) %>% 
  left_join(., variable_data, by = "name") %>%
  arrange(., mean)




marker_rf <- 
  boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  dplyr::filter(name %in% marker_rf) %>% 
  dplyr::select(name, everything())


ggplot(temp_data,aes(x = factor(Compound.name, Compound.name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF", width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

ggsave(filename = "marker_rf.pdf", width = 7, height = 7)

colnames(marker_rf)[-1] <-
  colnames(marker_rf)[-1] %>% 
  paste(., "importance", sep = "_") 

marker_rf <-
  marker_rf %>%
  left_join(variable_data, by = "name")

write.csv(marker_rf, "marker_rf_final.csv", row.names = FALSE)

####parameter tunning
sample_data_dis_x_rf <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  dplyr::select(., one_of(marker_rf$name)) %>% 
  as.matrix()

sample_data_dis_x_rf <- 
  apply(sample_data_dis_x_rf, 2, as.numeric)

sample_data_val_x_rf <- 
  sample_data_val_x %>% 
  as.data.frame() %>% 
  dplyr::select(., one_of(marker_rf$name)) %>% 
  as.matrix()

sample_data_val_x_rf <- 
  apply(sample_data_val_x_rf, 2, as.numeric)

set.seed(210)
fgl.res <- randomForest::tuneRF(sample_data_dis_x_rf, 
                  sample_data_dis_y[,1], 
                  mtryStart = 1,
                  stepFactor = 2, 
                  trace = TRUE, 
                  plot = TRUE)

fgl.res %>% 
  as.data.frame() %>% 
  ggplot(aes(mtry, y = OOBError)) +
  geom_line(colour = "#155F83FF") +
  geom_point(size = 4, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

ggsave(filename = "mtry_vs_error.pdf", width = 7, height = 7)

rf_regression <-
  randomForest(x = sample_data_dis_x_rf, 
               y = sample_data_dis_y[,1], 
               replace = TRUE, 
               importance = TRUE,
               proximity = TRUE, 
               mtry = 4)


##validate in validation dataset
###construct dataset for lasso regression
sample_data_val_x_rf <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_rf$name))

###use validation dataset for validation
predicted_y <-
  predict(
    object = rf_regression,
    newdata = sample_data_val_x_rf
    # type = "response"
  )

library(pROC)

roc <- 
roc(sample_data_val_y[,1], predicted_y)
pROC::auc(roc)

plot(roc)

temp_data <-
  data.frame(x = sample_data_val_y, 
             y = predicted_y, 
             stringsAsFactors = FALSE)

basicplot <- ggplot(temp_data, 
                    aes(d = x, m = y)) + 
  geom_roc()


##bootstrap
##should be 1000 for true 
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  patient_id <- 
  temp_phenotype_data$Patient_ID %>% 
    unique() 
  
  preterm_id <- 
  temp_phenotype_data %>% 
    dplyr::filter(diff_day >= 7) %>% 
    pull(Patient_ID) %>% 
    unique()
  
  normal_id <- 
    temp_phenotype_data %>% 
    dplyr::filter(diff_day < 7) %>% 
    pull(Patient_ID) %>% 
    unique()
  
  dis_id <-   
    c(sample(preterm_id, round(length(preterm_id)/2)),
      sample(normal_id, round(length(normal_id)/2)))
  
  val_id <- 
    setdiff(patient_id, dis_id)
  
  dis_index <- which(temp_phenotype_data$Patient_ID %in% dis_id)
  val_index <- which(temp_phenotype_data$Patient_ID %in% val_id)
  
  rf_regression_temp <-
    randomForest(x = temp_expression_data[dis_index,] %>% dplyr::select(marker_rf$name), 
                 y = sample_data_y[dis_index,1], 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE,
                 mtry = 4)
  

  temp_predict_y <- 
    predict(
      object = rf_regression_temp,
      newdata = as.matrix(temp_expression_data[val_index,] %>% dplyr::select(marker_rf$name))
    )

  predict_y[[i]] <- 
    temp_predict_y
  y[[i]] <- sample_data_y[val_index,1]
}


save(predict_y, file = "predict_y")
save(y, file = "y")

prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)


roc <- 
  roc(unlist(y), unlist(predict_y))
pROC::auc(roc)

plot(roc)







