##first set the work directory to project folder
setwd(r4projects::get_project_wd())
rm(list = ls())
load("data_analysis20200108/data_preparation_for_analysis/metabolite_tags")
setwd("data_analysis20200108/prediction/features/")

##load dataa
load("sample_data_dis")
load("sample_data_val")
load("sample_data_dis_x")
load("sample_data_val_x")
load("peak_table")
peak_tags <- peak_table %>% 
  select(name, mz, rt)
save(peak_tags, file = "peak_tags")
rm("peak_table")
##############################################################################
####random forest
#############################################################################

setwd("RF/GA_prediction/")
library(randomForest)
##use boruta method
library(Boruta)
library(tidyverse)


marker <- readr::read_csv("marker_rf_final.csv")


marker <- 
  marker %>%
  dplyr::select(name) %>% 
  dplyr::left_join(peak_tags, by = "name") %>% 
  dplyr::left_join(metabolite_tags, by = c("name"))

# xlsx::write.xlsx(marker, "marker_paper.xlsx")


sample_data_dis_y <- 
  sample_data_dis %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

set.seed(200)
boruta_test <- 
  Boruta(x = sample_data_dis_x,
         y = sample_data_dis_y, 
         doTrace = 3, 
         holdHistory = TRUE)

# plot(boruta_test)

# boruta_test

marker_rf <- 
  boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>% 
  names() %>% 
  sort()

marker_rf

###discard the peaks with bad peak shapes
is_table <- peak_tags %>% 
  filter(name %in% marker_rf) %>% 
  select(name, mz, rt)

is_table_pos <- 
  is_table %>% 
  filter(stringr::str_detect(name, "POS"))

is_table_neg <- 
  is_table %>% 
  filter(stringr::str_detect(name, "NEG"))

xlsx::write.xlsx2(as.data.frame(is_table_pos), 
                  file = "peak_shape/is_table_pos.xlsx",
                  row.names = FALSE)

xlsx::write.xlsx2(as.data.frame(is_table_neg), 
                  file = "peak_shape/is_table_neg.xlsx",
                  row.names = FALSE)

peak_data_pos <- 
  metflow2::extractPeaks(path = "./peak_shape/POS/", ppm = 15, 
                         threads = 4,
                         is.table = "is_table_pos.xlsx", 
                         rt.tolerance = 60)


for(i in 1:nrow(is_table_pos)){
  cat(i, " ")
  plot <- metflow2::showPeak(object = peak_data_pos, 
                             peak.index = i, 
                             alpha = 0.5, 
                             interactive = FALSE)
  
  ggsave(plot, filename = file.path("./peak_shape/POS", paste(is_table_pos$name[i], ".png", sep = "")), 
         width = 8, height = 6)
}


is_table_pos <- readxl::read_xlsx("./peak_shape/POS/is_table_pos.xlsx")
marker_rf_pos <-
  is_table_pos %>% 
  dplyr::filter(note == "G") %>% 
  pull(name)

peak_data_neg <- 
  metflow2::extractPeaks(path = "./peak_shape/NEG/", 
                         ppm = 15, 
                         threads = 4, 
                         is.table = "is_table_neg.xlsx")


for(i in 1:nrow(is_table_neg)){
  cat(i, " ")
  plot <- metflow2::showPeak(object = peak_data_neg, 
                             peak.index = i, 
                             alpha = 0.5, 
                             interactive = FALSE)
  
  ggsave(plot, filename = file.path("./peak_shape/NEG", paste(is_table_neg$name[i], ".png", sep = "")), 
         width = 8, height = 6)
}


no_infi <- function(x) all(!is.infinite(x))


is_table_neg <- readxl::read_xlsx("./peak_shape/NEG/is_table_neg.xlsx")
marker_rf_neg <-
  is_table_neg %>% 
  dplyr::filter(note == "G") %>% 
  pull(name)

marker_rf <- c(marker_rf_pos, marker_rf_neg)


temp_data <- 
boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  # left_join(., peak_tags, by = "name") %>% 
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
  filter(name %in% marker_rf) %>% 
  dplyr::select(name, everything())


rm(boruta_test)

  ggplot(temp_data,aes(x = factor(name, name), y = mean)) +
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

# marker_rf <- 
#   marker_rf %>% 
#   left_join(peak_tags, by = "name")

write.csv(marker_rf, "marker_rf.csv", row.names = FALSE)

# marker_rf <- readr::read_csv("marker_rf.csv")
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
fgl.res <- tuneRF(sample_data_dis_x_rf, 
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
               mtry = 2)

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

plot(
  sample_data_val_y[,1],
  predicted_y
)

abline(0, 1)

prediction_self <- 
  predict(object = rf_regression, 
          newx = as.matrix(sample_data_dis_x_rf))

plot(sample_data_dis_y[,1], prediction_self)
abline(0, 1)


data.frame(measured = sample_data_dis_y[,1], 
           predicted = prediction_self,
           stringsAsFactors = FALSE) %>% 
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15))

ggsave("discovery_data_measured_vs_predicted.pdf", width = 7, height = 7)


###from here we can see that why should us linear regression to correct prediction
data.frame(measured = sample_data_dis_y[,1],
           predicted = prediction_self,
           stringsAsFactors = FALSE) %>% 
  mutate(diff = predicted - measured) %>% 
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>% 
  ggplot(aes(measured, diff)) +
  geom_segment(aes(x = measured, xend = measured, 
                   y = 0, yend = diff, colour = colour), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "GA_week (measured)", y = "GA_week (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_error.pdf", width = 7, height = 7)



linear_regression <- 
  lm(formula = sample_data_dis_y[,1] ~ prediction_self)

linear_regression1 <- 
  lm(formula = prediction_self ~ sample_data_dis_y[,1])


data.frame(measured = sample_data_dis_y[,1], 
           predicted = prediction_self,
           stringsAsFactors = FALSE) %>% 
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_abline(intercept = coef(linear_regression1)[1], 
              slope = coef(linear_regression1)[2], 
              linetype = 2, colour = "red") +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15))

ggsave("discovery_data_measured_vs_predicted2.pdf", width = 7, height = 7)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y + coef(linear_regression)[1]


data.frame("measured" = sample_data_val_y[,1], 
           'predicted' = predicted_y2,
           stringsAsFactors = FALSE) %>%  
  ggplot(aes(x = measured, predicted)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15))


ggsave("measured_vs_predicted_val.pdf", width = 7, height = 7)


abs(sample_data_val_y[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~sample_data_val_y[,1]))
cor.test(predicted_y2, sample_data_val_y[,1])
write.table("Information", "information.txt")
cat("RMSE for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(abs(sample_data_val_y[,1] - predicted_y2) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(summary(lm(formula = predicted_y2~sample_data_val_y[,1]))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)

##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_rf), 
           size = nrow(sample_data_dis_x_rf), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_rf), dis_index)
  
  rf_regression_temp <-
    randomForest(x = sample_data_dis_x_rf[dis_index,], 
                 y = sample_data_dis_y[dis_index,1], 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE,
                 mtry = 4)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = rf_regression_temp,
      newx = as.matrix(sample_data_dis_x_rf[dis_index,])
    )
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_dis_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = rf_regression_temp,
      newdata = as.matrix(sample_data_dis_x_rf[val_index,])
    )
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    sample_data_dis_y[val_index,1]
  
}

save(predict_y, file = "predict_y")
save(y, file = "y")

prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)

temp <- 
  prediction %>% 
  arrange(., y) %>% 
  group_by(., y) %>% 
  dplyr::summarise(mean = mean(predict), sd = sd(predict))


abs(temp$y - temp$mean) %>% 
  mean()

summary(lm(formula = temp$y~temp$mean))
cor.test(temp$y, temp$mean)

cat("RMSE for internal validation dataset:\n", file = "information.txt", append = TRUE)

cat(abs(temp$y - temp$mean) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for internal validation dataset:\n", file = "information.txt", append = TRUE)

cat(summary(lm(formula = temp$y~temp$mean))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)

temp %>% 
  mutate(ymax = mean + sd, ymin = mean - sd) %>% 
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_dis.pdf", width = 7, height = 7)



###############################################################################
######for other prediction (time to due)
setwd("../time_to_due_prediction/")

sample_data_dis$Patient_ID

date_info_dis <- 
  sample_data_dis[,c("Patient_ID", "Sample_ID", "GA", "Visit","day_0", "Date.Acquired", "DD", "EDD")] %>% 
  mutate(begin_date = as.Date(day_0,"%m/%d/%Y"), 
         acquired_date = as.Date(Date.Acquired,"%m/%d/%y"),
         due_date = as.Date(DD, "%Y-%m-%d"),
         expected_due_date = as.Date(EDD, "%Y-%m-%d")
  ) %>% 
  dplyr::select(-c(day_0, Date.Acquired, DD, EDD))


expected_date_remained_dis <-
  as.numeric(date_info_dis$due_date - date_info_dis$acquired_date, units = "weeks") %>% 
  as.matrix()


date_info_val <- 
  sample_data_val[,c("Patient_ID", "Sample_ID", "GA", "Visit","day_0", "Date.Acquired", "DD", "EDD")] %>% 
  mutate(begin_date = as.Date(day_0,"%m/%d/%Y"), 
         acquired_date = as.Date(Date.Acquired,"%m/%d/%y"),
         due_date = as.Date(DD, "%Y-%m-%d"),
         expected_due_date = as.Date(EDD, "%Y-%m-%d")
  ) %>% 
  dplyr::select(-c(day_0, Date.Acquired, DD, EDD))

?difftime

expected_date_remained_val <-
  as.numeric(date_info_val$due_date - date_info_val$acquired_date, units = "weeks") %>% 
  as.matrix()

library(randomForest)

##use boruta method
library(Boruta)
set.seed(220)
boruta_test <- 
  Boruta(x = sample_data_dis_x,
         y = expected_date_remained_dis, 
         doTrace = 3, 
         holdHistory = TRUE)

# plot(boruta_test)

boruta_test

marker_rf <- 
  boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>% 
  names() %>% 
  sort()

marker_rf

###discard the peaks with bad peak shapes
is_table <- peak_tags %>% 
  filter(name %in% marker_rf) %>% 
  select(name, mz, rt)

is_table_pos <- 
  is_table %>% 
  filter(stringr::str_detect(name, "POS"))

is_table_neg <- 
  is_table %>% 
  filter(stringr::str_detect(name, "NEG"))

xlsx::write.xlsx2(as.data.frame(is_table_pos), 
                  file = "peak_shape/is_table_pos.xlsx",
                  row.names = FALSE)

xlsx::write.xlsx2(as.data.frame(is_table_neg), 
                  file = "peak_shape/is_table_neg.xlsx",
                  row.names = FALSE)

peak_data_pos <- 
  metflow2::extractPeaks(path = "./peak_shape/POS/", ppm = 15, 
                         threads = 4,
                         is.table = "is_table_pos.xlsx", 
                         rt.tolerance = 60)


for(i in 1:nrow(is_table_pos)){
  cat(i, " ")
  plot <- metflow2::showPeak(object = peak_data_pos, 
                             peak.index = i, 
                             alpha = 0.5, 
                             interactive = FALSE)
  
  ggsave(plot, filename = file.path("./peak_shape/POS", paste(is_table_pos$name[i], ".png", sep = "")), 
         width = 8, height = 6)
}




peak_data_neg <- 
  metflow2::extractPeaks(path = "./peak_shape/NEG/", 
                         ppm = 15, 
                         threads = 4, 
                         is.table = "is_table_neg.xlsx")


for(i in 1:nrow(is_table_neg)){
  cat(i, " ")
  plot <- metflow2::showPeak(object = peak_data_neg, 
                             peak.index = i, 
                             alpha = 0.5, 
                             interactive = FALSE)
  
  ggsave(plot, filename = file.path("./peak_shape/NEG", paste(is_table_neg$name[i], ".png", sep = "")), 
         width = 8, height = 6)
}


is_table_pos <- readxl::read_xlsx("./peak_shape/POS/is_table_pos.xlsx")
marker_rf_pos <-
  is_table_pos %>% 
  filter(note == "G") %>% 
  pull(name)


is_table_neg <- readxl::read_xlsx("./peak_shape/NEG/is_table_neg.xlsx")
marker_rf_neg <-
  is_table_neg %>% 
  filter(note == "G") %>% 
  pull(name)

marker_rf <- c(marker_rf_pos, marker_rf_neg)

no_infi <- function(x) all(!is.infinite(x))

boruta_test$ImpHistory[1:length(marker_rf),] %>% 
  t() %>% 
  as_tibble() %>% 
  select_if(., no_infi) %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  left_join(., peak_tags, by = "name") %>% 
  arrange(., mean) %>% 
  ggplot(.,aes(x = factor(name, name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF", width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

ggsave("marker_rf.pdf", width = 7, height = 7)

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
  filter(name %in% marker_rf) %>% 
  dplyr::select(name, everything())
colnames(marker_rf)[-1] <-
  colnames(marker_rf)[-1] %>% 
  paste(., "importance", sep = "_") 

marker_rf <-
  marker_rf %>%
  left_join(peak_tags, by = "name")

write.csv(marker_rf, "marker_rf.csv", row.names = FALSE)

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

set.seed(234)
fgl.res <- tuneRF(sample_data_dis_x_rf, 
                  expected_date_remained_dis[,1], 
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

ggsave("mtry_vs_error.pdf", width = 7, height = 7)

set.seed(235)
rf_regression <-
  randomForest(x = sample_data_dis_x_rf, 
               y = expected_date_remained_dis[,1], 
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


prediction_self <- 
  predict(object = rf_regression, 
          newx = as.matrix(sample_data_dis_x_rf))


###from here we can see that why should us linear regression to correct prediction
data.frame(measured = expected_date_remained_dis[,1],
           predicted = prediction_self,
           stringsAsFactors = FALSE) %>% 
  mutate(diff = predicted - measured) %>% 
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>% 
  ggplot(aes(measured, diff)) +
  geom_segment(aes(x = measured, xend = measured, 
                   y = 0, yend = diff, colour = colour), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "Time to due (measured)", y = "Predicted error (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_error.pdf", width = 7, height = 7)


linear_regression <- 
  lm(formula = expected_date_remained_dis[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y + coef(linear_regression)[1]

plot(expected_date_remained_val[,1], predicted_y)
abline(0, 1)

plot(expected_date_remained_val[,1], predicted_y2)
abline(0, 1)

abs(expected_date_remained_val[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~expected_date_remained_val[,1]))

write.table("Information", "information.txt")
cat("RMSE for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(abs(expected_date_remained_val[,1] - predicted_y2) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(summary(lm(formula = predicted_y2~expected_date_remained_val[,1]))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)


data.frame("measured" = expected_date_remained_val[,1], 
           'predicted' = predicted_y2,
           stringsAsFactors = FALSE) %>%  
  ggplot(aes(x = measured, predicted)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "Time to due (weeks, measured)", y = "Time to due (week, predicted)") +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, 25)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15))

ggsave('measured_vs_predicted_val.pdf', width = 7, height = 7)

##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_rf), 
           size = nrow(sample_data_dis_x_rf), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_rf), dis_index)
  
  rf_regression_temp <-
    randomForest(x = sample_data_dis_x_rf[dis_index,], 
                 y = expected_date_remained_dis[dis_index,1], 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE,
                 mtry = 4)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = rf_regression_temp,
      newx = as.matrix(sample_data_dis_x_rf[dis_index,])
    )
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = expected_date_remained_dis[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = rf_regression_temp,
      newdata = as.matrix(sample_data_dis_x_rf[val_index,])
    )
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    expected_date_remained_dis[val_index,1]
  
}

save(predict_y, file = "predict_y")
save(y, file = "y")


prediction <-
  data.frame(y = unlist(y),
             predict = unlist(predict_y),
             stringsAsFactors = FALSE)

temp <- 
  prediction %>% 
  arrange(., y) %>% 
  group_by(., y) %>% 
  dplyr::summarise(mean = mean(predict), sd = sd(predict))

plot(temp$y, temp$mean)


abline(0,1)

abs(temp$y - temp$mean) %>% 
  mean()
summary(lm(formula = temp$y~temp$mean))

cat("RMSE for internal validation dataset:\n", file = "information.txt", append = TRUE)

cat(abs(temp$y - temp$mean) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for internal validation dataset:\n", file = "information.txt", append = TRUE)

cat(summary(lm(formula = temp$y~temp$mean))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)

temp %>% 
  mutate(ymax = mean + sd, ymin = mean - sd) %>% 
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, 25)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_dis.pdf", width = 7, height = 7)




setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/prediction/features/")
load("metabolite_tags")
load("peak_tags")
setwd("RF/GA_prediction/")
marker <- readr::read_csv("marker_rf_final.csv")

marker <- 
marker %>% 
  left_join(metabolite_tags, by = "name")


marker$Compound.name
