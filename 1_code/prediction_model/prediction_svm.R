##first set the work directory to project folder
setwd(r4projects::get_project_wd())
setwd("data_analysis20191015/prediction/identification_table/")
##load dataa
rm(list = ls())
load("sample_data_dis")
load("sample_data_val")
load("sample_data_dis_x")
load("sample_data_val_x")
load("metabolite_tags")
####SVM regression
library(e1071)
##feature selection
library(caret)
library(penalizedSVM)
library(tidyverse)

sample_data_dis_y <- 
  sample_data_dis %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

setwd("SVM/GA_prediction/")
set.seed(300)
svm_test <- 
  caret::rfe(x = sample_data_dis_x, 
             y = sample_data_dis_y[,1], 
             sizes = seq(5,500,20),
             rfeControl	= rfeControl(functions = caretFuncs, 
                                     method = "cv",
                                     number = 7),
             method = "svmRadial")

typeof(svm_test)

svm_test$variables %>% 
  head()


svm_test$results %>% 
  head()


svm_test$bestSubset

marker_svm <-
  svm_test$optVariables

marker_svm <- 
  data.frame(marker_svm, stringsAsFactors = FALSE) %>% 
  left_join(., metabolite_tags, by = c("marker_svm" = "name")) %>% 
  dplyr::rename(name = marker_svm)


write.csv(marker_svm, file = "marker_svm.csv", row.names = FALSE)

####parameter tunning
sample_data_dis_x_svm <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  dplyr::select(., one_of(marker_svm$name)) %>% 
  as.matrix()

set.seed(310)
rbf_tune <- 
  tune.svm(x = sample_data_dis_x_svm, 
           y = sample_data_dis_y[,1], 
           kernel = "radial",
           gamma = 2^(-10:-3), 
           cost = 2^(-5:4),
           tunecontrol = tune.control(sampling = "cross", 
                                      cross = 7, 
                                      best.model = TRUE, 
                                      performances = TRUE)
  )


pdf(file = "parameter_optimization.pdf", width = 7, height = 7)
plot(x = rbf_tune, xlab = "Gamma", main = "", 
     cex.lab = 1.5, cex.axis = 1.3)
par(xpd = FALSE)
points(x = rbf_tune$best.parameters[1,1], 
       y = rbf_tune$best.parameters[1,2], 
       col = "#FFA319FF", cex = 1.5, pch = 19)

abline(v = rbf_tune$best.parameters[1,1],
       col = "#FFA319FF")

abline(h = rbf_tune$best.parameters[1,2],
       col = "#FFA319FF")
dev.off()

rbf_tune$best.performance

rbf_tune$best.parameters

##validate in validation dataset
sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_x_svm <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_svm$name))

sample_data_val_x_svm <- 
  apply(sample_data_val_x_svm, 2, as.numeric)

sample_data_dis_x_svm <- 
  sample_data_dis_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_svm$name))

sample_data_dis_x_svm <- 
  apply(sample_data_dis_x_svm, 2, as.numeric)


svm_regression <-
  svm(x = sample_data_dis_x_svm,
      y = sample_data_dis_y[, 1],
      scale = FALSE, 
      kernel = "radial", 
      cost = rbf_tune$best.parameters[1,2],
      gamma = rbf_tune$best.parameters[1,1]
  )

predicted_y <- 
  predict(
    object = svm_regression,
    newdata = sample_data_val_x_svm
  )


prediction_self <- 
  predict(object = svm_regression, 
          newdata = as.matrix(sample_data_dis_x_svm)
          )

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
  geom_smooth(colour = "black", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "GA (weeks, measured)", y = "Predicted error (weeks, predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


ggsave("measured_vs_predicted_error.pdf", width = 7, height = 7)

linear_regression <- 
  lm(formula = sample_data_dis_y[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y + coef(linear_regression)[1]

plot(sample_data_val_y[,1], predicted_y)
abline(0, 1)

plot(sample_data_val_y[,1], predicted_y2)
abline(0, 1)

abs(sample_data_val_y[,1] - predicted_y2) %>% 
  mean()
summary(lm(formula = predicted_y2~sample_data_val_y[,1]))

write.table("Information", "information.txt")
cat("RMSE for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(abs(sample_data_val_y[,1] - predicted_y2) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(summary(lm(formula = predicted_y2~sample_data_val_y[,1]))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)

data.frame(
  measured = sample_data_val_y[, 1],
  predicted = predicted_y2,
  stringsAsFactors = FALSE
) %>%
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_val.pdf", width = 7, height = 7)


##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)
for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_svm), 
           size = nrow(sample_data_dis_x_svm), 
           replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_svm), dis_index)
  
  svm_regression_temp <-
    svm(x = sample_data_dis_x_svm[dis_index,],
        y = sample_data_dis_y[dis_index, 1],
        scale = FALSE, kernel = "radial", 
        cost = rbf_tune$best.parameters[1,2],
        gamma = rbf_tune$best.parameters[1,1]
    )
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = svm_regression_temp,
      newx = as.matrix(sample_data_dis_x_svm[dis_index,])
    )
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_dis_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = svm_regression_temp,
      newdata = sample_data_dis_x_svm[val_index,]
    )
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    sample_data_dis_y[val_index,1]
}

save(predicted_y, file = "predicted_y")
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
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_dis.pdf", width = 7, height = 7)
















#################################################################################
####predicted time to due date
#################################################################################
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

set.seed(320)

svm_test <- 
  caret::rfe(x = sample_data_dis_x, 
             y = expected_date_remained_dis[,1], 
             sizes = seq(5,500,20),
             rfeControl	= rfeControl(functions = caretFuncs, 
                                     method = "cv",
                                     number = 7),
             method = "svmRadial")

typeof(svm_test)

svm_test$variables %>% 
  head()


svm_test$results %>% 
  head()

svm_test$bestSubset

marker_svm <-
  svm_test$optVariables

marker_svm <- 
  data.frame(marker_svm, stringsAsFactors = FALSE) %>% 
  left_join(., metabolite_tags, by = c("marker_svm" = "name")) %>% 
  dplyr::rename(name = marker_svm)


write.csv(marker_svm, file = "marker_svm.csv", row.names = FALSE)

####parameter tunning
sample_data_dis_x_svm <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  dplyr::select(., one_of(marker_svm$name)) %>% 
  as.matrix()

set.seed(330)

rbf_tune <- 
  tune.svm(x = sample_data_dis_x_svm, 
           y = expected_date_remained_dis[,1], 
           kernel = "radial",
           gamma = 2^(-10:-3), 
           cost = 2^(-5:4),
           tunecontrol = tune.control(sampling = "cross", 
                                      cross = 7, 
                                      best.model = TRUE, 
                                      performances = TRUE)
  )


pdf(file = "parameter_optimization.pdf", width = 7, height = 7)
plot(x = rbf_tune, xlab = "Gamma", main = "", 
     cex.lab = 1.5, cex.axis = 1.3)
par(xpd = FALSE)
points(x = rbf_tune$best.parameters[1,1], 
       y = rbf_tune$best.parameters[1,2], 
       col = "#FFA319FF", cex = 1.5, pch = 19)

abline(v = rbf_tune$best.parameters[1,1],
       col = "#FFA319FF")

abline(h = rbf_tune$best.parameters[1,2],
       col = "#FFA319FF")
dev.off()

rbf_tune$best.performance

rbf_tune$best.parameters

##validate in validation dataset
sample_data_val_x_svm <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_svm$name))

sample_data_val_x_svm <- 
  apply(sample_data_val_x_svm, 2, as.numeric)

sample_data_dis_x_svm <- 
  sample_data_dis_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_svm$name))

sample_data_dis_x_svm <- 
  apply(sample_data_dis_x_svm, 2, as.numeric)


svm_regression <-
  svm(x = sample_data_dis_x_svm,
      y = expected_date_remained_dis[, 1],
      scale = FALSE, 
      kernel = "radial", 
      cost = rbf_tune$best.parameters[1,2],
      gamma = rbf_tune$best.parameters[1,1]
  )

predicted_y <- 
  predict(
    object = svm_regression,
    newdata = sample_data_val_x_svm
  )


prediction_self <- 
  predict(object = svm_regression, 
          newdata = as.matrix(sample_data_dis_x_svm)
  )


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
  geom_smooth(colour = "black", fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "GA (weeks, measured)", y = "Predicted error (weeks, predicted - measured)") +
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

data.frame(
  measured = expected_date_remained_val[, 1],
  predicted = predicted_y2,
  stringsAsFactors = FALSE
) %>%
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_val.pdf", width = 7, height = 7)


##bootstrap
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x_svm), 
           size = nrow(sample_data_dis_x_svm), 
           replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_svm), dis_index)
  
  svm_regression_temp <-
    svm(x = sample_data_dis_x_svm[dis_index,],
        y = expected_date_remained_dis[dis_index, 1],
        scale = FALSE, kernel = "radial", 
        cost = rbf_tune$best.parameters[1,2],
        gamma = rbf_tune$best.parameters[1,1]
    )
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = svm_regression_temp,
      newx = as.matrix(sample_data_dis_x_svm[dis_index,])
    )
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = expected_date_remained_dis[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = svm_regression_temp,
      newdata = sample_data_dis_x_svm[val_index,]
    )
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  y[[i]] <- 
    expected_date_remained_dis[val_index,1]
}

save(predicted_y, file = "predicted_y")
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
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_dis.pdf", width = 7, height = 7)





