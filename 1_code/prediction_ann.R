##first set the work directory to project folder
sxtTools::setwd_project()
setwd("data_analysis20191015/prediction/identification_table/")
##load dataa
rm(list = ls())
load("sample_data_dis")
load("sample_data_val")
load("sample_data_dis_x")
load("sample_data_val_x")
load("metabolite_tags")

####ANN in R
library(neuralnet)
library(nnet)
library(caret)
library(tidyverse)

setwd("ANN/GA_prediction/")

##feature selection
sample_data_dis_y <- 
  sample_data_dis %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

set.seed(400)
model.nn <- train(x = sample_data_dis_x,
                  y = sample_data_dis_y[,1],
                  trControl = trainControl(method = "cv", number = 7),
                  method = "nnet")

print(model.nn)

plot(model.nn)
plot(varImp(model.nn))

importance <- 
  varImp(model.nn)

importance <- 
  importance$importance %>% 
  rownames_to_column(., var = "name") %>% 
  arrange(desc(Overall))


importance[1:80,] %>% 
  left_join(., metabolite_tags, by = "name") %>% 
  arrange(Overall) %>% 
  ggplot(aes(x = factor(Compound.name, levels = Compound.name), y = Overall)) +
  geom_segment(aes(x = factor(Compound.name, levels = Compound.name), 
                   xend = factor(Compound.name, levels = Compound.name), 
                   y = 55, yend = Overall),colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  labs(x = "", y = "Importance") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 8), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  coord_flip()

ggsave(filename = "marker_ann.pdf", width = 7, height = 7)

marker_ann <- 
  importance[1:80,]


marker_ann <- 
  left_join(marker_ann, metabolite_tags, by = "name")

write.csv(marker_ann, "marker_ann.csv", row.names = FALSE)

####parameter tunning
sample_data_dis_x_ann <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  dplyr::select(., one_of(marker_ann$name)) %>% 
  as.matrix()

sample_data_dis_x_y_ann <- 
  data.frame(sample_data_dis_x_ann, 
             y = sample_data_dis_y,
             stringsAsFactors = FALSE)


set.seed(410)
ann_regression <- 
  neuralnet(formula = GA~.,
            data = sample_data_dis_x_y_ann, 
            hidden =  c(10, 10), 
            # algorithm = "backprop",
            act.fct = "logistic",
            linear.output = TRUE)

ann_regression20 <- 
  neuralnet(formula = GA~.,
            data = sample_data_dis_x_y_ann, 
            hidden = 20, 
            # algorithm = "backprop",
            act.fct = "logistic",
            linear.output = TRUE)


# plot(ann_regression)
# plot(ann_regression8)

library(NeuralNetTools)

# plot <- 
#   garson(ann_regression8)
# 
# plot + 
#   coord_flip()


neuralnet::compute(x = ann_regression, 
                   covariate = sample_data_dis_x_ann)$net.result %>% 
  plot(sample_data_dis_y)


##validate in validation dataset
sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_x_ann <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_ann$name))

sample_data_val_x_ann <- 
  apply(sample_data_val_x_ann, 2, as.numeric)


###use validation dataset for validation
set.seed(420)
ann_regression <- 
  neuralnet(formula = GA~.,
            data = sample_data_dis_x_y_ann, 
            hidden =  c(10, 10), 
            act.fct = "logistic",
            linear.output = TRUE)


predicted_y <-
  neuralnet::compute(ann_regression,
                     sample_data_val_x_ann)$net.result

prediction_self <- 
  neuralnet::compute(ann_regression,
                     sample_data_dis_x_ann)$net.result

###from here we can see that why should us linear regression to correct prediction
data.frame(measured = sample_data_dis_y[,1],
           predicted = prediction_self[,1],
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

abs(sample_data_val_y[,1] - predicted_y) %>% 
  mean()

summary(lm(formula = predicted_y~sample_data_val_y[,1]))


write.table("Information", "information.txt")
cat("RMSE for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(abs(sample_data_val_y[,1] - predicted_y) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(summary(lm(formula = predicted_y~sample_data_val_y[,1]))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)


data.frame(
  measured = sample_data_val_y[, 1],
  predicted = predicted_y,
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

predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x), 
           size = nrow(sample_data_dis_x), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x), dis_index)
  
  ann_regression_temp <-
    neuralnet(formula = GA~.,
              data = data.frame(sample_data_dis_x_ann,
                                sample_data_dis_y)[dis_index,],
              hidden = c(10, 10),
              act.fct = "logistic",
              linear.output = TRUE)
  
  predict_y[[i]] <-
    neuralnet::compute(ann_regression_temp,
                       sample_data_dis_x_ann[val_index,])$net.result
  
  # predict_y[[i]] <- 
  #   predict(object = ann_regression_temp, 
  #           newdata = sample_data_x_ann[val_index,])[,1]
  
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
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
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




####parameter tunning
sample_data_dis_x_ann <- 
  sample_data_dis_x %>% 
  as.data.frame() %>% 
  dplyr::select(., one_of(marker_ann$name)) %>% 
  as.matrix()

sample_data_dis_x_y_ann <- 
  data.frame(sample_data_dis_x_ann, 
             y = expected_date_remained_dis,
             stringsAsFactors = FALSE)

set.seed(430)
ann_regression <- 
  neuralnet(formula = y~.,
            data = sample_data_dis_x_y_ann, 
            hidden =  c(10, 10), 
            act.fct = "logistic",
            linear.output = TRUE)

ann_regression20 <- 
  neuralnet(formula = y~.,
            data = sample_data_dis_x_y_ann, 
            hidden = 20, 
            act.fct = "logistic",
            linear.output = TRUE)


# plot(ann_regression)
# plot(ann_regression8)

library(NeuralNetTools)

# plot <- 
#   garson(ann_regression8)
# 
# plot + 
#   coord_flip()


neuralnet::compute(x = ann_regression, 
                   covariate = sample_data_dis_x_ann)$net.result %>% 
  plot(expected_date_remained_dis)


##validate in validation dataset
sample_data_val_x_ann <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_ann$name))

sample_data_val_x_ann <- 
  apply(sample_data_val_x_ann, 2, as.numeric)


###use validation dataset for validation
set.seed(440)
ann_regression <- 
  neuralnet(formula = y~.,
            data = sample_data_dis_x_y_ann, 
            hidden =  c(10, 10), 
            act.fct = "logistic",
            linear.output = TRUE)


predicted_y <-
  neuralnet::compute(ann_regression,
                     sample_data_val_x_ann)$net.result

prediction_self <- 
  neuralnet::compute(ann_regression,
                     sample_data_dis_x_ann)$net.result

###from here we can see that why should us linear regression to correct prediction
data.frame(measured = expected_date_remained_dis[,1],
           predicted = prediction_self[,1],
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

abs(expected_date_remained_val[,1] - predicted_y) %>% 
  mean()

summary(lm(formula = predicted_y~expected_date_remained_val[,1]))

write.table("Information", "information.txt")
cat("RMSE for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(abs(expected_date_remained_val[,1] - predicted_y) %>% 
      mean(), file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for external validation dataset:\n", file = "information.txt", append = TRUE)
cat(summary(lm(formula = predicted_y~expected_date_remained_val[,1]))$adj.r.squared, 
    file = "information.txt", append = TRUE)
cat("\n", file = "information.txt", append = TRUE)


data.frame(
  measured = expected_date_remained_val[, 1],
  predicted = predicted_y,
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

predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_dis_x), 
           size = nrow(sample_data_dis_x), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x), dis_index)
  
  ann_regression_temp <-
    neuralnet(formula = y~.,
              data = data.frame(sample_data_dis_x_ann,
                                y = expected_date_remained_dis)[dis_index,],
              hidden = c(10, 10),
              act.fct = "logistic",
              linear.output = TRUE)
  
  predict_y[[i]] <-
    neuralnet::compute(ann_regression_temp,
                       sample_data_dis_x_ann[val_index,])$net.result
  
  # predict_y[[i]] <- 
  #   predict(object = ann_regression_temp, 
  #           newdata = sample_data_x_ann[val_index,])[,1]
  
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
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_dis.pdf", width = 7, height = 7)




temp_data <- 
  data.frame(model = c("Lasso", "RF", "SVM", "ANN"),
             rmse_i = c(3.35, 3.27, 4.04, 3.89),
             r2_i = c(0.69, 0.74, 0.62, 0.66),
             rmse_e = c(3.59, 3.38, 4.16, 5.87),
             r2_e = c(0.63, 0.70, 0.57, 0.33),
             stringsAsFactors = FALSE)

plot1 <- 
  temp_data %>% 
  select(model, rmse_i, rmse_e) %>% 
  tidyr::pivot_longer(., -model, "Dataset") %>% 
  mutate(model = factor(model, levels = rev(temp_data$model))) %>% 
  ggplot(aes(model, value)) +
  geom_bar(aes(fill = Dataset), colour = "white",
           stat = "identity", position = position_dodge(), width = 0.5) +
  scale_colour_manual(values = c("rmse_i" = "#FFA319FF",
                                 "rmse_e" = "#155F83FF")) +
  scale_fill_manual(values = c("rmse_i" = "#FFA319FF",
                                 "rmse_e" = "#155F83FF")) +
  scale_y_reverse() +
  labs(x = "", y = "Root mean squared error (RMSE)") +
  coord_flip() +
  theme_bw() +
  geom_text(aes(label = value, fill = Dataset), hjust = -1,
            position = position_dodge(width = 0.5), size = 5, colour = "white") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), legend.position = "left", 
        panel.border = element_blank(), panel.grid = element_blank(),
        axis.line.x = element_line(colour = "black"),
        plot.margin = margin(0,0, 0, 0),)



plot1

plot2 <- 
temp_data %>% 
  select(model, r2_i, r2_e) %>% 
  tidyr::pivot_longer(., -model, "Dataset") %>% 
  mutate(model = factor(model, levels = rev(temp_data$model))) %>% 
  ggplot(aes(model, value)) +
  geom_bar(aes(fill = Dataset), colour = "white", 
           stat = "identity", position = position_dodge(), width = 0.5) +
  scale_colour_manual(values = c("r2_i" = "#FFA3197F",
                                 "r2_e" = "#155F837F")) +
  scale_fill_manual(values = c("r2_i" = "#FFA3197F",
                               "r2_e" = "#155F837F")) +
  labs(x = "", y = "Adjusted R2 (RMSE)") +
  coord_flip() +
  theme_bw() +
  geom_text(aes(label = value, fill = Dataset), hjust = 1,
            position = position_dodge(width = 0.5), size = 5, colour = "white") +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 13), legend.position = "right", 
        panel.border = element_blank(), panel.grid = element_blank(),
        axis.line.x = element_line(colour = "black"), 
        plot.margin = margin(0,0, 0, 0),
        axis.ticks.y = element_blank())

library(patchwork)


plot1 + plot2
