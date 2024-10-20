##first set the work directory to project folder
rm(list = ls())
source("R/tools.R")
setwd("data_analysis20191015/prediction/identification_table/")
##load dataa
load("sample_data_dis")
load("sample_data_val")
load("sample_data_dis_x")
load("sample_data_val_x")
load("metabolite_tags")
################################################################################
###################################LASSO regression#############################
################################################################################
##using lasso regressio to select variable
##
setwd("lasso/GA_prediction/")
###Predict GA first
library(glmnet)
sample_data_dis_y <- 
  sample_data_dis %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

set.seed(100)
lasso_regression2 <-
  cv.glmnet(
    x = sample_data_dis_x,
    y = sample_data_dis_y,
    family = "gaussian",
    type.measure = "mae",
    nfolds = 7,
    nlambda = 100,
    alpha = 1,
    standardize = FALSE
  )

plotLambdaVScoefficients(object = lasso_regression2$glmnet.fit)
ggsave(filename = "lambda_vs_coefficients.pdf", width = 7, height = 7)

plotLambdaVSerror(object = lasso_regression2)
ggsave(filename = "lambda_vs_error.pdf", width = 7, height = 7)

##construct the best model
best_lasso <-
  glmnet(sample_data_dis_x, 
         sample_data_dis_y,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se)

marker_lasso <- 
  which(lasso_regression2$lambda == lasso_regression2$lambda.1se) %>% 
  `[`(lasso_regression2$glmnet.fit$beta, ,.) %>% 
  # lasso_regression2$glmnet.fit$beta[,47] %>% 
  tibble(name = rownames(lasso_regression2$glmnet.fit$beta),
         coef = .) %>% 
  filter(coef != 0) 

##markers
marker_lasso %>% 
  left_join(., metabolite_tags, by = c("name")) %>% 
  arrange(coef) %>%
  # arrange(., desc(coef)) %>% 
  ggplot(aes(x = factor(Compound.name, Compound.name), y = coef)) +
  labs(x = "", y = "Coefficents") +
  geom_segment(aes(x = factor(Compound.name, Compound.name), 
                   xend = factor(Compound.name, Compound.name),
                   y = 0, yend = coef), 
               colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

ggsave(filename = "marker_lasso.pdf", width = 7, height = 7)

marker_lasso <- 
  left_join(marker_lasso, metabolite_tags, by = "name")

write.csv(marker_lasso, "marker_lasso.csv", row.names = FALSE)

##validate in validation dataset
sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_x_lasso <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_val_x_lasso <- 
  apply(sample_data_val_x_lasso, 2, as.numeric)

sample_data_dis_x_lasso <- 
  sample_data_dis_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_dis_x_lasso <- 
  apply(sample_data_dis_x_lasso, 2, as.numeric)

###use validation dataset for validation
best_lasso2 <-
  glmnet(as.matrix(sample_data_dis_x_lasso), 
         sample_data_dis_y,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

predicted_y <-
  predict(
    object = best_lasso2,
    newx = as.matrix(sample_data_val_x_lasso),
    s = lasso_regression2$lambda.1se
    # type = "response"
  )

prediction_self <- 
  predict(object = best_lasso2, 
          newx = as.matrix(sample_data_dis_x_lasso), 
          s = lasso_regression2$lambda.1se)

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

linear_regression <- 
  lm(formula = sample_data_dis_y[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y[,1] + coef(linear_regression)[1]

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
    sample(1:nrow(sample_data_dis_x_lasso), 
           size = nrow(sample_data_dis_x_lasso), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_lasso), dis_index)
  
  lasso_regression_best <- 
    glmnet(as.matrix(sample_data_dis_x_lasso[dis_index,]), 
           sample_data_dis_y[dis_index,],
           family = "gaussian", 
           alpha = 1, 
           lambda = lasso_regression2$lambda.1se)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[dis_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_dis_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[val_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
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
  dplyr::arrange(., y) %>% 
  dplyr::group_by(., y) %>% 
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
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
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


###construct dataset for lasso regression
set.seed(130)
lasso_regression2 <-
  cv.glmnet(
    x = sample_data_dis_x,
    y = expected_date_remained_dis,
    family = "gaussian",
    type.measure = "mae",
    nfolds = 7,
    nlambda = 100,
    alpha = 1,
    standardize = FALSE
  )

plotLambdaVScoefficients(object = lasso_regression2$glmnet.fit)
ggsave(filename = "lambda_vs_coefficients.pdf", width = 7, height = 7)

plotLambdaVSerror(object = lasso_regression2)
ggsave(filename = "lambda_vs_error.pdf", width = 7, height = 7)


##construct the best model
best_lasso <-
  glmnet(sample_data_dis_x, 
         expected_date_remained_dis,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se)

marker_lasso <- 
  which(lasso_regression2$lambda == lasso_regression2$lambda.1se) %>% 
  `[`(lasso_regression2$glmnet.fit$beta, ,.) %>% 
  # lasso_regression2$glmnet.fit$beta[,47] %>% 
  tibble(name = rownames(lasso_regression2$glmnet.fit$beta),
         coef = .) %>% 
  filter(coef != 0) 

##markers
marker_lasso %>% 
  left_join(., metabolite_tags, by = c("name")) %>% 
  arrange(coef) %>%
  # arrange(., desc(coef)) %>% 
  ggplot(aes(x = factor(Compound.name, Compound.name), y = coef)) +
  labs(x = "", y = "Coefficents") +
  geom_segment(aes(x = factor(Compound.name, Compound.name), 
                   xend = factor(Compound.name, Compound.name),
                   y = 0, yend = coef), 
               colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))

marker_lasso <- 
  left_join(marker_lasso, metabolite_tags, by = "name")

write.csv(marker_lasso, "marker_lasso.csv", row.names = FALSE)

##validate in validation dataset
sample_data_val_x_lasso <- 
  sample_data_val_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_val_x_lasso <- 
  apply(sample_data_val_x_lasso, 2, as.numeric)

sample_data_dis_x_lasso <- 
  sample_data_dis_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name))

sample_data_dis_x_lasso <- 
  apply(sample_data_dis_x_lasso, 2, as.numeric)

###use validation dataset for validation
best_lasso2 <-
  glmnet(as.matrix(sample_data_dis_x_lasso), 
         expected_date_remained_dis,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

predicted_y <-
  predict(
    object = best_lasso2,
    newx = as.matrix(sample_data_val_x_lasso),
    s = lasso_regression2$lambda.1se
    # type = "response"
  )

prediction_self <- 
  predict(object = best_lasso2, 
          newx = as.matrix(sample_data_dis_x_lasso), 
          s = lasso_regression2$lambda.1se)

plot(expected_date_remained_dis[,1], prediction_self)
abline(0, 1)

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
  labs(x = "Time to due (measured)", y = "Predicted error (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave("measured_vs_predicted_error.pdf", width = 7, height = 7)

linear_regression <- 
  lm(formula = expected_date_remained_dis[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y[,1] + coef(linear_regression)[1]

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
  labs(x = "Time to due (weeks, measured)", y = "Time to due (weeks, predicted)") +
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
    sample(1:nrow(sample_data_dis_x_lasso), 
           size = nrow(sample_data_dis_x_lasso), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_dis_x_lasso), dis_index)
  
  lasso_regression_best <- 
    glmnet(as.matrix(sample_data_dis_x_lasso[dis_index,]), 
           expected_date_remained_dis[dis_index,],
           family = "gaussian", 
           alpha = 1, 
           lambda = lasso_regression2$lambda.1se)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[dis_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = expected_date_remained_dis[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_dis_x_lasso[val_index,]),
      s = lasso_regression2$lambda.1se
    )[,1]
  
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
  dplyr::arrange(., y) %>% 
  dplyr::group_by(., y) %>% 
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
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  # scale_x_continuous(limits = c(10, 42)) +
  # scale_y_continuous(limits = c(10, 42)) +
  labs(x = "Time to due (weeks, measured)", 
       y = "Time to due (weeks, predicted)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))








