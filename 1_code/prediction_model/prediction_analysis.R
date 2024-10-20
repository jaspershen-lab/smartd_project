library(tidyverse)
library(plyr)
library(igraph)
library(dplyr)
setwd("data_analysis20190828/prediction/identification_table/")

###This is batch data set
smartd_rplc_batch1 <- 
  readr::read_csv("smartd_rplc_batch1.csv")

##remove the tags of metabolites
sample_data <- 
  smartd_rplc_batch1 %>% 
  dplyr::select(., -(name:Database))

###metabolite tags
metabolite_tags <- 
  smartd_rplc_batch1 %>% 
  dplyr::select(., name:Database)

##if there are NAs or not
sum(is.na(sample_data))

##what samples have NAs
which(is.na(sample_data), arr.ind = TRUE)[,2] %>% 
  unique()

colnames(sample_data)[c(146, 147, 148)]

###"SFU65" "SFU74" "SFU91" may have no positive data, so remove them from the dataset
sample_data <- 
  sample_data %>% 
  dplyr::select(., -c(SFU65, SFU74, SFU91))

###patient and sampel information
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

####log 10 and scale
sample_data <- 
  log(sample_data, 10) %>% 
  as.data.frame(sample_data)

sample_data <-
  apply(sample_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

sample_data <- 
  sample_data %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

colnames(sample_data)[-1] <- 
  smartd_rplc_batch1$name

###add the patient information to the sample_data
sample_data <- 
  right_join(x = sample_info_191021[,c(1:5)], sample_data, 
             by = "Sample_ID")

##remove GA is NA
sample_data <- 
  sample_data %>% 
  filter(!is.na(GA))

####cluster
# library(pheatmap)
# sample_data %>% 
#   dplyr::select(., -(subject_id:GA_week)) %>% 
#   pheatmap(., show_rownames = FALSE, 
#            show_colnames = FALSE)
#   
####FCM cluster
library(e1071)

####biomarker discovery and prediction
sample_data$Patient_ID %>% 
  unique() %>% 
  length()

##16 subjects in total
###################################LASSO regression
##using lasso regressio to select variables
library(glmnet)


###construct dataset for lasso regression
sample_data_x <- 
  sample_data %>% 
  dplyr::select(-(Patient_ID:GA)) %>% 
  as.matrix()


sample_data_y <- 
  sample_data %>% 
  dplyr::select(GA) %>% 
  as.matrix()

###use all the samples to do lasso regression
lasso_regression <-  
  glmnet(x = sample_data_x,
         y = sample_data_y, 
         family = "gaussian", 
         nlambda = 100,
         alpha = 1,
         standardize = FALSE)

##each row is a model, df is the variable number in this model and %Dev is R2 and lambda is the 
## penality for the model, bigger lambda, less variables in this model
print(lasso_regression)

coef(lasso_regression, 
     s = c(lasso_regression$lambda[40],0.1)) %>% 
  head()

###the coefficients of each variable change under different lambda
plot(lasso_regression, xvar = "lambda", label = TRUE)

lasso_regression$a0

###the coefficients of each variable change under different lambda
beta <- 
  lasso_regression$beta

lasso_regression$df
lasso_regression$lambda
lasso_regression$dev.ratio

plotLambdaVSdeviation(lasso_regression)

plotLambdaVScoefficients(object = lasso_regression)

grid::convertUnit(unit(13, "pt"), "mm", valueOnly=TRUE)

# plot(lasso_regression, xvar = "lambda", label = TRUE)
# plot(rev(pull(beta, colnames(beta)[484])), type = "l")
###we should use cross validation to get the best model and get the markers

lasso_regression2 <-
  cv.glmnet(
    x = sample_data_x,
    y = sample_data_y,
    family = "gaussian",
    type.measure = "mae",
    nfolds = 7,
    nlambda = 100,
    alpha = 1,
    standardize = FALSE
  )
        
lasso_regression2$lambda
lasso_regression2$cvm
lasso_regression2$cvsd
lasso_regression2$cvup
lasso_regression2$cvlo
lasso_regression2$nzero
lasso_regression2$name
lasso_regression2$name
lasso_regression2$glmnet.fit
lasso_regression2$lambda.min
lasso_regression2$lambda.1se

plotLambdaVScoefficients(object = lasso_regression2$glmnet.fit)

plotLambdaVSerror(object = lasso_regression2)

# plot(lasso_regression2)

# c(lasso_regression2$lambda.min,
#   lasso_regression2$lambda.1se)

##construct the best model
best_lasso <-
  glmnet(sample_data_x, 
         sample_data_y,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se)

marker_lasso <- 
  lasso_regression2$glmnet.fit$beta[,48] %>% 
  tibble(name = rownames(lasso_regression2$glmnet.fit$beta),
            coef = .) %>% 
  filter(coef != 0) 

predicted_y <-
  predict(
    object = best_lasso,
    newx = sample_data_x,
    s = lasso_regression2$lambda.min
    # type = "response"
  )

plot(
sample_data_y[,1],
predicted_y[,1]
)

abline(0, 1)

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


####
##validate in batch 2 dataset
smartd_rplc_batch2 <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS_NEG/ms1_data_rplc_batch2.csv")

##remove QC and blank form batch 2 dataset
smartd_rplc_batch2 <-
  smartd_rplc_batch2 %>% 
  select(-contains("blk")) %>% 
  select(-contains("QC"))


sample_removed <-
  which(is.na(smartd_rplc_batch2), arr.ind = TRUE)[,2] %>% 
  unique() %>% 
  `[`(colnames(smartd_rplc_batch2), .)

##These samples may have positive or negative samples but don't have, so we should remove them
smartd_rplc_batch2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(-sample_removed)

##remove the tags of metabolites
sample_data2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(-(name:rt))

###metabolite tags
metabolite_tags2 <- 
  smartd_rplc_batch2 %>% 
  dplyr::select(., name:rt)

####log 10 and scale
sample_data2 <- 
  log(sample_data2, 10) %>% 
  as.data.frame()

sample_data2 <-
  apply(sample_data2, 1, function(x){
    (x - mean(x))/sd(x)
  })

sample_data2 <- 
  sample_data2 %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_ID")

colnames(sample_data2)[-1] <- 
  smartd_rplc_batch2$name

####add GA information
sample_data2 <- 
  right_join(x = sample_info_191021[,c(1:5)], sample_data2, 
             by = "Sample_ID")

###remove GA is 0 and modify the GA
sample_data2 <- 
  sample_data2 %>% 
  filter(!is.na(GA))

###construct dataset for lasso regression
sample_data2_x <- 
  sample_data2 %>% 
  dplyr::select(-c(Patient_ID, GA)) %>% 
  as.matrix()

sample_data2_y <- 
  sample_data2 %>% 
  dplyr::select(GA) %>% 
  as.matrix()

##read the batch 1 and batch 2 match information
rplc_pos_batch1_batch2_matched_info <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS/rplc_pos_batch1_batch2_matched_info.csv")

rplc_pos_batch1_batch2_matched_info <- 
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS/rplc_pos_batch1_batch2_matched_info.csv")

rplc_matched_info <- rbind(rplc_pos_batch1_batch2_matched_info[,-1],
                           rplc_neg_batch1_batch2_matched_info)

rplc_matched_info %>% 
  filter(batch1 %in% marker_lasso$name) %>% 
  pull(batch2)

marker_lasso <- 
  marker_lasso %>% 
  dplyr::rename("name1" = name) %>% 
  mutate(., name2 = 
           match(name1, rplc_matched_info$batch1) %>% 
           `[`(rplc_matched_info$batch2, .)
  )


marker_lasso <- 
  marker_lasso %>% 
  filter(!is.na(name2))

####we should use the bootstrap method to validate our result
dim(sample_data_x)

sample_data_x_lasso <- 
  sample_data_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name1))


# sample_data_x_lasso <- 
#   apply(sample_data_x_lasso, 2, as.numeric)

###use validation dataset (batch 2) for validation
sample_data2_x_lasso <- 
  sample_data2_x %>% 
  as_tibble() %>% 
  dplyr::select(one_of(marker_lasso$name2))

colnames(sample_data2_x_lasso) <-
  marker_lasso$name1

sample_data2_x_lasso <- 
  apply(sample_data2_x_lasso, 2, as.numeric)

best_lasso2 <-
  glmnet(as.matrix(sample_data_x_lasso), 
         sample_data_y,
         family = "gaussian", 
         alpha = 1, 
         lambda = lasso_regression2$lambda.1se, 
         standardize = FALSE)

predicted_y <-
  predict(
    object = best_lasso2,
    newx = as.matrix(sample_data2_x_lasso),
    s = lasso_regression2$lambda.1se
    # type = "response"
  )

plot(
  sample_data2_y[,1],
  predicted_y[,1]
)

abline(0, 1)


prediction_self <- 
  predict(object = best_lasso2, 
          newx = as.matrix(sample_data_x_lasso), 
          s = lasso_regression2$lambda.1se)

plot(prediction_self, sample_data_y[,1])
abline(0, 1)

linear_regression <- 
  lm(formula = sample_data_y[,1] ~ prediction_self)

predicted_y2 <- 
  coef(linear_regression)[2] * predicted_y[,1] + coef(linear_regression)[1]


abs(sample_data2_y[,1] - predicted_y2) %>% 
  mean()

summary(lm(formula = predicted_y2~sample_data2_y[,1]))



plot(predicted_y[,1], sample_data2_y[,1])
abline(0, 1)

plot(predicted_y[,1], predicted_y2)
abline(0, 1)

predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)

for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_x_lasso), 
           size = nrow(sample_data_x_lasso), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_x_lasso), dis_index)
  
  # lasso_regression_temp <-
  #   cv.glmnet(
  #     x = as.matrix(sample_data_x_lasso[dis_index,]),
  #     y = sample_data_y[dis_index,],
  #     family = "gaussian",
  #     type.measure = "mse",
  #     nfolds = 7,
  #     nlambda = 50,
  #     alpha = 1,
  #     standardize = FALSE
  #   )
  # 
  # best_lambda <- lasso_regression_temp$lambda.1se
  
  lasso_regression_best <- 
    glmnet(as.matrix(sample_data_x_lasso[dis_index,]), 
           sample_data_y[dis_index,],
           family = "gaussian", 
           alpha = 1, 
           lambda = lasso_regression2$lambda.1se)
  
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_x_lasso[dis_index,]),
      s = best_lambda
      # type = "response"
    )[,1]
  
  
  # plot(prediction_self, sample_data_y[dis_index,1])
  # abline(0, 1)
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = lasso_regression_best,
      newx = as.matrix(sample_data_x_lasso[val_index,]),
      s = best_lambda
      # type = "response"
    )[,1]
  
  temp_predict_y2 <- 
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  # plot(sample_data_y[val_index,1], temp_predict_y, xlim = c(10, 40), ylim = c(10,40))
  # abline(0, 1)
  # 
  # plot(sample_data_y[val_index,1], temp_predict_y2, xlim = c(10, 40), ylim = c(10,40))
  # abline(0, 1)
  # 
  # plot(temp_predict_y, temp_predict_y2, xlim = c(10, 40), ylim = c(10,40))
  
  # sum(abs(sample_data_y[val_index,1] - temp_predict_y))
  # sum(abs(sample_data_y[val_index,1] - temp_predict_y2))
  
  # predict(
  #   object = linear_regression, 
  #   newx = temp_predict_y
  #   # type = "response"
  # )
  
  predict_y[[i]] <- 
    temp_predict_y2
  
  # predict_y[[i]] <- 
  #   predict(
  #     object = lasso_regression_best,
  #     newx = as.matrix(sample_data_x_lasso[val_index,]),
  #     s = best_lambda
  #     # type = "response"
  #   )[,1]
    
  y[[i]] <- 
    sample_data_y[val_index,1]
    
}

plot(unlist(y), unlist(predict_y))
abline(0,1)

(unlist(y) - unlist(predict_y))^2 %>% 
  mean()

abs(unlist(y) - unlist(predict_y)) %>% 
  mean()

summary(lm(formula = unlist(predict_y)~unlist(y)))

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


(temp$y - temp$mean)^2 %>% 
  mean()

abs(temp$y - temp$mean) %>% 
  mean()

summary(lm(formula = temp$y~temp$mean))

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


###try marker in validation dataset (batch2)







#####use other 
##random forest
###feature selection
library(randomForest)
##use boruta method
library(Boruta)
boruta_test <- 
  Boruta(x = sample_data_x,
         y = sample_data_y, 
         doTrace = 3, 
         holdHistory = TRUE)

plot(boruta_test)

boruta_test

marker_rf <- 
  boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>% 
  names() %>% 
  sort()

boruta_test$ImpHistory[1:16,] %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::transmute(., mean = apply(., 1, mean),
                   sd = apply(., 1, sd),
                   ymax = mean + sd, 
                   ymin = mean - sd) %>% 
  mutate(name = colnames(boruta_test$ImpHistory)) %>% 
  filter(name %in% marker_rf) %>% 
  left_join(., metabolite_tags, by = "name") %>% 
  arrange(., mean) %>% 
  ggplot(.,aes(x = factor(Compound.name, Compound.name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF", width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))
  

# plot(boruta_test)

####parameter tunning
sample_data_x_rf <- 
  sample_data_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_rf)) %>% 
  as.matrix()

fgl.res <- tuneRF(sample_data_x_rf, 
                  sample_data_y[,1], 
                  stepFactor = 1.5, 
                  trace = TRUE, plot = TRUE)


fgl.res %>% 
  as.data.frame() %>% 
  ggplot(aes(mtry, y = OOBError)) +
  geom_line(colour = "#155F83FF") +
  geom_point(size = 4, colour = "#FFA319FF") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))


rf_regression <-
  randomForest(x = sample_data_x_rf, 
               y = sample_data_y[,1], 
               replace = TRUE, 
               importance = TRUE,
               proximity = TRUE)




predict_y <- vector(mode = "list", length = 1000)
y <- vector(mode = "list", length = 1000)
for(i in 1:1000){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_x_rf), 
           size = nrow(sample_data_x_rf), 
           replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_x_rf), dis_index)
  
  rf_regression_temp <-
    randomForest(x = sample_data_x_rf[dis_index,], 
                 y = sample_data_y[dis_index,1], 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE)
  
  predict_y[[i]] <- 
    predict(
      object = rf_regression_temp,
      newdata = sample_data_x_rf[val_index,]
      # type = "response"
    )
  
  y[[i]] <- 
    sample_data_y[val_index,1]
}


(unlist(y) - unlist(predict_y))^2 %>% 
  mean()

abs(unlist(y) - unlist(predict_y)) %>% 
  mean()

summary(lm(formula = unlist(predict_y)~unlist(y)))


data.frame(marker_rf) %>% 
left_join(., metabolite_tags, by = c("marker_rf" = "name"))

plot(unlist(y), unlist(predict_y))


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


rf_regression$call
rf_regression$type
rf_regression$predicted
rf_regression$mse
rf_regression$rsq

rf_regression$importance %>% 
  as.data.frame() %>% 
  rownames_to_column(., "name") %>% 
  arrange(.,desc(`%IncMSE`)) %>% 
  `[`(.,1:length(marker1),) %>% 
  pull(., name)
  
varImpPlot(x = rf_regression)
plot(rf_regression)


sample_data %>% 
  select(., subject_id, GA_week, name = marker_rf[30]) %>% 
  ggplot(., aes(GA_week, name, colour = subject_id)) +
  geom_point() +
  geom_line()


sample_data_x_rf %>% 
  pheatmap::pheatmap()



####SVM regression
library(e1071)

##feature selection
library(caret)
library(penalizedSVM)

svm_test <- 
  caret::rfe(x = sample_data_x, 
           y = sample_data_y[,1], 
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


####parameter tunning
sample_data_x_svm <- 
  sample_data_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_svm)) %>% 
  as.matrix()



rbf_tune <- 
  tune.svm(x = sample_data_x_svm, 
           y = sample_data_y[,1], 
           kernel = "radial",
           gamma = 2^(-10:-3), 
           cost = 2^(-5:4),
           tunecontrol = tune.control(sampling = "cross", 
                                      cross = 7, 
                                      best.model = TRUE, 
                                      performances = TRUE)
           )


plot(x = rbf_tune, xlab = "Gamma", main = "", 
     cex.lab = 1.5, cex.axis = 1.3)

rbf_tune$best.performance

par(xpd = FALSE)
points(x = rbf_tune$best.parameters[1,1], 
       y = rbf_tune$best.parameters[1,2], 
       col = "#FFA319FF", cex = 1.5, pch = 19)

abline(v = rbf_tune$best.parameters[1,1],
       col = "#FFA319FF")

abline(h = rbf_tune$best.parameters[1,2],
       col = "#FFA319FF")

rbf_tune$best.parameters


predict_y <- vector(mode = "list", length = 1000)
y <- vector(mode = "list", length = 1000)
for(i in 1:1000){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_x_svm), 
           size = nrow(sample_data_x_svm), 
           replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_x_svm), dis_index)
  
  svm_regression_temp <-
    svm(x = sample_data_x_svm[dis_index,],
        y = sample_data_y[dis_index, 1],
        scale = FALSE, kernel = "radial", 
        cost = rbf_tune$best.parameters[1,2],
        gamma = rbf_tune$best.parameters[1,1]
        )
  
  predict_y[[i]] <- 
    predict(
      object = svm_regression_temp,
      newdata = sample_data_x_svm[val_index,]
      # type = "response"
    )
  
  y[[i]] <- 
    sample_data_y[val_index,1]
}


(unlist(y) - unlist(predict_y))^2 %>% 
  mean()

abs(unlist(y) - unlist(predict_y)) %>% 
  mean()

summary(lm(formula = unlist(predict_y)~unlist(y)))


data.frame(marker_rf) %>% 
  left_join(., metabolite_tags, by = c("marker_rf" = "name"))

plot(unlist(y), unlist(predict_y))


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


####ANN in R
library(neuralnet)
library(nnet)
library(caret)


##feature selection
model.nn <- train(x = sample_data_x,
                  y = sample_data_y[,1],
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


importance[1:50,] %>% 
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
  
  
marker_nnet <- 
  importance[1:50,]$name
  


####parameter tunning
sample_data_x_nnet <- 
  sample_data_x %>% 
  as.data.frame() %>% 
  select(., one_of(marker_nnet)) %>% 
  as.matrix()

sample_data_x_y_nnet <- 
  data.frame(sample_data_x_nnet, y = sample_data_y,
             stringsAsFactors = FALSE)


nnet_regression <- 
  neuralnet(formula = GA_week~.,
            data = sample_data_x_y_nnet, 
            hidden =  c(5, 3), 
            act.fct = "logistic",
            linear.output = TRUE)

nnet_regression8 <- 
  neuralnet(formula = GA_week~.,
            data = sample_data_x_y_nnet, 
            hidden = 8, 
            act.fct = "logistic",
            linear.output = TRUE)


plot(nnet_regression)

library(NeuralNetTools)

plot <- 
  garson(nnet_regression8)

plot + 
  coord_flip()


neuralnet::compute(nnet_regression, sample_data_x_nnet)$net.result %>% 
  plot(sample_data_y)


predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)
for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_x), size = nrow(sample_data_x), replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_x), dis_index)
  
  nnet_regression_temp <-
    neuralnet(formula = GA_week~.,
              data = data.frame(sample_data_x_nnet,
                                sample_data_y)[dis_index,],
              hidden = c(40, 20),
              act.fct = "logistic",
              linear.output = TRUE)
  
  # nnet_regression_temp <-
  #   train(
  #     x = sample_data_x_nnet[dis_index, ],
  #     y = sample_data_y[dis_index, 1],
  #     trControl = trainControl(method = "cv", number = 7),
  #     method = "nnet"
  #   )
    
  predict_y[[i]] <-
    neuralnet::compute(nnet_regression_temp,
                       sample_data_x_nnet[val_index,])$net.result
  
  # predict_y[[i]] <- 
  #   predict(object = nnet_regression_temp, 
  #           newdata = sample_data_x_nnet[val_index,])[,1]

  y[[i]] <- 
    sample_data_y[val_index,1]
  
}

plot(unlist(y), unlist(predict_y))
abline(0,1)

(unlist(y) - unlist(predict_y))^2 %>% 
  mean()

summary(lm(formula = unlist(predict_y)~unlist(y)))

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
