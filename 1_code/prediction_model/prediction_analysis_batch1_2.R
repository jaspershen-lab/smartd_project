library(tidyverse)
library(plyr)
library(igraph)

setwd("data_analysis20191015/prediction/identification_table/")

smartd_rplc <- 
  readr::read_csv("smartd_rplc.csv")

sample_data <- 
  smartd_rplc %>% 
  dplyr::select(., -(name:Database))

metabolite_tags <- 
  smartd_rplc %>% 
  dplyr::select(., name:Database)

sum(is.na(sample_data))

which(is.na(sample_data), arr.ind = TRUE)[,2] %>% 
  unique()

colnames(sample_data)[c(3, 4, 5, 338, 339, 340, 341, 342, 343, 344, 345)]

###"X100", "X101", "X102", "X114", "X115", "X151", "X152", "X159", "X85", "X92", "SFU65" 
##may have no positive data, so remove them from the dataset
sample_data <- 
  sample_data %>% 
  dplyr::select(., -c(X100, X101, X102, X114, X115, X151, X152, X159, X85, X92, SFU65))

###patient and sampel information
sfu1_148 <- 
  readr::read_csv("E:/project/smartD/patient information/SFU1-148_GA.csv")

patient_info <- 
  readr::read_csv("E:/project/smartD/data_analysis20190828/patient_info/patient_info.csv")

smartD_all346urine <- 
  readr::read_csv("E:/project/smartD/data_analysis20190828/patient_info/SmartD_all346urine.csv")

sfu1_148 <- 
  sfu1_148 %>% 
  mutate(subject_id = as.character(subject_id), visit = as.character(visit)) %>% 
  arrange(subject_id)

patient_info <- 
  patient_info %>% 
  mutate(Patient_ID = as.character(Patient_ID), Visit = as.character(Visit)) %>% 
  arrange(Patient_ID)

match(colnames(sample_data), sfu1_148$sample_id)

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
  smartd_rplc$name

# sample_data <- 
#   inner_join(x = sfu1_148[,-1],
#              y = sample_data, 
#              by = c("sample_id" = "Sample_ID"))

###sample_data contains batch 1 and batch 2 data, we should split them
load("E:/project/smartD/smartD_batch1_2/RPLC_xcms/POS/data_cleaning/smartd_rplc_pos_5")
sample_info <- 
  smartd_rplc_neg_5@sample.info

sample_data2 <-
  inner_join(x = sample_data,
             y = sample_info[, c("sample.name", "GA")],
             by = c("Sample_ID" = "sample.name"))

sample_data2 <- 
  sample_info %>% 
  plyr::dlply(., .variables = .(batch)) %>% 
  lapply(function(x)x$sample.name) %>% 
  lapply(function(x) sample_data2[sample_data2$Sample_ID %in% x ,])

##
sample_data_batch1 <- sample_data2[[1]]
sample_data_batch2 <- sample_data2[[2]]

####cluster
# library(pheatmap)
# sample_data %>% 
#   select(., -(subject_id:GA_week)) %>% 
#   pheatmap(., show_rownames = FALSE, 
#            show_colnames = FALSE)
  
####FCM cluster
# library(e1071)
####biomarker discovery and prediction
##use batch 1 data as discovery data
sample_data <- sample_data_batch1
# sample_data$subject_id %>% 
#   unique() %>% 
#   length()

##16 subjects in total
##using lasso regressio to select variables
library(glmnet)
###construct dataset for lasso regression
sample_data_x <- 
sample_data %>% 
  dplyr::select(-c(Sample_ID, GA)) %>% 
  as.matrix()

sample_data_y <- 
  sample_data %>% 
  dplyr::select(GA) %>% 
  as.matrix()

lasso_regression <-  
  glmnet(x = sample_data_x, y = sample_data_y, 
         family = "gaussian", 
         nlambda = 100, alpha = 1, 
         standardize = FALSE)

##each row is a model, DF is free degree, variable number is DF+1, and %Dev is R2. 
print(lasso_regression)


##get the coefficients of variables using one specific lambda value
coef(lasso_regression, 
     s = c(lasso_regression$lambda[40],0.1)) %>% 
  head()



plot(lasso_regression, xvar = "lambda", label = TRUE)


lasso_regression <-
  glmnet(
    x = sample_data_x,
    y = sample_data_y,
    family = "gaussian",
    nlambda = 100,
    alpha = 1,
    standardize = FALSE
  )

lasso_regression$a0

beta <- 
  lasso_regression$beta

lasso_regression$df
##lambdas for 100 models 
lasso_regression$lambda
##R2 for 100 models
lasso_regression$dev.ratio

###log lambda vs R2
data.frame(lambda = lasso_regression$lambda,
           dev.ratio = lasso_regression$dev.ratio,
           stringsAsFactors = FALSE) %>% 
  ggplot(aes(log(lambda), dev.ratio)) +
  labs(x = "log(lambda)", y = "Deviation ratio (%)") +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

###beta is the coefficients for ecah variable in different models
beta <-
  lasso_regression$beta %>%
  as.matrix() %>%
  t() %>%
  as_tibble() %>%
  mutate(lambda = lasso_regression$lambda) %>%
  gather(., key = "feature", value = "coef", -lambda)


library(ggsci)

beta %>% 
  ggplot(., aes(log(lambda), coef)) +
  geom_line(aes(colour = feature), show.legend = FALSE) +
  scale_x_continuous(position = "bottom",
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = c(-3,-2,-1,0,1), 
                                         name = "",
                                         labels = c(128, 107, 64, 21, 1))) +
  scale_colour_manual(values = colorRampPalette(pal_uchicago()(5))(600)) +
  labs(x = "Log Lambda", y = "Coefficients") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13)
    # plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
  ) 
  
grid::convertUnit(unit(13, "pt"), "mm", valueOnly=TRUE)

# plot(lasso_regression, xvar = "lambda", label = TRUE)
# plot(rev(pull(beta, colnames(beta)[484])), type = "l")
lasso_regression
# plot(lasso_regression, xvar = "lambda", label = TRUE)

##set alpha as 0, so there are 21 markers
lasso_regression2 <-
  cv.glmnet(
    x = sample_data_x,
    y = sample_data_y,
    family = "gaussian",
    type.measure = "mae",
    nfolds = 7,
    nlambda = 100,
    alpha = 1,
    standardize = FALSE,
    parallel = FALSE
  )
        
print(lasso_regression2)

lasso_regression2$lambda
lasso_regression2$cvm
lasso_regression2$cvsd
lasso_regression2$cvup
lasso_regression2$cvlo
lasso_regression2$nzero
lasso_regression2$name
lasso_regression2$glmnet.fit
lasso_regression2$lambda.min
lasso_regression2$lambda.1se


plot(lasso_regression2, xvar = "lambda", label = TRUE)

beta <-
  lasso_regression2$glmnet.fit$beta %>% 
  as.matrix() %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(lambda = lasso_regression2$lambda) %>% 
  gather(., key = "feature", value = "coef", -lambda)


plot(lasso_regression2$glmnet.fit, xvar = "lambda", label = TRUE)

beta %>% 
  ggplot(., aes(log(lambda), coef)) +
  geom_line(aes(colour = feature), show.legend = FALSE) +
  scale_x_continuous(position = "bottom",
                     sec.axis = sec_axis(trans = ~., 
                                         name = "", 
                                         breaks = c(-3, -2, -1, 0, 1),
                                         labels = c(129, 113, 62, 19, 3))) +
  scale_colour_manual(values = colorRampPalette(pal_uchicago()(5))(600)) +
  labs(x = "Log lambda", y = "Coefficients") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13)
    # plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
  ) 


cvm <-
  data.frame(
    lambda = lasso_regression2$lambda,
    df = lasso_regression2$glmnet.fit$df,
    cvm = lasso_regression2$cvm,
    cvup = lasso_regression2$cvup,
    cvlo = lasso_regression2$cvlo,
    stringsAsFactors = FALSE
  )

cvm %>% 
  ggplot(., aes(log(lambda), cvm)) +
  geom_vline(xintercept = log(c(lasso_regression2$lambda.min,
                                lasso_regression2$lambda.1se)),
             linetype = 2) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  scale_x_continuous(position = "bottom",
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = log(cvm$lambda)[seq(1, 100, by = 7)],
                                         labels = cvm$df[seq(1, 100, by = 7)],
                                         name = "")
                     ) +
  labs(x = "Log lambda", y = "Mean absolute error") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13)
    # plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
  )

plot(lasso_regression2)

# c(lasso_regression2$lambda.min,
#   lasso_regression2$lambda.1se)

##construct the best model
best_lasso <-
  glmnet(sample_data_x, 
         sample_data_y,
         family = "gaussian", 
         alpha = 1,
         standardize = FALSE,
         lambda = lasso_regression2$lambda.1se)


##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se) 

lasso_regression2$glmnet.fit$beta[,35] %>% 
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

####we should use the bootstrap method to validate our result
dim(sample_data_x)
predict_y <- vector(mode = "list", length = 100)
y <- vector(mode = "list", length = 100)
for(i in 1:100){
  cat(i, " ")
  dis_index <- 
    sample(1:nrow(sample_data_x), 
           size = nrow(sample_data_x), 
           replace = TRUE) %>% 
    unique() %>% 
    sort()
  
  val_index <-
    setdiff(1:nrow(sample_data_x), dis_index)
  
  lasso_regression_temp <-
    cv.glmnet(
      x = sample_data_x[dis_index,],
      y = sample_data_y[dis_index,],
      family = "gaussian",
      type.measure = "mse",
      nfolds = 7,
      nlambda = 100,
      alpha = 1,
      standardize = FALSE
    )
  
  best_lambda <- lasso_regression_temp$lambda.1se
  
  lasso_regression_best <- 
    glmnet(sample_data_x[dis_index,], 
           sample_data_y[dis_index,],
           family = "gaussian", 
           alpha = 1, 
           lambda = best_lambda)
  
  ##construct a new linear model to correct prediction and real value 
  prediction_self <- 
    predict(
      object = lasso_regression_best,
      newx = sample_data_x[dis_index,],
      s = best_lambda
      # type = "response"
    )[,1]
    
  
  plot(prediction_self, sample_data_y[dis_index,1])
  abline(0, 1)
  ###construct a new liner regression model to correct them
  linear_regression <- 
    lm(formula = sample_data_y[dis_index,1] ~ prediction_self)
  
  temp_predict_y <- 
    predict(
      object = lasso_regression_best,
      newx = sample_data_x[val_index,],
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
  
  # plot(sample_data_y[val_index,1],  predict_y[[i]])
  # 
  # plot(sample_data_y[val_index,1],  temp_predict_y)
  # abline(0, 1)
  
  
  # predict_y[[i]] <- 
  #   predict(
  #     object = liner_regression, 
  #     newx = predict_y[[i]]
  #     # type = "response"
  #   )
  
  # predict_y[[i]] <- 
  #   predict(
  #     object = lasso_regression_best,
  #     newx = sample_data_x[val_index,],
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

##markers
##markers
which(lasso_regression2$lambda == lasso_regression2$lambda.1se)
marker_lasso <-
lasso_regression2$glmnet.fit$beta[,21] %>% 
  tibble(name = rownames(lasso_regression2$glmnet.fit$beta),
         coef = .) %>%  
  filter(coef != 0)


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
               mtry = 7,
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
  select(., Sample_ID, GA, name = marker_rf[30]) %>% 
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
  dplyr::select(., one_of(marker_svm)) %>% 
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
  dplyr::select(., one_of(marker_nnet)) %>% 
  as.matrix()

sample_data_x_y_nnet <- 
  data.frame(sample_data_x_nnet, y = sample_data_y,
             stringsAsFactors = FALSE)


nnet_regression <- 
  neuralnet(formula = GA~.,
            data = sample_data_x_y_nnet, 
            hidden =  c(5, 3), 
            act.fct = "logistic",
            linear.output = TRUE)

nnet_regression8 <- 
  neuralnet(formula = GA~.,
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
