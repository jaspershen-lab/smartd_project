####for RF
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/")
###normal data
library(randomForest)

##load dataa
load("../../../sample_data_dis")
load("../../../sample_data_val")
load("../../../sample_data_dis_x")
load("../../../sample_data_val_x")
load("../../../metabolite_tags")

##############################################################################
####random forest
#############################################################################
library(randomForest)
##use boruta method
library(Boruta)
library(tidyverse)

sample_data_dis_y <- 
  sample_data_dis %>% 
  dplyr::select(GA) %>% 
  as.matrix()

sample_data_val_y <- 
  sample_data_val %>% 
  dplyr::select(GA) %>% 
  as.matrix()

marker_rf <- 
  readr::read_csv("marker_rf_final.csv")


info <-
  readxl::read_xlsx("E:/project/smartD/patient_information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

info <- 
  info %>% 
  filter(!is.na(`C/S`))

info <-
  info %>% 
  filter(`C/S` == "N")

idx_dis <- 
  which(sample_data_dis$Patient_ID %in% info$ID)


idx_val <- 
  which(sample_data_val$Patient_ID %in% info$ID)


sample_data_dis <- 
  sample_data_dis[idx_dis,]

sample_data_dis_x <- 
  sample_data_dis_x[idx_dis,]

sample_data_val <- 
  sample_data_val[idx_val,]

sample_data_val_x <- 
  sample_data_val_x[idx_val,]

unique(sample_data_dis$Patient_ID)

unique(sample_data_val$Patient_ID)


##remove CS
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


rf_regression <-
  randomForest(x = sample_data_dis_x_rf, 
               y = expected_date_remained_dis[,1], 
               replace = TRUE, 
               importance = TRUE,
               proximity = TRUE, 
               mtry = 2)


##validate in validation dataset
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

real_rmse <- 
  abs(expected_date_remained_val[,1] - predicted_y) %>% 
  mean()

real_r2 <- 
  summary(lm(formula = predicted_y~expected_date_remained_val[,1]))$adj.r.squared

fake_rmse <- vector(mode = "numeric", length = 100)
fake_r2 <- vector(mode = "numeric", length = 100)

for(i in 1:100){
  cat(i, " ")
  fake_dis_y <- 
    expected_date_remained_dis[,1][sample(1:nrow(expected_date_remained_dis))]
  fake_val_y <- 
    expected_date_remained_val[,1][sample(1:nrow(expected_date_remained_val))]
  
  rf_regression <-
    randomForest(x = sample_data_dis_x, 
                 y = fake_dis_y, 
                 replace = TRUE, 
                 importance = TRUE,
                 proximity = TRUE)
  ###use validation dataset for validation
  predicted_y <-
    predict(
      object = rf_regression,
      newdata = sample_data_val_x
      # type = "response"
    )
  
  fake_rmse[i] <- 
    abs(fake_val_y - predicted_y) %>% 
    mean()
  fake_r2[i] <- 
    summary(lm(formula = predicted_y~fake_val_y))$adj.r.squared
}

save(fake_rmse, file = "fake_rmse")
save(fake_r2, file = "fake_r2")

plot1 <- 
fake_rmse %>% 
  as_tibble() %>% 
ggplot(aes(value)) +
  # geom_histogram(aes(, y = ..density..), binwidth = 0.03,
  #                colour= "white",
  #                fill = "#8A9045FF") +
  geom_line(aes(x = value), colour = "#8A9045FF", size = 1, stat = "density") +
  geom_rug(aes(x = value), colour = "grey") +
  scale_x_continuous(limits = c(3, 7.3)) +
  labs(x = "RMSE", y = "Density") +
  theme_bw() +
  geom_vline(xintercept = real_rmse, linetype = 2, 
             colour = "#C16622FF", size = 1) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        plot.title = element_text(size = 15, hjust = 0.5)) +
  labs(title = "Permutation test (RMSE)", x = "RMSE",
       ylab = "Density")
plot1
library(grid)
vie <-
  viewport(x = 0.4, y = 0.5, width = 0.4, height = 0.5)

plot2 <- 
  fake_rmse %>% 
  as_tibble() %>% 
  ggplot(aes(value)) +
  # geom_histogram(aes(, y = ..density..), binwidth = 0.03,
  #                colour= "white",
  #                fill = "#8A9045FF") +
  geom_line(aes(x = value), colour = "#8A9045FF", size = 1, stat = "density") +
  geom_rug(aes(x = value), colour = "grey") +
  # scale_x_continuous(limits = c(3, 7.3)) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        plot.title = element_text(size = 15, hjust = 0.5)) +
  labs(x = "",
       y = "")

print(plot2, vp = vie)


plot1 <- 
fake_r2 %>% 
  as_tibble() %>% 
  ggplot(aes(value)) +
  # geom_histogram(aes(, y = ..density..), binwidth = 0.001,
  #                # colour= "white",
  #                fill = "#155F83FF") +
  geom_line(aes(x = value), colour = "#155F83FF", size = 1, stat = "density") +
  geom_rug(aes(x = value), colour = "grey") +
  scale_x_continuous(limits = c(-0.01, 0.8)) +
  theme_bw() +
  geom_vline(xintercept = real_r2, linetype = 2, 
             colour = "#C16622FF", size = 1) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15, hjust = 0.5)) +
  labs(title = "Permutation test (Adjusted R2)", 
       # title = expression(Permutation~test~(Adjusted~R^2)), 
       # x = expression(Adjusted~R^2),
       x = "Adjusted R2",
       y = "Density")

plot1
library(grid)
vie <-
  viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.5)

plot2 <- 
  fake_r2 %>% 
  as_tibble() %>% 
  ggplot(aes(value)) +
  # geom_histogram(aes(, y = ..density..), binwidth = 0.001,
  #                # colour= "white",
  #                fill = "#155F83FF") +
  geom_line(aes(x = value), colour = "#155F83FF", size = 1, stat = "density") +
  geom_rug(aes(x = value), colour = "grey") +
  scale_x_continuous(limits = c(-0.03, 0.05)) +
  theme_bw() +
  # geom_vline(xintercept = real_r2, linetype = 2, 
  #            colour = "#C16622FF", size = 1) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 15, hjust = 0.5)) +
  labs(x = "",
       y = "")

print(plot2, vp = vie)

###fake_rmse is a normal distributation
plot(density(fake_rmse))
mean(fake_rmse)
sd(fake_rmse)
real_rmse
##calculate P(x < 3.96)
pnorm(q = real_rmse, mean = mean(fake_rmse), sd = sd(fake_rmse))

plot(density(fake_r2))
mean(fake_r2)
sd(fake_r2)
real_r2

1 - pnorm(q = real_r2, mean = mean(fake_r2), sd = sd(fake_r2))



