setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("data_analysis20200108/prediction/metabolites/")

load("sample_data_dis")
load("sample_data_val")

setwd(r4projects::get_project_wd())
setwd("patient information/")
library(tidyverse)
sample_info <- readr::read_csv("sample_info_191021.csv")

temp_data <- 
  data.frame(subject_id = sample_info$Patient_ID, 
             sample_id = sample_info$Sample_ID,
             location = "UCSF",
             stringsAsFactors = FALSE)


temp_data <-
  temp_data %>% 
  mutate(batch = 
             case_when(subject_id %in% sample_data_dis$Patient_ID ~ 1,
                       subject_id %in% sample_data_val$Patient_ID ~ 2),
         TRUE ~ NA
         )

##remove data

temp_data <-
  temp_data %>% 
  dplyr::mutate(
    cleaning = case_when(
      sample_id %in% sample_data_dis$Sample_ID ~ "Remain", 
      sample_id %in% sample_data_val$Sample_ID ~ "Remain", 
      TRUE ~ "Removal"))


##discovery and validation
temp_data <-
  temp_data %>% 
  dplyr::mutate(
    prediction = case_when(
      sample_id %in% sample_data_dis$Sample_ID ~ "Discovery", 
      sample_id %in% sample_data_val$Sample_ID ~ "Validation", 
      TRUE ~ "NA")
    )


sample_information <- 
  temp_data %>% 
  dplyr::left_join(sample_info[,c(2,5)], by = c("sample_id" = "Sample_ID"))

xlsx::write.xlsx(sample_information, "sample_information.xlsx")


temp_data$Freq <- 1

temp_data2 <-  temp_data[,-c(1,2,7)]

ord2 <- list(NULL, NULL, NULL, nrow(temp_data2) : 1)

# ord <- list(NULL, with(data = tit3d, expr = order(Sex, Survived)), NULL)

temp_data2$prediction <- 
  factor(temp_data2$prediction, levels = c("Validation", "Discovery", "NA"))

alluvial(
  temp_data2,
  freq = temp_data$Freq,
  col = ifelse(temp_data$batch == "1", "#FFA3194C", "#155F834C"),
  # border = ifelse(temp_data$batch == "1", "#FFA3197F", "#155F837F"),
  border = NA,
  # hide = temp_data$prediction == "NA",
  cex = 1, 
  alpha = 0.5, 
  layer = temp_data$batch == "1"
  # ordering = ord2
  # gap.width = 1,  
)









library(alluvial)
tit <- as.data.frame(Titanic, stringsAsFactors = FALSE)
head(tit)

alluvial(
  tit[, 1:4],
  freq = tit$Freq,
  col = ifelse(tit$Survived == "Yes", "orange", "grey"),
  border = ifelse(tit$Survived == "Yes", "orange", "grey"),
  hide = tit$Freq == 0,
  cex = 0.7
)


library(alluvial)



tit <- as.data.frame(Titanic)

# 2d
tit2d <- aggregate( Freq ~ Class + Survived, data=tit, sum)
alluvial( tit2d[,1:2], freq=tit2d$Freq, xw=0.0, alpha=0.8,
          gap.width=0.1, col= "steelblue", border="white",
          layer = tit2d$Survived != "Yes" )

alluvial(tit2d[,1:2], 
         freq=tit2d$Freq, 
          hide=tit2d$Freq < 150,
          xw = 0.0, alpha=0.8,
          gap.width = 0.1, 
         col= "steelblue", border="white"
          # layer = tit2d$Survived != "Yes" 
         )

# 3d
tit3d <- aggregate( Freq ~ Class + Sex + Survived, data=tit, sum)

alluvial(tit3d[,1:3], freq=tit3d$Freq, alpha=1, xw=0.2,
         col=ifelse( tit3d$Survived == "No", "red", "gray"),
         layer = tit3d$Sex != "Female",
         border="white")


# 4d
alluvial( tit[,1:4], freq=tit$Freq, border=NA,
          hide = tit$Freq < quantile(tit$Freq, .50),
          col=ifelse( tit$Class == "3rd" & tit$Sex == "Male", "red", "gray") )

# 3d example with custom ordering
# Reorder "Sex" axis according to survival status
ord <- list(NULL, with(tit3d, order(Sex, Survived)), NULL)
alluvial(tit3d[,1:3], freq=tit3d$Freq, alpha=1, xw=0.2,
         col=ifelse( tit3d$Survived == "No", "red", "gray"),
         layer = tit3d$Sex != "Female",
         border="white", ordering = ord)

# Possible blocks options
for (blocks in c(TRUE, FALSE, "bookends")) {
  
  # Elaborate alluvial diagram from main examples file
  # temp_data <-  tit[, 1:4]
  alluvial( tit[, 1:4], freq = tit$Freq, border = NA,
            hide = tit$Freq < quantile(tit$Freq, .50),
            col = ifelse( tit$Class == "3rd" & tit$Sex == "Male",
                          "red", "gray" ),
            blocks = TRUE )
}


# Data returned
x <- alluvial( tit2d[,1:2], freq=tit2d$Freq, xw=0.0, alpha=0.8,
               gap.width=0.1, col= "steelblue", border="white",
               layer = tit2d$Survived != "Yes" )
points( rep(1, 16), x$endpoints[[1]], col="green")
points( rep(2, 16), x$endpoints[[2]], col="blue")



#






library(webr)
library(moonBook)
library(ggplot2)
PieDonut(acs,aes(pies=Dx,donuts=smoking))
