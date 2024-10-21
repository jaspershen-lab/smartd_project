setwd("E:/project/smartD/patient_information/")
library(tidyverse)
library(ggplot2)

###samnple ID and urine ID
##the first column is sample ID, and the second column is Urine ID. Urine is unique, one person one visit has one 
##Urine sample
sample_id_198_urine <- 
  readxl::read_xlsx(path = "198urine_SampleID_original IDs_RPLC.xlsx")

head(sample_id_198_urine)

##urine ID and patient ID, one patient may have multiple urine samples
##The first column (Urine ID) is urine ID, the second column (PTID) is patient ID.
smartD_urine_info <- 
  readr::read_csv(file = "SmartD Urine IDs-part1-confirmed.csv")

head(smartD_urine_info)

####how many samples for each patient
smartD_urine_info %>% 
  group_by(PTID) %>% 
  summarise(number = n()) %>% 
  arrange(number) %>% 
  mutate(PTID2 = factor(PTID, levels = PTID)) %>% 
  ggplot(aes(PTID2, number)) +
  geom_point() +
  scale_x_discrete() +
  theme_bw() +
  coord_flip()


###20 patients in total
smartD_urine_info$PTID %>% unique %>% 
  length()


smartD_urine_info %>% 
  group_by(PTID) %>% 
  summarise(n())


###patient ID and clinical information
##The first column (participant_is) is the patient ID. The second column (redcap_event_name) is the Visit
final_info <- 
  readr::read_csv(file = "FINALSMARTDiaphragm2_DATA_2019-07-01_1615.csv")

###final_info has 51 samples
final_info %>% 
  group_by(participant_id) %>% 
  summarise(n()) %>% 
  dim()

####unify the name of different files
colnames(sample_id_198_urine) <- c("Sample_ID", "Urine_ID")
colnames(smartD_urine_info)[1:2] <- c("Urine_ID", "Patient_ID")
###just make them all upper
smartD_urine_info$Patient_ID <- stringr::str_to_upper(smartD_urine_info$Patient_ID)
colnames(final_info)[1] <- c("Patient_ID")
final_info$Patient_ID <- stringr::str_to_upper(final_info$Patient_ID)

anti_join(sample_id_198_urine, smartD_urine_info, by = "Urine_ID") %>% 
  pull("Sample_ID")

# combine them according to Urine_ID
sample_id_198_urine$Urine_ID
smartD_urine_info$Urine_ID

setdiff(sample_id_198_urine$Urine_ID, smartD_urine_info$Urine_ID)
setdiff(smartD_urine_info$Urine_ID, sample_id_198_urine$Urine_ID)

##use the full_join
patient_info <- 
  dplyr::full_join(sample_id_198_urine, smartD_urine_info,
                   by = "Urine_ID")

patient_info <- 
  patient_info %>% 
  select(Patient_ID, everything()) %>% 
  arrange(Patient_ID, Visit)
  
##remove undefied columns
patient_info <-
  patient_info %>%
  select(., -matches("X[0-9]{1,2}")) %>% 
  arrange(Patient_ID, `Visit GA`) 
  
  
patient_info %>% 
  group_by(Patient_ID) %>% 
  summarise(n()) %>% 
  dim()
  
final_info %>% 
  group_by(Patient_ID) %>% 
  summarise(n()) %>% 
  dim()

# final_info <-
# final_info %>% 
#   arrange(Patient_ID) %>% 
#   group_by(Patient_ID) %>% 
#   filter(redcap_event_name == "study_visit_01_arm_1") %>% 
#   ungroup()
  

final_info$redcap_event_name <- 
  final_info$redcap_event_name %>% 
  stringr::str_replace_all("study_visit_|_arm_1", "") %>% 
  as.numeric()


# final_info <- 
patient_info <-
patient_info %>% 
  left_join(final_info, by = c("Patient_ID", "Visit" = "redcap_event_name"))


patient_info1 <- 
  patient_info %>% 
  filter(!is.na(Patient_ID), !is.na(`Visit GA`))


colour <- rep(NA, nrow(patient_info1))
ga <- round(patient_info1$`Visit GA`)
colour[ga < 10] <- "0-10"
colour[ga >= 10 & ga < 20] <- "10-20"
colour[ga >= 20 & ga < 30] <- "20-30"
colour[ga >= 30 & ga <= 40] <- "30-40"
patient_info1$`Visit GA` <- ga

patient_info1 <-
  patient_info1 %>%
mutate(colour = colour)

temp <- 
  patient_info1 %>% 
  group_by(Patient_ID) %>% 
  summarise(point = n()) %>% 
  arrange(point)

patient_info1$Patient_ID <- factor(patient_info1$Patient_ID, 
                                   levels = temp$Patient_ID)
(
plot1 <- 
patient_info1 %>% 
  ggplot(aes(y = Patient_ID, x = round(`Visit GA`))) +
  geom_vline(xintercept = seq(20, 40, by = 10), colour = "grey") +
  scale_x_continuous(breaks = c(10, 20, 30, 40), labels = c(10, 20, 30, 40)) +
  geom_point(aes(colour = colour), size = 3) +
  # guides(colour = guide_legend(title = "GA vist (week)",
  #                              title.position = "top", title.hjust = 0.5)) +
  ggsci::scale_color_jama() +
  theme_bw() +
  # coord_flip() +
  labs(x = "GA vist (weeks)", y ="Individual") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position = "none")
)


(plot2 <- 
  ggplot(patient_info1, aes(x = Patient_ID)) +
  geom_bar(colour = "skyblue", fill = "skyblue") +
    labs(x = "", y = "Tie point number") +
  theme_bw() +
    # geom_text(aes(x = Patient_ID, y = temp$point, label = temp$point)) +
    annotate(geom = "text", x = rep(1:20), 
             y = temp$point/2, label = temp$point, colour = "white", size = 5) +
    coord_flip() +
    # scale_x_discrete(label = "")+
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  )


temp2 <- 
patient_info1 %>% 
  group_by(`Visit GA`) %>% 
  summarise(People = n())

(plot3 <- 
    ggplot(temp2, aes(x = `Visit GA`, y = People)) +
    geom_bar(stat = "identity", colour = "#ED00007F", fill = "#ED00007F") +
    labs(x = "", y = "Panticipant number") +
    theme_bw() +
    annotate(geom = "text", x = temp2$`Visit GA`, 
             y = temp2$People/2, 
             label = temp2$People, colour = "white", size = 5) +
    # scale_x_discrete(label = "")+
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
)





library(gridExtra)


# grid.arrange(plot1, plot2, nrow = 1)

library(customLayout)

# lay1 <- lay_new(mat = matrix(1:2, ncol = 2), widths = c(3,1))
lay <- lay_new(mat = matrix(1:4, ncol = 2), widths = c(3,1), heights = c(1,3))
lay_show(lay)
# lay <- lay_bind_row(x = lay2, y = lay1, heights = c(1,3))

lay_grid(grobs = list(plot3, plot1, plot2, plot2), lay = lay)




patient_info %>% 
  filter(!is.na(Patient_ID), !is.na(`Visit GA`)) %>% 
  group_by(Patient_ID) %>% 
  summarise(number = n()) %>% 
  ungroup() %>% 
  group_by(number) %>% 
  summarise(number2 = n()) %>% 
  ggplot(aes(x = factor(number), y = number2)) +
  geom_bar(stat = "identity") +
  theme_bw()
  




write.csv(patient_info, "patient_info.csv", row.names = FALSE)

###patient information of batch 1
setwd("E:/project/smartD/patient_information")
sfu1_148_ga <- 
  readr::read_csv("SFU1-148_GA.csv")

head(sfu1_148_ga)



#####check the information of patients
setwd("E:/project/smartD/patient_information")
info <- readr::read_csv("SmartD_all346urine.csv")
load("E:/project/smartD/smartD_batch1/HILIC/POS/data_analysis/smartd_hilic_pos_batch1_5")

sample.info <- 
  smartd_hilic_pos_batch1_5@sample.info

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

match(sample.info$sample.name, info$sample_id_HILIC)

temp1 <- 
left_join(x = sample.info[,c("sample.name", "GA")], 
          y = info[,c("sample_id_HILIC", "Visit.GA")],
          by = c("sample.name" = "sample_id_HILIC"))

plot(temp1$GA, temp1$Visit.GA)






load("E:/project/smartD/smartD_batch2/HILIC/POS/data_analysis/smartd_hilic_pos_batch2_5")

sample.info <- 
  smartd_hilic_pos_batch2_5@sample.info

sample.info <- 
  sample.info %>% 
  filter(GA != 0)

match(sample.info$sample.name, info$sample_id_HILIC)


temp2 <- 
  left_join(x = sample.info[,c("sample.name", "GA")], 
            y = info[,c("sample_id_HILIC", "Visit.GA")],
            by = c("sample.name" = "sample_id_HILIC"))

plot(temp2$GA, temp2$Visit.GA)













