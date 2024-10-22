#to avoind source
no_exist_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(tidyverse)

setwd("3_data_analysis/data_analysis20200108/urine_metabolome/sample_information/")

all346_urine_info <-
  readr::read_csv("SmartD_all346urine.csv")

clinic_varibale_summary <-
  readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx",
                    sheet = 1
  )

### we must get the samples with 'P' coresponding to which sample
all346_urine_info$sample_id_RPLC %>%
  stringr::str_sort(numeric = TRUE) %>%
  grep("SFU_B", ., value = TRUE)

sample_id_RPLC <-
  all346_urine_info$sample_id_RPLC

sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] <-
  sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] %>%
  gsub(pattern = "SFU_B", replacement = "", x = .) %>%
  paste("X", ., sep = "")

all346_urine_info$sample_id_RPLC <-
  sample_id_RPLC

all346_info <-
  all346_urine_info

all346_info <-
  all346_info %>%
  select(
    Patient_ID = PTID, Sample_ID = sample_id_RPLC,
    Visit, Date.Acquired, GA = Visit.GA
  ) %>%
  arrange(Patient_ID, Visit)

## some sample's GA are NA, should be removed
# all346_info %>%
#   dplyr::filter(is.na(GA)) %>%
#   pull(Sample_ID)
#
# patient_info <-
#   readr::read_csv("patient_info.csv")
#
# patient_info$Sample_ID <-
#   patient_info$Sample_ID %>%
#   gsub(pattern = "SFU_B", replacement = "X", x = .)
#
# patient_info %>%
#   select(Sample_ID, `Visit GA`) %>%
#   arrange(Sample_ID)
#
# all346_info %>%
#   select(Sample_ID, GA) %>%
#   dplyr::filter(stringr::str_starts(Sample_ID, "X")) %>%
#   arrange(Sample_ID)
#
#
# intersect(all346_info$Sample_ID, patient_info$Sample_ID)

### modify GA
GA <- all346_info$GA

GA <-
  (GA - trunc(GA)) * 10 / 7 + trunc(GA)

all346_info$GA <-
  GA

## patient_info
final_info <-
  readr::read_csv("FINALSMARTDiaphragm2_DATA_2019-07-01_1615.csv")

### change the Patient ID
all346_info$Patient_ID[-grep(pattern = "SF", all346_info$Patient_ID)] <-
  all346_info$Patient_ID[-grep(pattern = "SF", all346_info$Patient_ID)] %>%
  paste("SF", ., sep = "")

participant_id <-
  final_info$participant_id

participant_id <-
  gsub("sf", "SF", participant_id)

participant_id[-grep("SF", participant_id)] <-
  paste("SF", participant_id[-grep("SF", participant_id)], sep = "")

final_info$participant_id <-
  participant_id

setdiff(all346_info$Patient_ID, final_info$participant_id)
setdiff(final_info$participant_id, all346_info$Patient_ID)

final_info$redcap_event_name <-
  gsub("study_visit_|_arm_1", "", final_info$redcap_event_name) %>%
  as.numeric()

sample_info_191021 <-
  left_join(all346_info, final_info, by = c(
    "Patient_ID" = "participant_id",
    "Visit" = "redcap_event_name"
  ))

# sample_info_191021 %>%
#   dplyr::filter(Sample_ID %in% name_batch1) %>%
#   pull(GA)
#
# sample_info_191021 %>%
#   dplyr::filter(Sample_ID %in% name_batch2) %>%
#   pull(GA)

clinic_varibale_summary$ID <-
  stringr::str_to_upper(clinic_varibale_summary$ID)

clinic_varibale_summary$ID[-grep("SF", clinic_varibale_summary$ID)] <-
  paste("SF", clinic_varibale_summary$ID[-grep("SF", clinic_varibale_summary$ID)], sep = "")

sample_info_191021$Patient_ID %in% clinic_varibale_summary$ID

sample_info_191021 <-
  sample_info_191021 %>%
  left_join(clinic_varibale_summary, by = c("Patient_ID" = "ID"))

## recalculate the GA
begin_date <- as.Date(sample_info_191021$EDD) - 280
GA2 <- (as.Date(sample_info_191021$Date.Acquired, format = "%m/%d/%y") - begin_date) / 7

GA <-
  data.frame(GA1 = sample_info_191021$GA, as.numeric(GA2), stringsAsFactors = FALSE)

GA <-
  apply(GA, 1, function(x) {
    if (is.na(x[1])) {
      return(NA)
    } else {
      return(x[2])
    }
  })

sample_info_191021$GA <- GA

sample_information <- sample_info_191021

##X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on
idx <- match(paste("X", 178:198, sep = ""), sample_information$Sample_ID)

sample_information$Sample_ID[idx] <-
  sample_information$Sample_ID[idx] %>% 
  stringr:::str_replace("X", "") %>% 
  as.numeric() %>% 
  `-`(177) %>% 
  paste("P", ., sep = "")

write.csv(sample_info_191021, "sample_info_191021.csv", row.names = FALSE)

save(sample_information, file = "sample_information")

# patient_info$Sample_ID
#
# info1 <-
#   sample_info_191021 %>%
#   dplyr::filter(stringr::str_starts(Sample_ID, "X")) %>%
#   select(Patient_ID:GA) %>%
#   arrange(Sample_ID)
#
# info2 <-
#   patient_info %>%
#   dplyr::filter(stringr::str_starts(Sample_ID, "SFU_B")) %>%
#   select(Patient_ID:`Visit GA`) %>%
#   arrange(Sample_ID) %>%
#   mutate(Sample_ID = stringr::str_replace(Sample_ID, "SFU_B", "X"))



##### the clinic information
sample_info_191021$Patient_ID %>%
  unique() %>%
  length()

sample_info_191021 %>%
  dplyr::filter(is.na(GA)) %>%
  dplyr::select(Sample_ID)

library(ggExtra)

temp_data <-
  sample_info_191021 %>%
  mutate(GA2 = case_when(
    is.na(GA) ~ 45,
    !is.na(GA) ~ GA
  )) %>%
  mutate(class = case_when(
    is.na(GA) ~ "PP",
    !is.na(GA) ~ "Normal"
  )) %>%
  mutate(Begin.Date = as.Date(EDD) - 280) %>%
  mutate(term.date = as.Date(DD) - Begin.Date) %>%
  mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  arrange(diff_day)


## which person is preterm
preterm21 <-
  temp_data %>%
  # mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::filter(diff_day >= 21) %>%
  dplyr::pull(Patient_ID) %>%
  unique()


preterm7 <-
  temp_data %>%
  # mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::filter(diff_day >= 7 & diff_day < 21) %>%
  dplyr::pull(Patient_ID) %>%
  unique()



text_colour <- unique(temp_data$Patient_ID)
text_colour <- case_when(
  text_colour %in% preterm21 ~ "#58593FFF",
  text_colour %in% preterm7 ~ "#800000FF",
  TRUE ~ "#BEBEBE"
)

text_face <- unique(temp_data$Patient_ID)
text_face <- case_when(
  text_face %in% preterm21 ~ "bold",
  text_face %in% preterm7 ~ "bold",
  TRUE ~ "plain"
)


(
  plot1 <-
    # temp_data %>%
    ggplot(
      data = temp_data,
      aes(x = GA2, y = factor(Patient_ID, level = unique(Patient_ID)))
    ) +
    geom_point(
      aes(
        x = GA2,
        y = factor(Patient_ID, level = unique(Patient_ID)),
        colour = class
      ),
      data = temp_data,
      show.legend = FALSE,
      size = 2
    ) +
    geom_hline(
      yintercept = match(preterm7, unique(temp_data$Patient_ID)),
      colour = "#800000FF",
      linetype = 1
    ) +
    geom_hline(
      yintercept = match(preterm21, unique(temp_data$Patient_ID)),
      colour = "#58593FFF",
      linetype = 1
    ) +
    scale_colour_manual(values = c(
      "PP" = "#155F83FF",
      "Normal" = "#FFA319FF"
      # "Preterm>7_days" = "#800000B2",
      # "Preterm>21_days" = "#350E20B2"
    )) +
    geom_rect(
      aes(
        xmin = 10,
        ymin = 0,
        xmax = 15,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_rect(
      aes(
        xmin = 20,
        ymin = 0,
        xmax = 25,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_rect(
      aes(
        xmin = 30,
        ymin = 0,
        xmax = 35,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_rect(
      aes(
        xmin = 40,
        ymin = 0,
        xmax = 45,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    labs(x = "Gestational age (weeks)", y = "Subject ID") +
    scale_x_continuous(
      breaks = c(10, 20, 30, 40, 45),
      labels = c(10, 20, 30, 40, "PP"),
      limits = c(10, 46)
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(colour = "#8A9045FF", linetype = 1),
      # plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      axis.title = element_text(size = 15),
      axis.text.y = element_text(
        size = 11,
        colour = text_colour,
        face = text_face
      ),
      axis.text.x = element_text(size = 13)
    )
)


term_date <-
  temp_data %>%
  select(Patient_ID, term.date) %>%
  mutate(term.date = as.numeric(term.date / 7)) %>%
  distinct()


plot1 <-
  plot1 +
  ggplot2::annotate(
    geom = "point", shape = 17,
    colour = "black", size = 2,
    x = term_date$term.date,
    y = term_date$Patient_ID
  )

(
  plot2 <-
    ggplot(temp_data, aes(x = GA2)) +
    geom_histogram(binwidth = 0.5, colour = "white", fill = "#8A9045B2") +
    labs(x = "GA (weeks)", y = "Sample number") +
    theme_bw() +
    scale_x_continuous(
      limits = c(10, 46),
      name = NULL, labels = NULL, breaks = NULL
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # plot.margin = margin(0,0,0,0),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13)
    )
)

(
  plot3 <-
    ggplot(temp_data, aes(x = factor(Patient_ID, levels = unique(Patient_ID)))) +
    geom_bar(width = 0.8, fill = "#8A9045B2") +
    labs(x = "GA (weeks)", y = "Sample number") +
    theme_bw() +
    scale_x_discrete(name = NULL, label = NULL, breaks = NULL) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    coord_flip() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # plot.margin = margin(0,0,0,0),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13)
    )
)


library(patchwork)

plot <-
  {
    plot2 + plot_spacer() + plot_layout(ncol = 2, widths = c(3, 1))
  } -
  {
    plot1 + plot3 + plot_layout(ncol = 2, widths = c(3, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 3))


plot

ggsave(plot,
       filename = "Sample_colection_distribution.pdf",
       width = 10, height = 7
)

ggsave(plot,
       filename = "Sample_colection_distribution2.pdf",
       width = 10, height = 7
)

plot



#### 20191030
setwd(r4projects::get_project_wd())
setwd("patient information/")
info <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
sample_info_191021 <- readr::read_csv("sample_info_191021.csv")
info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste("SF", ., sep = "")


info <-
  info %>%
  dplyr::filter(ID %in% sample_info_191021$Patient_ID)

## ID:participant ID
## Maternal DOB:participant ID
load("/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_dis")
load("/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_val")


dis_id <-
  sample_data_dis$Patient_ID %>%
  unique()

val_id <-
  sample_data_val$Patient_ID %>%
  unique()

dis_info <-
  info %>%
  dplyr::filter(ID %in% dis_id)

val_info <-
  info %>%
  dplyr::filter(ID %in% val_id)

## mother age
library(tidyverse)

dis_info <-
  dis_info %>% 
  dplyr:::mutate(Age = as.numeric(Age))

dis_info %>%
  dplyr::select(Age) %>%
  mutate(Age = as.numeric(Age)) %>%
  summarise(mean = mean(Age), sd = sd(Age))

val_info <-
  val_info %>% 
  dplyr:::mutate(Age = as.numeric(Age))

val_info %>%
  dplyr::select(Age) %>%
  mutate(Age = as.numeric(Age)) %>%
  summarise(mean = mean(Age), sd = sd(Age))


## BMI
dis_info %>%
  dplyr::select(Ht, Wt, BMI) %>%
  mutate(
    ht = trans_ht(Ht),
    wt = trans_wt(Wt),
    bmi = wt / (ht / 100)^2
  ) %>%
  summarise(
    mean.bmi = mean(bmi),
    sd.bmi = sd(bmi)
  )


dis_info <- 
  dis_info %>% 
  dplyr::mutate(
    ht = trans_ht(Ht),
    wt = trans_wt(Wt),
    bmi = wt / (ht / 100)^2
  )

val_info %>%
  dplyr::select(Ht, Wt, BMI) %>%
  mutate(
    ht = trans_ht(Ht),
    wt = trans_wt(Wt),
    bmi = wt / (ht / 100)^2
  ) %>%
  summarise(
    mean.bmi = mean(bmi),
    sd.bmi = sd(bmi)
  )

val_info <- 
  val_info %>% 
  dplyr::mutate(
    ht = trans_ht(Ht),
    wt = trans_wt(Wt),
    bmi = wt / (ht / 100)^2
  )


## GA
dis_info2 <-
  sample_data_dis %>%
  dplyr::filter(!is.na(day_0)) %>%
  select(-contains("_POS")) %>%
  select(-contains("_NEG")) %>%
  left_join(dis_info, ., by = c("ID" = "Patient_ID")) %>%
  mutate(
    day_0 = as.Date(day_0, "%m/%d/%Y"),
    edd = as.Date(EDD.x),
    begin = edd - 280
  ) %>%
  transmute(ga = as.Date(DD.x) - begin) %>%
  summarise(mean(as.numeric(ga)), sd(as.numeric(ga)))


dis_info2

val_info2 <-
  sample_data_val %>%
  dplyr::filter(!is.na(day_0)) %>%
  select(-contains("_POS")) %>%
  select(-contains("_NEG")) %>%
  left_join(val_info, ., by = c("ID" = "Patient_ID")) %>%
  mutate(
    day_0 = as.Date(day_0, "%m/%d/%Y"),
    edd = as.Date(EDD.x),
    begin = edd - 280
  ) %>%
  transmute(ga = as.Date(DD.x) - begin) %>%
  summarise(mean(as.numeric(ga)), sd(as.numeric(ga)))

val_info2

dis_info %>%
  select(`Birth wt`) %>%
  dplyr::filter(!is.na(`Birth wt`)) %>%
  pull(`Birth wt`) %>%
  stringr::str_split(";") %>%
  unlist() %>%
  as.numeric() %>%
  as_tibble() %>%
  summarise(mean = mean(value), sd = sd(value))


val_info %>%
  select(`Birth wt`) %>%
  dplyr::filter(!is.na(`Birth wt`)) %>%
  pull(`Birth wt`) %>%
  stringr::str_replace("/", ";") %>%
  stringr::str_split(";") %>%
  unlist() %>%
  as.numeric() %>%
  as_tibble() %>%
  summarise(mean = mean(value), sd = sd(value))



## gender of child
dis_info$Sex %>%
  as_tibble() %>%
  dplyr::filter(!is.na(value))


dis_info








####20200610
#### 20191030
setwd(r4projects::get_project_wd())
setwd("patient information/")
info <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
sample_info_191021 <- readr::read_csv("sample_info_191021.csv")
info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste("SF", ., sep = "")


info <-
  info %>%
  dplyr::filter(ID %in% sample_info_191021$Patient_ID)

## ID:participant ID
## Maternal DOB:participant ID
load("/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_dis")
load("/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_val")


dis_id <-
  sample_data_dis$Patient_ID %>%
  unique()

val_id <-
  sample_data_val$Patient_ID %>%
  unique()

dis_info <-
  info %>%
  dplyr::filter(ID %in% dis_id)

val_info <-
  info %>%
  dplyr::filter(ID %in% val_id)

dis_info$class <- "Training"
val_info$class <- "Testing"

patient_info <- rbind(dis_info, val_info)

## mother age
library(tidyverse)

patient_info <-
  patient_info %>% 
  dplyr:::mutate(Age = as.numeric(Age))

## BMI

patient_info <- 
  patient_info %>% 
  dplyr::mutate(
    ht = trans_ht(Ht),
    wt = trans_wt(Wt),
    bmi = wt / (ht / 100)^2
  )

## GA

##ethnicity
patient_info <-
  patient_info %>% 
  dplyr::mutate(
    ethinic =  case_when(
      ethinic == "1" ~ "White",
      ethinic == "2" ~ "Black",
      ethinic == "3" ~ "Latina",
      ethinic == "4" ~ "Pacific Islander",
      ethinic == "5" ~ "Asian",
      ethinic == "4 (Asian)" ~ "Asian",
      ethinic == "Afr Am" ~ "Black",
      TRUE ~ ethinic
    )
  )



patient_info <-
  patient_info %>% 
  dplyr::mutate(
    parity = stringr::str_extract(Parity, "G[0-9]{1}") %>%
      stringr::str_replace("G", "") %>%
      as.numeric()
  )




xlsx::write.xlsx(patient_info, "patient_info.xlsx")


