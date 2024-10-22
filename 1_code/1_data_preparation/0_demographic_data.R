library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(tidyverse)

###read data

all346_urine_info <-
  readr::read_csv("2_data/patient_information/SmartD_all346urine.csv")

clinic_varibale_summary <-
  readxl::read_xlsx(
    "2_data/patient_information/SmartD_ClinicalVariables_PartiallySummarized.xlsx",
    sheet = 1
  )

## patient_info
final_info <-
  readr::read_csv("2_data/patient_information/FINALSMARTDiaphragm2_DATA_2019-07-01_1615.csv")

dir.create("3_data_analysis/1_data_preparation/0_demographic_data",
           recursive = TRUE)

setwd("3_data_analysis/1_data_preparation/0_demographic_data")

### we must get the samples with 'P' coresponding to which sample
all346_urine_info$sample_id_RPLC %>%
  stringr::str_sort(numeric = TRUE) %>%
  grep("SFU_B", ., value = TRUE)

sample_id_RPLC <-
  all346_urine_info$sample_id_RPLC

sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] <-
  sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] %>%
  gsub(pattern = "SFU_B",
       replacement = "",
       x = .) %>%
  paste("X", ., sep = "")

all346_urine_info$sample_id_RPLC <-
  sample_id_RPLC

all346_info <-
  all346_urine_info

all346_info <-
  all346_info %>%
  dplyr::select(
    Patient_ID = PTID,
    Sample_ID = sample_id_RPLC,
    Visit,
    Date.Acquired,
    GA = Visit.GA
  ) %>%
  arrange(Patient_ID, Visit)

### modify GA
GA <- all346_info$GA

GA <-
  (GA - trunc(GA)) * 10 / 7 + trunc(GA)

all346_info$GA <-
  GA

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
  left_join(
    all346_info,
    final_info,
    by = c("Patient_ID" = "participant_id", "Visit" = "redcap_event_name")
  )

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
  data.frame(GA1 = sample_info_191021$GA,
             as.numeric(GA2),
             stringsAsFactors = FALSE)

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

library(plyr)

temp <-
  sample_information %>%
  plyr::dlply(.variables = .(Patient_ID))

sample_information <-
  sample_information

sample_information <-
  sample_information %>%
  dplyr::rename(sample_id = Sample_ID, subject_id = Patient_ID, ) %>%
  dplyr::select(subject_id, sample_id, Visit, dplyr::everything()) %>%
  dplyr::arrange(subject_id, Visit)

sample_information$GA
sample_information %>% 
  dplyr::filter(is.na(GA)) %>%
  dplyr::select(sample_id)

sample_information$post_partum <-
  sample_information$GA

sample_information$post_partum[is.na(sample_information$post_partum)] <- "YES"
sample_information$post_partum[sample_information$post_partum != "YES"] <- "NO"

save(sample_information, file = "sample_information.rda")

###demographic data
demographic_data <-
  clinic_varibale_summary %>%
  dplyr::filter(ID %in% sample_information$subject_id) %>%
  dplyr::rename(subject_id = ID)

demographic_data <-
  demographic_data %>%
  dplyr::rename(mother_age = Age) %>%
  dplyr::mutate(mother_age = as.numeric(mother_age))

demographic_data <-
  demographic_data %>%
  dplyr::rename(
    mother_dob = "Maternal DOB",
    mother_status = "Marital Status",
    mother_ethnicity = "Ethinic Group",
    mother_employment = "Employed",
    mother_housing = "Housing",
    mother_language = "Language",
    mother_insurance = "Insurance",
    mother_height = Ht,
    mother_weight = Wt,
    mother_bmi = BMI,
    mother_parity = "Parity",
    mother_mental_health = PMH,
    mother_past_surgical_history = "PSH",
    child_weight = "Birth wt",
    child_sex = "Sex",
    mother_estimated_due_date = EDD,
    mother_delivery_date = DD,
    mother_induction = Induction
  )

unique(sample_information$subject_id)

library(plyr)

temp <-
  sample_information %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x <-
      x %>%
      dplyr::select(
        subject_id,
        enrollment_date,
        EDD,
        DD,
        height_cm,
        weight_kg,
        term_birth,
        ptb,
        abortion,
        "BMI",
        "Parity",
        Induction,
        Sex,
        "Birth wt",
        term_birth
      ) %>%
      dplyr::rename(
        mother_enrollment_date = enrollment_date,
        mother_estimated_due_date = EDD,
        mother_delivery_date = DD,
        mother_height = height_cm,
        mother_weight = weight_kg,
        mother_bmi = BMI,
        mother_term_birth = term_birth,
        mother_aboration = abortion,
        mother_parity = Parity,
        mother_induction = Induction,
        child_sex = Sex,
        child_weight = "Birth wt",
        mother_preterm_birth = ptb,
        mother_term_birth = term_birth
      ) %>%
      dplyr::arrange(mother_enrollment_date)
    x[1, , drop  = FALSE]
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


dim(temp)
dim(demographic_data)

sort(temp$subject_id) == sort(demographic_data$subject_id)

temp <-
  temp %>%
  dplyr::arrange(subject_id)

demographic_data <-
  demographic_data %>%
  dplyr::arrange(subject_id)

colnames(temp)

intersect(colnames(temp), colnames(demographic_data))

###mother_height
demographic_data$mother_height
temp$mother_height[28] <- 160.02
demographic_data$mother_height[28]

demographic_data$mother_height <-
  temp$mother_height

###mother_height
demographic_data$mother_weight
temp$mother_weight
demographic_data$mother_weight[28]
temp$mother_weight[28] <- 59.87

demographic_data$mother_weight <-
  temp$mother_weight

##mother_parity
demographic_data$mother_parity <-
  temp$mother_parity

###child_sex
demographic_data$child_sex <-
  temp$child_sex

###mother_estimated_due_date
temp$mother_estimated_due_date ==
  demographic_data$mother_estimated_due_date

###mother_delivery_date
temp$mother_delivery_date ==
  demographic_data$mother_delivery_date

demographic_data <-
  demographic_data %>%
  dplyr::left_join(temp %>% dplyr::select(
    -c(
      mother_height,
      mother_weight,
      mother_bmi,
      mother_parity,
      child_sex,
      child_weight,
      mother_estimated_due_date,
      mother_delivery_date,
      mother_induction
    )
  ), by = "subject_id")


demographic_data$mother_parity <-
  demographic_data$mother_parity %>%
  stringr::str_extract("G[0-9]{1}") %>%
  stringr::str_replace("G", "") %>%
  as.numeric()

demographic_data <-
  demographic_data %>%
  dplyr::mutate(
    mother_ethnicity = case_when(
      mother_ethnicity == "1" ~ "White",
      mother_ethnicity == "2" ~ "Black",
      mother_ethnicity == "3" ~ "Latina",
      mother_ethnicity == "4" ~ "Pacific Islander",
      mother_ethnicity == "5" ~ "Asian",
      mother_ethnicity == "4 (Asian)" ~ "Asian",
      mother_ethnicity == "Afr Am" ~ "Black",
      TRUE ~ mother_ethnicity
    )
  )

demographic_data$mother_ethnicity[demographic_data$mother_ethnicity == "Latina"] <- "White"
demographic_data$mother_ethnicity[demographic_data$mother_ethnicity == "Caucasian"] <- "White"

demographic_data$mother_bmi <-
  demographic_data$mother_weight / (demographic_data$mother_height / 100) ^
  2

demographic_data <-
  demographic_data %>%
  dplyr::mutate(
    child_sex = case_when(
      child_sex == "F, F" ~ "Female_Female",
      child_sex == "M,M" ~ "Male_Male",
      child_sex == "M,M" ~ "Male_Male",
      child_sex == "M / F" ~ "Male_Female",
      child_sex == "M" ~ "Male",
      child_sex == "F" ~ "Female",
      TRUE ~ child_sex
    )
  )

demographic_data <-
  demographic_data %>%
  dplyr::mutate(
    mother_induction = case_when(
      mother_induction == "Y" ~ "YES",
      mother_induction == "Yes" ~ "YES",
      mother_induction == "N" ~ "NO",
      TRUE ~ mother_induction
    )
  )

demographic_data$child_weight <-
  demographic_data$child_weight %>%
  stringr::str_split(pattern = "; | /") %>%
  lapply(function(x) {
    paste(as.numeric(x), collapse = "{}")
  }) %>%
  unlist()

demographic_data$child_weight[demographic_data$child_weight == "NA"] <- NA

demographic_data$mother_delivery_weeks <-
  as.numeric((
    as.Date(demographic_data$mother_delivery_date) - (as.Date(
      demographic_data$mother_estimated_due_date
    ) - 280)
  ) / 7)

cbind(demographic_data$mother_delivery_weeks,
      demographic_data$mother_preterm_birth)

demographic_data %>%
  dplyr::select(mother_delivery_weeks,
                mother_preterm_birth,
                mother_term_birth)

demographic_data <-
  demographic_data %>%
  dplyr::mutate(
    mother_preterm_birth =
      case_when(
        mother_delivery_weeks < 37 ~ "YES",
        mother_delivery_weeks >= 37 ~ "NO"
      )
  )

save(demographic_data, file = "demographic_data.rda")

colnames(sample_information)
colnames(demographic_data)

intersect(colnames(sample_information),
          colnames(demographic_data))

new_columns <-
setdiff(colnames(demographic_data),
        colnames(sample_information))

sample_information <-
  sample_information %>%
  dplyr::left_join(demographic_data %>% dplyr::select(subject_id, new_columns),
                   by = "subject_id")

save(sample_information, file = "sample_information.rda")
