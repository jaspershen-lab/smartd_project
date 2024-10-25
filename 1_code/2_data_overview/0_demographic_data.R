library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

###read data
load("3_data_analysis/1_data_preparation/0_demographic_data/demographic_data.rda")

load("3_data_analysis/1_data_preparation/0_demographic_data/demographic_data.rda")

dir.create("3_data_analysis/2_data_overview/0_demographic_data",
           recursive = TRUE)

setwd("3_data_analysis/2_data_overview/0_demographic_data")

library(ggExtra)

## we have
colnames(demographic_data)

demographic_data <-
  demographic_data %>%
  dplyr::arrange(desc(mother_delivery_weeks))

child_weight <-
  demographic_data$child_weight %>%
  stringr::str_split("\\{\\}") %>%
  lapply(function(x) {
    mean(as.numeric(x))
  }) %>%
  unlist()

## delivery weeks
mother_delivery_weeks <-
  demographic_data$mother_delivery_weeks

#------------------------------------------------------------------------------
## mother_age
mother_age <-
  demographic_data$mother_age %>%
  as.numeric()

#------------------------------------------------------------------------------
## mother_ethnicity
mother_ethnicity <-
  demographic_data$mother_ethnicity

#------------------------------------------------------------------------------
## mother_bmi
mother_bmi <-
  demographic_data$mother_bmi

#------------------------------------------------------------------------------
## mother_parity
mother_parity <- demographic_data$mother_parity

## mother_induction
mother_induction <-
  demographic_data$mother_induction

mother_induction[is.na(mother_induction)] <- "Unknown"

## child_sex and twins
child_sex <- demographic_data$child_sex

child_sex[is.na(child_sex)] <- "Unknown"

# child_weight
child_weight <-
  demographic_data$child_weight %>%
  stringr::str_split("\\{\\}") %>%
  lapply(function(x) {
    mean(as.numeric(x)) / 1000
  }) %>%
  unlist()

###deliver weeks
mother_delivery_weeks <-
  demographic_data$mother_delivery_weeks


library(gghalves)

####continuous data distribution
plot_continuous_data <-
  data.frame(
    subject_id = demographic_data$subject_id,
    child_weight,
    mother_bmi,
    mother_delivery_weeks,
    mother_age,
    mother_parity
  ) %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "class",
                      values_to = "value") %>%
  dplyr::mutate(class = factor(
    class,
    levels = c(
      "mother_age",
      "mother_bmi",
      "mother_parity",
      "mother_delivery_weeks",
      "child_weight"
    )
  )) %>%
  ggplot(aes(class, value)) +
  geom_half_boxplot(side = "l") +
  geom_half_violin(side = "r") +
  geom_half_dotplot(side = "r") +
  facet_wrap(~ class, scales = "free", nrow = 1) +
  base_theme +
  labs(x = "")

plot_continuous_data

ggsave(
  plot_continuous_data,
  filename = "plot_continuous_data.pdf",
  width = 15,
  height = 5
)


###categorical data distribution

plot_discrete_data <-
  data.frame(subject_id = demographic_data$subject_id,
             mother_ethnicity,
             mother_induction,
             child_sex) %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "class",
                      values_to = "value") %>%
  dplyr::mutate(class = factor(
    class,
    levels = c("mother_ethnicity", "mother_induction", "child_sex")
  )) %>%
  ggplot(aes(class)) +
  geom_bar(stat = "count", aes(fill = value)) +
  scale_fill_manual(values = c(ethnicity_color, child_sex_color, induction_color)) +
  facet_wrap(~ class, scales = "free", nrow = 1) +
  base_theme +
  labs(x = "")

plot_discrete_data

ggsave(
  plot_discrete_data,
  filename = "plot_discrete_data.pdf",
  width = 6,
  height = 5
)

###circos plot
par(mar = c(2, 2, 2, 2))
par(xpd = TRUE)
library(circlize)
set.seed(999)

####colors for text
## which person is preterm
preterm1 <-
  demographic_data %>%
  dplyr::filter(mother_delivery_weeks < 37) %>%
  dplyr::pull(subject_id) %>%
  unique()

preterm2 <-
  demographic_data %>%
  dplyr::filter(mother_delivery_weeks >= 37 &
                  mother_delivery_weeks < 39) %>%
  dplyr::pull(subject_id) %>%
  unique()

text_colour <- unique(demographic_data$subject_id)

text_colour <- case_when(
  text_colour %in% preterm1 ~ "#58593FFF",
  text_colour %in% preterm2 ~ "#800000FF",
  TRUE ~ "#BEBEBE"
)

circos.clear()

## 36 participants
df <- data.frame(
  factors = demographic_data$subject_id,
  x = 1,
  y = 1,
  demographic_data
) %>%
  dplyr::mutate(factors = factor(factors, levels = factors))

circos.par(
  "track.height" = 0.1,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0.5, 35), 90),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.5, 1.5))

###Mother demographic data
### mother_age
range(mother_age)

circos.track(
  factors = df$factors,
  x = df$x,
  y = df$y,
  ylim = range(mother_age),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = round(c(min(mother_age), round((
        min(mother_age) + max(mother_age)
      ) / 2, 2), max(mother_age)), 3),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.text(
      x = mean(xlim),
      y = max(mother_age) + 13,
      labels = df$subject_id[i],
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.8,
      col = text_colour[i]
      # adj = c(0, 0.5)
    )
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = mother_age[i],
      col = "#8A9045FF",
      bg.border = "black"
    )
  }
)

##mother BMI
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = range(mother_bmi),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = round(c(min(mother_bmi), round((
        min(mother_bmi) + max(mother_bmi)
      ) / 2, 2), max(mother_bmi)), 3),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = mother_bmi[i],
      col = "#8A9045FF",
      bg.border = "black"
    )
  }
)

##mother_delivery_weeks
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = range(mother_delivery_weeks),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = round(c(min(mother_delivery_weeks), round((
        min(mother_delivery_weeks) + max(mother_delivery_weeks)
      ) / 2, 2), max(mother_delivery_weeks)), 3),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.points(
      x = mean(xlim),
      y = mother_delivery_weeks[i],
      cex = 1,
      pch = 19,
      col = "#8A9045FF"
    )
    circos.lines(
      x = mean(xlim),
      y =  mother_delivery_weeks[i],
      pch = 16,
      cex = 8,
      type = "h",
      col = "#8A9045FF",
      lwd = 2
    )
  }
)

## mother_ethnicity
library(RColorBrewer)
# display.brewer.all()
# display.brewer.pal(n = 7, name = "Set3")

temp_col <- mother_ethnicity
temp_col <-
  case_when(
    temp_col == "Asian" ~ ethnicity_color["Asian"],
    temp_col == "Black" ~ ethnicity_color["Black"],
    temp_col == "White" ~ ethnicity_color["White"],
    temp_col == "Pacific Islander" ~ ethnicity_color["Pacific Islander"],
    temp_col == "Other" ~ ethnicity_color["Other"]
  )

temp_col <- unname(temp_col)

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)

## mother_parity
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = range(mother_parity),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = round(c(min(mother_parity), round((
        min(mother_parity) + max(mother_parity)
      ) / 2, 2), max(mother_parity)), 3),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = mother_parity[i],
      col = "#8A9045FF",
      bg.border = "black"
    )
  }
)

####mother induction
temp_col <- mother_induction

temp_col <-
  case_when(
    temp_col == "YES" ~ induction_color["YES"],
    temp_col == "NO" ~ induction_color["NO"],
    temp_col == "Unknown" ~ induction_color["Unknown"]
  )

temp_col <-
  unname(temp_col)

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)


#### child_sex
temp_col <- child_sex
temp_col <-
  case_when(
    temp_col == "Male" ~ child_sex_color["Male"],
    temp_col == "Female" ~ child_sex_color["Female"],
    temp_col == "Unknown" ~ child_sex_color["Unknown"],
    temp_col == "Male_Female" ~ child_sex_color["Male_Female"],
    temp_col == "Female_Female" ~ child_sex_color["Female_Female"]
  )

temp_col <- unname(temp_col)

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)

## child_weight

range(child_weight, na.rm = TRUE)

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = range(child_weight, na.rm = TRUE),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = round(c(min(child_weight, na.rm = TRUE), round((
        min(child_weight, na.rm = TRUE) + max(child_weight, na.rm = TRUE)
      ) / 2, 2), max(child_weight, na.rm = TRUE)), 3),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = child_weight[i],
      col = "#8A9045FF",
      bg.border = "black"
    )
  }
)


circos.clear()


#
#
# #### 20191030
# setwd(r4projects::get_project_wd())
# setwd("patient information/")
# demographic_data <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
# demographic_data <- readr::read_csv("demographic_data.csv")
# demographic_data$ID <-
#   demographic_data$ID %>%
#   stringr::str_replace("sf", "") %>%
#   paste("SF", ., sep = "")
#
#
# demographic_data <-
#   demographic_data %>%
#   dplyr::filter(ID %in% demographic_data$subject_id)
#
# ## ID:participant ID
# ## Maternal DOB:participant ID
# load(
#   "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_dis"
# )
# load(
#   "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_val"
# )
#
#
# dis_id <-
#   sample_data_dis$subject_id %>%
#   unique()
#
# val_id <-
#   sample_data_val$subject_id %>%
#   unique()
#
# dis_info <-
#   demographic_data %>%
#   dplyr::filter(ID %in% dis_id)
#
# val_info <-
#   demographic_data %>%
#   dplyr::filter(ID %in% val_id)
#
# ## mother mother_age
# library(tidyverse)
#
# dis_info <-
#   dis_info %>%
#   dplyr:::mutate(mother_age = as.numeric(mother_age))
#
# dis_info %>%
#   dplyr::select(mother_age) %>%
#   mutate(mother_age = as.numeric(mother_age)) %>%
#   summarise(mean = mean(mother_age), sd = sd(mother_age))
#
# val_info <-
#   val_info %>%
#   dplyr:::mutate(mother_age = as.numeric(mother_age))
#
# val_info %>%
#   dplyr::select(mother_age) %>%
#   mutate(mother_age = as.numeric(mother_age)) %>%
#   summarise(mean = mean(mother_age), sd = sd(mother_age))
#
#
# ## mother_bmi
# dis_info %>%
#   dplyr::select(Ht, Wt, mother_bmi) %>%
#   mutate(ht = trans_ht(Ht),
#          wt = trans_wt(Wt),
#          mother_bmi = wt / (ht / 100) ^ 2) %>%
#   summarise(mean.mother_bmi = mean(mother_bmi),
#             sd.mother_bmi = sd(mother_bmi))
#
#
# dis_info <-
#   dis_info %>%
#   dplyr::mutate(ht = trans_ht(Ht),
#                 wt = trans_wt(Wt),
#                 mother_bmi = wt / (ht / 100) ^ 2)
#
# val_info %>%
#   dplyr::select(Ht, Wt, mother_bmi) %>%
#   mutate(ht = trans_ht(Ht),
#          wt = trans_wt(Wt),
#          mother_bmi = wt / (ht / 100) ^ 2) %>%
#   summarise(mean.mother_bmi = mean(mother_bmi),
#             sd.mother_bmi = sd(mother_bmi))
#
# val_info <-
#   val_info %>%
#   dplyr::mutate(ht = trans_ht(Ht),
#                 wt = trans_wt(Wt),
#                 mother_bmi = wt / (ht / 100) ^ 2)
#
#
# ## GA
# dis_info2 <-
#   sample_data_dis %>%
#   dplyr::filter(!is.na(day_0)) %>%
#   select(-contains("_POS")) %>%
#   select(-contains("_NEG")) %>%
#   left_join(dis_info, ., by = c("ID" = "subject_id")) %>%
#   mutate(
#     day_0 = as.Date(day_0, "%m/%d/%Y"),
#     edd = as.Date(EDD.x),
#     begin = edd - 280
#   ) %>%
#   transmute(ga = as.Date(DD.x) - begin) %>%
#   summarise(mean(as.numeric(ga)), sd(as.numeric(ga)))
#
#
# dis_info2
#
# val_info2 <-
#   sample_data_val %>%
#   dplyr::filter(!is.na(day_0)) %>%
#   select(-contains("_POS")) %>%
#   select(-contains("_NEG")) %>%
#   left_join(val_info, ., by = c("ID" = "subject_id")) %>%
#   mutate(
#     day_0 = as.Date(day_0, "%m/%d/%Y"),
#     edd = as.Date(EDD.x),
#     begin = edd - 280
#   ) %>%
#   transmute(ga = as.Date(DD.x) - begin) %>%
#   summarise(mean(as.numeric(ga)), sd(as.numeric(ga)))
#
# val_info2
#
# dis_info %>%
#   select(`Birth wt`) %>%
#   dplyr::filter(!is.na(`Birth wt`)) %>%
#   pull(`Birth wt`) %>%
#   stringr::str_split(";") %>%
#   unlist() %>%
#   as.numeric() %>%
#   as_tibble() %>%
#   summarise(mean = mean(value), sd = sd(value))
#
#
# val_info %>%
#   select(`Birth wt`) %>%
#   dplyr::filter(!is.na(`Birth wt`)) %>%
#   pull(`Birth wt`) %>%
#   stringr::str_replace("/", ";") %>%
#   stringr::str_split(";") %>%
#   unlist() %>%
#   as.numeric() %>%
#   as_tibble() %>%
#   summarise(mean = mean(value), sd = sd(value))
#
#
#
# ## gender of child
# dis_info$child_sex %>%
#   as_tibble() %>%
#   dplyr::filter(!is.na(value))
#
#
# dis_info
#
#
#
#
#
#
#
#
# ####20200610
# #### 20191030
# setwd(r4projects::get_project_wd())
# setwd("patient information/")
# demographic_data <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
# demographic_data <- readr::read_csv("demographic_data.csv")
# demographic_data$ID <-
#   demographic_data$ID %>%
#   stringr::str_replace("sf", "") %>%
#   paste("SF", ., sep = "")
#
#
# demographic_data <-
#   demographic_data %>%
#   dplyr::filter(ID %in% demographic_data$subject_id)
#
# ## ID:participant ID
# ## Maternal DOB:participant ID
# load(
#   "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_dis"
# )
# load(
#   "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_val"
# )
#
#
# dis_id <-
#   sample_data_dis$subject_id %>%
#   unique()
#
# val_id <-
#   sample_data_val$subject_id %>%
#   unique()
#
# dis_info <-
#   demographic_data %>%
#   dplyr::filter(ID %in% dis_id)
#
# val_info <-
#   demographic_data %>%
#   dplyr::filter(ID %in% val_id)
#
# dis_info$class <- "Training"
# val_info$class <- "Testing"
#
# patient_info <- rbind(dis_info, val_info)
#
# ## mother mother_age
# library(tidyverse)
#
# patient_info <-
#   patient_info %>%
#   dplyr:::mutate(mother_age = as.numeric(mother_age))
#
# ## mother_bmi
#
# patient_info <-
#   patient_info %>%
#   dplyr::mutate(ht = trans_ht(Ht),
#                 wt = trans_wt(Wt),
#                 mother_bmi = wt / (ht / 100) ^ 2)
#
# ## GA
#
# ##ethnicity
# patient_info <-
#   patient_info %>%
#   dplyr::mutate(
#     mother_ethnicity =  case_when(
#       mother_ethnicity == "1" ~ "White",
#       mother_ethnicity == "2" ~ "Black",
#       mother_ethnicity == "3" ~ "Latina",
#       mother_ethnicity == "4" ~ "Pacific Islander",
#       mother_ethnicity == "5" ~ "Asian",
#       mother_ethnicity == "4 (Asian)" ~ "Asian",
#       mother_ethnicity == "Afr Am" ~ "Black",
#       TRUE ~ mother_ethnicity
#     )
#   )
#
#
#
# patient_info <-
#   patient_info %>%
#   dplyr::mutate(
#     mother_parity = stringr::str_extract(mother_parity, "G[0-9]{1}") %>%
#       stringr::str_replace("G", "") %>%
#       as.numeric()
#   )
#
#
#
#
# xlsx::write.xlsx(patient_info, "patient_info.xlsx")
