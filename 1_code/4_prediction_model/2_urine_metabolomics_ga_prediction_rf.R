##avoid source
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
source("1_code/100_tools.R")

###load data
load(
  "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/metabolites/urine_metabolomics_data.rda"
)

##first set the work directory to project folder
dir.create(
  "3_data_analysis/4_prediction/2_urine_metabolomics_ga_prediction_rf",
  recursive = TRUE
)
setwd("3_data_analysis/4_prediction/2_urine_metabolomics_ga_prediction_rf")

###remove qc samples and postpartum samples
urine_metabolomics_data@sample_info$GA

urine_metabolomics_data <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(GA))

data_discovery <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(batch == 1)

data_validation <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(batch == 2)

##############################################################################
####random forest
#############################################################################
library(randomForest)
##use boruta method
library(Boruta)
library(tidyverse)

###discovery dataset
discovery_y <-
  data_discovery %>%
  extract_sample_info() %>%
  dplyr::select(GA) %>%
  as.matrix()

discovery_x <-
  data_discovery %>%
  extract_expression_data() %>%
  as.matrix()

discovery_x <-
  log(discovery_x, 10) %>%
  as.data.frame()

discovery_x <-
  apply(discovery_x, 1, function(x) {
    (x - mean(x)) / sd(x)
  })

###validation dataset
validation_y <-
  data_validation %>%
  extract_sample_info() %>%
  dplyr::select(GA) %>%
  as.matrix()

validation_x <-
  data_validation %>%
  extract_expression_data() %>%
  as.matrix()

validation_x <-
  log(validation_x, 10) %>%
  as.data.frame()

validation_x <-
  apply(validation_x, 1, function(x) {
    (x - mean(x)) / sd(x)
  })

set.seed(200)

marker_rf <-
  vector(mode = "list", length = 100)

set.seed(200)

# for(i in 1:100) {
#   cat(i, " ")
#   boruta_test <-
#     Boruta(x = sample_data_dis_x,
#            y = sample_data_dis_y,
#            doTrace = 3,
#            holdHistory = TRUE)
#   marker_rf[[i]] <-
#     boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>%
#     names() %>%
#     sort()
# }
#
# save(marker_rf, file = "marker_rf")
# load("marker_rf")
#
# marker_rf_new <-
#   unlist(marker_rf) %>%
#   data.frame(name = ., stringsAsFactors = FALSE) %>%
#   dplyr::group_by(name) %>%
#   dplyr::summarise(n = n()) %>%
#   arrange(desc(n))
#
# print(marker_rf_new, n = 100)
#
# marker_rf <-
#   marker_rf_new %>%
#   filter(n >= 50) %>%
#   pull(name)
#
# ###discard the peaks with bad peak shapes
# is_table <- metabolite_tags %>%
#   filter(name %in% marker_rf) %>%
#   select(name, mz, rt)
#
# is_table_pos <-
#   is_table %>%
#   filter(stringr::str_detect(name, "POS"))
#
# is_table_neg <-
#   is_table %>%
#   filter(stringr::str_detect(name, "NEG"))
#
#
# xlsx::write.xlsx2(as.data.frame(is_table_pos),
#                   file = "peak_shape/is_table_pos.xlsx",
#                   row.names = FALSE)
#
# xlsx::write.xlsx2(as.data.frame(is_table_neg),
#                   file = "peak_shape/is_table_neg.xlsx",
#                   row.names = FALSE)

peak_data_pos <-
  metflow2::extractPeaks(
    path = "./peak_shape/POS/",
    ppm = 15,
    threads = 4,
    is.table = "is_table_pos.xlsx"
  )


for (i in 1:nrow(is_table_pos)) {
  cat(i, " ")
  plot <- metflow2::showPeak(
    object = peak_data_pos,
    peak.index = i,
    alpha = 0.5,
    interactive = FALSE
  )
  
  ggsave(
    plot,
    filename = file.path(
      "./peak_shape/POS",
      paste(is_table_pos$name[i], ".png", sep = "")
    ),
    width = 8,
    height = 6
  )
}


is_table_pos <- readxl::read_xlsx("./peak_shape/POS/is_table_pos.xlsx")

marker_rf_pos <-
  is_table_pos %>%
  filter(peak_shape == "G") %>%
  pull(name)


peak_data_neg <-
  metflow2::extractPeaks(
    path = "./peak_shape/NEG/",
    ppm = 15,
    threads = 4,
    is.table = "is_table_neg.xlsx"
  )

for (i in 1:nrow(is_table_neg)) {
  cat(i, " ")
  plot <- metflow2::showPeak(
    object = peak_data_neg,
    peak.index = i,
    alpha = 0.5,
    interactive = FALSE
  )
  
  ggsave(
    plot,
    filename = file.path(
      "./peak_shape/NEG",
      paste(is_table_neg$name[i], ".png", sep = "")
    ),
    width = 8,
    height = 6
  )
}


is_table_neg <- readxl::read_xlsx("./peak_shape/NEG/is_table_neg.xlsx")
marker_rf_neg <-
  is_table_neg %>%
  filter(peak_shape == "G") %>%
  pull(name)

marker_rf <- c(marker_rf_pos, marker_rf_neg)

marker_rf <-
  metabolite_tags %>%
  filter(name %in% marker_rf)

readr::write_csv(marker_rf, "marker_rf.csv")
###check identifiction

load("./annotation_confirmation/massbankDatabase0.0.1")
load("./annotation_confirmation/metlinDatabase0.0.1")
load("./annotation_confirmation/monaDatabase0.0.1")
load("./annotation_confirmation/msDatabase_rplc0.0.1")
load("./annotation_confirmation/nistDatabase0.0.1")
load("./annotation_confirmation/orbitrapDatabase0.0.1")
load("./annotation_confirmation/result.pRPLC.nce50")
load("./annotation_confirmation/result.pRPLC.nce25")
load("./annotation_confirmation/result.nRPLC.nce50")
load("./annotation_confirmation/result.nRPLC.nce25")


marker_rf <-
  filter(metabolite_tags, name %in% marker_rf)

###POS
database_name <-
  marker_rf %>%
  filter(stringr::str_detect(name, "POS")) %>%
  pull(Database) %>%
  unique()

for (i in 1:length(database_name)) {
  cat(i, " ")
  temp_database <- database_name[i]
  idx <- lapply(result.pRPLC.nce25, function(x)
    x@database) %>%
    unlist() %>%
    unname() %>%
    `==`(temp_database) %>%
    which()
  temp_name <- marker_rf %>%
    filter(Database == temp_database) %>%
    filter(stringr::str_detect(name, "POS")) %>%
    pull(name)
  
  if (length(temp_name) == 0) {
    next()
  } else{
    metID::ms2plot(
      object = result.pRPLC.nce25[[idx]],
      which.peak = temp_name,
      database = get(temp_database),
      path = file.path(
        "annotation_confirmation",
        paste("POS", temp_database, "NCE25", sep = "_")
      ),
      show.plot = FALSE
    )
    
    metID::ms2plot(
      object = result.pRPLC.nce50[[idx]],
      which.peak = temp_name,
      database = get(temp_database),
      path =  file.path(
        "annotation_confirmation",
        paste("POS", temp_database, "NCE50", sep = "_")
      ),
      ,
      show.plot = FALSE
    )
  }
  
}

###neg
database_name <-
  marker_rf %>%
  filter(stringr::str_detect(name, "NEG")) %>%
  pull(Database) %>%
  unique()

for (i in 1:length(database_name)) {
  cat(i, " ")
  temp_database <- database_name[i]
  idx <- lapply(result.nRPLC.nce25, function(x)
    x@database) %>%
    unlist() %>%
    unname() %>%
    `==`(temp_database) %>%
    which()
  temp_name <- marker_rf %>%
    filter(Database == temp_database) %>%
    filter(stringr::str_detect(name, "NEG")) %>%
    pull(name)
  
  if (length(temp_name) == 0) {
    next()
  } else{
    metID::ms2plot(
      object = result.nRPLC.nce25[[idx]],
      which.peak = temp_name,
      database = get(temp_database),
      path = file.path(
        "annotation_confirmation",
        paste("NEG", temp_database, "NCE25", sep = "_")
      ),
      show.plot = FALSE
    )
    
    metID::ms2plot(
      object = result.nRPLC.nce50[[idx]],
      which.peak = temp_name,
      database = get(temp_database),
      path =  file.path(
        "annotation_confirmation",
        paste("NEG", temp_database, "NCE50", sep = "_")
      ),
      ,
      show.plot = FALSE
    )
  }
  
}

marker_rf <- readr::read_csv("marker_rf.csv")

marker_rf <-
  marker_rf %>%
  filter(MS2_match == "G") %>%
  pull(name)

temp_data <-
  boruta_test$ImpHistory[1:length(marker_rf), ] %>%
  t() %>%
  as_tibble() %>%
  select_if(., no_infi) %>%
  dplyr::transmute(
    .,
    mean = apply(., 1, mean),
    sd = apply(., 1, sd),
    ymax = mean + sd,
    ymin = mean - sd
  ) %>%
  mutate(name = colnames(boruta_test$ImpHistory)) %>%
  filter(name %in% marker_rf) %>%
  left_join(., metabolite_tags, by = "name") %>%
  arrange(., mean)


no_infi <- function(x)
  all(!is.infinite(x))

marker_rf <-
  boruta_test$ImpHistory[1:length(marker_rf), ] %>%
  t() %>%
  as_tibble() %>%
  select_if(., no_infi) %>%
  dplyr::transmute(
    .,
    mean = apply(., 1, mean),
    sd = apply(., 1, sd),
    ymax = mean + sd,
    ymin = mean - sd
  ) %>%
  mutate(name = colnames(boruta_test$ImpHistory)) %>%
  filter(name %in% marker_rf) %>%
  dplyr::select(name, everything())


# marker_rf <- readr::read_csv("marker_rf_final.csv")
#
#
# temp_data <- marker_rf
# colnames(temp_data) <-
#   stringr::str_replace(colnames(temp_data), "_importance", "")
#
# temp_data <-
#   temp_data %>%
#   arrange(mean)


text_colour <- temp_data$super_class
text_colour[is.na(text_colour)] <- "Unknown"

text_colour

colour <<-
  c(
    # "Clinical information" = "#FF6F00FF",
    # "Alkaloids and derivatives" = "#ADE2D0FF",
    # "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    # "Organic acids and derivatives" = "#8A4198FF",
    # "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    # "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

text_colour <-
  colour[match(text_colour, names(colour))] %>%
  unname()

ggplot(temp_data, aes(x = factor(Compound.name, Compound.name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                colour = "#155F83FF",
                width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  theme_bw() +
  coord_flip() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    # panel.grid.major.y = element_line(colour = text_colour),
    axis.text.y = element_text(size = 10, colour = text_colour)
  )

ggsave(filename = "marker_rf.pdf",
       width = 7,
       height = 7)


colour <-
  c(
    # "Clinical information" = "#FF6F00FF",
    # "Alkaloids and derivatives" = "#ADE2D0FF",
    # "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    # "Organic acids and derivatives" = "#8A4198FF",
    # "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    # "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "white"
  )

text_colour <-
  colour[match(text_colour, names(colour))] %>%
  unname()

ggplot(temp_data, aes(x = factor(Compound.name, Compound.name), y = mean)) +
  labs(x = "", y = "Importance") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                colour = "#155F83FF",
                width = 0) +
  geom_point(size = 2, colour = "#FFA319FF") +
  ggdark::dark_theme_classic() +
  coord_flip() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    # panel.grid.major.y = element_line(colour = text_colour),
    # panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, colour = text_colour),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(
  filename = "marker_rf_dark.png",
  width = 7,
  height = 7,
  bg = "transparent"
)

colnames(marker_rf)[-1] <-
  colnames(marker_rf)[-1] %>%
  paste(., "importance", sep = "_")

marker_rf <-
  marker_rf %>%
  left_join(metabolite_tags, by = "name")

write.csv(marker_rf, "marker_rf_final.csv", row.names = FALSE)

marker_rf <- readr::read_csv("marker_rf_final.csv")

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

set.seed(210)
fgl.res <- randomForest::tuneRF(
  sample_data_dis_x_rf,
  sample_data_dis_y[, 1],
  mtryStart = 1,
  stepFactor = 2,
  trace = TRUE,
  plot = TRUE
)

fgl.res %>%
  as.data.frame() %>%
  ggplot(aes(mtry, y = OOBError)) +
  geom_line(colour = "#155F83FF") +
  geom_point(size = 4, colour = "#FFA319FF") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )

ggsave(filename = "mtry_vs_error.pdf",
       width = 7,
       height = 7)

rf_regression <-
  randomForest(
    x = sample_data_dis_x_rf,
    y = sample_data_dis_y[, 1],
    replace = TRUE,
    importance = TRUE,
    proximity = TRUE,
    mtry = 4
  )

##validate in validation dataset
###construct dataset for lasso regression
sample_data_val_x_rf <-
  sample_data_val_x %>%
  as_tibble() %>%
  dplyr::select(one_of(marker_rf$name))

###use validation dataset for validation
predicted_y <-
  predict(object = rf_regression, newdata = sample_data_val_x_rf)

plot(sample_data_val_y[, 1], predicted_y)

abline(0, 1)

prediction_self <-
  predict(object = rf_regression, newx = as.matrix(sample_data_dis_x_rf))

plot(sample_data_dis_y[, 1], prediction_self)
abline(0, 1)

discovery_data_measured_vs_predicted <-
  data.frame(
    measured = sample_data_dis_y[, 1],
    predicted = prediction_self,
    stringsAsFactors = FALSE
  ) %>%
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 15)
  )
discovery_data_measured_vs_predicted
# save(discovery_data_measured_vs_predicted, file = "discovery_data_measured_vs_predicted")
# ggsave(plot = discovery_data_measured_vs_predicted,
#         filename = "discovery_data_measured_vs_predicted.pdf",
#        width = 7, height = 7)


####dark
discovery_data_measured_vs_predicted +
  ggdark::dark_theme_bw() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    color = "white"
  ) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "grey"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 15)
  )

ggsave(
  filename = "discovery_data_measured_vs_predicted_dark.png",
  bg = "transparent",
  width = 7,
  height = 7
)

###from here we can see that why should us linear regression to correct prediction
measured_vs_predicted_error <-
  data.frame(
    measured = sample_data_dis_y[, 1],
    predicted = prediction_self,
    stringsAsFactors = FALSE
  ) %>%
  mutate(diff = predicted - measured) %>%
  mutate(colour = ifelse(diff > 0, "pos", "neg")) %>%
  ggplot(aes(measured, diff)) +
  geom_segment(aes(
    x = measured,
    xend = measured,
    y = 0,
    yend = diff,
    colour = colour
  ),
  show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2, aes(colour = colour), show.legend = FALSE) +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  scale_colour_manual(values = c("pos" = "#800000FF", "neg" = "#155F83FF")) +
  theme_bw() +
  labs(x = "GA_week (measured)", y = "GA_week (predicted - measured)") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

measured_vs_predicted_error
# save(measured_vs_predicted_error, file = "measured_vs_predicted_error")
# ggsave(filename = "measured_vs_predicted_error.pdf",
#        plot = measured_vs_predicted_error,
#        width = 7, height = 7)



linear_regression <-
  lm(formula = sample_data_dis_y[, 1] ~ prediction_self)

linear_regression1 <-
  lm(formula = prediction_self ~ sample_data_dis_y[, 1])

discovery_data_measured_vs_predicted2 <-
  data.frame(
    measured = sample_data_dis_y[, 1],
    predicted = prediction_self,
    stringsAsFactors = FALSE
  ) %>%
  ggplot(aes(measured, predicted)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  geom_abline(
    intercept = coef(linear_regression1)[1],
    slope = coef(linear_regression1)[2],
    linetype = 2,
    colour = "red"
  ) +
  labs(x = "GA (weeks, measured)", y = "GA (weeks, predicted)") +
  scale_x_continuous(limits = c(10, 42)) +
  scale_y_continuous(limits = c(10, 42)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 15)
  )
discovery_data_measured_vs_predicted2
# save(discovery_data_measured_vs_predicted2, file = "discovery_data_measured_vs_predicted2")
# ggsave("discovery_data_measured_vs_predicted2.pdf", width = 7, height = 7)

predicted_y2 <-
  coef(linear_regression)[2] * predicted_y + coef(linear_regression)[1]

measured_vs_predicted_val <-
  data.frame(
    "measured" = sample_data_val_y[, 1],
    'predicted' = predicted_y2,
    stringsAsFactors = FALSE
  ) %>%
  ggplot(aes(x = measured, predicted)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    color = "black"
  ) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 15)
  )

measured_vs_predicted_val
# save(measured_vs_predicted_val, file = "measured_vs_predicted_val")
# ggsave(filename = "measured_vs_predicted_val.pdf",
#        plot = measured_vs_predicted_val,
#        width = 7, height = 7)


temp_data <-
  data.frame(
    ID = sample_data_val$Patient_ID,
    "measured" = sample_data_val_y[, 1],
    'predicted' = predicted_y2,
    stringsAsFactors = FALSE
  )

predicted_result <- temp_data
save(predicted_result, file = "predicted_result")

load("predicted_result")

rmse_r2 <-
  temp_data %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(rmse = mean(abs(measured - predicted)), r2 = summary(lm(formula = predicted ~
                                                                             measured))$adj.r.squared)

RMSE_R2_for_each_person <-
  rmse_r2 %>%
  # arrange(r2) %>%
  # mutate(ID = factor(ID, levels = unique(ID))) %>%
  ggplot() +
  geom_segment(
    aes(
      x = ID,
      y = 0,
      xend = ID,
      yend = rmse
    ),
    linetype = 2,
    colour = "#FF6F00FF"
  ) +
  geom_point(aes(x = ID, y = rmse), colour = "#FF6F00FF") +
  geom_segment(aes(
    x = ID,
    y = 0,
    xend = ID,
    yend = -r2
  ),
  colour = "#5A9599FF",
  linetype = 2) +
  geom_point(aes(x = ID, y = -r2), colour = "#5A9599FF") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  coord_flip() +
  labs(x = "") +
  scale_y_continuous(limits = c(-1.3, 5)) +
  geom_text(aes(
    x = ID,
    y = rmse,
    label = round(rmse, 2)
  ),
  hjust = -0.5,
  ,
  colour = "grey") +
  geom_text(aes(
    x = ID,
    y = -r2,
    label = round(r2, 2)
  ),
  hjust = 1.3,
  colour = "grey") +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))
RMSE_R2_for_each_person
# save(RMSE_R2_for_each_person, file = "RMSE_R2_for_each_person")
# ggsave(filename = "RMSE_R2_for_each_person.pdf",
#        plot = RMSE_R2_for_each_person,
#        width = 7, height = 7)



r2_rmse <-
  rmse_r2 %>%
  tidyr::pivot_longer(cols = -ID,
                      names_to = "class",
                      values_to = "value") %>%
  ggplot(aes(x = calss, y = value)) +
  geom_violin(aes(x = class, y = value, colour = class)) +
  geom_boxplot(aes(x = class, y = value, colour = class), width = 0.2) +
  geom_jitter(aes(x = class, y = value, colour = class)) +
  scale_colour_manual(values = c("r2" = "#5A9599FF", "rmse" = "#FF6F00FF")) +
  facet_wrap(. ~ class, scales = "free") +
  theme_bw() +
  labs(y = "", x =  "") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_blank()
  )

# save(r2_rmse, file = "r2_rmse")
# ggsave(filename = "r2_rmse.pdf", plot = r2_rmse, width = 7, height = 7)

measured_vs_predicted_val_for_each_person <-
  data.frame(
    ID = sample_data_val$Patient_ID,
    "measured" = sample_data_val_y[, 1],
    'predicted' = predicted_y2,
    stringsAsFactors = FALSE
  ) %>%
  ggplot(aes(x = measured, predicted)) +
  ggplot2::facet_wrap(. ~ ID) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2,
             colour = "#FFA319FF",
             method = "loess") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 15),
    strip.background = element_blank(),
    strip.text = element_text(size = 13)
  )


measured_vs_predicted_val_for_each_person

rmse_r2

measured_vs_predicted_val_for_each_person <-
  measured_vs_predicted_val_for_each_person +
  geom_text(
    mapping = aes(
      x = -Inf,
      y = Inf,
      label = paste(
        paste("RMSE: ", round(rmse, 2), "\n", sep = ""),
        paste("R2: ", round(r2, 2), sep = "")
      )
    ),
    data = rmse_r2,
    hjust = -0.1,
    vjust = 1,
    size = 4
  )
measured_vs_predicted_val_for_each_person
save(measured_vs_predicted_val_for_each_person, file = "measured_vs_predicted_val_for_each_person")
ggsave(
  "measured_vs_predicted_val_for_each_person.pdf",
  width = 9,
  height = 7
)

##dark
measured_vs_predicted_val_for_each_person <-
  data.frame(
    ID = sample_data_val$Patient_ID,
    "measured" = sample_data_val_y[, 1],
    'predicted' = predicted_y2,
    stringsAsFactors = FALSE
  ) %>%
  ggplot(aes(x = measured, predicted)) +
  ggplot2::facet_wrap(. ~ ID) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_point(size = 2,
             colour = "#FFA319FF",
             method = "loess") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 15),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


measured_vs_predicted_val_for_each_person

rmse_r2

measured_vs_predicted_val_for_each_person <-
  measured_vs_predicted_val_for_each_person +
  geom_text(
    mapping = aes(
      x = -Inf,
      y = Inf,
      label = paste(
        paste("RMSE: ", round(rmse, 2), "\n", sep = ""),
        paste("R2: ", round(r2, 2), sep = "")
      )
    ),
    data = rmse_r2,
    hjust = -0.1,
    vjust = 1,
    size = 4
  )



ggsave(
  plot = measured_vs_predicted_val_for_each_person,
  filename = "measured_vs_predicted_val_for_each_person_dark.png",
  width = 9,
  height = 7,
  bg = "transparent"
)



library(waffle)

test <- marker_rf$class
test[is.na(test)] <- "Other"

test <- table(test)
test1 <- as.numeric(test)
names(test1) <- names(test)
waffle(sort(test) * 3,
       colors = ggsci::pal_futurama()(12),
       rows = 8)



abs(sample_data_val_y[, 1] - predicted_y2) %>%
  mean()

# abs(predicted_result$measured - predicted_result$predicted) %>%
#   mean()
#
# summary(lm(formula = predicted_result$predicted~predicted_result$measured))

# cor.test(predicted_result$measured, predicted_result$predicted)

## T1 0-13, T2 14 -26, T3 27-40
predicted_result %>%
  dplyr::filter(measured >= 0 & measured < 14) %>%
  mutate(x = measured - predicted) %>%
  pull(x) %>%
  abs() %>%
  mean()

predicted_result %>%
  dplyr::filter(measured >= 14 & measured < 26) %>%
  mutate(x = measured - predicted) %>%
  pull(x) %>%
  abs() %>%
  mean()

predicted_result %>%
  dplyr::filter(measured >= 26 & measured < 40) %>%
  mutate(x = measured - predicted) %>%
  pull(x) %>%
  abs() %>%
  mean()

abs(predicted_result$measured - predicted_result$predicted) %>%
  mean()


summary(lm(formula = predicted_y2 ~ sample_data_val_y[, 1]))

write.table("Information", "information.txt")
cat("RMSE for external validation dataset:\n",
    file = "information.txt",
    append = TRUE)
cat(
  abs(sample_data_val_y[, 1] - predicted_y2) %>%
    mean(),
  file = "information.txt",
  append = TRUE
)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for external validation dataset:\n",
    file = "information.txt",
    append = TRUE)
cat(summary(lm(formula = predicted_y2 ~ sample_data_val_y[, 1]))$adj.r.squared,
    file = "information.txt",
    append = TRUE)
cat("\n", file = "information.txt", append = TRUE)

##bootstrap
predict_y <- vector(mode = "list", length = 1000)
y <- vector(mode = "list", length = 1000)

for (i in 1:1000) {
  cat(i, " ")
  patient_id <-
    sample_data_dis$Patient_ID %>%
    unique()
  
  dis_id <-
    sample(patient_id,
           size = length(patient_id),
           replace = TRUE) %>%
    unique()
  
  val_id <-
    setdiff(patient_id, dis_id)
  
  dis_index <- which(sample_data_dis$Patient_ID %in% dis_id)
  val_index <- which(sample_data_dis$Patient_ID %in% val_id)
  
  # dis_index <-
  #   sample(1:nrow(sample_data_dis_x_rf),
  #          size = nrow(sample_data_dis_x_rf), replace = TRUE) %>%
  #   unique() %>%
  #   sort()
  #
  # val_index <-
  #   setdiff(1:nrow(sample_data_dis_x_rf), dis_index)
  
  rf_regression_temp <-
    randomForest(
      x = sample_data_dis_x_rf[dis_index, ],
      y = sample_data_dis_y[dis_index, 1],
      replace = TRUE,
      importance = TRUE,
      proximity = TRUE,
      mtry = 4
    )
  
  ##construct a new linear model to correct prediction and real value
  prediction_self <-
    predict(object = rf_regression_temp, newx = as.matrix(sample_data_dis_x_rf[dis_index, ]))
  
  ###construct a new liner regression model to correct them
  linear_regression <-
    lm(formula = sample_data_dis_y[dis_index, 1] ~ prediction_self)
  
  temp_predict_y <-
    predict(object = rf_regression_temp,
            newdata = as.matrix(sample_data_dis_x_rf[val_index, ]))
  
  temp_predict_y2 <-
    coef(linear_regression)[2] * temp_predict_y + coef(linear_regression)[1]
  
  predict_y[[i]] <-
    temp_predict_y2
  
  y[[i]] <-
    sample_data_dis_y[val_index, 1]
  names(predict_y[[i]]) <-
    names(y[[i]]) <-
    sample_data_dis$Sample_ID[val_index]
  
}

save(predict_y, file = "predict_y")
save(y, file = "y")


temp_sample_id <-
  sample_data_dis %>%
  dplyr::filter(Patient_ID == "SF1562") %>%
  pull(Sample_ID)



temp_y1 <-
  data.frame(
    sample_id = names(unlist(y)[names(unlist(y)) %in% temp_sample_id]),
    y1 = unlist(y)[names(unlist(y)) %in% temp_sample_id],
    stringsAsFactors = FALSE
  )


temp_y2 <- data.frame(
  sample_id = names(unlist(predict_y)[names(unlist(y)) %in% temp_sample_id]),
  y2 = unlist(predict_y)[names(unlist(predict_y)) %in% temp_sample_id],
  stringsAsFactors = FALSE
)

temp_y <-
  temp_y1 %>%
  dplyr::left_join(temp_y2, by = "sample_id")


temp_y <-
  temp_y %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    y1_mean = mean(y1),
    y2_mean = mean(y2),
    y2_sd = sd(y2),
    n = dplyr::n()
  )


temp_y %>%
  arrange(y1_mean) %>%
  ggplot(aes(x = y1_mean, y = y2_mean)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  labs(x = "Gestational Age (week, Actual)", y = "Gestational Age (week, Predicted)") +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_errorbar(
    aes(ymin = y2_mean - y2_sd, ymax = y2_mean + y2_sd),
    data = temp_y,
    colour = "#155F83FF"
  ) +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  scale_x_continuous(limits = c(17.5, 39)) +
  scale_y_continuous(limits = c(17.5, 39)) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  theme_bw()


abs(temp_y$y1_mean - temp_y$y2_mean) %>%
  mean()

summary(lm(formula = temp_y$y1_mean ~ temp_y$y2_mean))







prediction <-
  data.frame(
    y = unlist(y),
    predict = unlist(predict_y),
    stringsAsFactors = FALSE
  )

temp <-
  prediction %>%
  arrange(., y) %>%
  group_by(., y) %>%
  dplyr::summarise(mean = mean(predict), sd = sd(predict))


abs(temp$y - temp$mean) %>%
  mean()

summary(lm(formula = temp$y ~ temp$mean))

cor.test(temp$y, temp$mean)



temp %>%
  dplyr::filter(y >= 0 & y < 14) %>%
  mutate(x = y - mean) %>%
  pull(x) %>%
  abs() %>%
  mean()

temp %>%
  dplyr::filter(y >= 14 & y < 26) %>%
  mutate(x = y - mean) %>%
  pull(x) %>%
  abs() %>%
  mean()


temp %>%
  dplyr::filter(y >= 26 & y < 40) %>%
  mutate(x = y - mean) %>%
  pull(x) %>%
  abs() %>%
  mean()




cat("RMSE for internal validation dataset:\n",
    file = "information.txt",
    append = TRUE)

cat(abs(temp$y - temp$mean) %>%
      mean(),
    file = "information.txt",
    append = TRUE)
cat("\n", file = "information.txt", append = TRUE)
cat("R2 for internal validation dataset:\n",
    file = "information.txt",
    append = TRUE)

cat(summary(lm(formula = temp$y ~ temp$mean))$adj.r.squared,
    file = "information.txt",
    append = TRUE)
cat("\n", file = "information.txt", append = TRUE)


measured_vs_predicted <-
  temp %>%
  mutate(ymax = mean + sd, ymin = mean - sd) %>%
  ggplot(aes(x = y, y = mean)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

save(measured_vs_predicted, file = "measured_vs_predicted")
ggsave(
  filename = "measured_vs_predicted_dis.pdf",
  plot = measured_vs_predicted,
  width = 7,
  height = 7
)



temp %>%
  mutate(ymax = mean + sd, ymin = mean - sd) %>%
  ggplot(aes(x = y, y = mean)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    color = "white"
  ) +
  labs(x = "GA_week (measured)", y = "GA_week (predicted)") +
  scale_x_continuous(limits = c(12, 40)) +
  scale_y_continuous(limits = c(12, 40)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), colour = "#155F83FF") +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "grey")
  )

ggsave(
  filename = "measured_vs_predicted_dark.png",
  bg = "transparent",
  width = 7,
  height = 7
)


###check the clinical information to affect the preidct
sample_info <-
  readr::read_csv("/Users/shenxt/projects/smartD/patient_information/sample_info_191021.csv")

##add clinical information
##due date
info <-
  readxl::read_xlsx(
    "/Users/shenxt/projects/smartD/patient_information/SmartD_ClinicalVariables_PartiallySummarized.xlsx"
  )
info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

info <-
  info %>%
  dplyr::filter(ID %in% sample_info$Patient_ID)

info$ID == unique(sample_info$Patient_ID)

due_date <-
  (as.Date(info$DD) - (as.Date(info$EDD) - 280)) / 7 %>%
  as.numeric()

##age
age <-
  info$Age %>%
  as.numeric()
age
#
##ethinic
ethinic <-
  info$`Ethinic Group`
#
ethinic <-
  case_when(
    ethinic == "1" ~ "White",
    ethinic == "Caucasian" ~ "White",
    ethinic == "2" ~ "Black",
    ethinic == "3" ~ "Latina",
    ethinic == "4" ~ "Pacific Islander",
    ethinic == "5" ~ "Asian",
    ethinic == "4 (Asian)" ~ "Asian",
    ethinic == "Asian" ~ "Asian",
    ethinic == "Afr Am" ~ "Black",
    TRUE ~ ethinic
  )
#
# ##BMI
ht <- trans_ht(info$Ht) / 100
wt <- trans_wt(info$Wt)
bmi <- wt / (ht ^ 2)
#
# #------------------------------------------------------------------------------
# ##parity
parity <- info$Parity
parity <-
  parity %>%
  stringr::str_extract("G[0-9]{1}") %>%
  stringr::str_replace("G", "") %>%
  as.numeric()
#
# ##sex and twins
# info$Sex
sex <- info$Sex
sex <-
  case_when(
    is.na(sex) ~ "NA",
    sex == "F, F" ~ 'F_F',
    sex == "M,M" ~ 'M_M',
    sex == "M,M" ~ 'M_M',
    sex == "M / F" ~ 'M_F',
    TRUE ~ sex
  )
# ##IVF
ivf <- rep(NA, nrow(info))
ivf[grep("ivf", stringr::str_to_lower(info$`Pregnancy dx`))] <- "YES"
ivf[grep("transfer", stringr::str_to_lower(info$Dating))] <- "YES"
ivf[is.na(ivf)] <- "NO"
#
##induction
induction <-
  info$Induction
#
induction <-
  case_when(
    is.na(induction) ~ "NA",
    induction == "Y" ~ "YES",
    induction == "Yes" ~ "YES",
    induction == "N" ~ "NO",
  )
#
#
birth_wt <-
  info$`Birth wt`
birth_wt[19] <- "3015;3170"
birth_wt <-
  sapply(birth_wt, function(x) {
    if (!is.na(x)) {
      x <- stringr::str_split(x, ";")[[1]] %>%
        as.numeric() %>%
        sum()
    } else{
      x
    }
    x
  })
#
birth_wt <- unname(birth_wt) %>%
  as.numeric()
#
clinical_data <- rbind(due_date, age, ethinic, bmi, parity, sex, ivf, induction, birth_wt)
#
colnames(clinical_data) <-
  info$ID


# ##due_date, age, bmi, birth_wt
# clinical_data <-
#   clinical_data %>%
#   # clinical_data[c("due_date", "age", "bmi", "birth_wt"),] %>%
#   apply(2, as.numeric)
#
rownames(clinical_data) <-
  c("due_date",
    "age",
    "ethinic",
    "bmi",
    "parity",
    "sex",
    "ivf",
    "induction",
    "birth_wt")


clinical_data <-
  t(clinical_data) %>%
  tibble::as_tibble()

rownames(clinical_data) <-
  info$ID

clinical_data <-
  clinical_data %>%
  rownames_to_column(var = "Patient_ID")

clinical_data <-
  clinical_data %>%
  mutate(
    class = case_when(
      Patient_ID %in% sample_data_dis$Patient_ID ~ "dis",
      Patient_ID %in% sample_data_val$Patient_ID ~ "val"
    )
  )



clinical_data %>%
  group_by(class, induction) %>%
  dplyr::summarise(n = n())


clinical_data %>%
  group_by(class) %>%
  dplyr::summarise(mean = mean(as.numeric(birth_wt), na.rm = TRUE),
                   sd = sd(as.numeric(birth_wt), na.rm = TRUE)) %>%
  view()


clinical_data %>%
  group_by(class) %>%
  dplyr::summarise(mean = mean(as.numeric(due_date) * 7),
                   sd = sd(as.numeric(due_date) * 7)) %>%
  view()


wilcox.test(as.numeric(clinical_data$age[clinical_data$class == "dis"]),
            as.numeric(clinical_data$age[clinical_data$class == "val"]))

wilcox.test(as.numeric(clinical_data$bmi[clinical_data$class == "dis"]),
            as.numeric(clinical_data$bmi[clinical_data$class == "val"]))

wilcox.test(as.numeric(clinical_data$due_date[clinical_data$class == "dis"]),
            as.numeric(clinical_data$due_date[clinical_data$class == "val"]))

wilcox.test(as.numeric(clinical_data$birth_wt[clinical_data$class == "dis"]),
            as.numeric(clinical_data$birth_wt[clinical_data$class == "val"]))


chisq.test(
  x = c(4, 12),
  y = c(3, 17),
  correct = TRUE,
  simulate.p.value = TRUE
)

chisq.test(
  x = c(5, 8),
  y = c(5, 5),
  correct = TRUE,
  simulate.p.value = TRUE
)

chisq.test(rbind(c(5, 8), c(5, 5)))

chisq.test(rbind(c(6, 10), c(5, 5)))

chisq.test(rbind(c(4, 12), c(3, 17)))

chisq.test(rbind(c(2, 6, 2, 0, 6, 0), c(4, 3, 7, 2, 3, 1)),
           correct = TRUE,
           simulate.p.value = TRUE)


sample_data_dis %>%
  select(Patient_ID, Sex) %>%
  distinct() %>%
  pull(Sex) %>%
  table()

sample_data_val %>%
  select(Patient_ID, Sex) %>%
  distinct() %>%
  pull(Sex) %>%
  table()


clinical_data %>%
  group_by(class) %>%
  dplyr::summarise(
    n1 = sum(parity == 0),
    n2 = sum(parity == 1),
    n3 = sum(parity >= 2)
  )

save(clinical_data, file = "clinical_data")


clinical_data <-
  rmse_r2 %>%
  left_join(clinical_data, by = c("ID" = "Patient_ID"))

save(clinical_data, file = "clinical_data")

xlsx::write.xlsx(clinical_data, "clinical_data.xlsx")

cor.test(clinical_data$r2, as.numeric(clinical_data$bmi))


cor.test(clinical_data$r2, as.numeric(clinical_data$due_date))

clinical_data %>%
  mutate(due_date = as.numeric(due_date)) %>%
  ggplot(aes(due_date, r2)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", , fill = "grey") +
  theme_bw() +
  labs(x = "Due date (days)", y = "Adjusted R2") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  geom_text(
    x = -Inf,
    y = Inf,
    label = "Correlation: 0.54\nP-value < 0.05",
    hjust = -1,
    vjust = 1,
    size = 5
  )


cor.test(clinical_data$r2, as.numeric(clinical_data$age))



clinical_data %>%
  mutate(due_date = as.numeric(age)) %>%
  ggplot(aes(due_date, r2)) +
  geom_point(size = 2, colour = "#FFA319FF") +
  geom_smooth(colour = "#8A9045FF", fill = "grey") +
  theme_bw() +
  labs(x = "Age (year)", y = "Adjusted R2") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  geom_text(
    x = -Inf,
    y = Inf,
    label = "Correlation: -0.46\nP-value < 0.05",
    hjust = -1,
    vjust = 1,
    size = 5
  )



cor.test(clinical_data$r2, as.numeric(clinical_data$parity))




clinical_data$rmse

clinical_data$r2



colnames(clinical_data)


plot1 <-
  clinical_data %>%
  ggplot(aes(x = ethinic, y = rmse)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "Ethinicity", y = "RMSE") +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )


plot1

ggsave(plot1,
       filename = "ethinicity_vs_RMSE.pdf",
       width = 7,
       height = 7)


plot2 <-
  clinical_data %>%
  ggplot(aes(x = ethinic, y = r2)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "Ethinicity", y = "R2") +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )
plot2
ggsave(plot2,
       filename = "ethinicity_vs_r2.pdf",
       width = 7,
       height = 7)


plot <-
  cowplot::plot_grid(
    plot1,
    plot2,
    align = "h",
    ncol = 2,
    labels = c("a", "b")
  )

ggsave(plot,
       filename = "ethinicity_vs_rmse_r2.pdf",
       width = 7,
       height = 3.5)

plot(clinical_data$age)

plot(clinical_data$bmi)

plot(clinical_data$parity)

plot(clinical_data$birth_wt)


colnames(clinical_data)

temp_data <-
  clinical_data %>%
  dplyr::select(rmse, r2, age, bmi, parity, birth_wt) %>%
  apply(2, as.numeric) %>%
  Hmisc::rcorr(type = "spearman")

temp_cor <- temp_data$r
temp_p <- temp_data$P


temp_cor <-
  temp_cor[-c(1, 2), -c(3, 4, 5, 6)] %>%
  as.data.frame() %>%
  rownames_to_column("clinical") %>%
  tidyr::pivot_longer(cols = -clinical,
                      names_to = "class",
                      values_to = "cor")

temp_p <-
  temp_p[-c(1, 2), -c(3, 4, 5, 6)] %>%
  as.data.frame() %>%
  rownames_to_column("clinical") %>%
  tidyr::pivot_longer(cols = -clinical,
                      names_to = "class",
                      values_to = "p")


temp_data <-
  temp_cor %>%
  left_join(temp_p, by = c("clinical", "class"))





temp_data %>%
  ggplot(aes(
    x = clinical,
    y = class,
    colour = cor,
    size = -log(p, 10)
  )) +
  geom_point() +
  # geom_text(aes(text = ))
  labs(x = "", y = "") +
  scale_colour_gradientn(colours = c(
    alpha("#84D7E1FF", 1),
    alpha("#84D7E1FF", 0.4),
    alpha("#C71000FF", 0.4),
    alpha("#C71000FF", 1)
  )) +
  scale_size_continuous(range = c(5, 20)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.position = "top"
  )

##dark
temp_data %>%
  ggplot(aes(
    x = clinical,
    y = class,
    colour = cor,
    size = -log(p, 10)
  )) +
  geom_point() +
  # guides(size = guide_legend(override.aes = list(color = "black"))) +
  # geom_text(aes(text = ))
  labs(x = "", y = "") +
  scale_colour_gradientn(colours = c(
    alpha("#84D7E1FF", 1),
    alpha("#84D7E1FF", 0.4),
    alpha("#C71000FF", 0.4),
    alpha("#C71000FF", 1)
  )) +
  scale_x_discrete(label = c("Age", "Birth weight", "BMI", 'Parity')) +
  scale_y_discrete(label = c("RMSE", "R2")) +
  scale_size_continuous(range = c(5, 20)) +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.position = "top",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "grey")
  )

ggsave(
  filename = "test.png",
  width = 10,
  height = 5,
  bg = "transparent"
)


temp_data[-c(1, 2), -c(3, 4, 5)] %>%
  pheatmap::pheatmap(
    color = colorRampPalette(c("#84D7E1FF", "white", "#C71000FF"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    display_numbers = TRUE
  )




clinical_data %>%
  ggplot(aes(x = rmse, y = r2)) +
  geom_point()





setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/prediction/metabolites/RF/GA_prediction/")
marker_rf <- readr::read_csv("marker_rf_final.csv")
test <-
  marker_rf %>%
  dplyr::mutate(super_class = case_when(is.na(super_class) ~ "Unknown", TRUE ~ super_class)) %>%
  dplyr::group_by(super_class) %>%
  dplyr::summarise(important_ratio = sum(mean_importance), n = n()) %>%
  mutate(x = important_ratio / sum(important_ratio))

df2 <-
  test %>%
  dplyr::select(super_class, x) %>%
  mutate(browser = super_class, share = x * 100) %>%
  select(browser, share) %>%
  as.data.frame()

library(webr)

PieDonut(
  df2,
  mapping = aes(browser, share),
  r0 = 0.4,
  start = 3 * pi / 2,
  labelpositionThreshold = 0.1
)


colour <<-
  c(
    # "Clinical information" = "#FF6F00FF",
    # "Alkaloids and derivatives" = "#ADE2D0FF",
    # "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    # "Organic acids and derivatives" = "#8A4198FF",
    # "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    # "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

df2 %>%
  dplyr::arrange(share) %>%
  dplyr::mutate(class = factor(browser, levels = browser)) %>%
  ggplot(aes(x = "", y = share, fill = class)) +
  geom_bar(width = 1,
           stat = "identity",
           position = position_stack()) +
  scale_fill_manual(values = colour) +
  # scale_x_continuous(limits = c(-20, 100)) +
  coord_polar(theta = "y") +
  theme_void()
