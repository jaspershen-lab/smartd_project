#to avoind source
no_exist_function()

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(tidyverse)

# load data
# load(
#   "3_data_analysis/1_urine_metabolomics_data/metabolites/urine_metabolomics_data.rda"
# )

load("3_data_analysis/1_data_preparation/0_demographic_data/sample_information.rda")

dir.create("3_data_analysis/2_data_overview/1_study_information",
           recursive = TRUE)
setwd("3_data_analysis/2_data_overview/1_study_information")

library(tidymass)

sample_info <-
  sample_information

###remove QC samples
# sample_info <-
#   sample_info %>%
#   dplyr::filter(class == "Subject")

##### the clinic information
###36 participants
sample_info$subject_id %>%
  unique() %>%
  length()

sample_info %>%
  dplyr::filter(is.na(GA)) %>%
  dplyr::select(sample_id)

library(ggExtra)

temp_data <-
  sample_info %>%
  dplyr::mutate(GA2 = case_when(is.na(GA) ~ 45, !is.na(GA) ~ GA)) %>%
  dplyr::mutate(class = case_when(is.na(GA) ~ "PP", !is.na(GA) ~ "Normal")) %>%
  dplyr::mutate(Begin.Date = as.Date(EDD) - 280) %>%
  dplyr::mutate(term.date = as.Date(DD) - Begin.Date) %>%
  dplyr::mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::arrange(diff_day)

## which person is preterm
preterm21 <-
  temp_data %>%
  # mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::filter(diff_day >= 21) %>%
  dplyr::pull(subject_id) %>%
  unique()

preterm7 <-
  temp_data %>%
  # mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::filter(diff_day >= 7 & diff_day < 21) %>%
  dplyr::pull(subject_id) %>%
  unique()

text_colour <- unique(temp_data$subject_id)

text_colour <- case_when(
  text_colour %in% preterm21 ~ "#58593FFF",
  text_colour %in% preterm7 ~ "#800000FF",
  TRUE ~ "#BEBEBE"
)

text_face <- unique(temp_data$subject_id)
text_face <- case_when(text_face %in% preterm21 ~ "plain",
                       text_face %in% preterm7 ~ "plain",
                       TRUE ~ "plain")

plot1 <-
  # temp_data %>%
  ggplot(data = temp_data, aes(x = GA2, y = factor(subject_id, level = unique(subject_id)))) +
  geom_point(
    aes(
      x = GA2,
      y = factor(subject_id, level = unique(subject_id)),
      colour = class
    ),
    data = temp_data,
    show.legend = FALSE,
    size = 2
  ) +
  scale_colour_manual(values = c("PP" = "#155F83FF", "Normal" = "#FFA319FF")) +
  labs(x = "Gestational age (GA, weeks)", y = "Subject ID") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(colour = "#8A9045FF", linetype = 1),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    axis.title = element_text(size = 13),
    axis.text.y = element_text(
      size = 12,
      colour = text_colour,
      face = text_face
    ),
    axis.text.x = element_text(size = 12)
  )

plot1

term_date <-
  temp_data %>%
  dplyr::select(subject_id, term.date) %>%
  mutate(term.date = as.numeric(term.date / 7)) %>%
  distinct()

plot1 <-
  plot1 +
  ggplot2::annotate(
    geom = "point",
    shape = 17,
    colour = "black",
    size = 2,
    x = term_date$term.date,
    y = term_date$subject_id
  )

plot1

plot2 <-
  ggplot(temp_data, aes(x = GA2)) +
  geom_histogram(
    binwidth = 0.5,
    colour = "#8DD3C7",
    fill = "#8DD3C7"
  ) +
  labs(x = "GA (weeks)", y = "Sample number") +
  theme_classic() +
  scale_x_continuous(# limits = c(10, 46),
    name = NULL,
    labels = NULL,
    breaks = NULL) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  theme(
    # panel.border = element_blank(),
    # axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )

plot2

plot3 <-
  ggplot(temp_data, aes(x = factor(subject_id, levels = unique(subject_id)))) +
  geom_bar(width = 0.8, fill = "#8DD3C7") +
  labs(x = "GA (weeks)", y = "Sample number") +
  theme_bw() +
  scale_x_discrete(name = NULL,
                   label = NULL,
                   breaks = NULL) +
  scale_y_continuous(
    expand = expand_scale(mult = c(0, .05)),
    breaks = c(0, 5, 10),
    labels = c(0, 5, 10)
  ) +
  coord_flip() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )

plot3

library(patchwork)

plot <-
  {
    plot2 + patchwork::plot_spacer() + plot_layout(ncol = 2, widths = c(4, 1))
  } -
  {
    plot1 + plot3 + plot_layout(ncol = 2, widths = c(4, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 4))


plot

ggsave(plot,
       filename = "sample_colection_distribution.pdf",
       width = 10,
       height = 7)
