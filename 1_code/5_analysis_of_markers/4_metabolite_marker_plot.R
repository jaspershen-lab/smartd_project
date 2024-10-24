##avoid source
##avoid source
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)

load(
  "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/metabolites/urine_metabolomics_data.rda"
)

urine_metabolomics_data <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "QC")

marker1 <- readr::read_csv(
  "3_data_analysis/4_prediction/2_urine_metabolomics_ga_prediction_rf/marker_rf_final.csv"
)
marker2 <- readr::read_csv(
  "3_data_analysis/4_prediction/3_urine_metabolomics_time_to_delivery_prediction_rf/remove_cs/marker_rf_final.csv"
)


cluster <- readr::read_csv(
  "3_data_analysis/5_analysis_of_markers/2_urine_metabolomics_clustering/fuzzy_c_means_cluster.csv"
)

dir.create("3_data_analysis/5_analysis_of_markers/4_metabolite_marker_plot",
           recursive = TRUE)
setwd("3_data_analysis/5_analysis_of_markers/4_metabolite_marker_plot")

marker <-
  rbind(marker1, marker2)

marker <-
  marker %>%
  dplyr::distinct(name, .keep_all = TRUE)

marker_data <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% marker$name) %>%
  extract_expression_data()

sample_info <-
  extract_sample_info(urine_metabolomics_data)

colnames(marker_data) == sample_info$sample_id

marker_data <-
  apply(marker_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

ga <- sample_info$ga_range

#
# marker_data[which(marker_data < -1.43, arr.ind = TRUE)] <- -1.43

####for each group
marker_data_mean <-
  data.frame(ga, t(marker_data), stringsAsFactors = FALSE) %>%
  split(f = ga) %>%
  lapply(function(x) {
    apply(x[, -1], 2, mean)
  }) %>%
  do.call(rbind, .)

marker_data_sem <-
  data.frame(ga, t(marker_data), stringsAsFactors = FALSE) %>%
  split(f = ga) %>%
  lapply(function(x) {
    apply(x[, -1], 2, function(y)
      sd(y) / (length(y) - 1))
  }) %>%
  do.call(rbind, .)

variable_info <-
  extract_variable_info(urine_metabolomics_data)

colnames(marker_data_mean) <-
  variable_info$Compound.name[match(colnames(marker_data_mean), variable_info$variable_id)]

colnames(marker_data_sem) <-
  variable_info$Compound.name[match(colnames(marker_data_sem), variable_info$variable_id)]

marker_data_mean <-
  t(marker_data_mean) %>%
  as.data.frame()

marker_data_sem <-
  t(marker_data_sem) %>%
  as.data.frame()

marker_data_mean <-
  marker_data_mean %>%
  tibble::rownames_to_column(var = "metabolite") %>%
  tidyr::pivot_longer(cols = -metabolite,
                      names_to = "ga",
                      values_to = "mean")

marker_data_sem <-
  marker_data_sem %>%
  tibble::rownames_to_column(var = "metabolite") %>%
  tidyr::pivot_longer(cols = -metabolite,
                      names_to = "ga",
                      values_to = "sem")

marker_data <-
  dplyr::left_join(marker_data_mean, marker_data_sem, by = c("metabolite", "ga"))

cluster1_metabolite <-
  cluster$metabolite[which(cluster$cluster == 1)]

cluster2_metabolite <-
  cluster$metabolite[which(cluster$cluster == 2)]

marker <- rbind(marker1, marker2)

marker_data <-
  marker_data %>%
  dplyr::left_join(cluster[, c("metabolite", "cluster")], by = "metabolite") %>%
  dplyr::left_join(marker[, c("Compound.name", "super_class")], by = c("metabolite" = "Compound.name"))


# marker_data$super_class[is.na(marker_data$super_class)] <- "Unknown"
#
# ga_level <- unique(marker_data$ga) %>% stringr::str_sort(numeric = TRUE)
#
# ga_level[15] <- "(11,13]"
#
# ga_level <- ga_level %>% stringr::str_sort(numeric = TRUE)
#
# ga_level[1] <- "[11,13]"

marker_data <-
  marker_data %>%
  dplyr::mutate(ga = factor(ga, levels = stringr::str_sort(unique(ga), numeric = TRUE)))

values <- c(
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

temp_data <-
  marker_data %>%
  dplyr::filter(cluster == 1) %>%
  dplyr::arrange(super_class) %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite)))

temp_data$super_class

plot1 <-
  temp_data %>%
  ggplot(aes(ga, mean, group = metabolite)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 14,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey",
    alpha = 0.7
  ) +
  geom_point(aes(color = super_class),
             shape = 16,
             show.legend = FALSE) +
  geom_line(aes(color = super_class), show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    color = "black",
    show.legend = FALSE,
    width = 0
  ) +
  scale_color_manual(values = values) +
  theme_bw() +
  facet_wrap(vars(metabolite)) +
  # scale_x_discrete(
  #   breaks = c("[11,13]", "(19,21]", "(29,31]", "PP"),
  #   labels = c("[11,13]", "(19,21]", "(29,31]", "PP")
  # ) +
  labs(x = "", y = "Scaled intensity") +
  theme(
    axis.title = element_text(size = 10),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = c("#FF6F00B2")),
    strip.text = element_text(size = 8, color = "white"),
    panel.grid = element_blank()
  )

plot1

ggsave(plot1,
       filename = "cluster1_plot.pdf",
       width = 10,
       height = 7)

temp_data <-
  marker_data %>%
  dplyr::filter(cluster == 2)  %>%
  dplyr::arrange(super_class) %>%
  dplyr::mutate(metabolite = factor(metabolite, levels = unique(metabolite)))

plot2 <-
  temp_data %>%
  ggplot(aes(ga, mean, group = metabolite)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 14,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey",
    alpha = 0.7
  ) +
  geom_point(aes(color = super_class),
             shape = 16,
             show.legend = FALSE) +
  geom_line(aes(color = super_class), show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    color = "black",
    show.legend = FALSE,
    width = 0
  ) +
  scale_color_manual(values = values) +
  theme_bw() +
  facet_wrap(vars(metabolite)) +
  # scale_x_discrete(
  #   breaks = c("[11,13]", "(19,21]", "(29,31]", "PP"),
  #   labels = c("[11,13]", "(19,21]", "(29,31]", "PP")
  # ) +
  labs(x = "", y = "Scaled intensity") +
  theme(
    axis.title = element_text(size = 10),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = c("#FF6F00B2")),
    strip.text = element_text(size = 8, color = "white"),
    panel.grid = element_blank()
  )

plot2

ggsave(plot2,
       filename = "cluster2_plot.pdf",
       width = 10,
       height = 5)
