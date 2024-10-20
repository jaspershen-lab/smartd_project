####for each group
sxtTools::setwd_project()
setwd("data_analysis20200108/")
rm(list = ls())
load("data_preparation_for_analysis/metabolite_table")
load("data_preparation_for_analysis/metabolite_tags")
load("data_preparation_for_analysis/peak_table")

info <-
  readxl::read_xlsx("/Users/shenxt/projects/smartD/patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

sample_info <-
  readr::read_csv("/Users/shenxt/projects/smartD/patient information/sample_info_191021.csv")

sxtTools::setwd_project()
marker1 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/GA_prediction/marker_rf_final.csv")
marker2 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/marker_rf_final.csv")

sxtTools::setwd_project()
setwd("data_analysis20200108/biological_analysis/")


marker <-
  rbind(marker1, marker2)


marker <-
  marker %>%
  distinct(name, .keep_all = TRUE)


marker_data <-
  metabolite_table %>%
  dplyr::filter(name %in% marker$name) %>%
  column_to_rownames("name") %>%
  dplyr::select(-c(mz, rt))

## rename P
str_sort(grep("P", colnames(marker_data), value = T), numeric = TRUE)
sample_info %>%
  dplyr::filter(is.na(GA)) %>%
  pull(Sample_ID) %>%
  str_sort(numeric = TRUE)

## So P1 is X178 and P2 is X179 and so on.
colnames(marker_data)[grep("P", colnames(marker_data))] <-
  str_replace(colnames(marker_data)[grep("P", colnames(marker_data))], "P", "") %>%
  as.numeric() %>%
  `+`(., 177) %>%
  paste("X", ., sep = "")

marker_data <-
  marker_data %>%
  select(one_of(sample_info$Sample_ID))

sum(is.na(marker_data))

marker_data <-
  apply(marker_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

ga <- sample_info$GA[
  match(colnames(marker_data), sample_info$Sample_ID)
]

ga[is.na(ga)] <- 50

range(ga)

ga <- cut_width(x = ga[ga != 0], width = 2)

ga <-
  as.character(ga)


ga[ga == "(49,51]"] <- "PP"


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
    apply(x[, -1], 2, function(y) sd(y)/(length(y) - 1))
  }) %>%
  do.call(rbind, .)


colnames(marker_data_mean) <-
  metabolite_tags$Compound.name[
    match(
      colnames(marker_data_mean),
      metabolite_tags$name
    )
  ]

colnames(marker_data_sem) <-
  metabolite_tags$Compound.name[
    match(
      colnames(marker_data_sem),
      metabolite_tags$name
    )
  ]

marker_data_mean <-
  t(marker_data_mean) %>% 
  as.data.frame()


marker_data_sem <-
  t(marker_data_sem) %>% 
  as.data.frame()

marker_data_mean <- 
marker_data_mean %>% 
  tibble::rownames_to_column(var = "metabolite") %>% 
  tidyr::pivot_longer(cols = -metabolite, names_to = "ga", values_to = "mean") 

marker_data_sem <- 
  marker_data_sem %>% 
  tibble::rownames_to_column(var = "metabolite") %>% 
  tidyr::pivot_longer(cols = -metabolite, names_to = "ga", values_to = "sem") 

marker_data <- 
  dplyr::left_join(marker_data_mean, marker_data_sem, by = c("metabolite", "ga"))

cluster <- readr::read_csv("fuzzy_c_means_cluster.csv")

cluster1_metabolite <-
  cluster$metabolite[which(cluster$cluster == 1)]

cluster2_metabolite <-
  cluster$metabolite[which(cluster$cluster == 2)]


marker <- rbind(marker1, marker2)

marker_data <- 
marker_data %>% 
  dplyr::left_join(cluster[,c("metabolite", "cluster")], by = "metabolite") %>% 
  dplyr::left_join(marker[,c("Compound.name", "super_class")], by = c("metabolite" = "Compound.name"))


marker_data$super_class[is.na(marker_data$super_class)] <- "Unknown"

ga_level <- unique(marker_data$ga) %>% stringr::str_sort(numeric = TRUE)

ga_level[15] <- "(11,13]"

ga_level <- ga_level %>% stringr::str_sort(numeric = TRUE)

ga_level[1] <- "[11,13]"

marker_data <- 
  marker_data %>% 
  dplyr::mutate(ga = factor(ga, levels = ga_level))


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
  dplyr::mutate(metabolite = factor(metabolite,
                                    levels = unique(metabolite)))

temp_data$super_class


plot1 <- 
temp_data %>% 
  ggplot(aes(ga, mean, group = metabolite)) +
  annotate("rect", 
           xmin= -Inf, 
           xmax=15, 
           ymin= -Inf, 
           ymax = Inf, 
           fill="grey", 
           alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem, 
                    color = super_class), show.legend = FALSE, 
                width = 0) +
  geom_point(aes(color = super_class), 
             shape = 16,
             show.legend = FALSE) +
  geom_line(aes(color = super_class), show.legend = FALSE) +
  scale_color_manual(values = values) +
  theme_bw() +
  facet_wrap(vars(metabolite)) +
  scale_x_discrete(breaks = c("[11,13]", "(19,21]", "(29,31]", "PP"), 
                     labels = c("[11,13]", "(19,21]", "(29,31]", "PP")) +
  labs(x = "", y = "Scaled intensity") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1, 
                                   size = 8),
        axis.text.y = element_text(size = 8), 
        strip.background = element_rect(fill = c("#FF6F00B2")),
        strip.text = element_text(size = 8, color = "white"),
        panel.grid = element_blank())

plot1 


ggsave(plot1, filename = "cluster1_plot.pdf", width = 7, height = 3.5)



temp_data <- 
  marker_data %>% 
  dplyr::filter(cluster == 2)  %>% 
  dplyr::arrange(super_class) %>% 
  dplyr::mutate(metabolite = factor(metabolite,
                                    levels = unique(metabolite)))

plot2 <- 
  temp_data %>% 
  ggplot(aes(ga, mean, group = metabolite)) +
  annotate("rect", 
           xmin= -Inf, 
           xmax=15, 
           ymin= -Inf, 
           ymax = Inf, 
           fill="grey", 
           alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem, 
                    color = super_class), show.legend = FALSE, 
                width = 0) +
  geom_point(aes(color = super_class), shape = 16, 
             show.legend = FALSE) +
  geom_line(aes(color = super_class), show.legend = FALSE) +
  scale_color_manual(values = values) +
  theme_bw() +
  facet_wrap(vars(metabolite)) +
  scale_x_discrete(breaks = c("[11,13]", "(19,21]", "(29,31]", "PP"), 
                   labels = c("[11,13]", "(19,21]", "(29,31]", "PP")) +
  labs(x = "", y = "Scaled intensity") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1, 
                                   size = 8),
        axis.text.y = element_text(size = 8), 
        strip.background = element_rect(fill = c("#FF6F00B2")),
        strip.text = element_text(size = 8, color = "white"),
        panel.grid = element_blank())

plot2



ggsave(plot2, filename = "cluster2_plot.pdf", width = 7, height = 5)





