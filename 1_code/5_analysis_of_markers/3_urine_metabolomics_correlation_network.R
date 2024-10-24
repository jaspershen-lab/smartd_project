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

dir.create(
  "3_data_analysis/5_analysis_of_markers/3_urine_metabolomics_correlation_network",
  recursive = TRUE
)
setwd(
  "3_data_analysis/5_analysis_of_markers/3_urine_metabolomics_correlation_network"
)

marker1$name
marker2$name

intersect(marker1$name, marker2$name)

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

marker_data <-
  data.frame(ga, t(marker_data), stringsAsFactors = FALSE) %>%
  split(f = ga) %>%
  lapply(function(x) {
    apply(x[, -1], 2, mean)
  }) %>%
  do.call(rbind, .)

range(marker_data)

marker_data[which(marker_data < -1.43, arr.ind = TRUE)] <- -1.43

variable_info <-
  extract_variable_info(urine_metabolomics_data)

colnames(marker_data) <-
  variable_info$Compound.name[match(colnames(marker_data), variable_info$variable_id)]

text_colour <- c("red", colorRampPalette(colors = c(
  alpha("#155F83FF", 1),
  alpha("#155F83FF", 0.4),
  alpha("#FFA319FF", 0.4),
  alpha("#FFA319FF", 1)
))(13))


## correlation network
## this networks should contain metabolite and clinical information
marker_data <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% marker$name) %>%
  extract_expression_data()

sum(is.na(marker_data))

marker_data <-
  apply(marker_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

subject_id <-
  sample_info$subject_id[match(colnames(marker_data), sample_info$sample_id)]

cor_met_data <-
  data.frame(subject_id = subject_id,
             t(marker_data),
             stringsAsFactors = FALSE)

cor_met_data <-
  cor_met_data %>%
  plyr::dlply(
    .variables = "subject_id",
    .fun = function(x) {
      x <- x[, -1]
      x <- apply(x, 2, mean)
      x
    }
  ) %>%
  do.call(rbind, .)

colnames(cor_met_data) <-
  variable_info$Compound.name[match(colnames(cor_met_data), variable_info$variable_id)]

cor_met_data <-
  cor_met_data %>%
  as.data.frame() %>%
  rownames_to_column(var = "subject_id")

### clinical information
info <-
  sample_info %>%
  dplyr::filter(subject_id %in% cor_met_data$subject_id) %>% 
  dplyr::distinct(subject_id, .keep_all = TRUE)

info <-
  info[match(cor_met_data$subject_id, info$subject_id), ]

cor_met_data$subject_id == info$subject_id

due_date <-
  info$mother_delivery_weeks

## age
age <-
  info$mother_age %>%
  as.numeric()
age
#
## ethinic
ethinic <-
  info$mother_ethnicity

# ##BMI
bmi <- info$mother_bmi

# #------------------------------------------------------------------------------
# ##parity
parity <- info$mother_parity

# ##sex and twins
# info$Sex
sex <- info$child_sex

sex <-
  case_when(
    is.na(sex) ~ "NA",
    sex == "Female_Female" ~ "F_F",
    sex == "Male_Male" ~ "M_M",
    sex == "Male_Female" ~ "M_F",
    TRUE ~ sex
  )

## induction
induction <-
  info$mother_induction
#
induction[is.na(induction)] <- "NA"
#
#
birth_wt <-
  info$child_weight
birth_wt

birth_wt <-
  sapply(birth_wt, function(x) {
    if (!is.na(x)) {
      x <- stringr::str_split(x, "\\{\\}")[[1]] %>%
        as.numeric() %>%
        sum()
    } else {
      x
    }
    x
  })
#
birth_wt <- unname(birth_wt) %>%
  as.numeric()
#
clinical_data <- rbind(# due_date,
  age, # ethinic,
  bmi, parity, # sex,
  # ivf,
  # induction,
  birth_wt)
#
colnames(clinical_data) <-
  info$subject_id

clinical_data <-
  t(clinical_data) %>%
  as.data.frame() %>%
  rownames_to_column(var = "subject_id")

clinical_data2 <-
  clinical_data %>%
  filter(!is.na(birth_wt))

clinical_data2$bmi <-
  (clinical_data2$bmi - mean(clinical_data2$bmi)) / sd(clinical_data2$bmi)

clinical_data2$birth_wt <-
  (clinical_data2$birth_wt - mean(clinical_data2$birth_wt)) / sd(clinical_data2$birth_wt)

cor_data <-
  clinical_data2 %>%
  left_join(cor_met_data, by = "subject_id")

cross_cor <-
  cor_data[, -1] %>%
  cor(., method = "spearman")

rownames(cross_cor) == colnames(cross_cor)

cross_cor[lower.tri(cross_cor)] <- 2

cross_cor <-
  as_tibble(cross_cor)

rownames(cross_cor) <-
  colnames(cross_cor)

test <-
  cross_cor %>%
  rownames_to_column(., var = "name1") %>%
  gather(., key = "name2", value = "Correlation", -name1) %>%
  distinct() %>%
  filter(., Correlation != 1 & Correlation != 2) %>%
  arrange(., desc(abs(Correlation)))

test <-
  test %>%
  dplyr::filter(stringr::str_detect(name1, "bmi"))

test_p <-
  pbapply::pbapply(test, 1, function(x) {
    peak1 <- as.character(x[1])
    peak2 <- as.character(x[2])
    int1 <-
      pull(cor_data, var = peak1)
    int2 <-
      pull(cor_data, var = peak2)
    p <- cor.test(int1,
                  int2,
                  alternative = "two.sided",
                  method = "spearman",
                  use = "na.or.complete")$p.value
    p
  })

test <- cbind(test, test_p, p.adjust = p.adjust(test_p, "fdr"))

write.csv(test, "bmi_cor.csv")

cross_cor <-
  cross_cor %>%
  rownames_to_column(., var = "name1") %>%
  gather(., key = "name2", value = "Correlation", -name1) %>%
  distinct() %>%
  filter(., Correlation != 1 & Correlation != 2) %>%
  arrange(., desc(abs(Correlation))) %>%
  filter(abs(Correlation) > 0.5)

cross_cor_p <-
  pbapply::pbapply(cross_cor, 1, function(x) {
    peak1 <- as.character(x[1])
    peak2 <- as.character(x[2])
    int1 <-
      pull(cor_data, var = peak1)
    int2 <-
      pull(cor_data, var = peak2)
    p <- cor.test(int1,
                  int2,
                  alternative = "two.sided",
                  method = "spearman",
                  use = "na.or.complete")$p.value
    p
  })

sum(cross_cor_p < 0.05)

### p adjustment
cross_cor_p2 <-
  p.adjust(cross_cor_p, method = "BH")

cross_cor <-
  cross_cor %>%
  mutate(P.adjust = cross_cor_p2)

cross_cor <-
  cross_cor %>%
  filter(P.adjust < 0.01 & abs(Correlation) > 0.5)

cross_cor <-
  cross_cor %>%
  dplyr::rename(from = name1, to = name2)

library(ggraph)
library(tidygraph)

cross_cor$P.adjust[which(cross_cor$P.adjust == 0)] <-
  min(cross_cor$P.adjust[cross_cor$P.adjust != 0])

range(-log(cross_cor$P.adjust, 10))

cross_cor2 <-
  cross_cor %>%
  mutate(abs.cor = abs(Correlation)) %>%
  arrange(desc(abs.cor)) %>%
  select(-abs.cor)

metabolite <-
  unique(c(cross_cor2$from, cross_cor2$to))

node_attr <-
  data.frame(name = metabolite, stringsAsFactors = FALSE) %>%
  left_join(variable_info, by = c("name" = "Compound.name"))

node_attr$node_class <- rep(NA, nrow(node_attr))

node_attr$node_class <-
  case_when(is.na(node_attr$Total.score) ~ "Clinical information",
            TRUE ~ "Metabolite")

node_attr$node_class[node_attr$node_class == "Metabolite"] <-
  node_attr$super_class[node_attr$node_class == "Metabolite"]

node_attr$node_class[is.na(node_attr$node_class)] <- "Unknown"

node_attr$Total.score[is.na(node_attr$Total.score)] <-
  mean(node_attr$Total.score, na.rm = TRUE)


node_attr$sort <-
  case_when(
    node_attr$node_class == "Clinical information" ~ 1,
    node_attr$node_class == "Lipids and lipid-like molecules" ~ 2,
    node_attr$node_class == "Nucleosides, nucleotides, and analogues" ~ 3,
    node_attr$node_class == "Organic oxygen compounds" ~ 4,
    node_attr$node_class == "Organoheterocyclic compounds" ~ 5,
    node_attr$node_class == "Unknown" ~ 6
  )

node_attr <-
  node_attr %>%
  arrange(sort)

node_attr$node_class <- factor(node_attr$node_class, levels = unique(node_attr$node_class))

cross_graph <-
  tidygraph::tbl_graph(nodes = node_attr,
                       edges = cross_cor2,
                       directed = FALSE)

## angle for label
metabolite <- igraph::V(cross_graph)$name
id <- 1:length(metabolite)
angle <- 360 * (id - 0.5) / length(metabolite)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

set_graph_style(family = "Arial")
extrafont::loadfonts()

plot <-
  ggraph(cross_graph, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(
    edge_colour = Correlation,
    edge_width = -log(P.adjust, 10)
  )) +
  scale_edge_colour_gradient2(low = "royalblue",
                              mid = "white",
                              high = "red") +
  scale_edge_width_continuous(range = c(0.8, 3.5)) +
  geom_node_point(aes(size = Total.score, colour = node_class)) +
  scale_colour_manual(
    values = c(
      "Clinical information" = "#C71000FF",
      # "Alkaloids and derivatives" = "#ADE2D0FF",
      # "Benzenoids" = "#C71000FF",
      "Lipids and lipid-like molecules" = "#FF6F00B2",
      "Nucleosides, nucleotides, and analogues" = "#C71000B2",
      # "Organic acids and derivatives" = "#8A4198FF",
      # "Organic nitrogen compounds" = "#5A9599FF",
      "Organic oxygen compounds" = "#008EA0B2",
      "Organoheterocyclic compounds" = "#8A4198B2",
      # "Organosulfur compounds" = "#3F4041FF",
      # "Phenylpropanoids and polyketides" = "#5A9599B2",
      "Unknown" = "#3F4041B2"
    )
  ) +
  # ggsci::scale_color_futurama(alpha = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous(range = c(3, 10),
                        guide = guide_legend(override.aes = list(colour = "black"))) +
  scale_linewidth_continuous(range = c(0.5, 3)) +
  geom_node_text(
    aes(
      x = x * 1.1,
      y = y * 1.1,
      label = name,
      colour = node_class
    ),
    angle = angle,
    hjust = hjust,
    # colour = "black",
    size = 3.5
  ) +
  # ggdark::dark_theme_void() +
  theme_void() +
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))

plot

ggsave(
  plot,
  filename = "correlation.pdf",
  width = 10,
  height = 7,
  bg = "transparent"
)

cor_data %>%
  filter(bmi < 3) %>%
  ggplot(aes(bmi, Pregnenolone)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

cor_data %>%
  filter(bmi < 3) %>%
  ggplot(aes(bmi, Progesterone)) +
  geom_point() +
  geom_smooth() +
  theme_bw()


plot(cor_data$bmi, cor_data$Pregnenolone)

plot(cor_data$bmi, cor_data$Progesterone)


# ## pathway enrichment analysis to get the insight of which pathways are altered in pregnancy
# setwd(r4projects::get_project_wd())
# setwd("data_analysis20200108/biological_analysis/pathway_enrichment")
# marker$Compound.name
# marker$KEGG.ID
# marker$HMDB.ID
# 
# load("hsa.kegg.pathway.rda")
# 
# path_result <-
#   enrichPathway(id = marker$KEGG.ID[!is.na(marker$KEGG.ID)],
#                 database = hsa.kegg.pathway,
#                 method = "hypergeometric")
# 
# path_result %>%
#   mutate(class = case_when(p.value < 0.05 ~ "Significant", TRUE ~ "No")) %>%
#   mutate(class = factor(class, levels = c("Significant", "No"))) %>%
#   ggplot(aes(x = Overlap, y = -log(p.value, 10))) +
#   geom_hline(yintercept = 1.3, linetype = 2) +
#   geom_point(aes(size = Pathway.length, colour = class)) +
#   scale_colour_manual(values = c(
#     "Significant" = "#C71000FF",
#     "No" = "#3F4041FF"
#   )) +
#   ggplot2::guides(colour = guide_legend(override.aes = list(size = 5))) +
#   theme_bw() +
#   ggrepel::geom_label_repel(mapping = aes(
#     x = Overlap,
#     y = -log(p.value, 10),
#     label = Pathway.name
#   )) +
#   theme(
#     axis.title = element_text(size = 15),
#     axis.text.x = element_text(size = 13),
#     axis.text.y = element_text(size = 13)
#   )
# 
# path_result %>%
#   mutate(class = case_when(p.value < 0.05 ~ "Significant", TRUE ~ "No")) %>%
#   mutate(class = factor(class, levels = c("Significant", "No"))) %>%
#   arrange(desc(p.value)) %>%
#   mutate(Pathway.name = factor(Pathway.name, levels = Pathway.name)) %>%
#   ggplot(aes(x = Pathway.name, y = -log(p.value, 10))) +
#   # geom_hline(yintercept = 1.3, linetype = 2) +
#   geom_bar(aes(fill = class), stat = "identity", show.legend = FALSE) +
#   scale_fill_manual(values = c(
#     "Significant" = "#C71000FF",
#     "No" = "#3F4041FF"
#   )) +
#   labs(x = "", y = "-log10(P-value)") +
#   # ggplot2::guides(colour = guide_legend(override.aes = list(size = 5))) +
#   theme_bw() +
#   # ggrepel::geom_label_repel(mapping = aes(x = Overlap,
#   #                                         y = -log(p.value, 10),
#   #                                         label = Pathway.name)) +
#   theme(
#     axis.title = element_text(size = 15),
#     axis.text.x = element_text(size = 13),
#     axis.text.y = element_text(size = 13)
#   ) +
#   coord_flip()
