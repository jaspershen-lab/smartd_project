##avoid source
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)

load("3_data_analysis/1_data_preparation/1_urine_metabolomics_data/metabolites/urine_metabolomics_data.rda")

urine_metabolomics_data <-
  urine_metabolomics_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(class != "QC")

marker1 <- readr::read_csv("3_data_analysis/4_prediction/2_urine_metabolomics_ga_prediction_rf/marker_rf_final.csv")
marker2 <- readr::read_csv(
  "3_data_analysis/4_prediction/3_urine_metabolomics_time_to_delivery_prediction_rf/remove_cs/marker_rf_final.csv"
)

dir.create("3_data_analysis/5_analysis_of_markers/1_urine_metabolomics/", recursive = TRUE)
setwd("3_data_analysis/5_analysis_of_markers/1_urine_metabolomics/")

marker1$name
marker2$name

intersect(marker1$name, marker2$name)

setdiff(marker1$name, marker2$name) %>%
  data.frame(name = ., stringsAsFactors = FALSE) %>%
  left_join(urine_metabolomics_data@variable_info, by = c("name" = "variable_id")) %>%
  pull(Compound.name)

setdiff(marker2$name, marker1$name) %>%
  data.frame(name = ., stringsAsFactors = FALSE) %>%
  left_join(urine_metabolomics_data@variable_info, by = c("name" = "variable_id")) %>%
  pull(Compound.name)

library(VennDiagram)

venn <- VennDiagram::venn.diagram(
  x = list("Mraker 1" = marker1$name, "Mraker 2" = marker2$name),
  filename = NULL,
  col = c("#C71000FF", "#84D7E1FF"),
  main.cex = 1.5,
  lwd = 2
)

grid.draw(venn)

setdiff(marker1$Compound.name, marker2$Compound.name)
setdiff(marker2$Compound.name, marker1$Compound.name)

library(waffle)

test <-
  dplyr::full_join(marker1[, c("name", "super_class")], marker2[, c("name", "super_class")], by = c("name", "super_class")) %>%
  pull(super_class)

test[is.na(test)] <- "Unknown"

test <- table(test)
test1 <- as.numeric(test)
names(test1) <- names(test)

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

plot <-
  waffle(test, colors = values, rows = 6, )

plot

ggsave(
  filename = "marker_super_class.pdf",
  bg = "transparent",
  width = 7,
  height = 4
)

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

## rename P
# str_sort(grep("P", colnames(marker_data), value = T), numeric = TRUE)

sample_info <-
  extract_sample_info(urine_metabolomics_data)

sample_info %>%
  dplyr::filter(is.na(GA)) %>%
  pull(subject_id) %>%
  str_sort(numeric = TRUE)

colnames(marker_data) == sample_info$sample_id

# ## So P1 is X178 and P2 is X179 and so on.
# colnames(marker_data)[grep("P", colnames(marker_data))] <-
#   str_replace(colnames(marker_data)[grep("P", colnames(marker_data))], "P", "") %>%
#   as.numeric() %>%
#   `+`(., 177) %>%
#   paste("X", ., sep = "")
# 
# marker_data <-
#   marker_data %>%
#   dplyr::select(one_of(sample_info$sample_id))

library(pheatmap)

sum(is.na(marker_data))

marker_data <-
  apply(marker_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

ga <- sample_info$GA[match(colnames(marker_data), sample_info$sample_id)]

ga[is.na(ga)] <- 50

range(ga)

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

# "ward.D", "ward.D2",
# "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

# "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
annotation_col <-
  marker %>%
  select(Compound.name, super_class) %>%
  mutate(super_class = case_when(is.na(super_class) ~ "Unknown", TRUE ~ super_class)) %>%
  column_to_rownames(var = "Compound.name")

anno_colors <-
  list(
    super_class = c(
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
  )

plot <-
  pheatmap(
    marker_data,
    color = colorRampPalette(c("#84D7E1FF", "white", "#C71000FF"))(100),
    border_color = "grey",
    clustering_method = "average",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    angle_col = 45,
    annotation_col = annotation_col,
    annotation_colors = anno_colors,
    col.axis = "red", cutree_rows = 2,
    # display_numbers = TRUE
    # legend_breaks = c(-1,0,1)
  )

library(ggplotify)
plot <- as.ggplot(plot)
plot
ggsave(plot,
       filename = "heatmap_marker.pdf",
       width = 10,
       height = 7)
##dark
plot <-
  pheatmap(
    marker_data,
    color = colorRampPalette(c("#84D7E1FF", "white", "#C71000FF"))(100),
    border_color = "grey",
    clustering_method = "mcquitty",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    angle_col = 45,
    annotation_col = annotation_col,
    annotation_colors = anno_colors,
    col.axis = "red"
    # display_numbers = TRUE
    # legend_breaks = c(-1,0,1)
  )
plot
rownames <- sort(rownames(marker_data))
rownames <- c("PP", rownames[-c(14)])

marker_data <-
  marker_data[rownames, ]

library(ggcorrplot)

corr <-
  cor(t(marker_data))

p.mat <- cor_pmat(t(marker_data))

library(corrplot)

plot <-
  ggcorrplot(
    corr = as.matrix(corr),
    legend.title = "Correlation",
    hc.order = FALSE,
    outline.col = "grey",
    type = "upper",
    ggtheme = ggplot2::theme_bw,
    colors = c(
      "#008EA0B2",
      alpha(colour = "#008EA0B2", alpha = 1),
      alpha(colour = "#FF6F00B2", alpha = 1),
      "#FF6F00B2"
    ),
    lab = TRUE,
    lab_col = "white",
    lab_size = 3,
    p.mat = p.mat,
    sig.level = 0.05,
    insig = "pch",
    pch.col = "red"
  )

text_colour <- c("red", colorRampPalette(colors = c(
  alpha("#155F83FF", 1),
  alpha("#155F83FF", 0.4),
  alpha("#FFA319FF", 0.4),
  alpha("#FFA319FF", 1)
))(13))

plot <-
  plot +
  theme(
    axis.text = element_text(colour = text_colour),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

distance <-
  as.matrix(dist(marker_data))[, 1][-1]

line_cols <-
  colorRampPalette(colors =
                     c("royalblue", "red"))(13)

p.value <- -log(p.mat[-1, 1], 10)

pp_data <-
  data.frame(
    x = rep(13, 13),
    y = rep(3, 13),
    xend = c(1:13) + 0.5,
    yend = c(1:13) - 0.5,
    curvature = seq(
      from = -0.14,
      to = 0.14,
      length.out = 13
    ),
    colour = line_cols,
    p.value = p.value,
    distance = distance,
    text_colour = text_colour[-1],
    stringsAsFactors = FALSE
  )

for (i in 1:nrow(pp_data)) {
  c(
    plot <-
      plot +
      annotate(
        geom = "curve",
        x = pp_data$x[i],
        y = pp_data$y[i],
        xend = pp_data$xend[i],
        yend = pp_data$yend[i],
        curvature = pp_data$curvature[i],
        colour = pp_data$colour[i],
        size = 3
      ) +
      annotate(
        geom = "point",
        x = i + 0.5,
        y = i - 0.5,
        colour = pp_data$text_colour[i],
        size = 3
      )
  )
}

plot <-
  plot +
  annotate(
    geom = "point",
    x = 13,
    y = 3,
    colour = "grey",
    size = 5
  )

corr_plot <- plot
corr_plot

# ##dark
# corr_plot <-
#   corr_plot +
#   ggdark::dark_mode() +
#   theme(
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA),
#     axis.text = element_text(colour = text_colour),
#     panel.grid.minor = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
#   ) +
#   labs(x = "", y = "")

ggsave(
  corr_plot,
  filename = "corr_plot.pdf",
  height = 7,
  width = 7,
  bg = "transparent"
)

## legend for corrplot
ggplot(pp_data) +
  geom_segment(aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend,
    curvature = curvature,
    colour = distance,
  ),
  size = 3) +
  scale_colour_gradientn(colours = pp_data$colour)


