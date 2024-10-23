##avoid source
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
source("1_code/100_tools.R")

##load data
##load data
load(
  "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/peaks/urine_metabolomics_data.rda"
)

dir.create(
  "3_data_analysis/3_biomakers/3_urine_metabolomics_k_means_consensus",
  recursive = TRUE
)
setwd("3_data_analysis/3_biomakers/3_urine_metabolomics_k_means_consensus")

##remove qc and blank samples
urine_metabolomics_data2 <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject")

expression_data <-
  extract_expression_data(urine_metabolomics_data2)

sample_info <-
  extract_sample_info(urine_metabolomics_data2)

variable_info <-
  extract_variable_info(urine_metabolomics_data2)

rownames(expression_data)
variable_info$variable_id

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

subject_data <-
  expression_data %>%
  dplyr::select(one_of(sample_info$sample_id))

#log transformation
subject_data <-
  log(subject_data + 1, 10)

subject_data[5, ] %>%
  as.numeric() %>%
  density() %>%
  plot()

dim(subject_data)
dim(variable_info)

rownames(subject_data) <- variable_info$variable_id

sample_info$GA[is.na(sample_info$GA)] <- 45


load("../1_urine_metabolomics_total_markers/sam_value.rda")

peak_name <-
  sam_value$name[which(sam_value$fdr < 0.05)]

temp_subject_data <-
  subject_data[peak_name, ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###k-mean
library(CancerSubtypes)

# result <-
#   CancerSubtypes::ExecuteCC(
#     clusterNum = 2,
#     d = as.matrix(temp_subject_data),
#     maxK = 6,
#     reps = 1000,
#     pItem = 0.8,
#     # pFeature = 0.8,
#     title = "k_means_consensus",
#     clusterAlg = "km",
#     distance = "euclidean",
#     plot = "png",
#     writeTable = TRUE
#   )
#
# save(result, file = "result.rda")
load("result.rda")

idx <- 3
sil = silhouette_SimilarityMatrix(result$originalResult[[idx]]$consensusClass,
                                  result$originalResult[[idx]]$consensusMatrix)
# sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)

sil_plot <-
  plot_silhouette(sil)

sil_plot

# ggsave(sil_plot,
#        file = "k_means_consensus/sil_plot3.pdf",
#        width = 7,
#        height = 7)

plot(result$originalResult[[idx]]$consensusTree)

name2 <- colnames(temp_subject_data)[result$originalResult[[idx]]$consensusTree$order]

temp_subject_data2 <- temp_subject_data[, name2]
cluster <- result$originalResult[[idx]]$consensusClass

temp_sample_info2 <-
  sample_info[match(colnames(temp_subject_data2), sample_info$sample_id), ]

cluster <- cluster[match(temp_sample_info2$sample_id, names(cluster))]

names(cluster) == colnames(temp_subject_data2)

###reorder cluster
cluster2 <- cluster[cluster == 2]
cluster3 <- cluster[cluster == 3]
cluster1 <- cluster[cluster == 1]
cluster <- c(rev(cluster2), rev(cluster3), rev(cluster1))

temp_sample_info2 <-
  temp_sample_info2[match(names(cluster), temp_sample_info2$sample_id), ]

temp_subject_data2 <-
  temp_subject_data2[, names(cluster)]

names(cluster) == temp_sample_info2$sample_id
names(cluster) == colnames(temp_subject_data2)

###complext heatamp
temp_data <- temp_subject_data2

ga <-
  temp_sample_info2 %>%
  dplyr::select(sample_id, GA) %>%
  dplyr::mutate(ga = ggplot2::cut_width(GA, width = 2, center = 3)) %>%
  dplyr::mutate(ga = stringr::str_replace(ga, "\\[", "(")) %>%
  pull(ga)

ga[ga == "(44,46]"] <- 'PP'

ga_level <- stringr::str_sort(unique(ga), numeric = TRUE)

library(ComplexHeatmap)

library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("#4292C6", "white", "red"))
col_fun(seq(-3, 3))
cluster_col_fun = colorRamp2(c(0, 2, 3), c("green", "white", "red"))

text_colour <- c(colorRampPalette(colors = c(
  alpha("#155F83FF", 1),
  alpha("#155F83FF", 0.4),
  alpha("#FFA319FF", 0.4),
  alpha("#FFA319FF", 1)
))(16), "red")

names(text_colour) <- ga_level

cluster_color <- ggsci::pal_d3()(10)[4:6]
names(cluster_color) <- c(1, 2, 3)

ha1 = HeatmapAnnotation(
  ga = factor(ga, levels = ga_level),
  cluster = factor(cluster, levels = as.character(c(1, 2, 3))),
  col = list(ga = text_colour, cluster = cluster_color),
  annotation_name_side = c("left")
)

ha2 = HeatmapAnnotation(
  "GA"  = anno_points(
    temp_sample_info2$GA[match(colnames(temp_subject_data2), temp_sample_info2$sample_id)],
    # ylim = c(0, 1),
    pch = 16,
    gp = gpar(col = c(
      rep(ggsci::pal_d3()(10)[5], length(cluster2)),
      rep(ggsci::pal_d3()(10)[6], length(cluster3)),
      rep(ggsci::pal_d3()(10)[4], length(cluster1))
    )),
    size = unit(2, "mm"),
    height = unit(4, "cm"),
    show_legend = c(TRUE, FALSE),
    axis_param = list(side = "left")
    # at = c(0, 0.5, 1),
    # labels = c("zero", "half", "one"))
  ),
  annotation_name_side = c("left")
)

temp_data <- temp_subject_data2
range((temp_data))

temp_data[temp_data > 3] <- 3
temp_data[temp_data < -3] <- -3

plot <-
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    border = FALSE,
    col = col_fun,
    name = "Int",
    clustering_method_rows = "ward.D",
    row_km = 2,
    top_annotation = ha1,
    bottom_annotation = ha2
  )

library(ggplotify)
plot <- as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "consensus_heatmap.pdf",
#        width = 12,
#        height = 7)


###test for two different clusters
index1 <- which(cluster == 1)
index2 <- which(cluster == 2)

data.frame(
  cluster = factor(cluster, levels = c(1, 2)),
  age = sample_info$mother_aboration,
  stringsAsFactors = FALSE
) %>%
  ggplot(aes(cluster, age)) +
  geom_boxplot() +
  geom_jitter()


t.test(sample_info$mother_aboration[index1],
       sample_info$mother_aboration[index2])

data.frame(
  cluster = factor(cluster, levels = c(1, 2)),
  bmi = sample_info$mother_bmi,
  stringsAsFactors = FALSE
) %>%
  ggplot(aes(cluster, bmi)) +
  geom_boxplot() +
  geom_jitter()

wilcox.test(sample_info$mother_bmi[index1], sample_info$mother_bmi[index2])

table(cluster, ethnicity2) %>%
  chisq.test()
