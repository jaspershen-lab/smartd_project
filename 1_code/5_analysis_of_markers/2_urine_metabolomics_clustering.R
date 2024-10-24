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

dir.create("3_data_analysis/5_analysis_of_markers/2_urine_metabolomics_clustering", recursive = TRUE)
setwd("3_data_analysis/5_analysis_of_markers/2_urine_metabolomics_clustering/")

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

#-------------------------------------------------------------------------------
### fuzzy C-means for metabolites
library(Mfuzz)
library(e1071)

marker_data <-
  t(marker_data)

mestimate <- function(df) {
  N <- dim(df)[[1]]
  D <- dim(df)[[2]]
  m.sj <- 1 + (1418 / N + 22.05) * D ^ (-2) + (12.33 / N + 0.243) * D ^
    (-0.0406 * log(N) - 0.1134)
  return(m.sj)
}

m <- mestimate(marker_data)
m

library(e1071)

# helper function for the within sum of squared error
sumsqr <- function(x, clusters) {
  sumsqr <- function(x)
    sum(scale(x, scale = FALSE) ^ 2)
  wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
  return(wss)
}

# get the wss for repeated clustering
iterate_fcm_WSS <- function(df, m) {
  totss <- numeric()
  for (i in 2:20) {
    FCMresults <- cmeans(df, centers = i, m = m)
    totss[i] <- sum(sumsqr(df, FCMresults$cluster))
  }
  return(totss)
}

wss_2to20 <- iterate_fcm_WSS(marker_data, m)

plot(1:20,
     wss_2to20[1:20],
     type = "b",
     xlab = "Number of Clusters",
     ylab = "wss")

data.frame(number = 1:20,
           error = wss_2to20[1:20],
           stringsAsFactors = FALSE) %>%
  ggplot(aes(number, error)) +
  geom_point(shape = 16, colour = "black") +
  geom_line(colour = "black") +
  labs(x = "Number of clusters", y = "WSS")

k <- 2
fcm_results <- cmeans(marker_data, centers = k, m = m)
# get the centroids into a long dataframe:
fcm_centroids <- fcm_results$centers
fcm_centroids_df <- data.frame(fcm_centroids, check.names = FALSE, check.rows = FALSE)
fcm_centroids_df$cluster <- row.names(fcm_centroids_df)

centroids_long <-
  fcm_centroids_df %>%
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "GA",
                      values_to = "value")

ga_level <- unique(centroids_long$GA)

# ga_level <-
#   c("[11,13]", ga_level[-c(15, 16)], "PP")

centroids_long$GA <- factor(centroids_long$GA, levels = ga_level)

ggplot(centroids_long,
       aes(
         x = GA,
         y = value,
         group = cluster,
         colour = as.factor(cluster)
       )) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  ggsci::scale_color_aaas() +
  labs(
    x = "",
    title = "",
    color = "Cluster",
    y = "Scaled intensity"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      vjust = 1,
      hjust = 1,
      colour = text_colour
    ),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

# start with the input data
fcm_plotting_df <- data.frame(marker_data, check.names = FALSE, check.rows = FALSE)

# add genes
fcm_plotting_df$metabolite <- row.names(fcm_plotting_df)

# bind cluster assinment
fcm_plotting_df$cluster <- fcm_results$cluster
# fetch the membership for each gene/top scoring cluster
fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster), function(row) {
  clust <- fcm_plotting_df$cluster[row]
  fcm_results$membership[row, clust]
})



k_to_plot <- 2

# subset the dataframe by the cluster and get it into long form
# using a little tidyr action
cluster_plot_df <-
  dplyr::filter(fcm_plotting_df, cluster == k_to_plot) %>%
  dplyr::select(everything(), membership, metabolite) %>%
  tidyr::pivot_longer(
    cols = -c(metabolite, cluster, membership),
    names_to = "GA",
    values_to = "value"
  )

# order the dataframe by score
cluster_plot_df <- cluster_plot_df[order(cluster_plot_df$membership), ]
# set the order by setting the factors using forcats
cluster_plot_df$metabolite <- forcats::fct_inorder(cluster_plot_df$metabolite)

cluster_plot_df$GA <- factor(cluster_plot_df$GA, levels = ga_level)

# subset the cores by cluster
core <- dplyr::filter(centroids_long, cluster == k_to_plot)

plot <-
  cluster_plot_df %>%
  dplyr::filter(membership > 0.8) %>%
  ggplot(aes(x = GA, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_line(aes(color = membership, group = metabolite), size = 1) +
  geom_point(
    aes(fill = membership, group = metabolite),
    size = 2,
    shape = 21,
    color = "black"
  ) +
  scale_fill_gradientn(colours = c("blue1", "red2")) +
  scale_color_gradientn(colours = c("blue1", "red2")) +
  # this adds the core
  geom_line(
    data = core,
    aes(GA, value, group = cluster),
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Scaled intensity", title = "") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      vjust = 1,
      hjust = 1,
      colour = text_colour
    ),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

plot

ggsave(
  filename = "cluster2.pdf",
  height = 7,
  width = 7,
  bg = "transparent"
)



k_to_plot <- 1

# subset the dataframe by the cluster and get it into long form
# using a little tidyr action
cluster_plot_df <-
  dplyr::filter(fcm_plotting_df, cluster == k_to_plot) %>%
  dplyr::select(everything(), membership, metabolite) %>%
  tidyr::pivot_longer(
    cols = -c(metabolite, cluster, membership),
    names_to = "GA",
    values_to = "value"
  )

# order the dataframe by score
cluster_plot_df <- cluster_plot_df[order(cluster_plot_df$membership), ]
# set the order by setting the factors using forcats
cluster_plot_df$metabolite <- forcats::fct_inorder(cluster_plot_df$metabolite)

cluster_plot_df$GA <- factor(cluster_plot_df$GA, levels = ga_level)

# subset the cores by cluster
core <- dplyr::filter(centroids_long, cluster == k_to_plot)

plot <-
  cluster_plot_df %>%
  dplyr::filter(membership > 0.8) %>%
  ggplot(aes(x = GA, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_line(aes(color = membership, group = metabolite), size = 1) +
  geom_point(
    aes(fill = membership, group = metabolite),
    size = 2,
    shape = 21,
    color = "black"
  ) +
  scale_fill_gradientn(colours = c("blue1", "red2")) +
  scale_color_gradientn(colours = c("blue1", "red2")) +
  # this adds the core
  geom_line(
    data = core,
    aes(GA, value, group = cluster),
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Scaled intensity", title = "") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      vjust = 1,
      hjust = 1,
      colour = text_colour
    ),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

plot

ggsave(
  filename = "cluster1.pdf",
  height = 7,
  width = 7,
  bg = "transparent"
)

cluster1_metabolite <-
  fcm_plotting_df$metabolite[which(fcm_plotting_df$cluster == 1)]

cluster2_metabolite <-
  fcm_plotting_df$metabolite[which(fcm_plotting_df$cluster == 2)]

write.csv(fcm_plotting_df, "fuzzy_c_means_cluster.csv", row.names = FALSE)


