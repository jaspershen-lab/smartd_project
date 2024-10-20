##avoid source 
no_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)
source("R/20200727/tools.R")

##load data
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/metabolites/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/metabolites/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/metabolites/variable_info")

setwd("data_analysis20200108/urine_metabolome/k_means_clustering/")

dim(expression_data)

rownames(expression_data)
variable_info$name

rownames(expression_data) == variable_info$name
colnames(expression_data) == sample_info$sample_id

##remove QC and blank samples
sample_info <-
  sample_info %>%
  dplyr::filter(!stringr::str_detect(sample_id, "QC")) 

expression_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

expression_data[5,] %>% 
  as.numeric() %>% 
  density() %>%
  plot()

subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

#log transformation
subject_data <-
  log(subject_data + 1, 10)

dim(subject_data)
dim(variable_info)

rownames(subject_data) <- variable_info$name

sample_info$GA[is.na(sample_info$GA)] <- 45

##ga range
dim(sample_info)
dim(subject_data)

sample_info$sample_id == colnames(subject_data)

ga_range <- 
  as.character(cut_width(x = sample_info$GA, width = 2, center = 3)) %>%
  stringr::str_replace("\\[", '(')

stringr::str_sort(unique(ga_range))

table(ga_range)

## for each person, just combine the samples in the sample ga range
sample_info$ga_range <-
  ga_range

sample_info$ga_range %>% unique() %>% sort()

sample_info$ga_range %>% table()

sample_info$ga_range[sample_info$ga_range == "(10,12]"] <- "(10,16]"
sample_info$ga_range[sample_info$ga_range == "(12,14]"] <- "(10,16]"
sample_info$ga_range[sample_info$ga_range == "(14,16]"] <- "(10,16]"

sample_info$ga_range[sample_info$ga_range == "(38,40]"] <- "(38,42]"
sample_info$ga_range[sample_info$ga_range == "(40,42]"] <- "(38,42]"

sample_info$ga_range[sample_info$ga_range == "(44,46]"] <- "PP"

###combine different samples in one ga range together
library(plyr)

subject_data <- 
  apply(subject_data, 1, function(x){
    (x) / sd(x)
  })

subject_data2 <-
  subject_data %>%
  data.frame(., ga_range = sample_info$ga_range, stringsAsFactors = FALSE) %>%
  mutate(ga_range = factor(ga_range,levels = sample_info$ga_range %>% unique() %>% stringr::str_sort(numeric = TRUE))) %>%
  plyr::dlply(.variables = .(ga_range))

# subject_data_mean <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, mean)
#   }) %>%
#   do.call(cbind, . )
# 
# subject_data_sd <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, sd)
#   }) %>%
#   do.call(cbind, . )
# 
# subject_data_sem <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, function(y){
#       sd(y)/(nrow(x) - 1)
#     })
#   }) %>%
#   do.call(cbind, . )
# 
# save(subject_data_mean, file = "subject_data_mean")
# save(subject_data_sd, file = "subject_data_sd")
# save(subject_data_sem, file = "subject_data_sem")

load("subject_data_mean")
load("subject_data_sd")
load("subject_data_sem")

subject_data2 <- 
  lapply(subject_data2, function(x){
    x <-
      x %>% 
      dplyr::select(-ga_range)
  })

##k-means analysis
#-------------------------------------------------------------------------------
## fuzzy C-means for metabolites
library(Mfuzz)
library(e1071)

# mestimate <- function(df) {
#   N <- dim(df)[[1]]
#   D <- dim(df)[[2]]
#   m.sj <- 1 + (1418 / N + 22.05) * D^(-2) + (12.33 / N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
#   return(m.sj)
# }
# 
# temp_data <- subject_data_mean
# 
# temp_data <-
# temp_data %>%
#   apply(1, function(x){
#     (x - mean(x))/sd(x)
#   }) %>%
#   t() %>%
#   as.data.frame()
# 
# m <- mestimate(temp_data)
# m
# 
# library(e1071)
# 
# # helper function for the within sum of squared error
# sumsqr <- function(x, clusters) {
#   sumsqr <- function(x) sum(scale(x, scale = FALSE)^2)
#   wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
#   return(wss)
# }
# 
# # get the wss for repeated clustering
# iterate_fcm_WSS <- function(df, m) {
#   totss <- numeric()
#   for (i in 2:20) {
#     FCMresults <- cmeans(df, centers = i, m = m)
#     totss[i] <- sum(sumsqr(df, FCMresults$cluster))
#   }
#   return(totss)
# }
# 
# wss_2to20 <- iterate_fcm_WSS(temp_data, m)
# 
# data.frame(number = 1:20, error = wss_2to20[1:20], stringsAsFactors = FALSE) %>%
#   ggplot(aes(number, error)) +
#   geom_point(shape = 16, colour = "black") +
#   geom_line(colour = "black") +
#   labs(x = "Number of clusters", y = "WSS") +
#   theme_bw()
# 
# text_colour <- c(
#   colorRampPalette(colors = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   ))(13),
#   "red"
# )
# 
# k <- 5
# fcm_results <- cmeans(temp_data,
#                       centers = k, m = m)
# 
# # get the centroids into a long dataframe:
# fcm_centroids <- fcm_results$centers
# fcm_centroids_df <- data.frame(fcm_centroids,
#                                check.names = FALSE,
#                                check.rows = FALSE)
# fcm_centroids_df$cluster <- row.names(fcm_centroids_df)
# 
# centroids_long <-
#   fcm_centroids_df %>%
#   tidyr::pivot_longer(
#     cols = -cluster,
#     names_to = "GA",
#     values_to = "value"
#   )
# 
# ga_level <- unique(centroids_long$GA)
# 
# centroids_long$GA <- factor(centroids_long$GA, levels = ga_level)
# 
# plot <-
# ggplot(
#   centroids_long,
#   aes(x = GA, y = value, group = cluster, colour = as.factor(cluster))
# ) +
#   geom_line(size = 1) +
#   # geom_point(size = 3, shape = 16) +
#   ggsci::scale_color_aaas() +
#   labs(
#     x = "",
#     title = "",
#     color = "Cluster",
#     y = "Scaled intensity"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 13),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 45,
#       vjust = 1,
#       hjust = 1, colour = text_colour
#     ),
#     axis.text.y = element_text(size = 12),
#     plot.title = element_text(size = 13, hjust = 0.5)
#   )
# 
# plot
# 
# ggsave(plot, filename = "main_plot.pdf", width = 7, height = 7)
# 
# # start with the input data
# fcm_plotting_df <- data.frame(temp_data,
#                               check.names = FALSE,
#                               check.rows = FALSE)
# 
# # add genes
# fcm_plotting_df$metabolite <- row.names(fcm_plotting_df)
# 
# # bind cluster assinment
# fcm_plotting_df$cluster <- fcm_results$cluster
# 
# # fetch the membership for each gene/top scoring cluster
# fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster), function(row) {
#   clust <- fcm_plotting_df$cluster[row]
#   fcm_results$membership[row, clust]
# })
# 
# 
# k_to_plot <- 5
# 
# # subset the dataframe by the cluster and get it into long form
# # using a little tidyr action
# cluster_plot_df <-
#   dplyr::filter(fcm_plotting_df, cluster == k_to_plot) %>%
#   dplyr::select(everything(), membership, metabolite) %>%
#   tidyr::pivot_longer(
#     cols = -c(metabolite, cluster, membership),
#     names_to = "GA", values_to = "value"
#   )
# 
# cluster_plot_df
# 
# # order the data frame by score
# cluster_plot_df <- cluster_plot_df[order(cluster_plot_df$membership), ]
# 
# # set the order by setting the factors using forcats
# cluster_plot_df$metabolite <- forcats::fct_inorder(cluster_plot_df$metabolite)
# 
# cluster_plot_df$GA <- factor(cluster_plot_df$GA, levels = ga_level)
# 
# #subset the cores by cluster
# core <- dplyr::filter(centroids_long, cluster == k_to_plot)
# 
# plot <-
# cluster_plot_df %>%
#   dplyr::filter(membership > 0.8) %>%
# ggplot(aes(x = GA, y = value)) +
#   geom_line(aes(group = metabolite),
#             size = 0.5,
#             color = "blue",
#             alpha = 0.8) +
#   # geom_point(aes(colour = membership, group = metabolite), size = 2) +
#   # scale_colour_gradientn(colours = c(alpha("red", 0.1),
#   #                                    alpha("red", 0.5),
#   #                                    "red")) +
#   # this adds the core
#   geom_line(data = core, aes(GA, value, group = cluster),
#             color = "black", inherit.aes = FALSE) +
#   labs(
#     x = "", y = "Scaled intensity",
#     title = ""
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = c(1,1),
#     legend.justification = c(1, 1),
#     axis.title = element_text(size = 13),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 45,
#       vjust = 1,
#       hjust = 1,
#       colour = text_colour
#     ),
#     axis.text.y = element_text(size = 12),
#     plot.title = element_text(size = 13, hjust = 0.5)
#   )
# 
# plot
# 
# ggsave(plot, 
#        filename = paste("cluster", k_to_plot, "_plot.pdf", sep = ""),
#        width = 7, height = 7)
# 
# 
# # cluster1_metabolite <-
# #   fcm_plotting_df$metabolite[which(fcm_plotting_df$cluster == 1)]
# # 
# # cluster2_metabolite <-
# #   fcm_plotting_df$metabolite[which(fcm_plotting_df$cluster == 2)]
# 
# 
# write.csv(fcm_plotting_df, "fuzzy_c_means_cluster.csv", row.names = FALSE)


fcm_plotting_df <- readr::read_csv("fuzzy_c_means_cluster.csv")

cluster1_metabolite <- 
fcm_plotting_df$metabolite[fcm_plotting_df$cluster == 1 & 
                             fcm_plotting_df$membership > 0.8]

cluster2_metabolite <- 
  fcm_plotting_df$metabolite[fcm_plotting_df$cluster == 2 & 
                               fcm_plotting_df$membership > 0.8]

cluster3_metabolite <- 
  fcm_plotting_df$metabolite[fcm_plotting_df$cluster == 3 & 
                               fcm_plotting_df$membership > 0.8]

cluster4_metabolite <- 
  fcm_plotting_df$metabolite[fcm_plotting_df$cluster == 4 & 
                               fcm_plotting_df$membership > 0.8]

cluster5_metabolite <- 
  fcm_plotting_df$metabolite[fcm_plotting_df$cluster == 5 & 
                               fcm_plotting_df$membership > 0.8]

Reduce(intersect, list(cluster1_metabolite, 
                       cluster2_metabolite,
                       cluster3_metabolite,
                       cluster4_metabolite,
                       cluster5_metabolite))

length(cluster1_metabolite)
length(cluster2_metabolite)
length(cluster3_metabolite)
length(cluster4_metabolite)
length(cluster5_metabolite)




###correlation network for different clusters
##cluster 1
dir.create("cluster1")

temp_data <- 
  subject_data_mean[cluster1_metabolite,] %>% 
  as.data.frame()

library(corrr)
cor <- corrr::correlate(x = t(temp_data),
                              method = "spearman")
cor <-
  cor %>%
  shave() %>%
  stretch() %>%
  dplyr::filter(!is.na(r))

p <-
  as.data.frame(t(cor)) %>%
  purrr::map(.f = function(x){
    cor.test(as.numeric(temp_data[x[1],]),
             as.numeric(temp_data[x[2],]), method = "spearman"
             )$p.value
  }) %>%
  unlist()

fdr <- p.adjust(p, method = "fdr")
fdr[fdr == 0] <- min(fdr[fdr != 0])

cor <- data.frame(cor, p, fdr, stringsAsFactors = FALSE) %>%
  dplyr::filter(fdr < 0.05)

save(cor, file = "cluster1/cor")

load("cluster1/cor")

library(igraph)
library(ggraph)
library(tidygraph)

sum(abs(cor$r) > 0.5 & cor$p < 0.01)

cor <- cor %>% 
  dplyr::filter(abs(r) > 0.8 & fdr < 0.01)

edge_data <- 
  cor %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE) %>% 
  dplyr::left_join(variable_info, by = c("node" = "name"))

node_data$super_class[is.na(node_data$super_class)] <- "Unknown"

node_data <- 
  node_data %>% 
  dplyr::arrange(super_class)

node_data$super_class <- 
  factor(node_data$super_class, levels = sort(unique(node_data$super_class)))

cor_graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

cor_graph <- tidygraph::as.igraph(x = cor_graph)
save(cor_graph, file = "cluster1/cor_graph")

  # degree1 <- igraph::degree(cor_graph, mode = "all", normalized = FALSE)
  # degree2 <- igraph::degree(cor_graph, mode = "all", normalized = TRUE)
  # 
  # betweenness1 <- igraph::betweenness(graph = cor_graph1, 
  #                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
  #                                     normalized = FALSE)
  # betweenness2 <- igraph::betweenness(graph = cor_graph1, 
  #                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
  #                                     normalized = TRUE)
  # closeness1 <- 
  #   igraph::closeness(graph = cor_graph1, mode = "all",
  #                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
  #                     normalized = FALSE)
  # closeness2 <- 
  #   igraph::closeness(graph = cor_graph1, mode = "all",
  #                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
  #                     normalized = TRUE)
  # 
  # importance <- 
  #   data.frame((degree2 - mean(degree2))/sd(degree2),
  #              (closeness2 - mean(closeness2))/sd(closeness2),
  #              (betweenness2 - mean(betweenness2))/sd(betweenness2),
  #              stringsAsFactors = FALSE
  #   ) %>% 
  #   apply(1, mean)
  # 
  # node_name <- igraph::vertex.attributes(cor_graph1)$node
  # 
  # node_info <-
  #   data.frame(
  #     node_name,
  #     degree1,
  #     degree2,
  #     betweenness1,
  #     betweenness2,
  #     closeness1,
  #     closeness2,
  #     importance,
  #     stringsAsFactors = FALSE
  #   ) %>% 
  #   dplyr::arrange(importance)
  # 
  # hub_gene1 <- node_info %>% 
  #   dplyr::filter(importance > quantile(importance, 0.75)) %>% 
  #   dplyr::pull(node_name)
  # 
  # hub_gene2 <- node_info %>% 
  #   dplyr::filter(degree1 > quantile(degree1, 0.75)) %>% 
  #   dplyr::pull(node_name)
  # 
  # hub_gene3 <- node_info %>% 
  #   dplyr::filter(betweenness1 > quantile(betweenness1, 0.75)) %>% 
  #   dplyr::pull(node_name)
  # 
  # hub_gene4 <- node_info %>% 
  #   dplyr::filter(closeness1 > quantile(closeness1, 0.75)) %>% 
  #   dplyr::pull(node_name)
  # 
  # hub_gene <- 
  #   reduce(list(hub_gene1 = hub_gene1, 
  #               hub_gene2 = hub_gene2,
  #               hub_gene3 = hub_gene3,
  #               hub_gene4 = hub_gene4), union)
  # 
  # 
  # temp_data <-
  #   node_info %>% 
  #   dplyr::filter(node_name %in% hub_gene) %>% 
  #   dplyr::select(node_name, importance, contains("1")) %>% 
  #   dplyr::arrange(importance)
  # 
  # plot <-
  #   temp_data %>% 
  #   dplyr::select(node_name, importance, contains("1")) %>% 
  #   dplyr::mutate(degree1 = (degree1 - mean(degree1))/sd(degree1)) %>% 
  #   dplyr::mutate(betweenness1 = (betweenness1 - mean(betweenness1))/sd(betweenness1)) %>% 
  #   dplyr::mutate(closeness1 = (closeness1 - mean(closeness1))/sd(closeness1)) %>% 
  #   tidyr::pivot_longer(-node_name, names_to = "class", values_to = "value") %>% 
  #   dplyr::mutate(node_name = factor(node_name, levels = node_info$node_name)) %>% 
  #   dplyr::mutate(class = 
  #                   case_when(
  #                     class == "importance" ~ "Importance",
  #                     class == "degree1" ~ "Degree",
  #                     class == "betweenness1" ~ "Betweenness",
  #                     class == "closeness1" ~ "Closeness"
  #                   )) %>% 
  #   dplyr::mutate(
  #     class = factor(class, levels = c(
  #       "Importance", "Degree", "Betweenness", "Closeness"
  #     ))
  #   ) %>% 
  #   ggplot(aes(value, node_name)) +
  #   geom_point(aes(color = class), shape = 16, show.legend = FALSE) +
  #   geom_segment(aes(x = -1, 
  #                    y = node_name, xend = value, yend = node_name,
  #                    color = class), show.legend = FALSE) +
  #   scale_x_continuous(expand = expansion(mult = c(-0,0.1))) +
  #   ggsci::scale_color_futurama() +
  #   labs(x = "", y = "") +
  #   facet_wrap(facets = "class", nrow = 1, scales = "free_x") +
  #   theme_bw() +
  #   theme()
  # 
  # ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.pdf", sep = ""),
  #        width = 8, height = 7)
  # ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.png", sep = ""),
  #        width = 8, height = 7)
  


###pathway enrichment
hmdb_id <- node_data$HMDB.ID
hmdb_id <- unique(hmdb_id[!is.na(hmdb_id)])

kegg_id <- 
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database", 
                      to = "KEGG", top = 1)
  })

kegg_id <- 
  kegg_id %>% 
  do.call(rbind, .)

kegg_id <- kegg_id$KEGG 
kegg_id <- kegg_id[!is.na(kegg_id)]

load("hsa_pathway")

path_result <- 
  enrichPathway(id = kegg_id, database = hsa_pathway)
  
path_result %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)


save(path_result, file = "cluster1/path_result")

colour <<- 
  c(
    "Clinical information" = "#FF6F00FF",
    "Alkaloids and derivatives" = "#ADE2D0FF",
    "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    "Organic acids and derivatives" = "#8A4198FF",
    "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

plot <-
  ggraph(cor_graph,
         layout = 'kk'
         # circular = TRUE
         ) +
         geom_edge_link(aes(edge_width = fdr,
                            color = cor),
                        alpha = 0.5,
                        show.legend = TRUE) +
           geom_node_point(aes(size = Degree,
                               fill = super_class),
                           alpha = 0.6,
                           shape = 21,
                           show.legend = TRUE) +
           w# geom_node_text(aes(label = node), repel = TRUE) +
           guides(fill = guide_legend(override.aes = list(size = 3))) +
           ggraph::scale_edge_color_gradient2(
             low = ggsci::pal_d3()(10)[1],
             mid = "white",
             high = ggsci::pal_d3()(10)[5],
             midpoint = 0.6
           ) +
           ggraph::scale_edge_width(range = c(0.05, 0.8)) +
           scale_size_continuous(range = c(1, 5)) +
           scale_fill_manual(values = colour) +
           ggraph::theme_graph() +
           theme(
             plot.background = element_rect(fill = "transparent", color = NA),
             panel.background = element_rect(fill = "transparent", color = NA),
             legend.position = "right",
             legend.background = element_rect(fill = "transparent", color = NA)
           )
         plot
# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "cluster1/cor_network.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)












##cluster 2
dir.create("cluster2")

temp_data <- 
  subject_data_mean[cluster2_metabolite,] %>% 
  as.data.frame()

library(corrr)
cor <- corrr::correlate(x = t(temp_data),
                        method = "spearman")
cor <-
  cor %>%
  shave() %>%
  stretch() %>%
  dplyr::filter(!is.na(r))

p <-
  as.data.frame(t(cor)) %>%
  purrr::map(.f = function(x){
    cor.test(as.numeric(temp_data[x[1],]),
             as.numeric(temp_data[x[2],]), method = "spearman"
    )$p.value
  }) %>%
  unlist()

fdr <- p.adjust(p, method = "fdr")
fdr[fdr == 0] <- min(fdr[fdr != 0])

cor <- data.frame(cor, p, fdr, stringsAsFactors = FALSE) %>%
  dplyr::filter(fdr < 0.05)

save(cor, file = "cluster2/cor")

load("cluster2/cor")

library(igraph)
library(ggraph)
library(tidygraph)

sum(abs(cor$r) > 0.5 & cor$p < 0.01)

cor <- cor %>% 
  dplyr::filter(abs(r) > 0.8 & fdr < 0.01)

edge_data <- 
  cor %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE) %>% 
  dplyr::left_join(variable_info, by = c("node" = "name"))

node_data$super_class[is.na(node_data$super_class)] <- "Unknown"

node_data <- 
  node_data %>% 
  dplyr::arrange(super_class)

node_data$super_class <- 
  factor(node_data$super_class, levels = sort(unique(node_data$super_class)))

cor_graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

cor_graph <- tidygraph::as.igraph(x = cor_graph)
save(cor_graph, file = "cluster2/cor_graph")

# degree1 <- igraph::degree(cor_graph, mode = "all", normalized = FALSE)
# degree2 <- igraph::degree(cor_graph, mode = "all", normalized = TRUE)
# 
# betweenness1 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = FALSE)
# betweenness2 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = TRUE)
# closeness1 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = FALSE)
# closeness2 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = TRUE)
# 
# importance <- 
#   data.frame((degree2 - mean(degree2))/sd(degree2),
#              (closeness2 - mean(closeness2))/sd(closeness2),
#              (betweenness2 - mean(betweenness2))/sd(betweenness2),
#              stringsAsFactors = FALSE
#   ) %>% 
#   apply(1, mean)
# 
# node_name <- igraph::vertex.attributes(cor_graph1)$node
# 
# node_info <-
#   data.frame(
#     node_name,
#     degree1,
#     degree2,
#     betweenness1,
#     betweenness2,
#     closeness1,
#     closeness2,
#     importance,
#     stringsAsFactors = FALSE
#   ) %>% 
#   dplyr::arrange(importance)
# 
# hub_gene1 <- node_info %>% 
#   dplyr::filter(importance > quantile(importance, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene2 <- node_info %>% 
#   dplyr::filter(degree1 > quantile(degree1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene3 <- node_info %>% 
#   dplyr::filter(betweenness1 > quantile(betweenness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene4 <- node_info %>% 
#   dplyr::filter(closeness1 > quantile(closeness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene <- 
#   reduce(list(hub_gene1 = hub_gene1, 
#               hub_gene2 = hub_gene2,
#               hub_gene3 = hub_gene3,
#               hub_gene4 = hub_gene4), union)
# 
# 
# temp_data <-
#   node_info %>% 
#   dplyr::filter(node_name %in% hub_gene) %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::arrange(importance)
# 
# plot <-
#   temp_data %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::mutate(degree1 = (degree1 - mean(degree1))/sd(degree1)) %>% 
#   dplyr::mutate(betweenness1 = (betweenness1 - mean(betweenness1))/sd(betweenness1)) %>% 
#   dplyr::mutate(closeness1 = (closeness1 - mean(closeness1))/sd(closeness1)) %>% 
#   tidyr::pivot_longer(-node_name, names_to = "class", values_to = "value") %>% 
#   dplyr::mutate(node_name = factor(node_name, levels = node_info$node_name)) %>% 
#   dplyr::mutate(class = 
#                   case_when(
#                     class == "importance" ~ "Importance",
#                     class == "degree1" ~ "Degree",
#                     class == "betweenness1" ~ "Betweenness",
#                     class == "closeness1" ~ "Closeness"
#                   )) %>% 
#   dplyr::mutate(
#     class = factor(class, levels = c(
#       "Importance", "Degree", "Betweenness", "Closeness"
#     ))
#   ) %>% 
#   ggplot(aes(value, node_name)) +
#   geom_point(aes(color = class), shape = 16, show.legend = FALSE) +
#   geom_segment(aes(x = -1, 
#                    y = node_name, xend = value, yend = node_name,
#                    color = class), show.legend = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(-0,0.1))) +
#   ggsci::scale_color_futurama() +
#   labs(x = "", y = "") +
#   facet_wrap(facets = "class", nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme()
# 
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.pdf", sep = ""),
#        width = 8, height = 7)
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.png", sep = ""),
#        width = 8, height = 7)



###pathway enrichment
hmdb_id <- node_data$HMDB.ID
hmdb_id <- unique(hmdb_id[!is.na(hmdb_id)])

kegg_id <- 
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database", 
                      to = "KEGG", top = 1)
  })

kegg_id <- 
  kegg_id %>% 
  do.call(rbind, .)

kegg_id <- kegg_id$KEGG 
kegg_id <- kegg_id[!is.na(kegg_id)]

load("hsa_pathway")

path_result <- 
  enrichPathway(id = kegg_id, database = hsa_pathway)

path_result %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)


save(path_result, file = "cluster2/path_result")

colour <<- 
  c(
    "Clinical information" = "#FF6F00FF",
    "Alkaloids and derivatives" = "#ADE2D0FF",
    "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    "Organic acids and derivatives" = "#8A4198FF",
    "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

plot <-
  ggraph(cor_graph,
         layout = 'kk'
         # circular = TRUE
  ) +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      fill = super_class),
                  alpha = 0.6,
                  shape = 21,
                  show.legend = TRUE) +
  # geom_node_text(aes(label = node), repel = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = ggsci::pal_d3()(10)[1],
    mid = "white",
    high = ggsci::pal_d3()(10)[5],
    midpoint = 0.6
  ) +
  ggraph::scale_edge_width(range = c(0.05, 0.8)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_manual(values = colour) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "cluster2/cor_network.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)









##cluster 3
dir.create("cluster3")

temp_data <- 
  subject_data_mean[cluster3_metabolite,] %>% 
  as.data.frame()

library(corrr)
cor <- corrr::correlate(x = t(temp_data),
                        method = "spearman")
cor <-
  cor %>%
  shave() %>%
  stretch() %>%
  dplyr::filter(!is.na(r))

p <-
  as.data.frame(t(cor)) %>%
  purrr::map(.f = function(x){
    cor.test(as.numeric(temp_data[x[1],]),
             as.numeric(temp_data[x[2],]), method = "spearman"
    )$p.value
  }) %>%
  unlist()

fdr <- p.adjust(p, method = "fdr")
fdr[fdr == 0] <- min(fdr[fdr != 0])

cor <- data.frame(cor, p, fdr, stringsAsFactors = FALSE) %>%
  dplyr::filter(fdr < 0.05)

save(cor, file = "cluster3/cor")

load("cluster3/cor")

library(igraph)
library(ggraph)
library(tidygraph)

sum(abs(cor$r) > 0.5 & cor$p < 0.01)

cor <- cor %>% 
  dplyr::filter(abs(r) > 0.8 & fdr < 0.01)

edge_data <- 
  cor %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE) %>% 
  dplyr::left_join(variable_info, by = c("node" = "name"))

node_data$super_class[is.na(node_data$super_class)] <- "Unknown"

node_data <- 
  node_data %>% 
  dplyr::arrange(super_class)

node_data$super_class <- 
  factor(node_data$super_class, levels = sort(unique(node_data$super_class)))

cor_graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

cor_graph <- tidygraph::as.igraph(x = cor_graph)
save(cor_graph, file = "cluster3/cor_graph")

# degree1 <- igraph::degree(cor_graph, mode = "all", normalized = FALSE)
# degree2 <- igraph::degree(cor_graph, mode = "all", normalized = TRUE)
# 
# betweenness1 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = FALSE)
# betweenness2 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = TRUE)
# closeness1 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = FALSE)
# closeness2 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = TRUE)
# 
# importance <- 
#   data.frame((degree2 - mean(degree2))/sd(degree2),
#              (closeness2 - mean(closeness2))/sd(closeness2),
#              (betweenness2 - mean(betweenness2))/sd(betweenness2),
#              stringsAsFactors = FALSE
#   ) %>% 
#   apply(1, mean)
# 
# node_name <- igraph::vertex.attributes(cor_graph1)$node
# 
# node_info <-
#   data.frame(
#     node_name,
#     degree1,
#     degree2,
#     betweenness1,
#     betweenness2,
#     closeness1,
#     closeness2,
#     importance,
#     stringsAsFactors = FALSE
#   ) %>% 
#   dplyr::arrange(importance)
# 
# hub_gene1 <- node_info %>% 
#   dplyr::filter(importance > quantile(importance, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene2 <- node_info %>% 
#   dplyr::filter(degree1 > quantile(degree1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene3 <- node_info %>% 
#   dplyr::filter(betweenness1 > quantile(betweenness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene4 <- node_info %>% 
#   dplyr::filter(closeness1 > quantile(closeness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene <- 
#   reduce(list(hub_gene1 = hub_gene1, 
#               hub_gene2 = hub_gene2,
#               hub_gene3 = hub_gene3,
#               hub_gene4 = hub_gene4), union)
# 
# 
# temp_data <-
#   node_info %>% 
#   dplyr::filter(node_name %in% hub_gene) %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::arrange(importance)
# 
# plot <-
#   temp_data %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::mutate(degree1 = (degree1 - mean(degree1))/sd(degree1)) %>% 
#   dplyr::mutate(betweenness1 = (betweenness1 - mean(betweenness1))/sd(betweenness1)) %>% 
#   dplyr::mutate(closeness1 = (closeness1 - mean(closeness1))/sd(closeness1)) %>% 
#   tidyr::pivot_longer(-node_name, names_to = "class", values_to = "value") %>% 
#   dplyr::mutate(node_name = factor(node_name, levels = node_info$node_name)) %>% 
#   dplyr::mutate(class = 
#                   case_when(
#                     class == "importance" ~ "Importance",
#                     class == "degree1" ~ "Degree",
#                     class == "betweenness1" ~ "Betweenness",
#                     class == "closeness1" ~ "Closeness"
#                   )) %>% 
#   dplyr::mutate(
#     class = factor(class, levels = c(
#       "Importance", "Degree", "Betweenness", "Closeness"
#     ))
#   ) %>% 
#   ggplot(aes(value, node_name)) +
#   geom_point(aes(color = class), shape = 16, show.legend = FALSE) +
#   geom_segment(aes(x = -1, 
#                    y = node_name, xend = value, yend = node_name,
#                    color = class), show.legend = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(-0,0.1))) +
#   ggsci::scale_color_futurama() +
#   labs(x = "", y = "") +
#   facet_wrap(facets = "class", nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme()
# 
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.pdf", sep = ""),
#        width = 8, height = 7)
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.png", sep = ""),
#        width = 8, height = 7)



###pathway enrichment
hmdb_id <- node_data$HMDB.ID
hmdb_id <- unique(hmdb_id[!is.na(hmdb_id)])

kegg_id <- 
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database", 
                      to = "KEGG", top = 1)
  })

kegg_id <- 
  kegg_id %>% 
  do.call(rbind, .)

kegg_id <- kegg_id$KEGG 
kegg_id <- kegg_id[!is.na(kegg_id)]

load("hsa_pathway")

path_result <- 
  enrichPathway(id = kegg_id, database = hsa_pathway)

path_result %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)


save(path_result, file = "cluster3/path_result")

colour <<- 
  c(
    "Clinical information" = "#FF6F00FF",
    "Alkaloids and derivatives" = "#ADE2D0FF",
    "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    "Organic acids and derivatives" = "#8A4198FF",
    "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

plot <-
  ggraph(cor_graph,
         layout = 'kk'
         # circular = TRUE
  ) +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      fill = super_class),
                  alpha = 0.6,
                  shape = 21,
                  show.legend = TRUE) +
  # geom_node_text(aes(label = node), repel = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = ggsci::pal_d3()(10)[1],
    mid = "white",
    high = ggsci::pal_d3()(10)[5],
    midpoint = 0.6
  ) +
  ggraph::scale_edge_width(range = c(0.05, 0.8)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_manual(values = colour) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "cluster3/cor_network.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)








##cluster 4
dir.create("cluster4")

temp_data <- 
  subject_data_mean[cluster4_metabolite,] %>% 
  as.data.frame()

library(corrr)
cor <- corrr::correlate(x = t(temp_data),
                        method = "spearman")
cor <-
  cor %>%
  shave() %>%
  stretch() %>%
  dplyr::filter(!is.na(r))

p <-
  as.data.frame(t(cor)) %>%
  purrr::map(.f = function(x){
    cor.test(as.numeric(temp_data[x[1],]),
             as.numeric(temp_data[x[2],]), method = "spearman"
    )$p.value
  }) %>%
  unlist()

fdr <- p.adjust(p, method = "fdr")
fdr[fdr == 0] <- min(fdr[fdr != 0])

cor <- data.frame(cor, p, fdr, stringsAsFactors = FALSE) %>%
  dplyr::filter(fdr < 0.05)

save(cor, file = "cluster4/cor")

load("cluster4/cor")

library(igraph)
library(ggraph)
library(tidygraph)

sum(abs(cor$r) > 0.5 & cor$p < 0.01)

cor <- cor %>% 
  dplyr::filter(abs(r) > 0.5 & fdr < 0.05)

edge_data <- 
  cor %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE) %>% 
  dplyr::left_join(variable_info, by = c("node" = "name"))

node_data$super_class[is.na(node_data$super_class)] <- "Unknown"

node_data <- 
  node_data %>% 
  dplyr::arrange(super_class)

node_data$super_class <- 
  factor(node_data$super_class, levels = sort(unique(node_data$super_class)))

cor_graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

cor_graph <- tidygraph::as.igraph(x = cor_graph)
save(cor_graph, file = "cluster4/cor_graph")

# degree1 <- igraph::degree(cor_graph, mode = "all", normalized = FALSE)
# degree2 <- igraph::degree(cor_graph, mode = "all", normalized = TRUE)
# 
# betweenness1 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = FALSE)
# betweenness2 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = TRUE)
# closeness1 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = FALSE)
# closeness2 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = TRUE)
# 
# importance <- 
#   data.frame((degree2 - mean(degree2))/sd(degree2),
#              (closeness2 - mean(closeness2))/sd(closeness2),
#              (betweenness2 - mean(betweenness2))/sd(betweenness2),
#              stringsAsFactors = FALSE
#   ) %>% 
#   apply(1, mean)
# 
# node_name <- igraph::vertex.attributes(cor_graph1)$node
# 
# node_info <-
#   data.frame(
#     node_name,
#     degree1,
#     degree2,
#     betweenness1,
#     betweenness2,
#     closeness1,
#     closeness2,
#     importance,
#     stringsAsFactors = FALSE
#   ) %>% 
#   dplyr::arrange(importance)
# 
# hub_gene1 <- node_info %>% 
#   dplyr::filter(importance > quantile(importance, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene2 <- node_info %>% 
#   dplyr::filter(degree1 > quantile(degree1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene3 <- node_info %>% 
#   dplyr::filter(betweenness1 > quantile(betweenness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene4 <- node_info %>% 
#   dplyr::filter(closeness1 > quantile(closeness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene <- 
#   reduce(list(hub_gene1 = hub_gene1, 
#               hub_gene2 = hub_gene2,
#               hub_gene3 = hub_gene3,
#               hub_gene4 = hub_gene4), union)
# 
# 
# temp_data <-
#   node_info %>% 
#   dplyr::filter(node_name %in% hub_gene) %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::arrange(importance)
# 
# plot <-
#   temp_data %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::mutate(degree1 = (degree1 - mean(degree1))/sd(degree1)) %>% 
#   dplyr::mutate(betweenness1 = (betweenness1 - mean(betweenness1))/sd(betweenness1)) %>% 
#   dplyr::mutate(closeness1 = (closeness1 - mean(closeness1))/sd(closeness1)) %>% 
#   tidyr::pivot_longer(-node_name, names_to = "class", values_to = "value") %>% 
#   dplyr::mutate(node_name = factor(node_name, levels = node_info$node_name)) %>% 
#   dplyr::mutate(class = 
#                   case_when(
#                     class == "importance" ~ "Importance",
#                     class == "degree1" ~ "Degree",
#                     class == "betweenness1" ~ "Betweenness",
#                     class == "closeness1" ~ "Closeness"
#                   )) %>% 
#   dplyr::mutate(
#     class = factor(class, levels = c(
#       "Importance", "Degree", "Betweenness", "Closeness"
#     ))
#   ) %>% 
#   ggplot(aes(value, node_name)) +
#   geom_point(aes(color = class), shape = 16, show.legend = FALSE) +
#   geom_segment(aes(x = -1, 
#                    y = node_name, xend = value, yend = node_name,
#                    color = class), show.legend = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(-0,0.1))) +
#   ggsci::scale_color_futurama() +
#   labs(x = "", y = "") +
#   facet_wrap(facets = "class", nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme()
# 
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.pdf", sep = ""),
#        width = 8, height = 7)
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.png", sep = ""),
#        width = 8, height = 7)



###pathway enrichment
hmdb_id <- node_data$HMDB.ID
hmdb_id <- unique(hmdb_id[!is.na(hmdb_id)])

kegg_id <- 
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database", 
                      to = "KEGG", top = 1)
  })

kegg_id <- 
  kegg_id %>% 
  do.call(rbind, .)

kegg_id <- kegg_id$KEGG 
kegg_id <- kegg_id[!is.na(kegg_id)]

load("hsa_pathway")

path_result <- 
  enrichPathway(id = kegg_id, database = hsa_pathway)

path_result %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_result, file = "cluster4/path_result")

colour <<- 
  c(
    "Clinical information" = "#FF6F00FF",
    "Alkaloids and derivatives" = "#ADE2D0FF",
    "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    "Organic acids and derivatives" = "#8A4198FF",
    "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

plot <-
  ggraph(cor_graph,
         layout = 'kk'
         # circular = TRUE
  ) +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      fill = super_class),
                  alpha = 0.6,
                  shape = 21,
                  show.legend = TRUE) +
  # geom_node_text(aes(label = node), repel = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = ggsci::pal_d3()(10)[1],
    mid = "white",
    high = ggsci::pal_d3()(10)[5],
    midpoint = 0.6
  ) +
  ggraph::scale_edge_width(range = c(0.05, 0.8)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_manual(values = colour) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "cluster4/cor_network.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)





##cluster 5
dir.create("cluster4")

temp_data <- 
  subject_data_mean[cluster4_metabolite,] %>% 
  as.data.frame()

library(corrr)
cor <- corrr::correlate(x = t(temp_data),
                        method = "spearman")
cor <-
  cor %>%
  shave() %>%
  stretch() %>%
  dplyr::filter(!is.na(r))

p <-
  as.data.frame(t(cor)) %>%
  purrr::map(.f = function(x){
    cor.test(as.numeric(temp_data[x[1],]),
             as.numeric(temp_data[x[2],]), method = "spearman"
    )$p.value
  }) %>%
  unlist()

fdr <- p.adjust(p, method = "fdr")
fdr[fdr == 0] <- min(fdr[fdr != 0])

cor <- data.frame(cor, p, fdr, stringsAsFactors = FALSE) %>%
  dplyr::filter(fdr < 0.05)

save(cor, file = "cluster4/cor")

load("cluster4/cor")

library(igraph)
library(ggraph)
library(tidygraph)

sum(abs(cor$r) > 0.5 & cor$p < 0.01)

cor <- cor %>% 
  dplyr::filter(abs(r) > 0.5 & fdr < 0.05)

edge_data <- 
  cor %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE) %>% 
  dplyr::left_join(variable_info, by = c("node" = "name"))

node_data$super_class[is.na(node_data$super_class)] <- "Unknown"

node_data <- 
  node_data %>% 
  dplyr::arrange(super_class)

node_data$super_class <- 
  factor(node_data$super_class, levels = sort(unique(node_data$super_class)))

cor_graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

cor_graph <- tidygraph::as.igraph(x = cor_graph)
save(cor_graph, file = "cluster4/cor_graph")

# degree1 <- igraph::degree(cor_graph, mode = "all", normalized = FALSE)
# degree2 <- igraph::degree(cor_graph, mode = "all", normalized = TRUE)
# 
# betweenness1 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = FALSE)
# betweenness2 <- igraph::betweenness(graph = cor_graph1, 
#                                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                                     normalized = TRUE)
# closeness1 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = FALSE)
# closeness2 <- 
#   igraph::closeness(graph = cor_graph1, mode = "all",
#                     weights = abs(igraph::edge_attr(graph = cor_graph1, name = "cor")), 
#                     normalized = TRUE)
# 
# importance <- 
#   data.frame((degree2 - mean(degree2))/sd(degree2),
#              (closeness2 - mean(closeness2))/sd(closeness2),
#              (betweenness2 - mean(betweenness2))/sd(betweenness2),
#              stringsAsFactors = FALSE
#   ) %>% 
#   apply(1, mean)
# 
# node_name <- igraph::vertex.attributes(cor_graph1)$node
# 
# node_info <-
#   data.frame(
#     node_name,
#     degree1,
#     degree2,
#     betweenness1,
#     betweenness2,
#     closeness1,
#     closeness2,
#     importance,
#     stringsAsFactors = FALSE
#   ) %>% 
#   dplyr::arrange(importance)
# 
# hub_gene1 <- node_info %>% 
#   dplyr::filter(importance > quantile(importance, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene2 <- node_info %>% 
#   dplyr::filter(degree1 > quantile(degree1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene3 <- node_info %>% 
#   dplyr::filter(betweenness1 > quantile(betweenness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene4 <- node_info %>% 
#   dplyr::filter(closeness1 > quantile(closeness1, 0.75)) %>% 
#   dplyr::pull(node_name)
# 
# hub_gene <- 
#   reduce(list(hub_gene1 = hub_gene1, 
#               hub_gene2 = hub_gene2,
#               hub_gene3 = hub_gene3,
#               hub_gene4 = hub_gene4), union)
# 
# 
# temp_data <-
#   node_info %>% 
#   dplyr::filter(node_name %in% hub_gene) %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::arrange(importance)
# 
# plot <-
#   temp_data %>% 
#   dplyr::select(node_name, importance, contains("1")) %>% 
#   dplyr::mutate(degree1 = (degree1 - mean(degree1))/sd(degree1)) %>% 
#   dplyr::mutate(betweenness1 = (betweenness1 - mean(betweenness1))/sd(betweenness1)) %>% 
#   dplyr::mutate(closeness1 = (closeness1 - mean(closeness1))/sd(closeness1)) %>% 
#   tidyr::pivot_longer(-node_name, names_to = "class", values_to = "value") %>% 
#   dplyr::mutate(node_name = factor(node_name, levels = node_info$node_name)) %>% 
#   dplyr::mutate(class = 
#                   case_when(
#                     class == "importance" ~ "Importance",
#                     class == "degree1" ~ "Degree",
#                     class == "betweenness1" ~ "Betweenness",
#                     class == "closeness1" ~ "Closeness"
#                   )) %>% 
#   dplyr::mutate(
#     class = factor(class, levels = c(
#       "Importance", "Degree", "Betweenness", "Closeness"
#     ))
#   ) %>% 
#   ggplot(aes(value, node_name)) +
#   geom_point(aes(color = class), shape = 16, show.legend = FALSE) +
#   geom_segment(aes(x = -1, 
#                    y = node_name, xend = value, yend = node_name,
#                    color = class), show.legend = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(-0,0.1))) +
#   ggsci::scale_color_futurama() +
#   labs(x = "", y = "") +
#   facet_wrap(facets = "class", nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme()
# 
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.pdf", sep = ""),
#        width = 8, height = 7)
# ggsave(plot, file = paste("trans_subnetwork",idx,"_hub_genes.png", sep = ""),
#        width = 8, height = 7)



###pathway enrichment
hmdb_id <- node_data$HMDB.ID
hmdb_id <- unique(hmdb_id[!is.na(hmdb_id)])

kegg_id <- 
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database", 
                      to = "KEGG", top = 1)
  })

kegg_id <- 
  kegg_id %>% 
  do.call(rbind, .)

kegg_id <- kegg_id$KEGG 
kegg_id <- kegg_id[!is.na(kegg_id)]

load("hsa_pathway")

path_result <- 
  enrichPathway(id = kegg_id, database = hsa_pathway)

path_result %>% 
  dplyr::filter(Overlap >= 3 & p.value.fdr < 0.05)

save(path_result, file = "cluster4/path_result")

colour <<- 
  c(
    "Clinical information" = "#FF6F00FF",
    "Alkaloids and derivatives" = "#ADE2D0FF",
    "Benzenoids" = "#C71000FF",
    "Lipids and lipid-like molecules" = "#FF6F00B2",
    "Nucleosides, nucleotides, and analogues" = "#C71000B2",
    "Organic acids and derivatives" = "#8A4198FF",
    "Organic nitrogen compounds" = "#5A9599FF",
    "Organic oxygen compounds" = "#008EA0B2",
    "Organoheterocyclic compounds" = "#8A4198B2",
    "Organosulfur compounds" = "#3F4041FF",
    "Phenylpropanoids and polyketides" = "#5A9599B2",
    "Unknown" = "#3F4041B2"
  )

plot <-
  ggraph(cor_graph,
         layout = 'kk'
         # circular = TRUE
  ) +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      fill = super_class),
                  alpha = 0.6,
                  shape = 21,
                  show.legend = TRUE) +
  # geom_node_text(aes(label = node), repel = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(
    low = ggsci::pal_d3()(10)[1],
    mid = "white",
    high = ggsci::pal_d3()(10)[5],
    midpoint = 0.6
  ) +
  ggraph::scale_edge_width(range = c(0.05, 0.8)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_manual(values = colour) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "cluster4/cor_network.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)

