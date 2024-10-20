##avoid source 
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
source("R/20200727/tools.R")


setwd("data_analysis20200108/urine_metabolome/DEG_analysis/marker_in_different_points/")

###28_30
load("(28,30]/piumet_output/Result/edge_data")
load("(28,30]/piumet_output/Result/node_data")

#remove peaks
node_data <- 
  node_data %>% 
  dplyr::filter(node_class != "m/z Peak")
edge_data <- 
edge_data %>% 
  dplyr::filter(from %in% node_data$node) %>% 
  dplyr::filter(to %in% node_data$node) 

graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = tidygraph::centrality_degree(mode = 'all'))


fill <- 
  c(
    "m/z Peak" = ggsci::pal_aaas()(10)[9],
    "Metabolite" = ggsci::pal_aaas()(10)[1],
    "Protein" = ggsci::pal_aaas(alpha = 0.1)(10)[2]
    # "Protein" = "tomato"
  )    

col <- 
  c(
    "m/z Peak" = ggsci::pal_aaas()(10)[9],
    "Metabolite" = ggsci::pal_aaas()(10)[1],
    "Protein" = ggsci::pal_aaas()(10)[2]
  )

shape = c(
  "m/z Peak" = 24,
  "Metabolite" = 22,
  # "Metabolite_others" = 22,
  "Protein" = 21
)

require(ggraph)

  plot <-
    ggraph(graph,
           layout = "kk") +
    geom_edge_link(aes(edge_width = V3),
                   alpha = 1,
                   color = "black",
                   show.legend = TRUE) +
    geom_node_point(aes(size = V2, 
                        fill = node_class,
                        shape = node_class),
                    alpha = 1, show.legend = TRUE) +
    scale_shape_manual(values = shape) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    # ggraph::geom_node_text(aes(label = node, 
    #                            color = node_class), 
    #                        repel = TRUE, size = 3) +
    ggraph::scale_edge_width(range = c(0.1,1)) +
    scale_size_continuous(range = c(2,5)) +
    scale_fill_manual(values = fill) +
    scale_color_manual(values = col) +
    ggraph::theme_graph() +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA))

  plot
  
