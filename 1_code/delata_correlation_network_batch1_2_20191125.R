#######construct delta correlation network
#####smartd_rplc_batch dataset is used to analysis for correlation networks analysis
#### step 1 remove the duplicated metabolites according to score
sxtTools::setwd_project()
setwd("data_analysis20191125/metabolite_table/correlation_network/delta_network/")
metabolite_table <- readr::read_csv("../cross_sectional_network/metabolite_table_cross.csv")
library(tidyverse)
###calculate correlation matrix
###delta network
####according to the GA information, we shoud class them to different peroid
metabolite_table %>% 
  ggplot(., aes(x = GA)) +
  geom_histogram(colour = "white", binwidth = 5, 
                 fill = "#8A9045CC", boundary = 10) +
  theme_bw() +
  labs(x = "GA (week)") +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15))

hist(metabolite_table$GA)
range(metabolite_table$GA)

library(plyr)

##for each patient, the  GA weeks for her samples
lapply(metabolite_table %>%
         plyr::dlply(., .(Patient_ID)),
       function(x)
         x$GA)

metabolite_table <- 
  metabolite_table %>% 
  mutate(GA_range = case_when(GA > 10 & GA <= 15 ~ {"10_15"},
                              GA > 15 & GA <= 20 ~ {"15_20"},
                              GA > 20 & GA <= 25 ~ {"20_25"},
                              GA > 25 & GA <= 30 ~ {"25_30"},
                              GA > 30 & GA <= 35 ~ {"30_35"},
                              GA > 35 & GA <= 40 ~ {"35_40"},
                              GA > 40 & GA <= 45 ~ {"40_45"}
                              ###10-15,15-20,20-25,25-30,30-35,35-40,40-45
  ))


lapply(
  plyr::dlply(metabolite_table, .variables = .(Patient_ID)),
  function(x) x$GA_range
)


lapply(
  plyr::dlply(metabolite_table, .variables = .(Patient_ID)),
  function(x) x$GA
)

####round the GA
metabolite_table_delta <- 
  metabolite_table %>% 
  mutate(GA_round = round(GA))

lapply(
  plyr::dlply(metabolite_table_delta, .variables = .(Patient_ID)),
  function(x) x$GA_round
) %>% 
  unlist() %>% 
  table() %>% 
  sort()

metabolite_table_delta <-
  metabolite_table_delta %>% 
  plyr::dlply(., .(Patient_ID)) %>% 
  lapply(., function(x){
    x <- x %>% arrange(., GA_range)
    x <- x[!duplicated(x$GA_range),,drop = FALSE]
  })

metabolite_table_delta <-
  metabolite_table_delta %>% 
  dplyr::bind_rows()

####upset to show the venn diagram
library(UpSetR)
temp_data <-
  lapply(
    plyr::dlply(metabolite_table_delta, .variables = .(Patient_ID)),
    function(x) x$GA_range
  )

total_week <-
  unique(unlist(temp_data)) %>%
  unname() %>%
  sort()

temp_data <-
  lapply(temp_data, function(x){
    sapply(total_week, function(y) {
      if(y %in% x){
        1
      } else{
        0
      }
    })
  })

temp_data <-
  temp_data %>%
  bind_cols()

# delta_temp_data <-
#   t(delta_temp_data) %>%
#   as.data.frame()
# colnames(delta_temp_data) <- total_week

rownames(temp_data) <- total_week

temp_data <-
  temp_data %>%
  rownames_to_column(., var = "GA_range")

upset(as.data.frame(temp_data),
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
      # sets = colnames(delta_temp_data)[-1],
      nintersects = NA, keep.order = FALSE,
      # matrix.color = c("grey","#155F83CC"),
      # main.bar.color = "red",
      mainbar.y.label = "Intercetion",
      sets.x.label = "Number of samples"
      # queries = list(
      #   list(
      #     query = intersects,
      #     params = list("1554", "1555", "1556", "1650", "1626"),
      #     active = TRUE,
      #     color = "#FFA319CC"
      #   )
      # )
)


#####
temp_data %>% 
  select(., c("SF1728", "SF1755", "SF1764", "SF1766", "SF1626")) %>% 
  apply(., 2, function(x){
    temp_data$GA_range[which(x == 1)]
  }) %>% 
  Reduce(intersect, .)

#####so we know that we should use the subjects: "SF1728", "SF1755", "SF1764", "SF1766", "SF1626", and time points
##"15_20" "20_25" "25_30" "30_35" "35_40" to construct the delta correlation network

metabolite_table_delta2 <- 
  metabolite_table_delta %>% 
  filter(., Patient_ID %in% c("SF1728", "SF1755", "SF1764", "SF1766", "SF1626")) %>% 
  filter(., GA_range %in% c("15_20", "20_25", "25_30", "30_35", "35_40")) %>% 
  # select(contains("_POS|subject_id")) %>% 
  plyr::dlply(., .(Patient_ID)) %>% 
  lapply(., function(x){
    x <- 
      x %>% 
      select(., matches("\\_POS|\\_NEG"))
  })


metabolite_table_delta2 <-
  # test <- 
  metabolite_table_delta2 %>% 
  lapply(., function(x){
    x1 <- x[-nrow(x),]
    x2 <- x[-1, ]
    x <- x2 - x1
  }) %>% 
  bind_rows()


delta_correlation_matrix <-
  metabolite_table_delta2 %>%
  as.matrix(.) %>%
  cor(., method = "spearman")


delta_correlation_data <- 
  delta_correlation_matrix 

rownames(delta_correlation_data) <- 
  colnames(delta_correlation_data)

##set lower tri as o
delta_correlation_data[lower.tri(delta_correlation_data)] <- 2

delta_correlation_data <- 
  delta_correlation_data %>% 
  as_tibble()

rownames(delta_correlation_data) <- 
  colnames(delta_correlation_data)

delta_correlation_data <- 
  delta_correlation_data %>% 
  rownames_to_column(., var = "Peak.name1") %>% 
  gather(., key = "Peak.name2", value = "Correlation", -Peak.name1) %>% 
  distinct() %>% 
  filter(., Correlation != 1 & Correlation != 2) %>% 
  arrange(., desc(abs(Correlation))) %>% 
  filter(abs(Correlation) > 0.6)


delta_correlation_p <-
  pbapply::pbapply(delta_correlation_data, 1, function(x) {
    peak1 <- as.character(x[1])
    peak2 <- as.character(x[2])
    int1 <-
      pull(metabolite_table_delta2, var = peak1)
    int2 <-
      pull(metabolite_table_delta2, var = peak2)
    p <-  cor.test(int1, int2, 
                   alternative = "two.sided",
                   method = "spearman")$p.value
    p
  })

sum(delta_correlation_p < 0.05)

###p adjustment
delta_correlation_p2 <- 
  p.adjust(delta_correlation_p, method = "BH")
sum(delta_correlation_p2 < 0.05)

delta_correlation_data <- 
  delta_correlation_data %>% 
  mutate(P.adjust = delta_correlation_p2) %>% 
  filter(P.adjust < 0.05 & abs(Correlation) > 0.6)

dim(delta_correlation_data)

###add compound name
delta_correlation_data <- 
  delta_correlation_data %>% 
  left_join(., smartd_rplc[,c("name", "Compound.name")], 
            by = c("Peak.name1" = "name")) %>% 
  dplyr::rename(Compound.name1 = Compound.name) %>% 
  left_join(., smartd_rplc[,c("name", "Compound.name")], 
            by = c("Peak.name2" = "name")) %>% 
  dplyr::rename(Compound.name2 = Compound.name)

######circos figure
library(circlize)
delta_correlation_data2 <- 
  delta_correlation_data[1:100,]
  # delta_correlation_data


ggplot(data = delta_correlation_data2, aes(x = Peak.name1, y = Peak.name2)) +
  geom_tile(aes(fill = Correlation)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw()

Metabolite <- 
  unique(c(delta_correlation_data2$Compound.name1, delta_correlation_data2$Compound.name2))
delta_temp_data <- data.frame(factor = "Metabolomics", 
                        x = 1:length(Metabolite), 
                        y = 1,
                        Metabolite = Metabolite)


library(circlize)
par(mar = c(2,2,2,2))
circos.par(start.degree = 75, clock.wise = TRUE, gap.after = 2.8)
circos.initialize(factors = delta_temp_data$factor, x = delta_temp_data$x)

circos.track(
  factors = delta_temp_data$factor,
  x = delta_temp_data$x,
  y = delta_temp_data$y - 1,
  ylim = c(0, 1),
  bg.border = NA,
  bg.col = NA,
  track.height = 0.4,
  # track.margin = c(0,0,0,0),
  panel.fun = function(x, y) {
    circos.points(x = x,
                  y = y,
                  pch = 19,
                  col = "black")
    circos.text(
      x = x,
      y = y + uy(2, "mm"),
      niceFacing = TRUE,
      facing = "clockwise",
      labels = delta_temp_data$Metabolite,
      adj = c(0, 0.5),
      col = "black",
      cex = 0.8
    )
  }
)

col_fun = colorRamp2(c(-1, 0, 1), c("#3C5488FF", "white", "#E64B35FF"))

for(i in 1:nrow(delta_correlation_data2)){
  cat(i, " ")
  point1 <- 
    match(delta_correlation_data[i,]$Compound.name1, delta_temp_data$Metabolite)
  point2 <- 
    match(delta_correlation_data[i,]$Compound.name2, delta_temp_data$Metabolite)
  correlation <- delta_correlation_data[i,]$Correlation
  col <- col_fun(x = correlation)
  circos.link(sector.index1 = "Metabolomics", point1 = point1, 
              sector.index2 = "Metabolomics", point2 = point2, 
              h = 1, col = col, lwd = 2)
}

circos.clear()


######network community analysis
library(igraph)

###construct igraph project of correlation data
delta_correlation_data3 <- 
  delta_correlation_data %>% 
  select(., from = Peak.name1, to = Peak.name2, 
         Correlation, P.adjust) %>% 
  as.data.frame()

node.attr <- smartd_rplc[,c(1:18)] %>% 
  filter(., name %in% (c(delta_correlation_data3$from, delta_correlation_data3$to)))

delta_correlation_data_igraph <- 
  graph_from_data_frame(d = delta_correlation_data3, 
                        directed = FALSE, 
                        vertices = node.attr
  )

# coords <- layout.davidson.harel(delta_correlation_data_igraph)
lay <- layout.auto(graph = delta_correlation_data_igraph)
plot(delta_correlation_data_igraph, 
     vertex.shape = "circle",
     vertex.label = NA,
     vertex.color = "#FFA319FF",
     vertex.size = vertex_attr(delta_correlation_data_igraph, "SS") * 10,
     # edge.color = "#767676FF",
     edge.color = ifelse(edge_attr(delta_correlation_data_igraph, "Correlation") > 0,
                         "#80000099", "#155F8399"),
     edge.width = abs(edge_attr(delta_correlation_data_igraph, "Correlation")) * 5,
     layout = lay
)

length(vertex_attr(delta_correlation_data_igraph, "name"))
dim(as_edgelist(delta_correlation_data_igraph))

##legend
# size_legend <-  
# quantile(x = vertex_attr(delta_correlation_data_igraph, "SS") * 10 %>% 
#              sort(), 
#            probs = c(0, 0.5, 1)) %>% 
#   round()
#   
# legend("topleft", 
#        legend = size_legend, 
#        bty = "n", 
#        pt.cex = size_legend/3,
#        pch = 19, 
#        # col = "black",
#        col = "#FFA319FF")



delta_subnetwork <- 
  cluster_edge_betweenness(graph = delta_correlation_data_igraph, 
                           weights = abs(edge_attr(delta_correlation_data_igraph, "Correlation")))

par(mar = c(5,5,4,2))

ggplot(
  data.frame(index = 1:length(delta_subnetwork$modularity),
             modu = delta_subnetwork$modularity, stringsAsFactors = FALSE),
  aes(index, modu) 
) +
  geom_vline(xintercept = which.max(delta_subnetwork$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line() +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))


membership(delta_subnetwork)
modularity(delta_subnetwork)
length(delta_subnetwork)
sizes(delta_subnetwork)
algorithm(delta_subnetwork)
merges(delta_subnetwork)
deltaing(communities = delta_subnetwork, graph = delta_correlation_data_igraph)
code_len(delta_subnetwork)
# is_hierarchical(delta_subnetwork)
# plot(as.dendrogram(delta_subnetwork))
# plot(as.hclust(delta_subnetwork))

which(sort(table(membership(communities = delta_subnetwork)), 
           decreasing = TRUE) > 2)
which(sort(table(membership(communities = delta_subnetwork)), 
           decreasing = TRUE) > 5)


temp_id <-
  sizes(delta_subnetwork) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 2) %>% 
  pull(., subnetwork_ID)



temp_subnetwork <-
lapply(temp_id, function(x){
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == x])
})


###node number
lapply(temp_subnetwork, function(x){
vertex_attr(graph = x, name = "name") %>% 
    length()
}) %>% 
  unlist()  %>% 
  sum()


###edge number
lapply(temp_subnetwork, function(x){
as_edgelist(x) %>% nrow()
}) %>% 
  unlist() %>% 
  sum()


temp_subnetwork <-
lapply(temp_id, function(x){
  subgraph(graph = cross_correlation_data_igraph, 
           v = names(membership(cross_subnetwork))[membership(cross_subnetwork) == x])
})


sizes(delta_subnetwork) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 2)

sizes(delta_subnetwork) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 5)

###there are 11 delta_subnetworks that with node bigger than 5.
##there are "5"  "1"  "4"  "7"  "2"  "10" "12" "6"  "18" "21" "20"

delta_subnetwork5 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 5])


save(delta_subnetwork5, file = "delta_subnetwork5")

lay <- layout.auto(graph = delta_subnetwork5)

delta_subnetwork5 %>% 
plot(., 
     vertex.label = NA,
     vertex.color = "#FFA319FF",
     vertex.size = vertex_attr(., "SS") * 10,
     # edge.color = "#767676FF",
     edge.color = ifelse(edge_attr(., "Correlation") > 0,
                         "#80000099", "#155F8399"),
     edge.width = abs(edge_attr(., "Correlation")) * 5,
     # edge.curved = 0.5,
     vertex.label.dist = 0.2,
     layout = lay
)


##############
##the interaction of subnetworks from cross and delta correlation network (vertexs)
# cross_subnetwork_member <- 
#   lapply(c("5", "1", "4", "7", "2","10","12","6","18","21","20"),
#          function(x){
#            vertex_attr(get(paste("cross_subnetwork", x, sep = "")), 
#                        name = "name")
#          })

# names(cross_subnetwork_member) <-
#   paste("Cross", 
#         c("5", "1", "4", "7", "2","10","12","6","18","21","20"),
#         sep = "_")
  


# delta_subnetwork_member <- 
#   lapply(c("5", "1", "4", "7", "2","10","12","6","18","21","20"),
#          function(x){
#            vertex_attr(get(paste("delta_subnetwork", x, sep= "")), 
#                        name = "name")
#          })

# names(delta_subnetwork_member) <-
#   paste("Delta",   c("5", "1", "4", "7", "2","10","12","6","18","21","20"), sep= "_")



####upset
# library(UpSetR)
# total_member <- 
#   unlist(c(cross_subnetwork_member, delta_subnetwork_member)) %>% 
#   unique()
# 
# data_for_upset1 <- 
# cross_subnetwork_member %>% 
#   lapply(., function(x){
#     sapply(total_member, function(y){
#       ifelse(y %in% x, 1, 0)
#     })
#   }) %>% 
#   bind_rows() %>% 
#   t() %>% 
#   as.data.frame()
# 
# 
# data_for_upset2 <- 
# delta_subnetwork_member %>% 
#   lapply(., function(x){
#     sapply(total_member, function(y){
#       ifelse(y %in% x, 1, 0)
#     })
#   }) %>% 
#   bind_rows() %>% 
#   t() %>% 
#   as.data.frame()
# 
# colnames(data_for_upset1) <-
#   colnames(data_for_upset2) <- 
#   total_member
# 
# data_for_upset <- 
#   rbind(data_for_upset1, data_for_upset2) %>% 
#   t() %>% 
#   as.data.frame()
#   
# 
# upset(data_for_upset,
#       group.by = "degree",
#       text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
#       sets = c(paste("Cross", c(5,26,53,9,33), sep = "_"),
#                paste("Delta", c(3,4,17,28,1), sep = "_")),
#       nintersects = NA, keep.order = TRUE,
#       # matrix.color = c("grey","#155F83CC"),
#       # main.bar.color = "red",
#       mainbar.y.label = "Intercetion",
#       sets.x.label = "Number of metabolites",
#       queries = list(
#         list(
#           query = intersects,
#           params = list("Delta_3", "Cross_5"),
#           active = TRUE,
#           color = "#FFA319CC"
#         )
#       )
# )


vertex_attr(cross_subnetwork11, name = "name")

cross_subnetwork11_member <- 
  vertex_attr(cross_subnetwork11, name = "name")
  
intersect(vertex_attr(cross_subnetwork11, name = "name"), 
          vertex_attr(delta_subnetwork5, name = "name"))


library(VennDiagram)

venn <- 
venn.diagram(x = list("Cross_11" = vertex_attr(cross_subnetwork11, name = "name"),
                      "Delta_5" = vertex_attr(delta_subnetwork5, name = "name")), 
             filename = NULL, 
             col = c("#8A9045CC", "#155F83CC"), lwd = 3,
             cex = 1.5, cat.cex = 1.5)

grid.draw(venn)






#####pathway enrichment for each delta_subnetwork
delta_subnetwork5_node_info <-
  vertex_attr(delta_subnetwork5) %>% 
  bind_cols()


write.csv(delta_subnetwork5_node_info,
          file = "delta_subnetwork5_node_info.csv",
          row.names = FALSE)

###delta_subnetwork5
delta_subnetwork5_ID_trans <-
  read.table(
    "delta_subnetwork5_ID_trans.txt",
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE
  )


delta_subnetwork5_delta_temp <- 
  delta_subnetwork5_node_info %>% 
  select(., name, Compound.name, HMDB.ID, KEGG.ID) %>% 
  left_join(., delta_subnetwork5_ID_trans, by = c("Compound.name" = "Query"))


# delta_subnetwork5_delta_temp <- 
delta_subnetwork5_delta_temp$KEGG[!is.na(delta_subnetwork5_delta_temp$KEGG) & delta_subnetwork5_delta_temp$KEGG == ''] <- NA

delta_subnetwork5_delta_temp$HMDB[!is.na(delta_subnetwork5_delta_temp$HMDB) & delta_subnetwork5_delta_temp$HMDB == ''] <- NA


delta_subnetwork5_delta_temp$KEGG

HMDB.KEGG <-
  apply(delta_subnetwork5_delta_temp, 1, function(x){
    HMDB1 <- x[3]
    KEGG1 <- x[4]
    HMDB2 <- x[6]
    KEGG2 <- x[8]
    
    if(is.na(HMDB1)){
      HMDB <- HMDB2
    }else{
      HMDB <- HMDB1
    }
    
    if(is.na(KEGG1)){
      KEGG <- KEGG2
    }else{
      KEGG <- KEGG1
    }
    
    c(HMDB, KEGG)
    
  }) %>% 
  t()


colnames(HMDB.KEGG) <- 
  c("HMDB.ID", "KEGG.ID")

cbind(delta_subnetwork5_node_info$HMDB.ID, HMDB.KEGG[,1])

delta_subnetwork5_node_info[,c("HMDB.ID", 'KEGG.ID')] <-
  HMDB.KEGG



write.csv(delta_subnetwork5_node_info, "delta_subnetwork5_node_info.csv", row.names = FALSE)

id <- delta_subnetwork5_node_info$KEGG.ID
id <- id[!is.na(id)]


delta_pathway5 <- 
  enrichPathway(id = id, database = hsa.kegg.pathway)

save(delta_pathway5, file = "delta_pathway5")


delta_pathway5 %>% 
  mutate(Marker = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(Overlap.frac * 100, -log(p.value, 10))) +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2, colour = "#FFA319FF") +
  geom_point(aes(size = Pathway.length, colour = Marker), shape = 16) +
  scale_colour_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  guides(colour = FALSE, 
         size = guide_legend(title = "Pathway size",
                             title.theme = element_text(size = 15),
                             label.theme = element_text(size = 13))) +
  labs(x = "Overlap (%)", y = "-log10(P value)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 13),
        legend.position = c(1, 0), 
        legend.justification = c(1,0),
        legend.background = element_blank()) +
  ggrepel::geom_label_repel(aes(label = Pathway.name), 
                            data = delta_pathway5 %>% filter(., p.value < 0.05))

delta_pathway5 %>% 
  mutate(log.p = -log(p.value, 10)) %>% 
  mutate(enrichment = ifelse(p.value < 0.05, "Yes", "No")) %>% 
  # `[`(., 1:15, ) %>% 
  # top_n(., 11, wt = "p.value") %>% 
  arrange(., log.p) %>% 
  ggplot(., aes(x = factor(Pathway.name, levels = Pathway.name), 
                y = log.p)) +  
geom_bar(stat = "identity", width = 0.8, 
         aes(fill = enrichment), show.legend = FALSE,
         # fill = "#FFA319FF"
         ) +
  scale_fill_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  labs(x = "", y = "-log10(P value)") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10))


###calculate the betweness and degree and all other
delta_betweenness5 <- betweenness(graph = delta_subnetwork5)
delta_degree5 <- igraph::degree(graph = delta_subnetwork5)
delta_closeness5 <- closeness(graph = delta_subnetwork5)

delta_importance5 <- 
  tibble(delta_betweenness5, delta_degree5, delta_closeness5) %>% 
  mutate_all(., function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  mutate(.,mean = apply(., 1 ,mean)) %>% 
  mutate(peak.name = vertex_attr(graph = delta_subnetwork5, name = "name")) %>% 
  left_join(., delta_subnetwork5_node_info[,c("name", "Compound.name", "KEGG.ID")], 
            by = c('peak.name' = "name")) %>%  
  select(peak.name, Compound.name, everything()) %>% 
  arrange(., desc(mean))


# grep("Cysteine and methionine metabolism", names(hsa.kegg.pathway)) %>% 
#   `[[`(hsa.kegg.pathway, .) %>% 
#   intersect(., delta_importance3[c(1:50),]$KEGG.ID)
# 
# c('C00021', 'C00022')


delta_importance5[c(1:30),] %>% 
  arrange(., mean) %>% 
  mutate(Compound.name = factor(Compound.name, levels = Compound.name)) %>% 
  # mutate(node = ifelse(is.element(KEGG.ID, c('C00021', 'C00022')), 'Yes', "No")) %>% 
  ggplot(.,aes(x = Compound.name, y = mean)) +
  # geom_point(size = 2, aes(colour = node)) +
  geom_point(size = 2, colour = "#155F83FF") +
  # scale_y_continuous(expand = c(0,0.1)) +
  labs(x = "", y = "Mean (Betweenness, Degree, Closeness)") +
  geom_segment(aes(x = Compound.name, y = 0, 
                   xend = Compound.name, yend = mean),
               colour = "#155F83FF") +
  scale_colour_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 10), legend.position = "none", 
        axis.text.x = element_text())


###betweenness, degree and closness
# test <- 
delta_importance5[c(1:30),] %>% 
  select(., Compound.name:mean) %>% 
  arrange(., mean) %>% 
  gather(., key = "Class", value = "Value", -Compound.name) %>% 
  mutate(Compound.name = factor(Compound.name, levels = Compound.name[1:30])) %>% 
  filter(., Class != "mean") %>% 
  ggplot(.,aes(x = Compound.name, 
               y = Value)) +
  geom_point(size = 2, 
             aes(colour = Class)
             # colour = "#155F83FF"
  ) +
  geom_segment(aes(x = Compound.name, y = 0, 
                   xend = Compound.name, 
                   yend = Value,
                   colour = Class)
               # colour = "#155F83FF"
  ) +
  scale_color_manual(values = c(
    "delta_betweenness5" = "#FFA319FF",
    "delta_closeness5" = "#8A9045FF",
    "delta_degree5" = "#C16622FF"
  )) +
  facet_grid(. ~ Class, 
             labeller = as_labeller(
               c("delta_betweenness5" = "Betweenness",
                 "delta_closeness5" = "Closeness",
                 "delta_degree5" = "Degree"
               )
             )
  ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 10), legend.position = "none", 
        axis.text.x = element_text(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 13)) +
  coord_flip()

###find the nodes with top 3 imporance
delta_importance5[1,]$peak.name

delta_temp_subgraph <- 
  neighbors(graph = delta_subnetwork5, 
            v = delta_importance5[1,]$peak.name) %>% 
  names() %>% 
  `c`(., delta_importance5[1,]$peak.name) %>% 
  unique() %>% 
  subgraph(delta_subnetwork5, v = .)


delta_temp_subgraph <- 
  set_vertex_attr(graph = delta_temp_subgraph, 
                  name = "label", 
                  value = ifelse(vertex_attr(delta_temp_subgraph, "name") == delta_importance5[1,]$peak.name, "Yes", "No"))

##get top 20 correlations
E(delta_temp_subgraph)
delta_temp_subgraph <- 
  as_long_data_frame(delta_temp_subgraph) %>% 
  select(., from_name, to_name, Correlation) %>%
  filter(., from_name == delta_importance5[1,]$peak.name | to_name == delta_importance5[1,]$peak.name) %>% 
  mutate(Correlation = abs(Correlation)) %>% 
  arrange(., desc(Correlation)) %>% 
  top_n(., 20)


delta_temp_subgraph <- 
  subgraph(graph = delta_subnetwork5, 
           v = unique(c(delta_temp_subgraph$from_name, delta_temp_subgraph$to_name)))


delta_temp_subgraph <- 
  set_vertex_attr(graph = delta_temp_subgraph, name = "label", 
                  value = ifelse(vertex_attr(delta_temp_subgraph, "name") == delta_importance5[1,]$peak.name, "Yes", "No"))

which(vertex_attr(delta_temp_subgraph, "name") == delta_importance5[1,]$peak.name)

lay <- layout_as_tree(delta_temp_subgraph, root = 21)

plot(delta_temp_subgraph,
     vertex.color = ifelse(vertex_attr(delta_temp_subgraph, "label")=="Yes", "#FFA31999", "#8A904599"),
     vertex.label = vertex_attr(delta_temp_subgraph, "Compound.name"),
     vertex.label.color = "black",
     vertex.size = vertex_attr(delta_temp_subgraph, "SS") * 20,
     vertex.label.dist = 1,
     
     # edge.color = "#767676FF",
     edge.color = ifelse(edge_attr(delta_temp_subgraph, "Correlation") > 0,
                         "#80000099", "#155F8399"),
     edge.width = abs(edge_attr(delta_temp_subgraph, "Correlation")) * 5,
     layout = -lay[, 2:1]
)









