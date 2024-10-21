#######construct delta correlation network
#####smartd_rplc_batch dataset is used to analysis for correlation networks analysis
#### step 1 remove the duplicated metabolites according to score
delta_temp <- 
  smartd_rplc_batch1 %>% 
  select(-(name:Database))

table(which(is.na(delta_temp), arr.ind = TRUE)[,2])
colnames(delta_temp)[c(146, 147, 148)]
sum(is.na(delta_temp))

##SFU65, SFU74, SFU91 have no positive data, so remove it
delta_temp <- 
  delta_temp %>% 
  select(., -c(SFU65, SFU74, SFU91))

###calculate correlation matrix
###delta network
sfu1_148 <- 
  readr::read_csv("E:/project/smartD/patient_information/SFU1-148_GA.csv")

patient_info <- 
  readr::read_csv("E:/project/smartD/data_analysis20190828/patient_info/patient_info.csv")

sfu1_148 <- 
  sfu1_148 %>% 
  mutate(subject_id = as.character(subject_id), visit = as.character(visit)) %>% 
  arrange(subject_id)

patient_info <- 
  patient_info %>% 
  mutate(Patient_ID = as.character(Patient_ID), Visit = as.character(Visit)) %>% 
  arrange(Patient_ID)

match(colnames(delta_temp), sfu1_148$sample_id)


delta_temp <- 
  log(delta_temp, 10) %>% 
  as.data.frame(delta_temp)

####log 10 and scale
delta_temp <-
  apply(delta_temp, 1, function(x){
    (x - mean(x))/sd(x)
  })

delta_temp <- 
  delta_temp %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_name")

colnames(delta_temp)[-1] <- 
  smartd_rplc_batch1$name

delta_temp <- 
  inner_join(x = sfu1_148[,-1], delta_temp, by = c("sample_id" = "Sample_name"))

####according to the GA information, we shoud class them to different peroid
delta_temp %>% 
  ggplot(., aes(x = GA_week)) +
  geom_histogram(colour = "white", binwidth = 5, 
                 fill = "#8A9045CC", boundary = 10) +
  theme_bw() +
  labs(x = "GA (week)") +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15))

hist(delta_temp$GA_week)
range(delta_temp$GA_week)

library(plyr)

lapply(delta_temp %>%
         plyr::dlply(., .(subject_id)),
       function(x)
         x$GA_week)

delta_temp <- 
  delta_temp %>% 
  mutate(GA_range = case_when(delta_temp$GA_week > 10 & delta_temp$GA_week <= 15 ~ {"10_15"},
                              delta_temp$GA_week > 15 & delta_temp$GA_week <= 20 ~ {"15_20"},
                              delta_temp$GA_week > 20 & delta_temp$GA_week <= 25 ~ {"20_25"},
                              delta_temp$GA_week > 25 & delta_temp$GA_week <= 30 ~ {"25_30"},
                              delta_temp$GA_week > 30 & delta_temp$GA_week <= 35 ~ {"30_35"},
                              delta_temp$GA_week > 35 & delta_temp$GA_week <= 40 ~ {"35_40"},
                              delta_temp$GA_week > 40 & delta_temp$GA_week <= 45 ~ {"40_45"}
                              ###10-15,15-20,20-25,25-30,30-35,35-40,40-45
  ))


lapply(
  plyr::dlply(delta_temp, .variables = .(subject_id)),
  function(x) x$GA_range
)


lapply(
  plyr::dlply(delta_temp, .variables = .(subject_id)),
  function(x) x$GA_week
)


####round the GA
delta_temp <- 
  delta_temp %>% 
  mutate(GA_round = round(GA_week))

lapply(
  plyr::dlply(delta_temp, .variables = .(subject_id)),
  function(x) x$GA_round
) %>% 
  unlist() %>% 
  table() %>% 
  sort()

delta_temp <-
  delta_temp %>% 
  plyr::dlply(., .(subject_id)) %>% 
  lapply(., function(x){
    x <- x %>% arrange(., GA_range)
    x <- x[!duplicated(x$GA_range),,drop = FALSE]
  })

delta_temp <-
  delta_temp %>% 
  dplyr::bind_rows()

####upset to show the venn diagram
library(UpSetR)
delta_temp_data <-
  lapply(
    plyr::dlply(delta_temp, .variables = .(subject_id)),
    function(x) x$GA_range
  )

total_week <-
  unique(unlist(delta_temp_data)) %>%
  unname() %>%
  sort()

delta_temp_data <-
  lapply(delta_temp_data, function(x){
    sapply(total_week, function(y) {
      if(y %in% x){
        1
      } else{
        0
      }
    })
  })

delta_temp_data <-
  delta_temp_data %>%
  bind_cols()

# delta_temp_data <-
#   t(delta_temp_data) %>%
#   as.data.frame()
# colnames(delta_temp_data) <- total_week

rownames(delta_temp_data) <- total_week

delta_temp_data <-
  delta_temp_data %>%
  rownames_to_column(., var = "GA_range")

upset(as.data.frame(delta_temp_data),
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
      # sets = colnames(delta_temp_data)[-1],
      nintersects = NA, keep.order = FALSE,
      # matrix.color = c("grey","#155F83CC"),
      # main.bar.color = "red",
      mainbar.y.label = "Intercetion",
      sets.x.label = "Number of samples",
      queries = list(
        list(
          query = intersects,
          params = list("1554", "1555", "1556", "1650", "1626"),
          active = TRUE,
          color = "#FFA319CC"
        )
      )
)


#####
delta_temp_data %>% 
  select(., c("1554", "1555", "1556", "1650", "1626")) %>% 
  apply(., 2, function(x){
    delta_temp_data$GA_range[which(x == 1)]
  }) %>% 
  Reduce(intersect, .)

#####so we know that we should use the subjects: 1554, 1555, 1556, 1650, 1626, and time points
##"15_20" "20_25" "25_30" "30_35" "35_40" to construct the delta correlation network

delta_temp2 <- 
  delta_temp %>% 
  filter(., subject_id %in% c("1554", "1555", "1556", "1650", "1626")) %>% 
  filter(., GA_range %in% c("15_20", "20_25", "25_30", "30_35", "35_40")) %>% 
  # select(contains("_POS|subject_id")) %>% 
  plyr::dlply(., .(subject_id)) %>% 
  lapply(., function(x){
    x <- 
      x %>% 
      select(., matches("\\_POS|\\_NEG"))
  })


delta_temp2 <-
  # test <- 
  delta_temp2 %>% 
  lapply(., function(x){
    x1 <- x[-nrow(x),]
    x2 <- x[-1, ]
    x <- x2 - x1
  }) %>% 
  bind_rows()


setwd("E:/project/smartD/data_analysis20190828/correlation_network_analysis/identification_table/POS_NEG/delta_netowork")

delta_correlation_matrix <-
  delta_temp2 %>%
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
      pull(delta_temp2, var = peak1)
    int2 <-
      pull(delta_temp2, var = peak2)
    p <-  cor.test(int1, int2, 
                   alternative = "two.sided",
                   method = "spearman")$p.value
    p
  })

###p adjustment
delta_correlation_p2 <- 
  p.adjust(delta_correlation_p, method = "BH")

delta_correlation_data <- 
  delta_correlation_data %>% 
  mutate(P.adjust = delta_correlation_p2) %>% 
  filter(P.adjust < 0.05 & abs(Correlation) > 0.6)

###add compound name
delta_correlation_data <- 
  delta_correlation_data %>% 
  left_join(., smartd_rplc_batch1[,c("name", "Compound.name")], 
            by = c("Peak.name1" = "name")) %>% 
  dplyr::rename(Compound.name1 = Compound.name) %>% 
  left_join(., smartd_rplc_batch1[,c("name", "Compound.name")], 
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

node.attr <- smartd_rplc_batch1[,c(1:18)] %>% 
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

###there are 26 delta_subnetworks that with node bigger than 5.
##there are "3"  "4"  "17" "28" "1"  "12" "22" "20" "11" "21" "19" "33" "8"  "16" "27" "18" "10" "29"
## "32" "15" "24" "41" "36" "9"  "14" "23"

delta_subnetwork3 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 3])

delta_subnetwork4 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 4])

delta_subnetwork17 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 17])

delta_subnetwork28 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 28])

delta_subnetwork1 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 1])

delta_subnetwork12 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 12])

delta_subnetwork22 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 22])

delta_subnetwork20 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 20])

delta_subnetwork11 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 11])

delta_subnetwork21 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 21])

delta_subnetwork19 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 19])

delta_subnetwork33 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 33])

delta_subnetwork8 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 8])

delta_subnetwork16 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 16])

delta_subnetwork27 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 27])

delta_subnetwork18 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 18])

delta_subnetwork10 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 10])

delta_subnetwork29 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 29])

delta_subnetwork32 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 32])

delta_subnetwork15 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 15])

delta_subnetwork24 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 24])

delta_subnetwork41 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 41])

delta_subnetwork36 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 36])

delta_subnetwork9 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 9])

delta_subnetwork14 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 14])

delta_subnetwork23 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 23])

delta_subnetwork3 <- 
  subgraph(graph = delta_correlation_data_igraph, 
           v = names(membership(delta_subnetwork))[membership(delta_subnetwork) == 3])

save(delta_subnetwork3, file = "delta_subnetwork3")
save(delta_subnetwork4, file = "delta_subnetwork4")
save(delta_subnetwork17, file = "delta_subnetwork17")
save(delta_subnetwork28, file = "delta_subnetwork28")
save(delta_subnetwork1, file = "delta_subnetwork1")
save(delta_subnetwork12, file = "delta_subnetwork12")
save(delta_subnetwork22, file = "delta_subnetwork22")
save(delta_subnetwork20, file = "delta_subnetwork20")
save(delta_subnetwork11, file = "delta_subnetwork11")
save(delta_subnetwork21, file = "delta_subnetwork21")
save(delta_subnetwork19, file = "delta_subnetwork19")
save(delta_subnetwork33, file = "delta_subnetwork33")
save(delta_subnetwork8, file = "delta_subnetwork8")
save(delta_subnetwork16, file = "delta_subnetwork16")
save(delta_subnetwork27, file = "delta_subnetwork27")
save(delta_subnetwork18, file = "delta_subnetwork18")
save(delta_subnetwork10, file = "delta_subnetwork10")
save(delta_subnetwork29, file = "delta_subnetwork29")
save(delta_subnetwork32, file = "delta_subnetwork32")
save(delta_subnetwork15, file = "delta_subnetwork15")
save(delta_subnetwork24, file = "delta_subnetwork24")
save(delta_subnetwork41, file = "delta_subnetwork41")
save(delta_subnetwork36, file = "delta_subnetwork36")
save(delta_subnetwork9, file = "delta_subnetwork9")
save(delta_subnetwork14, file = "delta_subnetwork14")
save(delta_subnetwork23, file = "delta_subnetwork23")


lay <- layout.auto(graph = delta_subnetwork3)

delta_subnetwork3 %>% 
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
cross_subnetwork_member <- 
  lapply(c("5", "26","53","9","33","37","43","41","55","71","29","31","49","84","96"),
         function(x){
           vertex_attr(get(paste("cross_subnetwork", x, sep= "")), 
                       name = "name")
         })

names(cross_subnetwork_member) <-
  paste("Cross", 
        c("5", "26","53","9","33","37","43","41","55","71","29","31","49","84","96"),
        sep = "_")
  


delta_subnetwork_member <- 
  lapply(c("3", "4", "17", "28", "1", "12", "22", "20", "11", "21", "19", "33",
           "8","16", "27", "18", "10", "29",
          "32", "15", "24", "41", "36", "9", "14", "23"),
         function(x){
           vertex_attr(get(paste("delta_subnetwork", x, sep= "")), 
                       name = "name")
         })

names(delta_subnetwork_member) <-
  paste("Delta",   c("3", "4", "17", "28", "1", "12", "22", "20", "11", "21", "19", "33",
           "8","16", "27", "18", "10", "29",
          "32", "15", "24", "41", "36", "9", "14", "23"), sep= "_")



####upset
library(UpSetR)
total_member <- 
  unlist(c(cross_subnetwork_member, delta_subnetwork_member)) %>% 
  unique()

data_for_upset1 <- 
cross_subnetwork_member %>% 
  lapply(., function(x){
    sapply(total_member, function(y){
      ifelse(y %in% x, 1, 0)
    })
  }) %>% 
  bind_rows() %>% 
  t() %>% 
  as.data.frame()


data_for_upset2 <- 
delta_subnetwork_member %>% 
  lapply(., function(x){
    sapply(total_member, function(y){
      ifelse(y %in% x, 1, 0)
    })
  }) %>% 
  bind_rows() %>% 
  t() %>% 
  as.data.frame()

colnames(data_for_upset1) <-
  colnames(data_for_upset2) <- 
  total_member

data_for_upset <- 
  rbind(data_for_upset1, data_for_upset2) %>% 
  t() %>% 
  as.data.frame()
  

upset(data_for_upset,
      group.by = "degree",
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
      sets = c(paste("Cross", c(5,26,53,9,33), sep = "_"),
               paste("Delta", c(3,4,17,28,1), sep = "_")),
      nintersects = NA, keep.order = TRUE,
      # matrix.color = c("grey","#155F83CC"),
      # main.bar.color = "red",
      mainbar.y.label = "Intercetion",
      sets.x.label = "Number of metabolites",
      queries = list(
        list(
          query = intersects,
          params = list("Delta_3", "Cross_5"),
          active = TRUE,
          color = "#FFA319CC"
        )
      )
)




vertex_attr(cross_subnetwork5, name = "name")

cross_subnetwork5_member <- 
  vertex_attr(cross_subnetwork5, name = "name")
  
intersect(vertex_attr(cross_subnetwork5, name = "name"), 
          vertex_attr(delta_subnetwork3, name = "name"))


library(VennDiagram)

venn <- 
venn.diagram(x = list("Cross_5" = vertex_attr(cross_subnetwork5, name = "name"),
                      "Delta_3" = vertex_attr(delta_subnetwork3, name = "name")), 
             filename = NULL, 
             col = c("#8A9045CC", "#155F83CC"), lwd = 3,
             cex = 1.5, cat.cex = 1.5)

grid.draw(venn)






#####pathway enrichment for each delta_subnetwork
delta_subnetwork3_node_info <-
  vertex_attr(delta_subnetwork3) %>% 
  bind_cols()


write.csv(delta_subnetwork3_node_info,
          file = "delta_subnetwork3_node_info.csv",
          row.names = FALSE)

###delta_subnetwork3
delta_subnetwork3_ID_trans <-
  read.table(
    "delta_subnetwork3_ID_trans.txt",
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE
  )


delta_subnetwork3_delta_temp <- 
  delta_subnetwork3_node_info %>% 
  select(., name, Compound.name, HMDB.ID, KEGG.ID) %>% 
  left_join(., delta_subnetwork3_ID_trans, by = c("Compound.name" = "Query"))


# delta_subnetwork3_delta_temp <- 
delta_subnetwork3_delta_temp$KEGG[!is.na(delta_subnetwork3_delta_temp$KEGG) & delta_subnetwork3_delta_temp$KEGG == ''] <- NA

delta_subnetwork3_delta_temp$HMDB[!is.na(delta_subnetwork3_delta_temp$HMDB) & delta_subnetwork3_delta_temp$HMDB == ''] <- NA


delta_subnetwork3_delta_temp$KEGG

HMDB.KEGG <-
  apply(delta_subnetwork3_delta_temp, 1, function(x){
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

cbind(delta_subnetwork3_node_info$HMDB.ID, HMDB.KEGG[,1])

delta_subnetwork3_node_info[,c("HMDB.ID", 'KEGG.ID')] <-
  HMDB.KEGG



write.csv(delta_subnetwork3_node_info, "delta_subnetwork3_node_info.csv", row.names = FALSE)

id <- delta_subnetwork3_node_info$KEGG.ID
id <- id[!is.na(id)]


delta_pathway3 <- 
  enrichPathway(id = id, database = hsa.kegg.pathway)

save(delta_pathway3, file = "delta_pathway3")


delta_pathway3 %>% 
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
                            data = delta_pathway3 %>% filter(., p.value < 0.05))

delta_pathway3 %>% 
  mutate(log.p = -log(p.value, 10)) %>% 
  mutate(enrichment = ifelse(p.value < 0.05, "Yes", "No")) %>% 
  `[`(., 1:15, ) %>% 
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
delta_betweenness3 <- betweenness(graph = delta_subnetwork3)
delta_degree3 <- igraph::degree(graph = delta_subnetwork3)
delta_closeness3 <- closeness(graph = delta_subnetwork3)

delta_importance3 <- 
  tibble(delta_betweenness3, delta_degree3, delta_closeness3) %>% 
  mutate_all(., function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  mutate(.,mean = apply(., 1 ,mean)) %>% 
  mutate(peak.name = vertex_attr(graph = delta_subnetwork3, name = "name")) %>% 
  left_join(., delta_subnetwork3_node_info[,c("name", "Compound.name", "KEGG.ID")], 
            by = c('peak.name' = "name")) %>%  
  select(peak.name, Compound.name, everything()) %>% 
  arrange(., desc(mean))


# grep("Cysteine and methionine metabolism", names(hsa.kegg.pathway)) %>% 
#   `[[`(hsa.kegg.pathway, .) %>% 
#   intersect(., delta_importance3[c(1:50),]$KEGG.ID)
# 
# c('C00021', 'C00022')


delta_importance3[c(1:30),] %>% 
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
delta_importance3[c(1:30),] %>% 
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
    "delta_betweenness3" = "#FFA319FF",
    "delta_closeness3" = "#8A9045FF",
    "delta_degree3" = "#C16622FF"
  )) +
  facet_grid(. ~ Class, 
             labeller = as_labeller(
               c("delta_betweenness3" = "Betweenness",
                 "delta_closeness3" = "Closeness",
                 "delta_degree3" = "Degree"
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
delta_importance3[1,]$peak.name

delta_temp_subgraph <- 
  neighbors(graph = delta_subnetwork3, 
            v = delta_importance3[1,]$peak.name) %>% 
  names() %>% 
  `c`(., delta_importance3[1,]$peak.name) %>% 
  unique() %>% 
  subgraph(delta_subnetwork3, v = .)


delta_temp_subgraph <- 
  set_vertex_attr(graph = delta_temp_subgraph, name = "label", 
                  value = ifelse(vertex_attr(delta_temp_subgraph, "name") == delta_importance3[1,]$peak.name, "Yes", "No"))

##get top 20 correlations
E(delta_temp_subgraph)
delta_temp_subgraph <- 
  as_long_data_frame(delta_temp_subgraph) %>% 
  select(., from_name, to_name, Correlation) %>%
  filter(., from_name == delta_importance3[1,]$peak.name | to_name == delta_importance3[1,]$peak.name) %>% 
  mutate(Correlation = abs(Correlation)) %>% 
  arrange(., desc(Correlation)) %>% 
  top_n(., 20)


delta_temp_subgraph <- 
  subgraph(graph = delta_subnetwork3, 
           v = unique(c(delta_temp_subgraph$from_name, delta_temp_subgraph$to_name)))


delta_temp_subgraph <- 
  set_vertex_attr(graph = delta_temp_subgraph, name = "label", 
                  value = ifelse(vertex_attr(delta_temp_subgraph, "name") == delta_importance3[1,]$peak.name, "Yes", "No"))

which(vertex_attr(delta_temp_subgraph, "name") == delta_importance3[1,]$peak.name)

lay <- layout_as_tree(delta_temp_subgraph, root = 13)

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









