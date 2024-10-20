setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/")
rm(list=ls())
load("data_preparation_for_analysis/metabolite_table")
load("data_preparation_for_analysis/metabolite_tags")
load("data_preparation_for_analysis/peak_table")

info <-
  readxl::read_xlsx("E:/project/smartD/patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

sample_info <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

setwd(r4projects::get_project_wd())
marker1 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/GA_prediction/marker_rf_final.csv")
marker2 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/marker_rf_final.csv")

marker1$name

marker2$name


intersect(marker1$name, marker2$name)

library(VennDiagram)

venn <- VennDiagram::venn.diagram(x = list("Mraker 1" = marker1$name, 
                                           "Mraker 2" = marker2$name), 
                                  filename = NULL, col = c("#C71000FF", "#84D7E1FF"), 
                                  main.cex = 1.5, lwd = 2)


grid.draw(venn)


setdiff(marker1$Compound.name, marker2$Compound.name)
setdiff(marker2$Compound.name, marker1$Compound.name)


library(waffle)

test <- c(marker1$super_class, marker2$super_class)
test[is.na(test)] <- "Unknown"

test <- table(test)
test1 <- as.numeric(test)
names(test1) <- names(test)

waffle(test, 
       colors = values, rows = 6)




values = c(
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




setwd(r4projects::get_project_wd())
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


##rename P
str_sort(grep("P", colnames(marker_data), value = T), numeric = TRUE)
sample_info %>% 
  filter(is.na(GA)) %>% 
  pull(Sample_ID) %>% 
  str_sort(numeric = TRUE)

##So P1 is X178 and P2 is X179 and so on.


colnames(marker_data)[grep("P", colnames(marker_data))] <-
  str_replace(colnames(marker_data)[grep("P", colnames(marker_data))],"P", "") %>% 
    as.numeric() %>% 
    `+`(., 177) %>% 
    paste("X", ., sep = "")

marker_data <-
  marker_data %>% 
  select(one_of(sample_info$Sample_ID))



library(pheatmap)

# marker_data <- 
#   log(marker_data, 10)

sum(is.na(marker_data))

marker_data <- 
  apply(marker_data, 1, function(x){
    (x - mean(x))/sd(x)
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


marker_data <-
  data.frame(ga, t(marker_data), stringsAsFactors = FALSE) %>% 
  split(f = ga) %>% 
    lapply(function(x){
      apply(x[,-1], 2, mean)
    }) %>% 
    do.call(rbind, .)


range(marker_data)


marker_data[which(marker_data < -1.5, arr.ind = TRUE)] <- -1.5


colnames(marker_data) <- 
metabolite_tags$Compound.name[
  match(colnames(marker_data), 
        metabolite_tags$name)  
]

# "ward.D", "ward.D2", 
# "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).


  #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"

annotation_col <- 
marker %>% 
  select(Compound.name, super_class) %>% 
  mutate(super_class = case_when(is.na(super_class) ~ "Unknown",
                                 TRUE ~ super_class)) %>% 
  column_to_rownames(var = 'Compound.name')

anno_colors <- 
  list(super_class = c(
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
  ))


pheatmap(marker_data, 
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

rownames <- sort(rownames(marker_data))
rownames <- c('PP',"[11,13]", rownames[-c(15,16)])

marker_data <-
  marker_data[rownames,]

library(ggcorrplot)

corr <- 
  cor(t(marker_data))
p.mat <- cor_pmat(t(marker_data))

library(corrplot)

plot <- 
ggcorrplot(corr = as.matrix(corr),
           legend.title = "Correlation",
           hc.order = FALSE,
           outline.col = "grey",
           type = "upper",
           ggtheme = ggplot2::theme_bw,
           colors = c("#008EA0B2", 
                      alpha(colour = "#008EA0B2", alpha = 0.4),
                      alpha(colour = "#FF6F00B2", alpha = 0.4),
                      "#FF6F00B2"),
           lab = TRUE,
           lab_col = "white",
           lab_size = 3,
           p.mat = p.mat, 
           sig.level = 0.05, 
           insig = "pch", 
           pch.col = "red")


text_colour <- c(
  "red",
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(15)
)
  


plot <- 
plot +
  theme(axis.text = element_text(colour = text_colour),
        panel.grid.minor  = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_blank()
        )


distance <- 
  as.matrix(dist(marker_data))[,1][-1]
  

line_cols <- 
  colorRampPalette(colors = 
                     c("royalblue", "red")
                     )(15)
  
p.value <- -log(p.mat[-1,1], 10)

pp_data <-
  data.frame(x = rep(13, 15), 
             y = rep(3, 15), 
             xend = c(1:15) + 0.5,
             yend = c(1:15) - 0.5,
             curvature = seq(from = -0.14, to = 0.14, length.out = 15),
             colour = line_cols,
             p.value = p.value,
             distance = distance,
             text_colour = text_colour[-1],
             stringsAsFactors = FALSE)


for(i in 1:nrow(pp_data)){
  c(
    plot <- 
    plot +
    annotate(geom = "curve", 
             x = pp_data$x[i], 
             y = pp_data$y[i], 
             xend = pp_data$xend[i], 
             yend = pp_data$yend[i],
             curvature  = pp_data$curvature[i],
             colour = pp_data$colour[i],
             size = 3
    ) +
      annotate(geom = "point", x = i + 0.5, y = i - 0.5,
               colour = pp_data$text_colour[i],
               size = 3)
  )
}



plot <-
  plot +
  annotate(geom = "point", x = 13, y = 3, colour = "grey", size = 5)


corr_plot <- plot
corr_plot


##legend for corrplot
ggplot(pp_data) +
  geom_segment(aes(x = x, 
                   y = y, 
                   xend = xend, 
                   yend = yend,
                   curvature  = curvature,
                   colour = distance,
                   ),size = 3) +
  scale_colour_gradientn(colours = pp_data$colour)








##markers in different stages
marker_data <- 
  metabolite_table %>% 
  dplyr::filter(name %in% marker$name) %>% 
  column_to_rownames("name") %>% 
  dplyr::select(-c(mz, rt))

##rename P
str_sort(grep("P", colnames(marker_data), value = T), numeric = TRUE)
sample_info %>% 
  filter(is.na(GA)) %>% 
  pull(Sample_ID) %>% 
  str_sort(numeric = TRUE)

##So P1 is X178 and P2 is X179 and so on.

colnames(marker_data)[grep("P", colnames(marker_data))] <-
  str_replace(colnames(marker_data)[grep("P", colnames(marker_data))],"P", "") %>% 
  as.numeric() %>% 
  `+`(., 177) %>% 
  paste("X", ., sep = "")

marker_data <-
  marker_data %>% 
  select(one_of(sample_info$Sample_ID))


sum(is.na(marker_data))

marker_data <- 
  apply(marker_data, 1, function(x){
    (x - mean(x))/sd(x)
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


rownames(marker_data) <- 
metabolite_tags$Compound.name[match(rownames(marker_data), metabolite_tags$name)]
  
marker_data2 <-
  data.frame(ga, t(marker_data), stringsAsFactors = FALSE) %>% 
  split(f = ga) %>% 
  lapply(function(x){
    x[,-1]
  })



#####
text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(15),
  "red"
)

box_colour <-
  text_colour

names(box_colour) <- ga_level

idx <- 1
test <- 
mapply(FUN = function(x,y){
  list(data.frame(GA = y, value = x[,idx], stringsAsFactors = FALSE))
},
x = marker_data2,
y = names(marker_data2))


test <- do.call(rbind, test)

ga_level <- unique(test$GA)
ga_level <- 
  c("[11,13]", ga_level[-c(15,16)], "PP")
test$GA <- factor(test$GA, levels = ga_level)



ggplot(data = test, aes(x = GA, y = value)) +
  geom_boxplot(aes(colour = GA), 
               outlier.shape = NA, show.legend = FALSE) +
  # geom_jitter(aes(colour = GA)) +
  scale_colour_manual(values = box_colour) +
  scale_y_continuous(limits = c(-1 ,1)) +
  theme_bw() +
  labs(x = "", y = "Scaled intensity", 
       title = rownames(marker_data)[idx]) +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1, 
                                   colour = text_colour),
        axis.text.y = element_text(size = 13), 
        plot.title = element_text(size = 15, hjust = 0.5))

rownames(marker_data)[idx]







#-------------------------------------------------------------------------------
###fuzzy C-means for metabolites
library(Mfuzz)
library(e1071)

setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/")
rm(list=ls())
load("data_preparation_for_analysis/metabolite_table")
load("data_preparation_for_analysis/metabolite_tags")
load("data_preparation_for_analysis/peak_table")

info <-
  readxl::read_xlsx("E:/project/smartD/patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

sample_info <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

setwd(r4projects::get_project_wd())
marker1 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/GA_prediction/marker_rf_final.csv")
marker2 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/marker_rf_final.csv")

setwd(r4projects::get_project_wd())
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


##rename P
str_sort(grep("P", colnames(marker_data), value = T), numeric = TRUE)
sample_info %>% 
  filter(is.na(GA)) %>% 
  pull(Sample_ID) %>% 
  str_sort(numeric = TRUE)

##So P1 is X178 and P2 is X179 and so on.


colnames(marker_data)[grep("P", colnames(marker_data))] <-
  str_replace(colnames(marker_data)[grep("P", colnames(marker_data))],"P", "") %>% 
  as.numeric() %>% 
  `+`(., 177) %>% 
  paste("X", ., sep = "")

marker_data <-
  marker_data %>% 
  select(one_of(sample_info$Sample_ID))


sum(is.na(marker_data))

marker_data <- 
  apply(marker_data, 1, function(x){
    (x - mean(x))/sd(x)
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


marker_data <-
  data.frame(ga, t(marker_data), stringsAsFactors = FALSE) %>% 
  split(f = ga) %>% 
  lapply(function(x){
    apply(x[,-1], 2, mean)
  }) %>% 
  do.call(rbind, .)


range(marker_data)


colnames(marker_data) <- 
  metabolite_tags$Compound.name[
    match(colnames(marker_data), 
          metabolite_tags$name)  
    ]

marker_data <- 
  t(marker_data)


mestimate<- function(df){
  N <-  dim(df)[[1]]
  D <- dim(df)[[2]]
  m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
  return(m.sj)
}

m <- mestimate(marker_data)
m


library(e1071)
#helper function for the within sum of squared error
sumsqr <- function(x, clusters){
  sumsqr <- function(x) sum(scale(x, scale = FALSE)^2)
  wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
  return(wss)
}

#get the wss for repeated clustering
iterate_fcm_WSS <- function(df,m){
  totss <- numeric()
  for (i in 2:20){
    FCMresults <- cmeans(df,centers=i,m=m)
    totss[i] <- sum(sumsqr(df,FCMresults$cluster))
  }
  return(totss)
}

wss_2to20 <- iterate_fcm_WSS(marker_data, m)
plot(1:20, wss_2to20[1:20], 
     type = "b",
     xlab = "Number of Clusters", 
     ylab = "wss")



text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(15),
  "red"
)

k <- 2
fcm_results <- cmeans(marker_data, 
                      centers = k, m = m)


#get the centroids into a long dataframe:
fcm_centroids <- fcm_results$centers
fcm_centroids_df <- data.frame(fcm_centroids, check.names = FALSE, check.rows = FALSE)
fcm_centroids_df$cluster <- row.names(fcm_centroids_df)

centroids_long <- 
  fcm_centroids_df %>% 
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "GA",
                      values_to = "value")
  
ga_level <- unique(centroids_long$GA)
ga_level <- 
  c("[11,13]", ga_level[-c(15,16)], "PP")

centroids_long$GA <- factor(centroids_long$GA, levels = ga_level)

ggplot(centroids_long, 
       aes(x = GA,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  ggsci::scale_color_aaas() +
  labs(x = "",
       title = "", 
       color = "Cluster", 
       y = "Scaled intensity") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1, colour = text_colour),
        axis.text.y = element_text(size = 13), 
        plot.title = element_text(size = 15, hjust = 0.5))



#start with the input data
fcm_plotting_df <- data.frame(marker_data, check.names = FALSE, check.rows = FALSE)

#add genes
fcm_plotting_df$metabolite <- row.names(fcm_plotting_df)

#bind cluster assinment
fcm_plotting_df$cluster <- fcm_results$cluster
#fetch the membership for each gene/top scoring cluster
fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster),function(row){
  clust <- fcm_plotting_df$cluster[row]
  fcm_results$membership[row,clust]
})




k_to_plot = 2

#subset the dataframe by the cluster and get it into long form
#using a little tidyr action
cluster_plot_df <- 
  dplyr::filter(fcm_plotting_df, cluster == k_to_plot) %>%
  dplyr::select(everything(),membership,metabolite) %>%
  tidyr::pivot_longer(cols = -c(metabolite, cluster, membership), 
                      names_to = "GA", values_to = "value")
  

#order the dataframe by score
cluster_plot_df <- cluster_plot_df[order(cluster_plot_df$membership),]
#set the order by setting the factors using forcats
cluster_plot_df$metabolite = forcats::fct_inorder(cluster_plot_df$metabolite)

cluster_plot_df$GA <- factor(cluster_plot_df$GA, levels = ga_level)

#subset the cores by cluster
core <- dplyr::filter(centroids_long, cluster == k_to_plot)

ggplot(cluster_plot_df, aes(x = GA, y = value)) + 
  geom_line(aes(colour = membership, group = metabolite), size = 1) +
  geom_point(aes(colour = membership, group = metabolite), size = 2) +
  scale_colour_gradientn(colours = c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core, aes(GA, value, group=cluster), color="black",inherit.aes = FALSE) +
  labs(x = "", y = "Scaled intensity", 
       title= "") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1,
                                   colour = text_colour),
        axis.text.y = element_text(size = 13), 
        plot.title = element_text(size = 15, hjust = 0.5))
  


cluster1_metabolite <-
fcm_plotting_df$metabolite[which(fcm_plotting_df$cluster == 1)]

cluster2_metabolite <-
  fcm_plotting_df$metabolite[which(fcm_plotting_df$cluster == 2)]


write.csv(fcm_plotting_df, "fuzzy_c_means_cluster.csv", row.names = FALSE)













##correlation network
##this networks should contain metabolite and clinical information
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/")
rm(list=ls())
load("data_preparation_for_analysis/metabolite_table")
load("data_preparation_for_analysis/metabolite_tags")
load("data_preparation_for_analysis/peak_table")

info <-
  readxl::read_xlsx("E:/project/smartD/patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
info <-
  info %>%
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>%
  mutate(ID = paste("SF", ID, sep = ""))

sample_info <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")

setwd(r4projects::get_project_wd())
marker1 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/GA_prediction/marker_rf_final.csv")
marker2 <- readr::read_csv("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/marker_rf_final.csv")


setwd(r4projects::get_project_wd())
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


##rename P
str_sort(grep("P", colnames(marker_data), value = T), numeric = TRUE)
sample_info %>% 
  filter(is.na(GA)) %>% 
  pull(Sample_ID) %>% 
  str_sort(numeric = TRUE)

##So P1 is X178 and P2 is X179 and so on.


colnames(marker_data)[grep("P", colnames(marker_data))] <-
  str_replace(colnames(marker_data)[grep("P", colnames(marker_data))],"P", "") %>% 
  as.numeric() %>% 
  `+`(., 177) %>% 
  paste("X", ., sep = "")

marker_data <-
  marker_data %>% 
  select(one_of(sample_info$Sample_ID))


sum(is.na(marker_data))

marker_data <- 
  apply(marker_data, 1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()


subject_id <- 
  sample_info$Patient_ID[match(colnames(marker_data), sample_info$Sample_ID)]


cor_met_data <-
  data.frame(subject_id = subject_id,
    t(marker_data),
    stringsAsFactors = FALSE)


cor_met_data <-
  cor_met_data %>% 
  plyr::dlply(.variables = 'subject_id', .fun = function(x){
    x <- x[,-1]
    x <- apply(x, 2, mean)
    x
  }) %>% 
    do.call(rbind, .)


colnames(cor_met_data) <-
metabolite_tags$Compound.name[match(colnames(cor_met_data), metabolite_tags$name)]

cor_met_data <-
  cor_met_data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "subject_id")
  

###clinical information
info <- 
  info %>% 
  filter(ID %in% cor_met_data$subject_id)

due_date <-
  (as.Date(info$DD) - (as.Date(info$EDD) - 280))/7 %>%
  as.numeric()

##age
age <-
  info$Age %>%
  as.numeric()
age
# 
##ethinic
ethinic <-
  info$`Ethinic Group`
# 
ethinic <-
  case_when(ethinic == "1" ~ "White",
            ethinic == "Caucasian" ~ "White",
            ethinic == "2" ~ "Black",
            ethinic == "3" ~ "Latina",
            ethinic == "4" ~ "Pacific Islander",
            ethinic == "5" ~ "Asian",
            ethinic == "4 (Asian)" ~ "Asian",
            ethinic == "Asian" ~ "Asian",
            ethinic == "Afr Am" ~ "Black",
            TRUE ~ ethinic
  )
# 
# ##BMI
ht <- trans_ht(info$Ht)/100
wt <- trans_wt(info$Wt)
bmi <- wt/(ht^2)
# 
# #------------------------------------------------------------------------------
# ##parity
parity <- info$Parity
parity <-
  parity %>%
  stringr::str_extract("G[0-9]{1}") %>%
  stringr::str_replace("G", "") %>%
  as.numeric()
# 
# ##sex and twins
# info$Sex
sex <- info$Sex
sex <-
  case_when(is.na(sex) ~ "NA",
            sex == "F, F" ~ 'F_F',
            sex == "M,M" ~ 'M_M',
            sex == "M,M" ~ 'M_M',
            sex == "M / F" ~ 'M_F',
            TRUE ~ sex
  )
# ##IVF
ivf <- rep(NA, nrow(info))
ivf[grep("ivf", stringr::str_to_lower(info$`Pregnancy dx`))] <- "YES"
ivf[grep("transfer", stringr::str_to_lower(info$Dating))] <- "YES"
ivf[is.na(ivf)] <- "NO"
# 
##induction
induction <-
  info$Induction
# 
induction <-
  case_when(is.na(induction) ~ "NA",
            induction == "Y" ~ "YES",
            induction == "Yes" ~ "YES",
            induction == "N" ~ "NO",
  )
# 
# 
birth_wt <-
  info$`Birth wt`
birth_wt[19] <- "3015;3170"
birth_wt <-
  sapply(birth_wt, function(x){
    if(!is.na(x)){
      x <- stringr::str_split(x, ";")[[1]] %>%
        as.numeric() %>%
        sum()
    }else{
      x
    }
    x
  })
# 
birth_wt <- unname(birth_wt) %>%
  as.numeric()
# 
clinical_data <- rbind(
  # due_date,
  age,
  # ethinic,
  bmi,
  parity,
  # sex,
  # ivf,
  # induction,
  birth_wt
)
# 
colnames(clinical_data) <-
  info$ID


clinical_data <-
  t(clinical_data) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "subject_id")


clinical_data2 <- 
  clinical_data %>% 
  filter(!is.na(birth_wt))

clinical_data2$bmi <-
  (clinical_data2$bmi - mean(clinical_data2$bmi))/sd(clinical_data2$bmi)

clinical_data2$birth_wt <-
  (clinical_data2$birth_wt - mean(clinical_data2$birth_wt))/sd(clinical_data2$birth_wt)


cor_data <-
  clinical_data2 %>% 
  left_join(cor_met_data, by = "subject_id")



cross_cor <- 
  cor_data[,-1] %>% 
  cor(., method = "spearman")


rownames(cross_cor) == colnames(cross_cor)


cross_cor[lower.tri(cross_cor)] <- 2

cross_cor <- 
  as_tibble(cross_cor)
  

rownames(cross_cor) <- 
  colnames(cross_cor)

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
    p <-  cor.test(int1, int2, alternative = "two.sided",
                   method = "spearman",  use = "na.or.complete")$p.value
    p
  })

sum(cross_cor_p < 0.05)



###p adjustment
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
  left_join(metabolite_tags, by = c("name" = "Compound.name"))

node_attr$node_class <- rep(NA, nrow(node_attr))

node_attr$node_class <- 
  case_when(is.na(node_attr$Total.score) ~ "Clinical information",
            TRUE ~ 'Metabolite')


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

node_attr$node_class <- factor(node_attr$node_class, 
                               levels = unique(node_attr$node_class))

cross_graph <-
  tidygraph::tbl_graph(nodes = node_attr, 
                       edges = cross_cor2, 
                       directed = FALSE)





##angle for label
metabolite <- igraph::V(cross_graph)$name
id = 1:length(metabolite)
angle <- 360 * (id - 0.5)/length(metabolite)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

set_graph_style(family = "Arial")
extrafont::loadfonts()
plot <- 
  ggraph(cross_graph, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(edge_colour = Correlation, 
                    edge_width = -log(P.adjust, 10))) +
  scale_edge_colour_gradient2(low = "royalblue", 
                              mid = "white",
                              high = "red") +
  scale_edge_width_continuous(range = c(0.8, 3.5)) +
  geom_node_point(aes(size = Total.score, colour = node_class)) +
  scale_colour_manual(values = c(    "Clinical information" = "#C71000FF",
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
                                     "Unknown" = "#3F4041B2")) +
  # ggsci::scale_color_futurama(alpha = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous(range = c(3,10),
                        guide = guide_legend(override.aes = list(colour = "black"))) +
  geom_node_text(aes(x = x * 1.1,
                     y = y * 1.1,
                     label = name, 
                     colour = node_class), 
                 angle = angle, 
                 hjust = hjust,
                 # colour = "black",
                 size = 3.5) +
  # ggdark::dark_theme_void() +
  theme_void() +
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))

plot


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

##pathway enrichment analysis to get the insight of which pathways are altered in pregnancy
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/biological_analysis/pathway_enrichment")
marker$Compound.name
marker$KEGG.ID
marker$HMDB.ID

load("hsa.kegg.pathway.rda")

path_result <- 
enrichPathway(id = marker$KEGG.ID[!is.na(marker$KEGG.ID)], 
              database = hsa.kegg.pathway, 
              method = "hypergeometric")




path_result %>% 
  mutate(class = case_when(
    p.value < 0.05 ~ "Significant",
    TRUE ~ "No"
  )) %>% 
  mutate(class = factor(class, levels = c("Significant", "No"))) %>% 
  ggplot(aes(x = Overlap, y = -log(p.value, 10) )) +
  geom_hline(yintercept = 1.3, linetype = 2) +
  geom_point(aes(size = Pathway.length, colour = class)) +
  scale_colour_manual(values = c("Significant" = "#C71000FF", 
                                 "No" = "#3F4041FF")) +
  ggplot2::guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  ggrepel::geom_label_repel(mapping = aes(x = Overlap,
                                          y = -log(p.value, 10), 
                                          label = Pathway.name)) +
    theme(axis.title = element_text(size = 15),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13))





path_result %>% 
  mutate(class = case_when(
    p.value < 0.05 ~ "Significant",
    TRUE ~ "No"
  )) %>% 
  mutate(class = factor(class, levels = c("Significant", "No"))) %>% 
  arrange(desc(p.value)) %>% 
  mutate(Pathway.name = factor(Pathway.name, levels = Pathway.name)) %>% 
  ggplot(aes(x = Pathway.name, y = -log(p.value, 10) )) +
  # geom_hline(yintercept = 1.3, linetype = 2) +
  geom_bar(aes(fill = class), stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("Significant" = "#C71000FF", 
                                 "No" = "#3F4041FF")) +
  labs(x = "", y = "-log10(P-value)") +
  # ggplot2::guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  # ggrepel::geom_label_repel(mapping = aes(x = Overlap,
  #                                         y = -log(p.value, 10), 
  #                                         label = Pathway.name)) +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  coord_flip()






