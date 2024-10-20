##avoid source 
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
source("R/20200727/tools.R")

##load data
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info")

setwd("data_analysis20200108/urine_metabolome/DEG_analysis/")
dim(expression_data)

rownames(expression_data)
variable_info$name

rownames(expression_data) == variable_info$name
colnames(expression_data) == sample_info$sample_id

###find marker which are change according to pregnancy
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

subject_data[5,] %>% 
  as.numeric() %>% 
  density() %>%
  plot()

dim(subject_data)
dim(variable_info)

rownames(subject_data) <- variable_info$name

sample_info$GA[is.na(sample_info$GA)] <- 45

###remove samples which are after birth
remove_idx <- which(sample_info$GA == 45)
# remove_name <- sample_info$sample_id[remove_idx]

##SAM analysi
# sam_test <-
# samr::SAM(x = subject_data[,-remove_idx],
#           y = (sample_info$GA)[-remove_idx],
#           resp.type = "Quantitative",
#           geneid = rownames(subject_data),
#           genenames = rownames(subject_data),
#           return.x = FALSE,
#           fdr.output = 0.05,
#           regression.method = "ranks",
#           random.seed = "123",
#           nperms = 1000)
# 
# save(sam_test, file = "sam_test")

load("sam_test")

gene_up <- sam_test$siggenes.table$genes.up
gene_down <- sam_test$siggenes.table$genes.lo

peak_marker <- c(gene_up[,1], gene_down[,1])

plot(sample_info$GA, subject_data[peak_marker[1],])

plot(sample_info$GA, expression_data[peak_marker[1],])

plot(sample_info$GA, expression_data[gene_down[4,1],])

###linear regression
subject_data <-
  subject_data[peak_marker,]

rownames(subject_data) <- peak_marker

colnames(subject_data) == sample_info$sample_id

ga <- sample_info$GA
batch <- sample_info$batch
bmi <- sample_info$bmi
age <- sample_info$Age
parity <- sample_info$parity
ethnicity <- sample_info$ethnicity

remove_idx <- which(ga == 45)

# lm_p_value <-
# pbapply::pbapply(subject_data[,-remove_idx], 1, function(x){
#   x <- as.numeric(x)
#   lm_reg <- lm(formula = ga[-remove_idx] ~ x + batch[-remove_idx] +
#                  bmi[-remove_idx] + age[-remove_idx] + parity[-remove_idx] +
#                  ethnicity[-remove_idx])
#   temp <- summary(lm_reg)$coefficients  %>% as.data.frame()
#   as.numeric(temp["x",4])
# })
# 
# 
# save(lm_p_value, file = 'lm_p_value')
load("lm_p_value")
lm_fdr <- p.adjust(lm_p_value, method = "fdr")

lm_fdr <- 
lm_fdr %>% 
  data.frame()

colnames(lm_fdr) <- "fdr"

lm_fdr <- 
  lm_fdr %>% 
  tibble::rownames_to_column(var = "name")

###volcano plot
sam_value <- rbind(gene_up, gene_down) %>% 
  as.data.frame() %>% 
  dplyr::select(`Gene ID`, `Score(d)`) %>% 
  dplyr::rename(name = `Gene ID`, score = `Score(d)`) %>% 
  dplyr::left_join(lm_fdr, by = "name")

plot <- 
volcano_plot(fc = as.numeric(sam_value$score),
             log2_fc = FALSE,
             p_value = sam_value$fdr, 
             p.cutoff = 0.05, 
             fc.cutoff = 1,
             size_range = c(0.1,5)) 

sum(lm_fdr < 0.05 & as.numeric(sam_value$score) > 0)
sum(lm_fdr < 0.05 & as.numeric(sam_value$score) < 0)


plot <- 
plot +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw()

plot
ggsave(plot, filename = "sam_lm_regression.pdf", width = 7, height = 7)

peak_name <- 
sam_value$name[which(sam_value$fdr < 0.05)]


temp_subject_data <- 
  subject_data[peak_name,] %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

###k-mean
library(CancerSubtypes)
result <-
  ExecuteCC(
    clusterNum = 2,
    d = as.matrix(temp_subject_data),
    maxK = 6,
    reps = 1000,
    pItem = 0.8,
    pFeature = 0.8,
    title = "k_means_consensus",
    clusterAlg = "km",
    distance = "euclidean",
    plot = "png",
    writeTable = TRUE
  )

save(result, file = "result")
load("result")

idx <- 3
sil=silhouette_SimilarityMatrix(result$originalResult[[idx]]$consensusClass, 
                                result$originalResult[[idx]]$consensusMatrix)
# sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)

sil_plot <- 
  plot_silhouette(sil)

sil_plot

ggsave(sil_plot, file = "k_means_consensus/sil_plot3.pdf", width = 7, height = 7)

plot(result$originalResult[[idx]]$consensusTree)

name2 <- colnames(temp_subject_data)[result$originalResult[[idx]]$consensusTree$order]

temp_subject_data2 <- temp_subject_data[,name2]
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
  temp_sample_info2[match(names(cluster), temp_sample_info2$sample_id),]

temp_subject_data2 <-
  temp_subject_data2[,names(cluster)]

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


text_colour <- c(
  colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(16),
  "red"
)

names(text_colour) <- ga_level

cluster_color <- ggsci::pal_d3()(10)[4:6]
names(cluster_color) <- c(1,2,3)

ha1 = HeatmapAnnotation(
  ga = factor(ga, levels = ga_level),
  cluster = factor(cluster, levels = as.character(c(1, 2, 3))),
  col = list(
    ga = text_colour,
    cluster = cluster_color
  ),
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
    axis_param = list(side = "left"
                      # at = c(0, 0.5, 1),
                      # labels = c("zero", "half", "one")
                      )
    ),
    annotation_name_side = c("left")
  )

temp_data <- temp_subject_data2
range((temp_data))

temp_data[temp_data > 3] <- 3
temp_data[temp_data < -3] <- -3

plot <-     
  Heatmap(temp_data, 
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
          bottom_annotation = ha2)

library(ggplotify)
plot <- as.ggplot(plot)
plot
ggsave(plot, filename = "consensus_heatmap.pdf", width = 12, height = 7)


###test for two different clusters
index1 <- which(cluster == 1)
index2 <- which(cluster == 2)

data.frame(cluster = factor(cluster, levels = c(1,2)),
           age, stringsAsFactors = FALSE) %>% 
  ggplot(aes(cluster, age)) +
  geom_boxplot() +
  geom_jitter()


t.test(age[index1], age[index2])

data.frame(cluster = factor(cluster, levels = c(1,2)),
           bmi, stringsAsFactors = FALSE) %>% 
  ggplot(aes(cluster, bmi)) +
  geom_boxplot() +
  geom_jitter()

wilcox.test(bmi[index1], bmi[index2])

table(cluster, ethnicity2) %>% 
  chisq.test()
