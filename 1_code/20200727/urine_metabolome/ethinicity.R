##avoid source 
no_function()

###we want to know which clinical information make samples are different
##-----------------------------------------------------------------------------
##RPLC pos and neg
###----------------------------------------------------------------------------
sxtTools::setwd_project()
rm(list=ls())
source("R/20200727/tools.R")

##peaks
load(
  "data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data"
)

load(
  "data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info"
)

load(
  "data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info"
)

setwd("data_analysis20200108/urine_metabolome/ethnicity/")

library(metflow2)
library(tidyverse)

subject_data <- 
  expression_data[,sample_info$class == "Subject"]

qc_data <- 
  expression_data[,sample_info$class == "QC"]

qc_data <- 
  qc_data %>% 
  select(-c(QC2.1, QC2.2, QC2.3))

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)

subject_data <- subject_data[remain_idx,]
qc_data <- qc_data[remain_idx,]

###log
subject_data <- 
  log(subject_data + 1, 2)

qc_data <- 
  log(qc_data + 1, 2)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )

qc_data <- 
  t(
    apply(qc_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  tibble::as_tibble()

ga <- sample_info$GA[match(colnames(subject_data), sample_info$sample_id)]
ga[is.na(ga)] <- 0
batch <- sample_info$batch[match(colnames(subject_data), sample_info$sample_id)]
subject_id <- sample_info$subject_id[match(colnames(subject_data), sample_info$sample_id)]
age <- sample_info$Age[match(colnames(subject_data), sample_info$sample_id)]
ethnicity <- sample_info$ethnicity[match(colnames(subject_data), sample_info$sample_id)]
parity <- sample_info$parity[match(colnames(subject_data), sample_info$sample_id)]
bmi <- sample_info$bmi[match(colnames(subject_data), sample_info$sample_id)]

subject_data2 <- 
  t(subject_data) %>% 
  tibble::as_tibble()

rownames(subject_data2)

qc_data2 <- 
  t(qc_data) %>% 
  tibble::as_tibble()

rownames(subject_data2) <- colnames(subject_data)


# ####PCA without QC only for subjects
# temp_data <- 
#   subject_data2
# 
# pca_object <- 
#   prcomp(x = temp_data)
# 
# library(ggfortify)
# 
# x <- pca_object$x
# 
# x <- x[,1:2]
# 
# subject_id <- sample_info$subject_id[match(colnames(subject_data), sample_info$sample_id)]
# age <- sample_info$Age[match(colnames(subject_data), sample_info$sample_id)]
# ethnicity <- sample_info$ethnicity[match(colnames(subject_data), sample_info$sample_id)]
# parity <- sample_info$parity[match(colnames(subject_data), sample_info$sample_id)]
# bmi <- sample_info$bmi[match(colnames(subject_data), sample_info$sample_id)]
# 
# 
# x <- data.frame(x, 
#                 GA = ga,
#                 subject_id = subject_id,
#                 age,
#                 ethnicity,
#                 parity,
#                 bmi,
#                 stringsAsFactors = FALSE)
# 
# plot <-
#   ggplot(x[x$GA!=0,], aes(PC1, PC2, fill = GA)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   )) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,1), legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
#            colour = "#C71000FF", size = 5) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))
# 
# 
# plot
# 
# ggsave(plot, filename = "pca_ga_plot.pdf",
#        width = 7, height = 7)
# 
# 
# ##pca plot subject_id
# color <-
#   colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(subject_id)))
# 
# plot <- 
#   ggplot(x, aes(PC1, PC2, fill = subject_id)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   scale_fill_manual(values = color) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = "right",
#         legend.background = element_blank()) +
#   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
#   #          colour = "#C71000FF", size = 5) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))
# 
# plot
# 
# ggsave(plot, filename = "pca_plot_subject.pdf", 
#        width = 9, height = 7)
# 
# 
# ##pca plot age
# plot <-
#   ggplot(x, aes(PC1, PC2, fill = age)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   guides(fill = guide_colourbar(title = "Age")) +
#   scale_fill_gradient(low = "blue", high = "red") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,1),
#         legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
#   #          colour = "#C71000FF", size = 5) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))
# 
# plot
# 
# ggsave(plot, filename = "pca_plot_age.pdf", 
#        width = 7, height = 7)
# 
# ##pca plot ethnicity
# color <- 
#   ggsci::pal_futurama()(7)[1:7]
# 
# names(color) <- sort(unique(ethnicity))
# 
# plot <- 
#   ggplot(x, aes(PC1, PC2, fill = ethnicity)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   scale_fill_manual(values = color) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(0,1),
#         legend.justification = c(0,1),
#         legend.background = element_blank()) +
#   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
#   #          colour = "#C71000FF", size = 5) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))
# 
# plot
# 
# ggsave(plot, filename = "pca_plot_ethnicity.pdf", 
#        width = 7, height = 7)
# 
# ##pca plot parity
# plot <- 
#   ggplot(x, aes(PC1, PC2, fill = parity)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   scale_fill_gradient(low = alpha("red", 0.3), 
#                       high = alpha("red", 1)) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,1),
#         legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
#   #          colour = "#C71000FF", size = 5) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))
# 
# plot
# 
# ggsave(plot, filename = "pca_plot_parity.pdf", 
#        width = 7, height = 7)
# 
# ##pca plot bmi
# plot <- 
#   ggplot(x, aes(PC1, PC2, fill = bmi)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   scale_fill_gradient(low = alpha("blue", 0.3), high = "blue") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,1),
#         legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
#   #          colour = "#C71000FF", size = 5) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))
# 
# plot
# 
# ggsave(plot, filename = "pca_plot_bmi.pdf", 
#        width = 7, height = 7)
# 
# 
# x2 <-
#   x %>%
#   dplyr::filter( GA == 0)
# 
# ###GA for each ethnicity
# plot <-
#   ggplot(x[x$GA!=0,], aes(PC1, PC2)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 3, shape = 21, 
#              aes(fill = GA), 
#              show.legend = FALSE) +
#   geom_point(mapping = aes(x = PC1, y = PC2),
#              data = x2, fill = "#C71000FF", 
#              shape = 21, size = 3) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   )) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         # legend.position = c(1,1), legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = "")) +
#   facet_wrap(facets = ~ethnicity, scales = "free", ncol = 3)
# 
# 
# plot
# 
# ggsave(plot, filename = "pca_ethnicity_ga_plot.pdf",
#        width = 7, height = 7)
# 
# 
# ##tsne analysis
# tsne_object <- Rtsne::Rtsne(
#   X = as.matrix(temp_data),
#   dims = 2,
#   perplexity = 30,
#   verbose = TRUE
# )
# 
# Y <- tsne_object$Y
# Y <-
#   data.frame(Y, 
#              ga,
#              subject_id,
#              age,
#              ethnicity,
#              parity,
#              bmi,
#              stringsAsFactors = FALSE)
# 
# ####t-sne for ga
# plot <- 
#   ggplot(Y[Y$ga!=0,], aes(X1, X2, fill = ga)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 5, shape = 21) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   )) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,1), legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   annotate(geom = "point", x = Y$X1[Y$ga==0], y = Y$X2[Y$ga == 0], 
#            fill = "#C71000FF", size = 5, shape = 21) +
#   labs(x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2")
# 
# plot
# 
# ggsave(plot, filename = "tsne_plot_ga.pdf", 
#        width = 7, height = 7)
# 
# ####t-sne for subject_id
# color <- 
#   colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(Y$subject_id)))
# 
# names(color) <- unique(Y$subject_id)
# 
# plot <- Y %>% 
#   ggplot(aes(X1, X2)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(aes(fill = subject_id), size = 5,
#              shape = 21, color = "black", alpha = 0.8) +
#   scale_fill_manual(values =color) +
#   labs(x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2") +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12),
#     strip.background = element_rect(fill = "#0099B47F"),
#     strip.text = element_text(color = "white", size = 13)
#   )
# 
# plot
# 
# ggsave(plot, filename = "tsne_plot_subject_id.pdf", width = 9, height = 7)
# 
# 
# 
# ####t-sne for age
# color <- 
#   colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(Y$subject_id)))
# 
# names(color) <- unique(Y$subject_id)
# 
# plot <- Y %>% 
#   ggplot(aes(X1, X2)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(aes(fill = age), size = 5,
#              shape = 21, color = "black", alpha = 0.8) +
#   scale_fill_gradient(low = "blue", high = "red") +
#   labs(x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2") +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12),
#     strip.background = element_rect(fill = "#0099B47F"),
#     strip.text = element_text(color = "white", size = 13)
#   )
# 
# plot
# 
# ggsave(plot, filename = "tsne_plot_age.pdf", width = 9, height = 7)
# 
# 
# ####t-sne for ethnicity
# color <- 
#   ggsci::pal_futurama()(7)[1:7]
# 
# names(color) <- sort(unique(ethnicity))
# 
# 
# library(ggfortify)
# 
# plot <- Y %>% 
#   ggplot(aes(X1, X2)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(aes(fill = ethnicity), size = 5,
#              shape = 21, color = "black", alpha = 0.8) +
#   scale_fill_manual(values = color) +
#   labs(x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2") +
#   theme_bw() +
#   theme(
#     legend.position = c(0,1),
#     legend.justification = c(0,1),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 12),
#     strip.background = element_rect(fill = "#0099B47F"),
#     strip.text = element_text(color = "white", size = 13)
#   )
# 
# plot
# 
# ggsave(plot, filename = "tsne_plot_ethnicity.pdf", width = 7, height = 7)
# 
# 
# 
# 
# 
# 
# ###GA for each ethnicity
# Y2 <-
#   Y %>%
#   dplyr::filter(ga == 0)
# 
# plot <-
#   ggplot(Y[Y$ga!=0,], aes(X1, X2)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(size = 3, shape = 21, 
#              aes(fill = ga), 
#              show.legend = FALSE) +
#   geom_point(mapping = aes(x = X1, y = X2),
#              data = Y2, fill = "#C71000FF", 
#              shape = 21, size = 3) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   )) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         # legend.position = c(1,1), legend.justification = c(1,1),
#         legend.background = element_blank()) +
#   labs(x = "t-SNE Dimension 1",
#        y = "t-SNE Dimension 2") +
#   facet_wrap(facets = ~ethnicity, scales = "free", ncol = 3)
# 
# 
# plot
# 
# ggsave(plot, filename = "pca_ethnicity_ga_plot.pdf",
#        width = 7, height = 7)


###consensus analysis
#to avoind source
##load data

#here we use the K-means consensus clustering
library(CancerSubtypes)

###remove QC cv > 0.3
temp_sample_info <- 
  sample_info %>% 
  dplyr::filter(!stringr::str_detect(sample_id, "QC"))

temp_subject_data <- 
  expression_data %>% 
  dplyr::select(temp_sample_info$sample_id)

temp_qc_data <- 
  expression_data %>% 
  dplyr::select(-temp_sample_info$sample_id)

qc_cv <-
  temp_qc_data %>% 
  apply(1, function(x){
    sd(x)/mean(x)
  })

# subject_cv <-
#   temp_subject_data %>% 
#   apply(1, function(x){
#     sd(x)/mean(x)
#   })

remain_idx <- 
  which(qc_cv < 0.3)

temp_subject_data <- 
  temp_subject_data[remain_idx,]

temp_qc_data <-
  temp_qc_data[remain_idx,]

##remain subject CV top 100
subject_cv <-
  temp_subject_data %>%
  apply(1, function(x){
    sd(x)/mean(x)
  })

remain_idx <- 
  which(subject_cv > quantile(subject_cv, probs = c(0.5)))

temp_subject_data <- temp_subject_data[remain_idx,]
temp_qc_data <- temp_qc_data[remain_idx,]


##log and scale
temp_subject_data <- 
  log(temp_subject_data + 1, 2)

temp_subject_data <- 
  temp_subject_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

# ###k-mean
# result <-
#   ExecuteCC(
#     clusterNum = 2,
#     d = as.matrix(temp_subject_data),
#     maxK = 6,
#     reps = 1000,
#     pItem = 0.8,
#     pFeature = 0.8,
#     title = "k_means_consensus",
#     clusterAlg = "km",
#     distance = "euclidean",
#     plot = "png",
#     writeTable = TRUE
#   )
# 
# save(result, file = "result")
load("result")

idx <- 2
sil=silhouette_SimilarityMatrix(result$originalResult[[idx]]$consensusClass, 
                                result$originalResult[[idx]]$consensusMatrix)
# sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)

sil_plot <- 
  plot_silhouette(sil)

sil_plot

ggsave(sil_plot, file = "k_means_consensus/sil_plot2.pdf", width = 7, height = 7)

plot(result$originalResult[[idx]]$consensusTree)

name2 <- colnames(temp_subject_data)[result$originalResult[[idx]]$consensusTree$order]

temp_subject_data2 <- temp_subject_data[,name2]
cluster <- result$originalResult[[idx]]$consensusClass

temp_sample_info2 <-
  sample_info[match(colnames(temp_subject_data2), sample_info$sample_id), ]

cluster <- cluster[match(temp_sample_info2$sample_id, names(cluster))]

names(cluster) == colnames(temp_subject_data2)

temp_sample_info2$GA[is.na(temp_sample_info2$GA)] <- 50

###complext heatamp
temp_data <- temp_subject_data2

ga <- 
  temp_sample_info2 %>% 
  dplyr::select(sample_id, GA) %>% 
  dplyr::mutate(ga = ggplot2::cut_width(GA, width = 2)) %>% 
  dplyr::mutate(ga = stringr::str_replace(ga, "\\[", "(")) %>% 
  pull(ga)

ga_level <- unique(ga) %>% stringr::str_sort(numeric = TRUE)

##other information
age <- 
  temp_sample_info2$Age

bmi <- 
  temp_sample_info2$bmi

sum(is.na(bmi))

parity <- 
  temp_sample_info2$parity

sum(is.na(parity))

ethnicity <- 
  temp_sample_info2$ethnicity

sum(is.na(ethnicity))

library(ComplexHeatmap)

library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("#4292C6", "white", "red"))
col_fun(seq(-3, 3))
cluster_col_fun = colorRamp2(c(0, 2, 3), c("green", "white", "red")) 

ethnicity_color <-
  ggsci::pal_futurama()(7)[1:7]

names(ethnicity_color) <- sort(unique(ethnicity))

ethnicity2 <- ethnicity

ethnicity2[ethnicity2 == "Caucasian" | ethnicity2 == "White" |
             ethnicity2 == "Asian" | ethnicity2 == "Other" |
             ethnicity2 == "Pacific Islander"] <- "Asian&White"

ethnicity2[ethnicity2 == "Latina" | ethnicity2 == "Black"] <- 
  "Black&Latina"

ethnicity_color2 <-
  ggsci::pal_futurama()(7)[c(7,2)]

names(ethnicity_color2) <- sort(unique(ethnicity2))

ha1 = HeatmapAnnotation(
  age = age,
  bmi = bmi,
  ethnicity = ethnicity,
  ethnicity2 = ethnicity2,
  parity = parity,
  cluster = factor(cluster, levels = as.character(c(1, 2))),
  #   scale_fill_gradient(low = alpha("red", 0.3), 
  #                       high = alpha("red", 1)) +
  col = list(
    ethnicity = ethnicity_color,
    ethnicity2 = ethnicity_color2,
    age = circlize::colorRamp2(
      breaks = c(min(age), max(age)),
      colors = c("blue", "red")),
    parity = circlize::colorRamp2(
      breaks = c(min(parity), max(parity)),
      colors = c("lightpink", alpha("red", 1))),
    bmi = circlize::colorRamp2(
      breaks = c(min(bmi), max(bmi)),
      colors = c("lightblue", alpha("blue", 1)))
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
          top_annotation = ha1)

library(ggplotify)
plot <- as.ggplot(plot)
plot
# ggsave(plot, filename = "heatmap.pdf", width = 12, height = 7)
# ggsave(plot, filename = "heatmap.png", width = 12, height = 7)

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





###find the different peaks in different time points
##avoid source 
no_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)
source("R/20200727/tools.R")

##load data
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info")

setwd("data_analysis20200108/urine_metabolome/ethnicity/")

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

# ###remove samples which are after birth
# remove_idx <- which(sample_info$GA == 45)
# # remove_name <- sample_info$sample_id[remove_idx]
# 
# #SAM analysi
# # sam_test <-
# # samr::SAM(x = subject_data[,-remove_idx],
# #           y = (sample_info$GA)[-remove_idx],
# #           resp.type = "Quantitative",
# #           geneid = rownames(subject_data),
# #           genenames = rownames(subject_data),
# #           return.x = FALSE,
# #           fdr.output = 0.05,
# #           regression.method = "ranks",
# #           random.seed = "123",
# #           nperms = 1000)
# # 
# # save(sam_test, file = "sam_test")
# 
# load("sam_test")
# 
# gene_up <- sam_test$siggenes.table$genes.up
# gene_down <- sam_test$siggenes.table$genes.lo
# 
# peak_marker <- c(gene_up[,1], gene_down[,1])
# 
# plot(sample_info$GA, subject_data[peak_marker[1],])
# 
# plot(sample_info$GA, expression_data[peak_marker[1],])
# 
# plot(sample_info$GA, expression_data[gene_down[4,1],])
# 
# ###linear regression
# subject_data <- 
#   subject_data[peak_marker,]
# 
# rownames(subject_data) <- peak_marker
# 
# colnames(subject_data) == sample_info$sample_id
# 
# ga <- sample_info$GA
# batch <- sample_info$batch
# bmi <- sample_info$bmi
# age <- sample_info$Age
# parity <- sample_info$parity
# ethnicity <- sample_info$ethnicity
# 
# remove_idx <- which(ga == 45)
# 
# lm_p_value <- 
# apply(subject_data[,-remove_idx], 1, function(x){
#   x <- as.numeric(x)
#   lm_reg <- lm(formula = ga[-remove_idx] ~ x + batch[-remove_idx] + 
#                  bmi[-remove_idx] + age[-remove_idx] + parity[-remove_idx] +
#                  ethnicity[-remove_idx])
#   temp <- summary(lm_reg)$coefficients  %>% as.data.frame()
#   as.numeric(temp["x",4])
# })
# 
# lm_fdr <- p.adjust(lm_p_value, method = "fdr")
# 
# ###correlation
# 
# # cor_value <-
# #   purrr::map(.x = as.data.frame(t(subject_data[,-remove_idx])), .f = function(x){
# #     temp1 <- cor.test(x, sample_info$GA[-remove_idx], method = "spearman")
# #     temp2 <- cor.test(sxt_rank(x), sample_info$GA[-remove_idx], method = "spearman")
# #     c(temp1$estimate, temp1$p.value, temp2$estimate, temp2$p.value)
# #   })
# # 
# # cor_value <-
# #   do.call(rbind, cor_value)
# # 
# # colnames(cor_value) <- c("correlation1", "p1", "correlation2", "p2")
# # 
# # cor_value <- data.frame(name = rownames(subject_data),
# #                         cor_value, stringsAsFactors = FALSE)
# # 
# # 
# # cor_value$p1_adj <- p.adjust(p = cor_value$p1, method = "fdr")
# # 
# # cor_value$p2_adj <- p.adjust(p = cor_value$p2, method = "fdr")
# # 
# # plot(cor_value$correlation1)
# # 
# # plot(cor_value$correlation2)
# # 
# # save(cor_value, file = "cor_value")
# 
# load("cor_value")
# 
# peak_marker <-
#   rbind(gene_up, gene_down) %>%
#   as.data.frame() %>%
#   dplyr::mutate(score = as.numeric(`Score(d)`)) %>%
#   dplyr::left_join(cor_value, by = c("Gene ID" = "name")) %>%
#   # dplyr::arrange(score) %>%
#   dplyr::mutate(index = 1:(nrow(gene_up) + nrow(gene_down))) %>%
#   dplyr::mutate(class = case_when(
#     score > 0 & correlation1 > 0 ~ "up",
#     score < 0 & correlation1 < 0 ~ "down",
#     TRUE ~ "no"
#   ))
# 
# 
# peak_marker$lm_fdr <- lm_fdr
# 
# save(peak_marker, file = "peak_marker")
# load("peak_marker")
# 
# 
# plot <- peak_marker %>% 
#   ggplot(aes(score, correlation1)) +
#   geom_hline(yintercept = 0, linetype = 2, color = "black") +
#   geom_vline(xintercept = 0, linetype = 2, color = "black") +
#   geom_point(aes(color = class, size = -log(p1_adj, 10)), 
#              show.legend = TRUE, alpha = 0.8) +
#   labs(x = "Score (SAM test)", y = "Correlation (Spearman)") +
#   scale_color_manual(values = c("up" = ggsci::pal_aaas()(10)[2],
#                                 "down" = ggsci::pal_aaas()(10)[1],
#                                 "no" = "#D9D9D9")) +
#   scale_size_continuous(range = c(0.05, 5)) +
#   guides(size = guide_legend(override.aes = list(color = "black"))) +
#   # scale_colour_brewer() +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         legend.position = c(1,0), 
#         legend.justification = c(1,0),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 13)) 
# 
# plot
# 
# ggsave(plot, filename = "dem_plot_light.pdf",
#        width = 7, height = 7, bg = "transparent")
# 
# ggsave(plot, filename = "dem_plot_light.png",
#        width = 7, height = 7, bg = "transparent")
# 
# table(peak_marker$class)
# 
# variable_info <- 
#   variable_info[match(rownames(subject_data), variable_info$name),]
# 
# 
# ###remove lm_fdr > 0.05
# remain_idx <- which(peak_marker$lm_fdr < 0.05)
# 
# peak_marker <- peak_marker[remain_idx,]
# variable_info <- variable_info[remain_idx,]
# subject_data <- subject_data[remain_idx,]

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

sample_info$ga_range[sample_info$ga_range == "(10,12]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(12,14]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(14,16]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(16,18]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(18,20]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(20,22]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(22,24]"] <- "(10,26]"
sample_info$ga_range[sample_info$ga_range == "(24,26]"] <- "(10,26]"

sample_info$ga_range[sample_info$ga_range == "(26,28]"] <- "(26,34]"
sample_info$ga_range[sample_info$ga_range == "(28,30]"] <- "(26,34]"
sample_info$ga_range[sample_info$ga_range == "(30,32]"] <- "(26,34]"
sample_info$ga_range[sample_info$ga_range == "(32,34]"] <- "(26,34]"

sample_info$ga_range[sample_info$ga_range == "(34,36]"] <- "(34,42]"
sample_info$ga_range[sample_info$ga_range == "(36,38]"] <- "(34,42]"
sample_info$ga_range[sample_info$ga_range == "(38,40]"] <- "(34,42]"
sample_info$ga_range[sample_info$ga_range == "(40,42]"] <- "(34,42]"

table(sample_info$ga_range)

###combine different samples in one ga range together
library(plyr)

ethnicity <- sample_info$ethnicity
ethnicity[ethnicity == "Black" | ethnicity == "Latina"] <- "Black&Latina"
ethnicity[ethnicity != "Black&Latina"] <- "Asian&White"

names(ethnicity) <- sample_info$sample_id


##total different peaks
black_idx <- which(ethnicity == "Black&Latina")
asian_idx <- which(ethnicity == "Asian&White")

# p_value <- 
# apply(subject_data, 1, function(x){
#   x <- as.numeric(x)
#   wilcox.test(x[asian_idx], x[black_idx])$p.value
# })
# 
# sum(p_value < 0.05)
# fdr <- p.adjust(p_value, method = "fdr")
# 
# fold_change <- apply(subject_data, 1, function(x){
#   x <- as.numeric(x)
#   mean(x[black_idx]) / mean(x[asian_idx])
# })
# 
# p_fc <- data.frame(fc = fold_change, fdr = fdr, stringsAsFactors = FALSE) %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "name") %>% 
#   dplyr::mutate(class = 
#                   case_when(
#                     fc > 1 & fdr < 0.05 ~ "Increase",
#                     fc < 1 & fdr < 0.05 ~ "Decrease",
#                     TRUE ~ "No"
#                   ))
# 
# save(p_fc, file = "p_fc")
load("p_fc")

plot <- 
volcano_plot(fc = p_fc$fc, p_value = p_fc$fdr, p.cutoff = 0.01, 
             fc.cutoff = 1.5)

plot <- 
plot +
  theme_bw() +
  geom_vline(xintercept = c(-0.58,0.58), linetype = 2) +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

plot
# ggsave(plot, filename = "volcano.pdf", width = 7, height = 7)




# remain_idx <- which(p_fc$class != "No")
# 
# temp_data <- 
#   subject_data[remain_idx,]

##find the different peaks in each time point
ga_range_level <- sort(unique(sample_info$ga_range))

# ga_range_level <- ga_range_level[-4]

p_fc_each_point <- 
purrr::map(ga_range_level[-4], function(x){
  temp_idx <- 
    which(sample_info$ga_range == x)
  
  temp_subject_data <- 
  subject_data[,temp_idx]
  
  temp_ethnicity <- ethnicity[match(colnames(temp_subject_data), names(ethnicity))]
  
  temp_asian_idx <- which(
    temp_ethnicity == "Asian&White"
  )  
  
  temp_black_idx <- which(
    temp_ethnicity == "Black&Latina"
  )  
  
  temp_p_fc <- 
  apply(temp_subject_data, 1, function(x){
    x <- as.numeric(x)
    temp_p <- wilcox.test(x[temp_asian_idx], 
                x[temp_black_idx])$p.value
    temp_fc <- mean(x[temp_black_idx])/mean(x[temp_asian_idx])
    c(temp_p, temp_fc)
  }) %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "name") %>% 
    dplyr::rename(p = V1, fc = V2)
  
  temp_fdr <- p.adjust(temp_p_fc$p, method = "fdr")
  temp_p_fc <- 
    data.frame(temp_p_fc, fdr = temp_fdr,
               stringsAsFactors = FALSE)
  temp_p_fc
  
}) 


temp_p_fc <-
  data.frame(name = p_fc$name,
             fc = p_fc$fc,
             fdr = p_fc$fdr,
             stringsAsFactors = FALSE) %>%
  as.data.frame()


names(p_fc_each_point) <- ga_range_level[-4]
p_fc_each_point <- c(p_fc_each_point, "(10,42]" = list(temp_p_fc))

marker_each_point <- 
  lapply(p_fc_each_point, function(x){
    x <- x %>% 
      dplyr::filter(fdr < 0.01)
    
    if(nrow(x) == 0){
      return(NULL)
    }
    
    x <- 
      x %>% 
      dplyr::mutate(
        class = case_when(
          fc > 1.5 ~ "increase",
          fc < 1/1.5 ~ "decrease",
          TRUE ~ "no"
        )
      ) %>% 
      dplyr::select(name, class)
    x
  })


####
#####a sankey 
marker_each_point %>% 
  lapply(nrow) %>% 
  unlist()

all_marker_name <- 
  lapply(marker_each_point, function(x){
    x$name
  }) %>% 
  unlist() %>% 
  unique()

library(ggalluvial)

temp_data <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      return(NULL)
    }
    x <- 
      data.frame(name = all_marker_name,
                 stringsAsFactors = FALSE) %>% 
      left_join(x, by = "name") %>% 
      dplyr::select(name, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
    
  })

temp_data <-
  purrr::map2(.x = temp_data, .y = names(temp_data), .f = function(x,y){
    if(is.null(x)){
      return(NULL)
    }
    data.frame(x, point = y, stringsAsFactors = FALSE)
  })

temp_data <- 
  do.call(rbind, temp_data)


temp_data %>% 
  dplyr::group_by(point, class) %>% 
  dplyr::summarise(n = n())


plot1 <- 
  temp_data %>% 
  dplyr::mutate(point = factor(point, levels = names(marker_each_point))) %>% 
  ggplot(aes(x = point, 
             y = freq,
             stratum = class, 
             alluvium = name,
             fill = class, 
             label = class)) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c(
    "increase" = ggsci::pal_futurama()(10)[2],
    "decrease" = ggsci::pal_futurama()(10)[3],
    "no" = "grey"
  )) +
  ggalluvial::geom_stratum(alpha = 1) +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot1


# ggsave(plot = plot1, filename = "sankey_plot.pdf", width = 10, height = 7)


###upset plot
library(UpSetR)

increase_marker <- 
purrr::map(marker_each_point, .f = function(x){
  x %>% 
    dplyr::filter(class == "increase") %>% 
    pull(name)
})

decrease_marker <- 
  purrr::map(marker_each_point, .f = function(x){
    x %>% 
      dplyr::filter(class == "decrease") %>% 
      pull(name)
  })


library(VennDiagram)

library(ggVennDiagram)

plot <- 
ggVennDiagram(increase_marker) +
  scale_fill_gradient(low="white",
                      high = ggsci::pal_futurama()(10)[2])

plot

# ggsave(plot, filename = "increase_marker_venn.pdf", width = 7, height = 7)



plot <- 
  ggVennDiagram(decrease_marker) +
  scale_fill_gradient(low="white",
                      high = ggsci::pal_futurama()(10)[3])
plot
# ggsave(plot, filename = "decrease_marker_venn.pdf", width = 7, height = 7)


sum(unlist(decrease_marker) %>% table() == 2)
sum(unlist(decrease_marker) %>% table() == 3)
sum(unlist(decrease_marker) %>% table() == 4)
sum(unlist(decrease_marker) %>% table() >= 3)

sum(unlist(increase_marker) %>% table() == 2)
sum(unlist(increase_marker) %>% table() == 3)
sum(unlist(increase_marker) %>% table() == 4)
sum(unlist(increase_marker) %>% table() >= 3)

##PIUMet analysis
dir.create("PIUMet_analysis")


marker_peak1 <- 
  unlist(increase_marker) %>% table()
marker_peak1 <- 
names(marker_peak1[marker_peak1 >= 3])

marker_peak2 <- 
  unlist(decrease_marker) %>% table()
marker_peak2 <- 
  names(marker_peak2[marker_peak2 >= 3])

marker_peak <- unique(c(marker_peak1, marker_peak2))

marker_peak <- 
data.frame(name = marker_peak, stringsAsFactors = FALSE) %>% 
  dplyr::left_join(p_fc, by = "name") %>% 
  dplyr::mutate(
    polarity = case_when(
      stringr::str_detect(name, "POS") ~ "positive",
      TRUE ~ "negative"
    )
  ) %>% 
  dplyr::left_join(variable_info[,c("name", "mz")], by = "name") %>% 
  dplyr::mutate(fdr = -log(fdr, 10)) %>% 
  dplyr::select(mz, polarity, fdr)


# write.table(marker_peak, file.path("PIUMet_analysis/","marker_peak.txt"), 
#             row.names = FALSE, quote = FALSE,
#             sep = "\t", col.names = FALSE)


names(marker_each_point)

lapply(marker_each_point, function(x){
  table(x$class)
})


##read PIUMet
# readPIUMet(path = "PIUMet_analysis/piumet_output/",
#            variable_info = variable_info,
#            marker_name = unique(c(marker_peak1, marker_peak2)),
#            size_range = c(1,5), text = FALSE)


load("PIUMet_analysis/piumet_output/Result/annotation_result")
annotation_result1 <- annotation_result

annotation_result1 <- 
  annotation_result1 %>% 
  dplyr::distinct() %>% 
  plyr::dlply(.variables = .(HMDB.ID)) %>% 
  lapply(function(x){
    x$name <- paste(x$name, collapse = ";")
    x <- 
      x %>% 
      dplyr::select(HMDB.ID, name,Metabolite.Name) %>% 
      dplyr::distinct()
    x
  }) %>% 
  do.call(rbind, .)


##from HMDB to KEGG
KEGG_ID <- 
  purrr::map(annotation_result1$HMDB.ID, .f = function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database",, to = "KEGG", top = 1)
  })

KEGG_ID <-
KEGG_ID %>% 
  do.call(rbind, .)

load("hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

KEGG_ID2 <-
  hmdb_data$KEGG.ID[match(KEGG_ID$`Human Metabolome Database`, hmdb_data$HMDB.ID)]


KEGG_ID <- data.frame(KEGG_ID, KEGG_ID2, stringsAsFactors = FALSE)

colnames(KEGG_ID) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
apply(KEGG_ID, 1, function(x){
  x <- as.character(x)
  if(is.na(x[2]) & is.na(x[3])){
    return(NA)
  }

  if(!is.na(x[2])){
    return(x[2])
  }

  if(!is.na(x[3])){
    return(x[3])
  }


})

annotation_result1$KEGG_ID <- KEGG_ID

load("../DEG_analysis/marker_in_different_points/annotation_result")
annotation_result2 <- annotation_result

stringr::str_split(annotation_result2$name, pattern = ";") %>% 
  unlist() %>% 
  unique() %>% 
  intersect(annotation_result1$name)


intersect(annotation_result1$Metabolite.Name, annotation_result2$Metabolite.Name)  


annotation_result1$Metabolite.Name



annotation_result1$name
##pathway enrichment
load("hsa_pathway")
kegg_id <- annotation_result1$KEGG_ID[!is.na(annotation_result1$KEGG_ID)]
path_result <- 
  enrichPathway(id = kegg_id, database = hsa_pathway)

load("hsa_disease_pathway")
kegg_id <- annotation_result1$KEGG_ID[!is.na(annotation_result1$KEGG_ID)]
path_disease_result <- 
  enrichPathway(id = kegg_id, database = hsa_disease_pathway)


###graph
load("PIUMet_analysis/piumet_output/Result/edge_data")
load("PIUMet_analysis/piumet_output/Result/node_data")

node_data <- 
node_data %>% 
  dplyr::left_join(annotation_result1[,c("name", "Metabolite.Name")], 
                   by = c("node" = "Metabolite.Name"))

node_data <- 
node_data %>% 
  dplyr::mutate(
    node_class = case_when(
      is.na(name) & node_class == "Metabolite" ~ "Hidden metabolite",
      node_class == "Protein" ~ "Hidden protein",
      TRUE ~ node_class
    )
  )


###calculate fold change for each node
fc <- 
lapply(node_data$name, function(x){
  if(is.na(x)){
    return(NA)
  }
  
  x <- stringr::str_split(x, ";")[[1]]
  
  x <-
  subject_data[x,] %>% 
    apply(2, mean)

  mean(x[which(ethnicity == "Black&Latina")])/
  mean(x[which(ethnicity == "Asian&White")])
  
}) %>% 
  unlist()

# fc[is.na(fc)] <- 1

node_data$fc <- fc

node_data <-
node_data %>% 
  dplyr::mutate(fc2 = 
                  dplyr::case_when(
                    node_class != "Metabolite" ~ 0,
                   TRUE ~ log(fc, 2)
                  )
                ) %>% 
  dplyr::mutate(
    node_class = case_when(
      is.na(fc) ~ node_class,
      fc > 1 & node_class == "Metabolite" ~ "Increase metabolite",
      fc < 1 & node_class == "Metabolite" ~ "Decrease metabolite"
    )
  )


node_data <- 
node_data %>% 
  dplyr::rename(peak_name = name)

graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = tidygraph::centrality_degree(mode = 'all'))

degree <- igraph::vertex_attr(graph)$Degree
node_class <- igraph::vertex_attr(graph)$node_class
degree[node_class == "m/z Peak"] <- 0.5

graph <- 
igraph::set_vertex_attr(graph = graph, name = "Degree", value = degree)


col <- 
  fill <- 
  c(
    "m/z Peak" = "gray87",
    "Hidden metabolite" = "gray87",
    "Metabolite" = "red",
    "Hidden protein" = "steelblue1"
  )

shape = c(
  "m/z Peak" = 24,
  "Hidden metabolite" = 22,
  "Increase metabolite" = 22,
  "Decrease metabolite" = 22,
  # "Metabolite_others" = 22,
  "Hidden protein" = 21
)


plot <-
  ggraph(graph,
         layout = "kk") +
  geom_edge_link(aes(edge_width = V3),
                 alpha = 1,
                 color = "black",
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree, 
                      fill = fc2,
                      shape = node_class),
                  alpha = 1, 
                  show.legend = TRUE) +
  scale_shape_manual(values = shape) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_width(range = c(0.1,0.3)) +
  scale_size_continuous(range = c(2,12)) +
  scale_fill_gradient2(low = ggsci::pal_futurama()(10)[3], 
                       mid = "white", 
                       high = ggsci::pal_futurama()(10)[2], 
                       midpoint = 0) +
  # scale_color_manual(values = col) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))
plot
  

ggsave(plot, filename = "PIUMet_analysis/graph_plot.pdf", width = 10, height = 7)


node_attr <- igraph::vertex_attr(graph = graph)
node_attr <- 
node_attr %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  dplyr::filter(node_class != "m/z Peak") %>% 
dplyr::filter(node_class != "Hidden metabolite") %>% 
  dplyr::arrange(Degree) %>% 
  dplyr::mutate(Degree = as.numeric(Degree)) %>% 
  dplyr::filter(Degree > 4) %>% 
  dplyr::mutate(node = factor(node, levels = node)) %>% 
  dplyr::mutate(node_class = case_when(
    stringr::str_detect(node_class, "metabolite") ~ "Metabolite",
    TRUE ~ node_class
  ))
  

fill <- 
  c(
    "m/z Peak" = "gray87",
    "Metabolite" = "red",
    "Hidden protein" = "steelblue1"
  )


plot <- 
node_attr %>%
  dplyr::mutate(fc2 = as.numeric(fc2)) %>% 
  dplyr::mutate(node_class = factor(node_class, levels = c("Metabolite", "Hidden protein"))) %>% 
ggplot(aes(x = Degree, y = node)) +
  geom_bar(stat = "identity", 
           aes(fill = fc2, color = node_class), show.legend = FALSE) +
  scale_fill_gradient2(low = ggsci::pal_futurama()(10)[3], 
                       mid = "white", 
                       high = ggsci::pal_futurama()(10)[2], 
                       midpoint = 0) +
  scale_color_manual(values = fill) +
  theme_bw() +
  labs(x = "Degree", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
  facet_wrap(~node_class, scales = "free") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 10))

plot

ggsave(plot, filename = "PIUMet_analysis/degree_of_graph.pdf", width = 10, height = 7) 


plot <- 
path_result %>% 
  dplyr::filter(Overlap >= 3) %>% 
  ggplot(aes(x = Overlap, y = -log(p.value.fdr))) +
  geom_point(aes(size = Overlap)) +
  theme_bw() +
  ggrepel::geom_label_repel(mapping = aes(Overlap, 
                                          -log(p.value.fdr), 
                                          label = Pathway.name)) +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 10))
  
ggsave(plot, filename = "PIUMet_analysis/pathway_enrichment_plot_text.pdf", 
       width = 7, height = 7)



###
path_result$Pathway.name[2]
path_result$Pathway.ID[2]

grep("hsa00140",names(hsa_pathway))

##pathway quantative
idx <- which(annotation_result2$KEGG_ID %in% hsa_pathway[[15]])

peak_name <- annotation_result2$name[idx] %>% 
  stringr::str_split(";") %>% 
  unlist() %>% 
  unique()

temp_data <- 
subject_data[peak_name,] %>% 
  apply(1, function(x){
    x/sd(x)
  }) %>%
  t() %>%
  as.data.frame() %>%
  apply(2, mean)

ga <- sample_info$GA

ga_range <- 
  as.character(cut_width(x = ga, width = 2, center = 3)) %>%
  stringr::str_replace("\\[", '(')

table(ga_range)

ga_range[ga_range == "(10,12]"] <- "(10,26]"
ga_range[ga_range == "(12,14]"] <- "(10,26]"
ga_range[ga_range == "(14,16]"] <- "(10,26]"
ga_range[ga_range == "(16,18]"] <- "(10,26]"
ga_range[ga_range == "(18,20]"] <- "(10,26]"
ga_range[ga_range == "(20,22]"] <- "(10,26]"
ga_range[ga_range == "(22,24]"] <- "(10,26]"
ga_range[ga_range == "(24,26]"] <- "(10,26]"

ga_range[ga_range == "(26,28]"] <- "(26,34]"
ga_range[ga_range == "(28,30]"] <- "(26,34]"
ga_range[ga_range == "(30,32]"] <- "(26,34]"
ga_range[ga_range == "(32,34]"] <- "(26,34]"

ga_range[ga_range == "(34,36]"] <- "(34,42]"
ga_range[ga_range == "(36,38]"] <- "(34,42]"
ga_range[ga_range == "(38,40]"] <- "(34,42]"
ga_range[ga_range == "(40,42]"] <- "(34,42]"

ga_range[ga_range == "(44,46]"] <- "PP"

table(ga_range)

ga_level <- sort(unique(ga_range))

ethnicity_color <-
  ggsci::pal_futurama()(7)[c(7,2)]

names(ethnicity_color) <- sort(unique(ethnicity))


temp_data <- 
data.frame(ga_range, ethnicity, value = temp_data, stringsAsFactors = FALSE) 

head(temp_data)

wilcox.test(
  temp_data %>% 
    dplyr::filter(ga_range == "(34,42]" & ethnicity == "Asian&White") %>% 
    pull(value),
  temp_data %>% 
    dplyr::filter(ga_range == "(34,42]" & ethnicity == "Black&Latina") %>% 
    pull(value)  
)



plot <- 
temp_data %>%
  # dplyr::mutate(ga_range = factor(ga_range, levels = ga_level))
  group_by(ethnicity,
           ga_range) %>%
  dplyr::mutate(mean = mean(value),
                sem = sd(value) / sqrt(nrow(temp_data))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ga_range = factor(ga_range, levels = ga_level)) %>%
  ggplot(aes(ga_range, mean)) +
  # geom_point(aes(color = ethnicity)) +
  geom_line(aes(group = ethnicity, color = ethnicity), show.legend = FALSE) +
  geom_point(aes(group = ethnicity, color = ethnicity), 
             shape = 16,
             size = 3,
             show.legend = FALSE) +
  geom_errorbar(
    mapping = aes(
      ymin = mean - sem,
      ymax = mean + sem,
      color = ethnicity
    ),
    width = 0, show.legend = FALSE
  ) +
  scale_color_manual(values = ethnicity_color) +
  labs(x = "", y = "PIUMet pathway intensity") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))

plot

ggsave(plot, filename = "PIUMet_analysis/Tryptophan metabolism_in_different_point.pdf", 
       width = 14, height = 7)


name <-
annotation_result1$Metabolite.Name[match(hsa_pathway[[31]], annotation_result1$KEGG_ID)]

name <- name[!is.na(name)]

super_class <- 
  node_data$super.class

super_class[!is.na(super_class)]

plot <-
  ggraph(graph,
         layout = "kk") +
  geom_edge_link(aes(edge_width = V3),
                 alpha = 1,
                 color = "black",
                 show.legend = TRUE) +
  geom_node_point(aes(size = fc2, 
                      fill = node_class,
                      shape = node_class),
                  alpha = 1, show.legend = TRUE) +
  ggraph::geom_node_text(aes(label = ifelse(node %in% name, node, NA), 
                             color = node_class), 
                         repel = TRUE, size = 3) +
  scale_shape_manual(values = shape) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_width(range = c(0.1,0.3)) +
  scale_size_continuous(range = c(4,12)) +
  scale_fill_manual(values = fill) +
  scale_color_manual(values = col) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))
plot


ggsave(plot, filename = "PIUMet_analysis/graph_plot_text.pdf", 
       width = 20, height = 20)









###correlation network
dim(annotation_result1)
##from peak to metabolite
temp_data <- 
annotation_result1$name %>% 
  purrr::map(
    .f = function(x){
      x <- stringr::str_split(x, ";")[[1]]
      subject_data_mean[x,,drop = FALSE] %>% 
        apply(2, mean)
    }
  ) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()


rownames(temp_data) <- annotation_result1$HMDB.ID

##calculate correlation
library(corrr)
cor_value <- corrr::correlate(x = t(temp_data), method = "spearman")
cor_value <-
  cor_value %>%
  shave() %>%
  stretch() %>%
  dplyr::filter(!is.na(r))

p <-
  as.data.frame(t(cor_value)) %>%
  purrr::map(.f = function(x){
    cor.test(as.numeric(temp_data[x[1],]),
             as.numeric(temp_data[x[2],]), method = "spearman"
             )$p.value
  }) %>%
  unlist()

fdr <- p.adjust(p, method = "fdr")
fdr[fdr == 0] <- min(fdr[fdr != 0])

cor_value <- data.frame(cor_value, p, fdr, stringsAsFactors = FALSE) %>%
  dplyr::filter(fdr < 0.05)

save(cor_value, file = "cor_value")

cor_value <- cor_value %>% 
  # dplyr::filter(abs(r) > quantile(abs(r), 0.75)) %>% 
  dplyr::filter(abs(r) > 0.5)


edge_data <- 
  cor_value %>% 
  dplyr::mutate(from = x, 
                to = y,
                fdr = -log(fdr, 10),
                cor = r,
                abs.cor = abs(r)) %>% 
  dplyr::select(from, to, cor, abs.cor, fdr) 


node_data <- data.frame(node = unique(c(edge_data$from, edge_data$to)),
                        stringsAsFactors = FALSE)

node_data <- 
node_data %>% 
  dplyr::left_join(annotation_result1, by = c("node" = "HMDB.ID")) %>% 
  dplyr::rename(peak_name = name)

library(igraph)
library(tidygraph)

graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


igraph::diversity(graph = trans_graph, weights = abs(edge_attr(trans_graph,
                                                               "cor")))

graph <- tidygraph::as.igraph(x = graph)

subnetworks <-
  igraph::cluster_edge_betweenness(graph = graph,
                                   weights = abs(edge_attr(graph,
                                                           "cor")))
save(subnetworks, file = "subnetworks")
load("subnetworks")


table(membership(subnetworks))


plot <- 
  ggplot(
    data.frame(index = 1:length(subnetworks$modularity),
               modu = subnetworks$modularity, stringsAsFactors = FALSE),
    aes(index, modu) 
  ) +
  geom_vline(xintercept = which.max(subnetworks$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line(colour = "black") +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12))

plot <-
  plot + 
  ggplot2::annotate(geom = "point", 
                    x = which.max(subnetworks$modularity),
                    y = max(subnetworks$modularity), 
                    size = 3, 
                    colour = "#FFA319FF") +
  annotate(geom = "text", 
           x = which.max(subnetworks$modularity),
           y = max(subnetworks$modularity), 
           label = paste("(",  which.max(subnetworks$modularity),
                         ",", 
                         max(subnetworks$modularity) %>% round(3),
                         ")"),
           size = 5,
           colour = "#FFA319FF"
  )

plot

ggsave(plot, filename = "modularity.pdf", width = 7, height = 7)
ggsave(plot, filename = "modularity.png", width = 7, height = 7)


table(membership(communities = subnetworks))
which(table(membership(communities = subnetworks)) >= 3)
##cluster 1 2 3 4 5 7 8
node <- vertex_attr(graph = graph, name = "node")

membership <- 
data.frame(node, membership = as.numeric(membership(communities = subnetworks)), 
           stringsAsFactors = FALSE) %>% 
  # dplyr::arrange(membership) %>% 
  # dplyr::filter(membership %in% c(1,2,3,4,5,7,8)) %>% 
  dplyr::mutate(membership = paste("Cluser", membership, sep = "")) %>% 
  dplyr::mutate(membership = case_when(
    membership %in% c(paste("Cluser", c(1,2,3,4,5,7,8), sep = "")) ~ membership,
    TRUE ~ "Other"
  ))

graph <- 
igraph::set_vertex_attr(graph = graph, name = "membership", value = membership$membership)


save(graph, file = "graph")
# angle <- 360 * (c(1:nrow(node_data)) - 0.5)/nrow(node_data)
# hjust <- ifelse(angle > 180, 1, 0)
# angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)


plot <-
  ggraph(graph,
         layout = 'fr') +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE
  ) +
  geom_node_point(aes(size = Degree,
                      fill = membership),
                  alpha = 1, 
                  shape = 21,
                  show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(low = "#3B4992FF", mid = "white",
                                     high = "#EE0000FF", midpoint = 0) +
  ggraph::scale_edge_width(range = c(0.3,1)) +
  geom_node_text(aes(label = Metabolite.Name), repel = TRUE) +
  scale_size_continuous(range = c(2, 10)) +
  ggsci::scale_fill_futurama() +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))
# coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-1.4,1.4))

plot

# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "PIUMet_analysis/cor_network_all.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)

ggsave(
  plot,
  filename = "PIUMet_analysis/cor_network_all.png",
  width = 8,
  height = 7,
  bg = "transparent"
)



###from metabolite to module
membership1 <- 
membership %>% 
  dplyr::filter(membership != "Other")

temp_data <- 
lapply(unique(membership1$membership), function(x){
  x <- membership1$node[membership1$membership == x]
  peak_name <- annotation_result1$name[match(x, annotation_result1$HMDB.ID)]
    peak_name <- stringr::str_split(peak_name, ";") %>%
      unlist() %>% 
      unique()
    subject_data_mean[peak_name,] %>% 
      apply(2, mean)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()


rownames(temp_data) <- unique(membership1$membership)

temp_data <- 
temp_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()


library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c(ggsci::pal_aaas()(10)[5], "white", ggsci::pal_aaas()(10)[6]))


range(temp_data)

temp_data[temp_data < -2.05] <- -2.05

plot <- 
  Heatmap(temp_data, 
          cluster_columns = FALSE, 
          cluster_rows = TRUE, 
          show_row_names = TRUE, 
          show_column_names = TRUE,
          border = TRUE, 
          col = col_fun)

plot

library(ggplotify)

plot <- as.ggplot(plot)


ggsave(
  plot,
  file = "PIUMet_analysis/module_heatmap.pdf",
  width = 7,
  height = 7,
  bg = "transparent"
)




####
##module quantative
##cluster 
temp_data <- 
  lapply(unique(membership1$membership), function(x){
    x <- membership1$node[membership1$membership == x]
    peak_name <- annotation_result1$name[match(x, annotation_result1$HMDB.ID)]
    peak_name <- stringr::str_split(peak_name, ";") %>%
      unlist() %>% 
      unique()
    
    subject_data[peak_name,] %>% 
      apply(1, function(x){
        x/sd(x)
      }) %>%
      t() %>%
      as.data.frame() %>%
      apply(2, mean)
    
  })

names(temp_data) <- unique(membership1$membership)

ga <- sample_info$GA

ga_range <- 
  as.character(cut_width(x = ga, width = 2, center = 3)) %>%
  stringr::str_replace("\\[", '(')

table(ga_range)

# ga_range[ga_range == "(10,12]"] <- "(10,26]"
# ga_range[ga_range == "(12,14]"] <- "(10,26]"
# ga_range[ga_range == "(14,16]"] <- "(10,26]"
# ga_range[ga_range == "(16,18]"] <- "(10,26]"
# ga_range[ga_range == "(18,20]"] <- "(10,26]"
# ga_range[ga_range == "(20,22]"] <- "(10,26]"
# ga_range[ga_range == "(22,24]"] <- "(10,26]"
# ga_range[ga_range == "(24,26]"] <- "(10,26]"
# 
# ga_range[ga_range == "(26,28]"] <- "(26,34]"
# ga_range[ga_range == "(28,30]"] <- "(26,34]"
# ga_range[ga_range == "(30,32]"] <- "(26,34]"
# ga_range[ga_range == "(32,34]"] <- "(26,34]"
# 
# ga_range[ga_range == "(34,36]"] <- "(34,42]"
# ga_range[ga_range == "(36,38]"] <- "(34,42]"
# ga_range[ga_range == "(38,40]"] <- "(34,42]"
# ga_range[ga_range == "(40,42]"] <- "(34,42]"
# 
# ga_range[ga_range == "(44,46]"] <- "PP"


ga_range[ga_range == "(10,12]"] <- "(10,16]"
ga_range[ga_range == "(12,14]"] <- "(10,16]"
ga_range[ga_range == "(14,16]"] <- "(10,16]"

ga_range[ga_range == "(38,40]"] <- "(38,42]"
ga_range[ga_range == "(40,42]"] <- "(38,42]"

ga_range[ga_range == "(44,46]"] <- "PP"


table(ga_range)

ga_level <- sort(unique(ga_range))

ethnicity_color <-
  ggsci::pal_futurama()(7)[c(7,2)]

names(ethnicity_color) <- sort(unique(ethnicity))


temp_data <- 
  purrr::map2(.x = temp_data, 
              .y = names(temp_data),
    .f = function(x, y){
      data.frame(ga_range, ethnicity,
                 value = x, 
                 cluster = y,
                 stringsAsFactors = FALSE)    
    }
  )


wilcox.test(
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "(10,26]" & ethnicity == "Asian&White") %>%
    pull(value),
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "(10,26]" & ethnicity == "Black&Latina") %>%
    pull(value)
)

wilcox.test(
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "(26,34]" & ethnicity == "Asian&White") %>%
    pull(value),
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "(26,34]" & ethnicity == "Black&Latina") %>%
    pull(value)
)

wilcox.test(
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "(34,42]" & ethnicity == "Asian&White") %>%
    pull(value),
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "(34,42]" & ethnicity == "Black&Latina") %>%
    pull(value)
)

wilcox.test(
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "PP" & ethnicity == "Asian&White") %>%
    pull(value),
  temp_data[[7]] %>%
    dplyr::filter(ga_range == "PP" & ethnicity == "Black&Latina") %>%
    pull(value)
)


temp_data <- 
temp_data %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

plot <- 
  temp_data %>%
  dplyr::filter(cluster %in% c("Cluser5", "Cluser8")) %>% 
  # dplyr::mutate(ga_range = factor(ga_range, levels = ga_level))
  group_by(ethnicity,
           ga_range,
           cluster) %>%
  dplyr::mutate(mean = mean(value),
                sem = sd(value) / sqrt(nrow(temp_data))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ga_range = factor(ga_range, levels = ga_level)) %>%
  ggplot(aes(ga_range, mean)) +
  # geom_point(aes(color = ethnicity)) +
  geom_line(aes(group = ethnicity, color = ethnicity), show.legend = FALSE) +
  # geom_point(aes(group = ethnicity, color = ethnicity), 
  #            shape = 16,
  #            size = 3,
  #            show.legend = FALSE) +
  geom_errorbar(
    mapping = aes(
      ymin = mean - sem,
      ymax = mean + sem,
      color = ethnicity
    ),
    width = 0, show.legend = FALSE
  ) +
  scale_color_manual(values = ethnicity_color) +
  labs(x = "", y = "PIUMet module intensity") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12)) + 
  facet_wrap(~cluster, ncol = 1)

plot

ggsave(plot, filename = "PIUMet_analysis/module_in_different_point.pdf", 
       width = 7, height = 7)


#####new grah
membership <- igraph::vertex_attr(graph = graph)$membership
membership[membership != "Cluser5"  & membership != "Cluser8"] <- "Other"

graph2 <- graph
graph2 <- 
igraph::set_vertex_attr(graph = graph2, name = "membership", value = membership)


plot <-
  ggraph(graph2,
         layout = 'fr') +
  geom_edge_link(aes(edge_width = fdr,
                     color = cor),
                 alpha = 0.5,
                 show.legend = TRUE
  ) +
  geom_node_point(aes(size = Degree,
                      fill = membership),
                  alpha = 1, 
                  shape = 21,
                  show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggraph::scale_edge_color_gradient2(low = "#3B4992FF", mid = "white",
                                     high = "#EE0000FF", midpoint = 0) +
  ggraph::scale_edge_width(range = c(0.1,1)) +
  geom_node_text(aes(label = ifelse(membership == "Cluser5" | membership == "Cluser8" ,
                                    Metabolite.Name,NA)), 
                 repel = TRUE) +
  scale_size_continuous(range = c(2, 10)) +
  ggsci::scale_fill_futurama() +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))
# coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-1.4,1.4))

plot

# extrafont::font_import()
extrafont::loadfonts()

ggsave(
  plot,
  filename = "PIUMet_analysis/cor_network_all2.pdf",
  width = 8,
  height = 7,
  bg = "transparent"
)





###
membership <- igraph::vertex_attr(graph2) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  dplyr::select(node, membership) %>% 
  dplyr::arrange(membership)

write.csv(membership, "membership.csv", row.names = FALSE)

cluster5 <- 
membership %>% 
  dplyr::filter(membership == "Cluser5") %>% 
  pull(node)

hmdb_data[match(cluster5, hmdb_data$Lab.ID),c(33,34,35,36)]$Class %>% table()


cluster8 <- 
  membership %>% 
  dplyr::filter(membership == "Cluser8") %>% 
  pull(node)

hmdb_data[match(cluster8, hmdb_data$Lab.ID),c(33,34,35,36)]$Class %>% table()





