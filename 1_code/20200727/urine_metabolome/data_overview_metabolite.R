#to avoind source
no_exist_function()
#-------------------------------------------------------------------------------
##RPLC pos and neg
###-------------------------------------------------------------------------------
sxtTools::setwd_project()
rm(list=ls())

##metabolites
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/metabolites/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/metabolites/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/metabolites/variable_info")

setwd("data_analysis20200108/urine_metabolome/data_overview/metabolites/")

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

subject_data2 <- 
  t(subject_data) %>% 
  tibble::as_tibble()

rownames(subject_data2)

qc_data2 <- 
  t(qc_data) %>% 
  tibble::as_tibble()

rownames(subject_data2) <- colnames(subject_data)

####subject and QC together
temp_subject <- 
  subject_data2 %>% 
  mutate(GA = ga,
         batch = batch)

temp_qc <- 
  qc_data2 %>% 
  mutate(GA = 0,
         batch = "QC")

rownames(temp_subject) <-
  colnames(subject_data)

rownames(temp_qc) <-
  colnames(qc_data)

colnames(temp_qc) <- colnames(temp_subject)

temp_data <- 
  rbind(temp_subject, temp_qc)

pca_object <- 
  prcomp(x = 
           temp_data %>% select(-c(GA, batch)))

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                batch = temp_data$batch,
                stringsAsFactors = FALSE)

x <- 
  x %>% 
  rownames_to_column(var = "name")


plot <- 
  ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(colour = batch), size = 5, shape = 16) +
  ggsci::scale_colour_futurama(alpha = 0.7) +
  # guides(colour = guide_colourbar(title = "Class", title.position = "top")) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_data_quality_plot.pdf", 
       width = 7, height = 7)


####PCA without PCA only for subjects
temp_data <- 
  subject_data2

pca_object <- 
  prcomp(x = temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = ga,
                stringsAsFactors = FALSE)

plot <- 
  ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 16) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#C71000FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))


plot


ggsave(plot, filename = "pca_plot.pdf", 
       width = 7, height = 7)


####for each participant
x <- 
  data.frame(sample_id = colnames(subject_data),
             x, stringsAsFactors = FALSE)

x <- 
  x %>% 
  left_join(sample_info[,c("sample_id", "subject_id")],
            by = c("sample_id"))

x2 <-
  x %>% 
  mutate(y = 0) %>% 
  dplyr::filter( GA == 0)


plot <- 
x %>% 
  dplyr::filter(GA != 0) %>% 
  mutate(y = 0) %>% 
  ggplot(aes(x = PC2, y = y, colour = GA)) +
  geom_point(shape = 16) +
  geom_point(mapping = aes(x = PC2, y = y),
             data = x2, color = "#C71000FF", shape = 16) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = 0, 
  #          colour = "#C71000FF") +
  ggrepel::geom_text_repel(aes(PC2, y = y, label = round(GA, 2)), 
                           size = 2.5) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  labs(y = "") +
  facet_wrap(~subject_id, ncol = 6) +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "top",
        legend.background = element_blank(),
        strip.background = element_rect(fill = "grey"), 
        strip.text = element_text(color = "white"),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank()
  )

plot

ggsave(plot, filename = "pca_for_each_person.pdf", width = 10, height = 8)

###tsne analysis
# tsne_object <- Rtsne::Rtsne(
#   X = as.matrix(temp_data),
#   dims = 2,
#   perplexity = 30,
#   verbose = TRUE
# )
# 
# Y <- tsne_object$Y
# Y <-
#   data.frame(Y, "GA" = sample_info$GA[sample_info$batch == 1],
#              stringsAsFactors = FALSE)
# 
# (
#   plot <- ggplot(Y, aes(X1, X2, colour = GA)) +
#     geom_point(size = 3) +
#     labs(x = "Dimension 1",
#          y = "Dimension 2") +
#     theme_bw() +
#     scale_colour_gradient(low = alpha("#155F83FF", 0.1),
#                           high = alpha("#FFA319FF", 1)) +
#     guides(colour = guide_colourbar(title = "GA (week)")) +
#     theme(
#       axis.title = element_text(size = 15),
#       axis.text = element_text(size = 12),
#       legend.title = element_text(size = 15),
#       legend.text = element_text(size = 12),
#       strip.background = element_rect(fill = "#0099B47F"),
#       strip.text = element_text(color = "white", size = 15)
#     )
# )
