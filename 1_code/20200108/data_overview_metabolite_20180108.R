#-------------------------------------------------------------------------------
##RPLC pos
#-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())
setwd("data_analysis20200108/data_overview/")
load("../data_cleaning/RPLC/POS/rplc_pos_6")
library(metflow2)
library(tidyverse)

sample_info <- rplc_pos_6@sample.info

subject_data <- metflow2::getData(object = rplc_pos_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(rplc_pos_6,
                             slot = "QC")

qc_data <- 
  qc_data %>% 
  select(-c(QC2.1, QC2.2, QC2.3, QC_0_25))

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)


subject_data <- subject_data[remain_idx,]
qc_data <- qc_data[remain_idx,]

sample_info <- 
  sample_info %>% 
  dplyr::filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))

###log
subject_data <- 
  log(subject_data)

qc_data <- 
  log(qc_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  tibble::as_tibble()


qc_data <- 
  t(
    apply(qc_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  tibble::as_tibble()

subject_data2 <- 
  t(subject_data) %>% 
  tibble::as_tibble()

rownames(subject_data2)

qc_data2 <- 
  t(qc_data) %>% 
  tibble::as_tibble()

####subject and QC together
temp_subject <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA,
         batch = sample_info$batch)

temp_qc <- 
  qc_data2 %>% 
  mutate(GA = 0,
         batch = "QC")

rownames(temp_subject) <-
  colnames(subject_data)

rownames(temp_qc) <-
  colnames(qc_data)

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

pca_data_quality_plot <- 
ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(colour = batch), size = 5) +
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
  
pca_data_quality_plot
save(pca_data_quality_plot, file = "pca_data_quality_plot")
ggsave(pca_data_quality_plot, filename = "pca_data_quality_plot.pdf", 
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
                sample_info,
                stringsAsFactors = FALSE)


pca_pos_plot <- 
ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
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
           colour = "#8A9045FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))


pca_pos_plot

save(pca_pos_plot, file = "pca_pos_plot")
ggsave(pca_pos_plot, filename = "pca_pos_plot.pdf", 
       width = 7, height = 7)




#-------------------------------------------------------------------------------
##RPLC pneg
#-------------------------------------------------------------------------------
load("../data_cleaning/RPLC/NEG/rplc_neg_6")
library(metflow2)
library(tidyverse)

sample_info <- rplc_neg_6@sample.info

subject_data <- metflow2::getData(object = rplc_neg_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(rplc_neg_6,
                             slot = "QC")

# qc_data <- 
#   qc_data %>% 
#   select(-c(QC2.1, QC2.2, QC2.3, QC_0_25))

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)


subject_data <- subject_data[remain_idx,]
qc_data <- qc_data[remain_idx,]

sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))

###log
subject_data <- 
  log(subject_data)

qc_data <- 
  log(qc_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  tibble::as_tibble()


qc_data <- 
  t(
    apply(qc_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>% 
  tibble::as_tibble()

subject_data2 <- 
  t(subject_data) %>% 
  tibble::as_tibble()

rownames(subject_data2)

qc_data2 <- 
  t(qc_data) %>% 
  tibble::as_tibble()

####subject and QC together
temp_subject <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA,
         batch = sample_info$batch)

temp_qc <- 
  qc_data2 %>% 
  mutate(GA = 0,
         batch = "QC")

rownames(temp_subject) <-
  colnames(subject_data)

rownames(temp_qc) <-
  colnames(qc_data)

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

pca_neg_data_quality_plot <- 
  ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(colour = batch), size = 5) +
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

pca_neg_data_quality_plot

save(pca_neg_data_quality_plot, file = "pca_neg_data_quality_plot")
ggsave(pca_neg_data_quality_plot, filename = "pca_neg_data_quality_plot.pdf", 
       width = 7, height = 7)



####PCA with out PCA only for subjects
temp_data <- 
  subject_data2

pca_object <- 
  prcomp(x = temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                sample_info,
                stringsAsFactors = FALSE)


pca_neg_plot <- 
  ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
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
           colour = "#8A9045FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))


pca_neg_plot

save(pca_neg_plot, file = "pca_neg_plot")
ggsave(pca_neg_plot, filename = "pca_neg_plot.pdf", 
       width = 7, height = 7)






###RPLC POS and NEG
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/features/")
load("../../data_cleaning/RPLC/POS/rplc_pos_6")
load("../../data_cleaning/RPLC/NEG/rplc_neg_6")

library(metflow2)
library(tidyverse)

sample_info_pos <- rplc_pos_6@sample.info

subject_data_pos <- metflow2::getData(object = rplc_pos_6, 
                                      slot = "Subject")

qc_data_pos <- metflow2::getData(rplc_pos_6,
                                 slot = "QC")
###remove peaks with large RSD
qc_rsd_pos <- apply(qc_data_pos, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_pos <- 
  which(qc_rsd_pos < 30)


subject_data_pos <- subject_data_pos[remain_idx_pos,]
qc_data_pos <- qc_data_pos[remain_idx_pos,]


sample_info_pos <- 
  sample_info_pos %>% 
  dplyr::filter(class == "Subject")

subject_data_pos <- 
  subject_data_pos %>% 
  dplyr::select(one_of(sample_info_pos$sample.name))

###log
subject_data_pos <- 
  log(subject_data_pos)

qc_data_pos <- 
  log(qc_data_pos)

sample_info_neg <- rplc_neg_6@sample.info

subject_data_neg <- metflow2::getData(object = rplc_neg_6, 
                                      slot = "Subject")

qc_data_neg <- metflow2::getData(rplc_neg_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_neg <- apply(qc_data_neg, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_neg <- 
  which(qc_rsd_neg < 30)


subject_data_neg <- subject_data_neg[remain_idx_neg,]
qc_data_neg <- qc_data_neg[remain_idx_neg,]


sample_info_neg <- 
  sample_info_neg %>% 
  dplyr::filter(class == "Subject")

subject_data_neg <- 
  subject_data_neg %>% 
  dplyr::select(one_of(sample_info_neg$sample.name))

###log
subject_data_neg <- 
  log(subject_data_neg)


##combine pos and neg
subject_data_pos <- 
  subject_data_pos %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


subject_data_neg <- 
  subject_data_neg %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


colnames(subject_data_pos) == colnames(subject_data_neg)

subject_data <- 
  rbind(subject_data_pos, subject_data_neg)


qc_data_pos <- 
  qc_data_pos %>% 
  select(one_of(intersect(colnames(qc_data_pos), colnames(qc_data_neg))))


qc_data_neg <- 
  qc_data_neg %>% 
  select(one_of(intersect(colnames(qc_data_pos), colnames(qc_data_neg))))


colnames(qc_data_pos) == colnames(qc_data_neg)

qc_data <- 
  rbind(qc_data_pos, qc_data_neg)

sample_info <- 
  sample_info_pos %>% 
  dplyr::filter(sample.name %in% colnames(subject_data))


colnames(subject_data) == sample_info$sample.name


subject_data <-
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>%
  as_tibble()


qc_data <-
  t(
    apply(qc_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  ) %>%
  as_tibble()


qc_data <- 
  qc_data %>% 
  select(-c(QC2.1,QC2.2, QC2.3))

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)


qc_data2 <- t(qc_data)
qc_data2 <- as_tibble(qc_data2)
rownames(qc_data2)


####subject and QC together
temp_subject <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA,
         batch = sample_info$batch)

temp_qc <- 
  qc_data2 %>% 
  mutate(GA = 0,
         batch = "QC")

rownames(temp_subject) <-
  colnames(subject_data)

rownames(temp_qc) <-
  colnames(qc_data)


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

pca_pos_neg_data_quality_plot <- 
  ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(colour = batch), size = 5) +
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
  # ggrepel::geom_text_repel(aes(x = PC1, PC2, label = name))

pca_pos_neg_data_quality_plot

save(pca_pos_neg_data_quality_plot, file = "pca_pos_neg_data_quality_plot")
ggsave(pca_pos_neg_data_quality_plot, filename = "pca_pos_neg_data_quality_plot.pdf", 
       width = 7, height = 7)


##PCA without QC
####PCA with out PCA only for subjects
temp_data <- 
  subject_data2

pca_object <- 
  prcomp(x = temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                sample_info,
                stringsAsFactors = FALSE)


pca_pos_neg_plot <- 
  ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
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


pca_pos_neg_plot

save(pca_pos_neg_plot, file = "pca_pos_neg_plot")
ggsave(pca_pos_neg_plot, filename = "pca_pos_neg_plot.pdf", 
       width = 7, height = 7)


#######################################################
##use metabolite to check the overview
#-------------------------------------------------------------------------------
##RPLC pos
#-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())
setwd("data_analysis20200108/data_overview/metabolites/")
load("../../data_cleaning/RPLC/POS/rplc_pos_6")
load("../../data_preparation_for_analysis/metabolite_tags")
library(metflow2)
library(tidyverse)

sample_info <- rplc_pos_6@sample.info

subject_data <- metflow2::getData(object = rplc_pos_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(rplc_pos_6,
                             slot = "QC")


rownames(subject_data) <-
  rownames(qc_data) <-
  rplc_pos_6@ms1.data[[1]]$name

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)

name <- rownames(subject_data)[remain_idx]

subject_data <- 
  subject_data[remain_idx,]

rownames(subject_data) <- name

sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))


###log
subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


subject_data <-
  subject_data %>% 
  as_tibble()

rownames(subject_data) <- name

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)
colnames(subject_data2) <- 
  rownames(subject_data)
subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA)) %>% 
  select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA)) %>% 
  select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA)) %>% 
  select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)





#####RPLC negative
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/metabolites/")
rm(list=ls())
load("../../data_cleaning/RPLC/NEG/rplc_neg_6")
load("../../data_preparation_for_analysis/metabolite_tags")
library(metflow2)
library(tidyverse)

sample_info <- rplc_neg_6@sample.info

subject_data <- metflow2::getData(object = rplc_neg_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(rplc_neg_6,
                             slot = "QC")


rownames(subject_data) <- rownames(qc_data) <-
  rplc_neg_6@ms1.data[[1]]$name

subject_data <- as.data.frame(subject_data)

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)


subject_data <- subject_data[remain_idx,]

# sample_info <-
#   sample_info %>%
#   filter(GA != 0)


sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))


###log
subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


# subject_data <-
#   subject_data %>% 
#   as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as.data.frame(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))
  # select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)





###RPLC POS and NEG
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/")
load("../data_cleaning/RPLC/POS/rplc_pos_6")
load("../data_cleaning/RPLC/NEG/rplc_neg_6")

library(metflow2)
library(tidyverse)

sample_info_pos <- rplc_pos_6@sample.info

subject_data_pos <- metflow2::getData(object = rplc_pos_6, 
                                      slot = "Subject")

qc_data_pos <- metflow2::getData(rplc_pos_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_pos <- apply(qc_data_pos, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_pos <- 
  which(qc_rsd_pos < 30)


subject_data_pos <- subject_data_pos[remain_idx_pos,]


sample_info_pos <- 
  sample_info_pos %>% 
  filter(class == "Subject")

subject_data_pos <- 
  subject_data_pos %>% 
  dplyr::select(one_of(sample_info_pos$sample.name))

###log
subject_data_pos <- 
  log(subject_data_pos)

# subject_data_pos <- 
#   t(
#     apply(subject_data_pos, 1, function(x){
#       (x - mean(x))/sd(x)
#     })
#   )
# 
# subject_data_pos <-
#   subject_data_pos %>% 
#   as_tibble()


sample_info_neg <- rplc_neg_6@sample.info

subject_data_neg <- metflow2::getData(object = rplc_neg_6, 
                                      slot = "Subject")

qc_data_neg <- metflow2::getData(rplc_neg_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_neg <- apply(qc_data_neg, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_neg <- 
  which(qc_rsd_neg < 30)


subject_data_neg <- subject_data_neg[remain_idx_neg,]


sample_info_neg <- 
  sample_info_neg %>% 
  filter(class == "Subject")

subject_data_neg <- 
  subject_data_neg %>% 
  dplyr::select(one_of(sample_info_neg$sample.name))

###log
subject_data_neg <- 
  log(subject_data_neg)


##combine pos and neg

subject_data_pos <- 
  subject_data_pos %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


subject_data_neg <- 
  subject_data_neg %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


colnames(subject_data_pos) == colnames(subject_data_neg)

subject_data <- 
  rbind(subject_data_pos, subject_data_neg)

sample_info <- 
  sample_info_pos %>% 
  filter(sample.name %in% colnames(subject_data))


colnames(subject_data) == sample_info$sample.name

subject_data <-
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )

subject_data <-
  subject_data %>%
  as_tibble()


subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA == 0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)










































#-------------------------------------------------------------------------------
##HILIC pos
#-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())
setwd("data_analysis20200108/data_overview/")
load("../data_cleaning/HILIC/POS/hilic_pos_6")
library(metflow2)
library(tidyverse)

sample_info <- hilic_pos_6@sample.info

subject_data <- metflow2::getData(object = hilic_pos_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(hilic_pos_6,
                             slot = "QC")


###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)


subject_data <- subject_data[remain_idx,]


sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))


###log
subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


subject_data <-
  subject_data %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)





#####HILIC negative
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/")
rm(list = ls())
load("../data_cleaning/HILIC/NEG/hilic_neg_6")
library(metflow2)
library(tidyverse)

sample_info <- hilic_neg_6@sample.info

subject_data <- metflow2::getData(object = hilic_neg_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(hilic_neg_6,
                             slot = "QC")


###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)


subject_data <- subject_data[remain_idx,]

# sample_info <-
#   sample_info %>%
#   filter(GA != 0)


sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))


###log
subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


subject_data <-
  subject_data %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)





###HILIC POS and NEG
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/")
rm(list = ls())
load("../data_cleaning/HILIC/POS/hilic_pos_6")
load("../data_cleaning/HILIC/NEG/hilic_neg_6")

library(metflow2)
library(tidyverse)

sample_info_pos <- hilic_pos_6@sample.info

subject_data_pos <- metflow2::getData(object = hilic_pos_6, 
                                      slot = "Subject")

qc_data_pos <- metflow2::getData(hilic_pos_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_pos <- apply(qc_data_pos, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_pos <- 
  which(qc_rsd_pos < 30)


subject_data_pos <- subject_data_pos[remain_idx_pos,]


sample_info_pos <- 
  sample_info_pos %>% 
  filter(class == "Subject")

subject_data_pos <- 
  subject_data_pos %>% 
  dplyr::select(one_of(sample_info_pos$sample.name))

###log
subject_data_pos <- 
  log(subject_data_pos)

# subject_data_pos <- 
#   t(
#     apply(subject_data_pos, 1, function(x){
#       (x - mean(x))/sd(x)
#     })
#   )
# 
# subject_data_pos <-
#   subject_data_pos %>% 
#   as_tibble()


sample_info_neg <- hilic_neg_6@sample.info

subject_data_neg <- metflow2::getData(object = hilic_neg_6, 
                                      slot = "Subject")

qc_data_neg <- metflow2::getData(hilic_neg_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_neg <- apply(qc_data_neg, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_neg <- 
  which(qc_rsd_neg < 30)


subject_data_neg <- subject_data_neg[remain_idx_neg,]


sample_info_neg <- 
  sample_info_neg %>% 
  filter(class == "Subject")

subject_data_neg <- 
  subject_data_neg %>% 
  dplyr::select(one_of(sample_info_neg$sample.name))

###log
subject_data_neg <- 
  log(subject_data_neg)


##combine pos and neg

subject_data_pos <- 
  subject_data_pos %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


subject_data_neg <- 
  subject_data_neg %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


colnames(subject_data_pos) == colnames(subject_data_neg)

subject_data <- 
  rbind(subject_data_pos, subject_data_neg)

sample_info <- 
  sample_info_pos %>% 
  filter(sample.name %in% colnames(subject_data))


colnames(subject_data) == sample_info$sample.name

subject_data <-
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )

subject_data <-
  subject_data %>%
  as_tibble()


subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA == 0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)







#######################################################
##use metabolite to check the overview
#-------------------------------------------------------------------------------
##HILIC pos
#-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())
setwd("data_analysis20200108/data_overview/metabolites/")
load("../../data_cleaning/HILIC/POS/hilic_pos_6")
load("../../data_preparation_for_analysis/metabolite_tags")
library(metflow2)
library(tidyverse)

sample_info <- hilic_pos_6@sample.info

subject_data <- metflow2::getData(object = hilic_pos_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(hilic_pos_6,
                             slot = "QC")


rownames(subject_data) <-
  rownames(qc_data) <-
  hilic_pos_6@ms1.data[[1]]$name

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)

name <- rownames(subject_data)[remain_idx]

subject_data <- 
  subject_data[remain_idx,]

rownames(subject_data) <- name

sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))


###log
subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


subject_data <-
  subject_data %>% 
  as_tibble()

rownames(subject_data) <- name

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)
colnames(subject_data2) <- 
  rownames(subject_data)
subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA)) %>% 
  select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA)) %>% 
  select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA)) %>% 
  select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)





#####HILIC negative
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/metabolites/")
rm(list=ls())
load("../../data_cleaning/HILIC/NEG/hilic_neg_6")
load("../../data_preparation_for_analysis/metabolite_tags")
library(metflow2)
library(tidyverse)

sample_info <- hilic_neg_6@sample.info

subject_data <- metflow2::getData(object = hilic_neg_6, 
                                  slot = "Subject")

qc_data <- metflow2::getData(hilic_neg_6,
                             slot = "QC")


rownames(subject_data) <- rownames(qc_data) <-
  hilic_neg_6@ms1.data[[1]]$name

subject_data <- as.data.frame(subject_data)

###remove peaks with large RSD
qc_rsd <- apply(qc_data, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx <- 
  which(qc_rsd < 30)


subject_data <- subject_data[remain_idx,]

# sample_info <-
#   sample_info %>%
#   filter(GA != 0)


sample_info <- 
  sample_info %>% 
  filter(class == "Subject")

subject_data <- 
  subject_data %>% 
  dplyr::select(one_of(sample_info$sample.name))


###log
subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


# subject_data <-
#   subject_data %>% 
#   as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as.data.frame(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))
# select(one_of(metabolite_tags$name))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.negition = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)





###HILIC POS and NEG
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/")
load("../data_cleaning/HILIC/POS/hilic_pos_6")
load("../data_cleaning/HILIC/NEG/hilic_neg_6")

library(metflow2)
library(tidyverse)

sample_info_pos <- hilic_pos_6@sample.info

subject_data_pos <- metflow2::getData(object = hilic_pos_6, 
                                      slot = "Subject")

qc_data_pos <- metflow2::getData(hilic_pos_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_pos <- apply(qc_data_pos, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_pos <- 
  which(qc_rsd_pos < 30)


subject_data_pos <- subject_data_pos[remain_idx_pos,]


sample_info_pos <- 
  sample_info_pos %>% 
  filter(class == "Subject")

subject_data_pos <- 
  subject_data_pos %>% 
  dplyr::select(one_of(sample_info_pos$sample.name))

###log
subject_data_pos <- 
  log(subject_data_pos)

# subject_data_pos <- 
#   t(
#     apply(subject_data_pos, 1, function(x){
#       (x - mean(x))/sd(x)
#     })
#   )
# 
# subject_data_pos <-
#   subject_data_pos %>% 
#   as_tibble()


sample_info_neg <- hilic_neg_6@sample.info

subject_data_neg <- metflow2::getData(object = hilic_neg_6, 
                                      slot = "Subject")

qc_data_neg <- metflow2::getData(hilic_neg_6,
                                 slot = "QC")


###remove peaks with large RSD
qc_rsd_neg <- apply(qc_data_neg, 1, function(x){
  sd(as.numeric(x))*100/mean(as.numeric(x))
})

remain_idx_neg <- 
  which(qc_rsd_neg < 30)


subject_data_neg <- subject_data_neg[remain_idx_neg,]


sample_info_neg <- 
  sample_info_neg %>% 
  filter(class == "Subject")

subject_data_neg <- 
  subject_data_neg %>% 
  dplyr::select(one_of(sample_info_neg$sample.name))

###log
subject_data_neg <- 
  log(subject_data_neg)


##combine pos and neg

subject_data_pos <- 
  subject_data_pos %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


subject_data_neg <- 
  subject_data_neg %>% 
  select(one_of(intersect(colnames(subject_data_pos), colnames(subject_data_neg))))


colnames(subject_data_pos) == colnames(subject_data_neg)

subject_data <- 
  rbind(subject_data_pos, subject_data_neg)

sample_info <- 
  sample_info_pos %>% 
  filter(sample.name %in% colnames(subject_data))


colnames(subject_data) == sample_info$sample.name

subject_data <-
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )

subject_data <-
  subject_data %>%
  as_tibble()


subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)


##batch1 
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 1) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 1], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) 
# annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC1[x$GA==0], colout = "black")






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


##batch2
##PCA analysis
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  filter(batch == 2) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA[sample_info$batch == 2], 
                stringsAsFactors = FALSE)

ggplot(x[x$GA!=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)



##batch1 and batch 2
temp_data <- 
  subject_data2 %>% 
  mutate(batch = sample_info$batch) %>% 
  select(-c(batch, GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA !=0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.7),
  #                       high = alpha("#FFA319FF", 0.7)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  annotate(geom = "point", x = x$PC1[x$GA == 0], y = x$PC2[x$GA == 0], 
           colour = "#8A9045FF", size = 5)


















