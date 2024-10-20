#to avoind source
no_exist_function()
#-------------------------------------------------------------------------------
# 
##RPLC pos
#-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())

##peaks
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info")

setwd("data_analysis20200108/urine_metabolome/data_overview/")

idx_pos <-
grep("POS", variable_info$name)

expression_data <- expression_data[idx_pos,]

variable_info <- variable_info[idx_pos,]

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

x <- 
  x %>% 
  rownames_to_column(var = "name")

library(ggfortify)

plot <- 
  autoplot(
    object = pca_object,
    data = x,
    fill = "batch",
    size = 5, 
    shape = 21,
    alpha = 1,
    frame = TRUE, 
    frame.type = 'norm',
    variance_percentage = TRUE,
    frame.colour  = "batch",
    scale = 0
    # loadings = TRUE, 
    # loadings.label = TRUE
  ) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  guides(fill = guide_legend(title = "Batch")) +
  ggsci::scale_fill_d3() +
  ggsci::scale_color_d3() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  ) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_data_quality_plot_pos.pdf",
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
  ggplot(x[x$GA!=0,], aes(PC1, PC2, fill = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  guides(fill = guide_colourbar(title = "GA (week)")) +
  scale_fill_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
           fill = "#C71000FF", size = 5, shape = 21) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))


plot

ggsave(plot, filename = "pca_pos_plot.pdf",
       width = 7, height = 7)


##RPLC neg
###-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())

##peaks
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/expression_data")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/sample_info")
load("data_analysis20200108/urine_metabolome/data_preparation_for_analysis/peaks/variable_info")

setwd("data_analysis20200108/urine_metabolome/data_overview/")

idx_neg <-
  grep("NEG", variable_info$name)

expression_data <- expression_data[idx_neg,]

variable_info <- variable_info[idx_neg,]

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

x <- 
  x %>% 
  rownames_to_column(var = "name")

library(ggfortify)

plot <- 
  autoplot(
    object = pca_object,
    data = x,
    fill = "batch",
    size = 5, 
    shape = 21,
    alpha = 1,
    frame = TRUE, 
    frame.type = 'norm',
    variance_percentage = TRUE,
    frame.colour  = "batch",
    scale = 0
    # loadings = TRUE, 
    # loadings.label = TRUE
  ) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  guides(fill = guide_legend(title = "Batch")) +
  ggsci::scale_fill_d3() +
  ggsci::scale_color_d3() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  ) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_data_quality_plot_neg.pdf",
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
  ggplot(x[x$GA!=0,], aes(PC1, PC2, fill = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  guides(fill = guide_colourbar(title = "GA (week)")) +
  scale_fill_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
           fill = "#C71000FF", size = 5, shape = 21) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))


plot

ggsave(plot, filename = "pca_neg_plot.pdf",
       width = 7, height = 7)




##-----------------------------------------------------------------------------
##RPLC pos and neg
###----------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())

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

setwd("data_analysis20200108/urine_metabolome/data_overview/")

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

####subject and QC together for data quality
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

library(ggfortify)

plot <- 
autoplot(
  object = pca_object,
  data = x,
  fill = "batch",
  size = 5, 
  shape = 21,
  alpha = 1,
  frame = TRUE, 
  frame.type = 'norm',
  variance_percentage = TRUE,
  frame.colour  = "batch",
  scale = 0
  # loadings = TRUE, 
  # loadings.label = TRUE
) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  guides(fill = guide_legend(title = "Batch")) +
  ggsci::scale_fill_d3() +
  ggsci::scale_color_d3() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  ) +
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

subject_id <- sample_info$subject_id[match(colnames(subject_data), sample_info$sample_id)]
age <- sample_info$Age[match(colnames(subject_data), sample_info$sample_id)]
ethnicity <- sample_info$ethnicity[match(colnames(subject_data), sample_info$sample_id)]
parity <- sample_info$parity[match(colnames(subject_data), sample_info$sample_id)]
bmi <- sample_info$bmi[match(colnames(subject_data), sample_info$sample_id)]


x <- data.frame(x, 
                GA = ga,
                subject_id = subject_id,
                age,
                ethnicity,
                parity,
                bmi,
                stringsAsFactors = FALSE)

plot <- 
  ggplot(x[x$GA!=0,], aes(PC1, PC2, fill = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  guides(fill = guide_colourbar(title = "GA (week)")) +
  scale_fill_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
           fill = "#C71000FF", size = 5, shape = 21) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))


plot

# plot2 <-
#   x %>% 
#   ggplot(aes(y = PC2, x = GA, fill = GA)) +
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
#         legend.background = element_blank())


ggsave(plot, filename = "pca_plot_ga.pdf", 
       width = 7, height = 7)


##pca plot subject_id
color <-
  colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(subject_id)))

plot <- 
  ggplot(x, aes(PC1, PC2, fill = subject_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.background = element_blank()) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
  #          colour = "#C71000FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_plot_subject.pdf", 
       width = 9, height = 7)


##pca plot age
plot <-
  ggplot(x, aes(PC1, PC2, fill = age)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  guides(fill = guide_colourbar(title = "Age")) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank()) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
  #          colour = "#C71000FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_plot_age.pdf", 
       width = 7, height = 7)

##pca plot ethinicty
color <- 
  ggsci::pal_futurama()(7)[1:7]

names(color) <- sort(unique(ethnicity))

plot <- 
  ggplot(x, aes(PC1, PC2, fill = ethnicity)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank()) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
  #          colour = "#C71000FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_plot_ethnicity.pdf", 
       width = 7, height = 7)

##pca plot parity
plot <- 
  ggplot(x, aes(PC1, PC2, fill = parity)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  scale_fill_gradient(low = alpha("red", 0.3), 
                        high = alpha("red", 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank()) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
  #          colour = "#C71000FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_plot_parity.pdf", 
       width = 7, height = 7)


##pca plot bmi
plot <- 
  ggplot(x, aes(PC1, PC2, fill = bmi)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  scale_fill_gradient(low = alpha("blue", 0.3), high = "blue") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank()) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0], 
  #          colour = "#C71000FF", size = 5) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = ""))

plot

ggsave(plot, filename = "pca_plot_bmi.pdf", 
       width = 7, height = 7)


####for each participant
x <-
  data.frame(sample_id = colnames(subject_data),
             x, stringsAsFactors = FALSE)

# x <-
#   x %>%
#   left_join(sample_info[,c("sample_id", "subject_id")],
#             by = c("sample_id"))

x2 <-
  x %>%
  mutate(y = 0) %>%
  dplyr::filter( GA == 0)


plot <-
x %>%
  dplyr::filter(GA != 0) %>%
  mutate(y = 0) %>%
  ggplot(aes(x = PC2, y = y, fill = GA)) +
  geom_point(shape = 21) +
  geom_point(mapping = aes(x = PC2, y = y),
             data = x2, fill = "#C71000FF", shape = 21) +
  # annotate(geom = "point", x = x$PC1[x$GA==0], y = 0,
  #          colour = "#C71000FF") +
  ggrepel::geom_text_repel(aes(PC2, y = y, label = round(GA, 2)),
                           size = 2.5) +
  scale_fill_gradientn(colours = c(
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


##tsne analysis
tsne_object <- Rtsne::Rtsne(
  X = as.matrix(temp_data),
  dims = 2,
  perplexity = 30,
  verbose = TRUE
)

Y <- tsne_object$Y
Y <-
  data.frame(Y, 
             ga,
             subject_id,
             age,
             ethnicity,
             parity,
             bmi,
             stringsAsFactors = FALSE)

####t-sne for ga
plot <- 
  ggplot(Y[Y$ga!=0,], aes(X1, X2, fill = ga)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  guides(fill = guide_colourbar(title = "GA (week)")) +
  scale_fill_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  annotate(geom = "point", x = Y$X1[Y$ga==0], y = Y$X2[Y$ga == 0], 
           fill = "#C71000FF", size = 5, shape = 21) +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2")

plot

ggsave(plot, filename = "tsne_plot_ga.pdf", 
       width = 7, height = 7)


####t-sne for subject_id
color <- 
colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(Y$subject_id)))

names(color) <- unique(Y$subject_id)

plot <- Y %>% 
  ggplot(aes(X1, X2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(fill = subject_id), size = 5,
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_manual(values =color) +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  )

plot

ggsave(plot, filename = "tsne_plot_subject_id.pdf", width = 9, height = 7)



####t-sne for age
color <- 
  colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(Y$subject_id)))

names(color) <- unique(Y$subject_id)

plot <- Y %>% 
  ggplot(aes(X1, X2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(fill = age), size = 5,
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  )

plot

ggsave(plot, filename = "tsne_plot_age.pdf", width = 9, height = 7)


####t-sne for ethnicity
color <- 
  ggsci::pal_futurama()(7)[1:7]

names(color) <- sort(unique(ethnicity))


library(ggfortify)

plot <- Y %>% 
  ggplot(aes(X1, X2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(fill = ethnicity), size = 5,
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_manual(values = color) +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") +
  theme_bw() +
  theme(
    legend.position = c(0,1),
    legend.justification = c(0,1),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  )

plot

ggsave(plot, filename = "tsne_plot_ethnicity.pdf", width = 7, height = 7)









