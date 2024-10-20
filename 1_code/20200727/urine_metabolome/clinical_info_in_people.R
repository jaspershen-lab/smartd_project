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

sample_info <- 
sample_info %>% 
mutate(GA2 = case_when(
  is.na(GA) ~ 45,
  !is.na(GA) ~ GA
)) %>%
  mutate(class = case_when(
    is.na(GA) ~ "PP",
    !is.na(GA) ~ "Normal"
  )) %>%
  mutate(Begin.Date = as.Date(EDD) - 280) %>%
  mutate(term.date = as.Date(DD) - Begin.Date) %>%
  mutate(diff_day = as.Date(EDD) - as.Date(DD)) 

ga <- sample_info$GA[match(colnames(subject_data), sample_info$sample_id)]
ga[is.na(ga)] <- 0
batch <- sample_info$batch[match(colnames(subject_data), sample_info$sample_id)]
subject_id <- sample_info$subject_id[match(colnames(subject_data), sample_info$sample_id)]
age <- sample_info$Age[match(colnames(subject_data), sample_info$sample_id)]
ethnicity <- sample_info$ethnicity[match(colnames(subject_data), sample_info$sample_id)]
parity <- sample_info$parity[match(colnames(subject_data), sample_info$sample_id)]
bmi <- sample_info$bmi[match(colnames(subject_data), sample_info$sample_id)]
diff_day <- sample_info$diff_day[match(colnames(subject_data), sample_info$sample_id)]

subject_data2 <- 
  t(subject_data) %>% 
  tibble::as_tibble()

rownames(subject_data2)

qc_data2 <- 
  t(qc_data) %>% 
  tibble::as_tibble()

rownames(subject_data2) <- colnames(subject_data)


####PCA without QC only for subjects
temp_data <-
  subject_data2

pca_object <-
  prcomp(x = temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

ga <- sample_info$GA[match(colnames(subject_data), sample_info$sample_id)]
ga[is.na(ga)] <- 0
batch <- sample_info$batch[match(colnames(subject_data), sample_info$sample_id)]
subject_id <- sample_info$subject_id[match(colnames(subject_data), sample_info$sample_id)]
age <- sample_info$Age[match(colnames(subject_data), sample_info$sample_id)]
ethnicity <- sample_info$ethnicity[match(colnames(subject_data), sample_info$sample_id)]
parity <- sample_info$parity[match(colnames(subject_data), sample_info$sample_id)]
bmi <- sample_info$bmi[match(colnames(subject_data), sample_info$sample_id)]
diff_day <- sample_info$diff_day[match(colnames(subject_data), sample_info$sample_id)]
sex <- sample_info$Sex[match(colnames(subject_data), sample_info$sample_id)]
sex[is.na(sex)] <- "Unknown"
sex[sex == "F, F"] <- "F,F"
sex[sex == "M / F"] <- "M,F"

x <- data.frame(x,
                GA = ga,
                subject_id = subject_id,
                age,
                ethnicity,
                parity,
                bmi,
                diff_day,
                sex,
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

ggsave(plot, filename = "pca_ga_plot.pdf",
       width = 7, height = 7)

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

##pca plot ethnicity
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



##pca plot diff_day
x$diff_day <- as.numeric(x$diff_day)
x <- 
x %>% 
  dplyr::mutate(preterm = 
                  case_when(
                    diff_day >= 21 ~ "Preterm 21 days",
                    diff_day >= 7 & diff_day < 21 ~ "Preterm 7 days",
                    TRUE ~ "Normal"
                  ))

color <- ggsci::pal_uchicago()(10)[c(1,3,4)]
names(color) <- c("Preterm 21 days", "Preterm 7 days", "Normal")
  
plot <-
  ggplot(x, aes(PC1, PC2, fill = preterm)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
scale_fill_manual(values = color) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank()
  ) +
  labs(
    x = paste("PC1 (", round(summary(pca_object)$importance[1, 1], 2), "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object)$importance[1, 2], 2), "%)", sep = "")
  )

plot

ggsave(plot, filename = "pca_plot_preterm.pdf",
       width = 7, height = 7)


##sex
plot <-
  ggplot(x, aes(PC1, PC2, fill = sex)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  ggsci::scale_fill_aaas() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank()
  ) +
  labs(
    x = paste("PC1 (", round(summary(pca_object)$importance[1, 1], 2), "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object)$importance[1, 2], 2), "%)", sep = "")
  )

plot

ggsave(plot, filename = "pca_plot_sex.pdf",
       width = 7, height = 7)


x2 <-
  x %>%
  dplyr::filter( GA == 0)

###GA for each ethnicity
plot <-
  ggplot(x[x$GA!=0,], aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 3, shape = 21,
             aes(fill = GA),
             show.legend = FALSE) +
  geom_point(mapping = aes(x = PC1, y = PC2),
             data = x2, fill = "#C71000FF",
             shape = 21, size = 3) +
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
        # legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  labs(x = paste("PC1 (", round(summary(pca_object)$importance[1,1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_object)$importance[1,2], 2), "%)", sep = "")) +
  facet_wrap(facets = ~ethnicity, scales = "free", ncol = 3)


plot

ggsave(plot, filename = "pca_ethnicity_ga_plot.pdf",
       width = 7, height = 7)


##tsne analysis
# tsne_object <- Rtsne::Rtsne(
#   X = as.matrix(temp_data),
#   dims = 2,
#   perplexity = 30,
#   verbose = TRUE
# )

save(tsne_object, file = "tsne_object")

Y <- tsne_object$Y

Y <-
  data.frame(Y,
             ga,
             subject_id,
             age,
             ethnicity,
             parity,
             bmi,
             diff_day,
             sex,
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


###GA for each ethnicity
Y2 <-
  Y %>%
  dplyr::filter(ga == 0)

plot <-
  ggplot(Y[Y$ga!=0,], aes(X1, X2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 3, shape = 21,
             aes(fill = ga),
             show.legend = FALSE) +
  geom_point(mapping = aes(x = X1, y = X2),
             data = Y2, fill = "#C71000FF",
             shape = 21, size = 3) +
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
        # legend.position = c(1,1), legend.justification = c(1,1),
        legend.background = element_blank()) +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") +
  facet_wrap(facets = ~ethnicity, scales = "free", ncol = 3)


plot

ggsave(plot, filename = "pca_ethnicity_ga_plot.pdf",
       width = 7, height = 7)



##pca plot diff_day
Y$diff_day <- as.numeric(Y$diff_day)
Y <-
  Y %>%
  dplyr::mutate(
    preterm =
      case_when(
        diff_day >= 21 ~ "Preterm 21 days",
        diff_day >= 7 &
          diff_day < 21 ~ "Preterm 7 days",
        TRUE ~ "Normal"
      )
  )

color <- ggsci::pal_uchicago()(9)[c(1,3,4)]
names(color) <- c("Preterm 21 days", "Preterm 7 days", "Normal")

plot <-
  ggplot(Y, aes(X1, X2, fill = preterm)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank()
  ) +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2")

plot

ggsave(plot, filename = "tsne_plot_preterm.pdf",
       width = 7, height = 7)



##sex
plot <-
  ggplot(Y, aes(X1, X2, fill = sex)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5, shape = 21) +
  ggsci::scale_fill_aaas() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank()
  ) +
  labs(x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") 

plot

ggsave(plot, filename = "tsne_plot_sex.pdf",
       width = 7, height = 7)
