##avoid source
no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
source("1_code/100_tools.R")

##load data
load(
  "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/peaks/urine_metabolomics_data.rda"
)

dir.create("3_data_analysis/3_biomakers/1_urine_metabolomics_total_markers",
           recursive = TRUE)
setwd("3_data_analysis/3_biomakers/1_urine_metabolomics_total_markers")

##remove qc and blank samples
urine_metabolomics_data2 <-
  urine_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject")

expression_data <-
  extract_expression_data(urine_metabolomics_data2)

sample_info <-
  extract_sample_info(urine_metabolomics_data2)

variable_info <-
  extract_variable_info(urine_metabolomics_data2)

dim(expression_data)

###find marker which are change according to pregnancy

#log transformation
subject_data <-
  log(expression_data + 1, 10)

subject_data[5, ] %>%
  as.numeric() %>%
  density() %>%
  plot()

expression_data[5, ] %>%
  as.numeric() %>%
  density() %>%
  plot()

dim(subject_data)
dim(variable_info)

sample_info$GA[is.na(sample_info$GA)] <- 45

###remove samples which are after birth
remove_idx <- which(sample_info$GA == 45)
# remove_name <- sample_info$sample_id[remove_idx]

#SAM analysi
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

peak_marker <- c(gene_up[, 1], gene_down[, 1])

plot(sample_info$GA, subject_data[peak_marker[1], ])

plot(sample_info$GA, expression_data[peak_marker[1], ])

plot(sample_info$GA, expression_data[gene_down[4, 1], ])

###linear regression
subject_data <-
  subject_data[peak_marker, ]

rownames(subject_data) <- peak_marker

colnames(subject_data) == sample_info$sample_id

ga <- sample_info$GA
batch <- sample_info$batch
# bmi <- sample_info$bmi
bmi <- sample_info$mother_bmi
# age <- sample_info$Age
age <- sample_info$mother_age
parity <- sample_info$mother_parity
ethnicity <- sample_info$mother_ethnicity
# ethnicity <- sample_info$ethnicity

remove_idx <- which(ga == 45)

# lm_p_value <-
#   apply(subject_data[, -remove_idx], 1, function(x) {
#     x <- as.numeric(x)
#     lm_reg <- lm(formula = ga[-remove_idx] ~ x + batch[-remove_idx] +
#                    bmi[-remove_idx] + age[-remove_idx] + parity[-remove_idx] +
#                    ethnicity[-remove_idx])
#     temp <- summary(lm_reg)$coefficients  %>% as.data.frame()
#     as.numeric(temp["x", 4])
#   })
#
# lm_fdr <- p.adjust(lm_p_value, method = "fdr")
# save(lm_fdr, file = "lm_fdr.rda")
load("lm_fdr.rda")

###correlation

# cor_value <-
#   purrr::map(.x = as.data.frame(t(subject_data[,-remove_idx])), .f = function(x){
#     temp1 <- cor.test(x, sample_info$GA[-remove_idx], method = "spearman")
#     temp2 <- cor.test(sxt_rank(x), sample_info$GA[-remove_idx], method = "spearman")
#     c(temp1$estimate, temp1$p.value, temp2$estimate, temp2$p.value)
#   })
#
# cor_value <-
#   do.call(rbind, cor_value)
#
# colnames(cor_value) <- c("correlation1", "p1", "correlation2", "p2")
#
# cor_value <- data.frame(name = rownames(subject_data),
#                         cor_value, stringsAsFactors = FALSE)
#
#
# cor_value$p1_adj <- p.adjust(p = cor_value$p1, method = "fdr")
#
# cor_value$p2_adj <- p.adjust(p = cor_value$p2, method = "fdr")
#
# plot(cor_value$correlation1)
#
# plot(cor_value$correlation2)
#
# save(cor_value, file = "cor_value")

load("cor_value")

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
load("peak_marker")

plot <- peak_marker %>%
  ggplot(aes(score, correlation1)) +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_point(aes(color = class, size = -log(p1_adj, 10)),
             show.legend = TRUE,
             alpha = 0.8) +
  labs(x = "Score (SAM test)", y = "Correlation (Spearman)") +
  scale_color_manual(values = c(
    "up" = ggsci::pal_aaas()(10)[2],
    "down" = ggsci::pal_aaas()(10)[1],
    "no" = "#D9D9D9"
  )) +
  scale_size_continuous(range = c(0.05, 5)) +
  guides(size = guide_legend(override.aes = list(color = "black"))) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  )

plot

ggsave(
  plot,
  filename = "dem_plot.pdf",
  width = 7,
  height = 7,
  bg = "transparent"
)

table(peak_marker$class)

variable_info <-
  variable_info[match(rownames(subject_data), variable_info$variable_id), ]


###volcano plot
lm_fdr <-
  lm_fdr %>%
  data.frame()

colnames(lm_fdr) <- "fdr"

lm_fdr <-
  lm_fdr %>%
  tibble::rownames_to_column(var = "name")

# sam_value <- rbind(gene_up, gene_down) %>%
#   as.data.frame() %>%
#   dplyr::select(`Gene ID`, `Score(d)`) %>%
#   dplyr::rename(name = `Gene ID`, score = `Score(d)`) %>%
#   dplyr::left_join(lm_fdr, by = "name")
# 
# save(sam_value, file = "sam_value.rda")
load("sam_value.rda")

plot <-
  volcano_plot(
    fc = as.numeric(sam_value$score),
    log2_fc = FALSE,
    p_value = sam_value$fdr,
    p.cutoff = 0.05,
    fc.cutoff = 1,
    size_range = c(0.1, 5)
  )
plot
sum(lm_fdr < 0.05 & as.numeric(sam_value$score) > 0)
sum(lm_fdr < 0.05 & as.numeric(sam_value$score) < 0)

plot <-
  plot +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw()

plot
ggsave(plot,
       filename = "volcano_plot.pdf",
       width = 7,
       height = 7)
