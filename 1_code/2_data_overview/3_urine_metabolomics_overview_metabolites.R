# ####only use the peaks in the study
# 
# #to avoind source
# no_exist_function()
# library(tidymass)
# library(tidyverse)
# library(r4projects)
# setwd(get_project_wd())
# rm(list = ls())
# source('1_code/100_tools.R')
# 
# ###quality assessment
# #-------------------------------------------------------------------------------
# ##metabolites
# load(
#   "3_data_analysis/1_data_preparation/1_urine_metabolomics_data/metabolites/urine_metabolomics_data.rda"
# )
# 
# dir.create(
#   "3_data_analysis/2_data_overview/2_urine_metabolomics_overview_metabolites",
#   recursive = TRUE
# )
# setwd("3_data_analysis/2_data_overview/2_urine_metabolomics_overview_metabolites")
# 
# ####positive mode
# urine_metabolomics_data_pos <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(stringr::str_detect(variable_id, "POS"))
# 
# ###remove metabolites with large RSD
# ####remove some qc samples
# urine_metabolomics_data_pos <-
#   urine_metabolomics_data_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% c("QC2.1", "QC2.2", "QC2.3"))
# 
# qc_sample_id <-
#   urine_metabolomics_data_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "QC") %>%
#   pull(sample_id)
# 
# ###add RSD
# urine_metabolomics_data_pos <-
#   urine_metabolomics_data_pos %>%
#   mutate_rsd(according_to_samples = qc_sample_id)
# 
# urine_metabolomics_data_pos <-
#   urine_metabolomics_data_pos %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(rsd < 30)
# 
# dim(urine_metabolomics_data_pos)
# 
# urine_metabolomics_data_pos@sample_info
# 
# #####batch effect
# pca_object <-
#   urine_metabolomics_data_pos %>%
#   scale_data() %>%
#   run_pca()
# 
# temp_data <-
#   urine_metabolomics_data_pos
# 
# temp_data@sample_info$batch[temp_data@sample_info$class == "QC"] <- "QC"
# 
# plot <-
#   massstat::pca_score_plot(object = temp_data,
#                            pca_object = pca_object,
#                            color_by = "batch") +
#   scale_fill_manual(values = batch_color)
# 
# plot
# 
# ggsave(plot,
#        filename = "pca_data_quality_plot_pos.pdf",
#        width = 8,
#        height = 7)
# 
# #####GA
# qc_id <-
#   urine_metabolomics_data_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "QC") %>%
#   pull(sample_id)
# 
# temp_data <-
#   urine_metabolomics_data_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% qc_id)
# 
# pca_object <-
#   temp_data %>%
#   scale_data() %>%
#   run_pca()
# 
# plot <-
#   massstat::pca_score_plot(
#     object = temp_data,
#     pca_object = pca_object,
#     color_by = "GA",
#     frame = FALSE
#   ) +
#   geom_point(size = 5, shape = 21) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   ),
#   na.value = "#C71000FF")
# 
# plot
# 
# ggsave(plot,
#        filename = "pca_ga_pos.pdf",
#        width = 8,
#        height = 7)
# 
# 
# ####negative mode
# urine_metabolomics_data_neg <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(stringr::str_detect(variable_id, "NEG"))
# 
# ###remove metabolites with large RSD
# ####remove some qc samples
# urine_metabolomics_data_neg <-
#   urine_metabolomics_data_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% c("QC2.1", "QC2.2", "QC2.3"))
# 
# qc_sample_id <-
#   urine_metabolomics_data_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "QC") %>%
#   pull(sample_id)
# 
# ###add RSD
# urine_metabolomics_data_neg <-
#   urine_metabolomics_data_neg %>%
#   mutate_rsd(according_to_samples = qc_sample_id)
# 
# urine_metabolomics_data_neg <-
#   urine_metabolomics_data_neg %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(rsd < 30)
# 
# dim(urine_metabolomics_data_neg)
# 
# urine_metabolomics_data_neg@sample_info
# 
# #####batch effect
# pca_object <-
#   urine_metabolomics_data_neg %>%
#   scale_data() %>%
#   run_pca()
# 
# temp_data <-
#   urine_metabolomics_data_neg
# 
# temp_data@sample_info$batch[temp_data@sample_info$class == "QC"] <- "QC"
# 
# plot <-
#   massstat::pca_score_plot(object = temp_data,
#                            pca_object = pca_object,
#                            color_by = "batch") +
#   scale_fill_manual(values = batch_color)
# 
# plot
# 
# ggsave(plot,
#        filename = "pca_data_quality_plot_neg.pdf",
#        width = 8,
#        height = 7)
# 
# #####GA
# qc_id <-
#   urine_metabolomics_data_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "QC") %>%
#   pull(sample_id)
# 
# temp_data <-
#   urine_metabolomics_data_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% qc_id)
# 
# pca_object <-
#   temp_data %>%
#   scale_data() %>%
#   run_pca()
# 
# plot <-
#   massstat::pca_score_plot(
#     object = temp_data,
#     pca_object = pca_object,
#     color_by = "GA",
#     frame = FALSE
#   ) +
#   geom_point(size = 5, shape = 21) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   ),
#   na.value = "#C71000FF")
# 
# plot
# 
# ggsave(plot,
#        filename = "pca_ga_neg.pdf",
#        width = 8,
#        height = 7)
# 
# 
# ####positive and negative mode
# ###remove metabolites with large RSD
# ####remove some qc samples
# urine_metabolomics_data <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% c("QC2.1", "QC2.2", "QC2.3"))
# 
# qc_sample_id <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "QC") %>%
#   pull(sample_id)
# 
# ###add RSD
# urine_metabolomics_data <-
#   urine_metabolomics_data %>%
#   mutate_rsd(according_to_samples = qc_sample_id)
# 
# urine_metabolomics_data <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(rsd < 30)
# 
# dim(urine_metabolomics_data)
# 
# urine_metabolomics_data@sample_info
# 
# #####batch effect
# pca_object <-
#   urine_metabolomics_data %>%
#   scale_data() %>%
#   run_pca()
# 
# temp_data <-
#   urine_metabolomics_data
# 
# temp_data@sample_info$batch[temp_data@sample_info$class == "QC"] <- "QC"
# 
# plot <-
#   massstat::pca_score_plot(object = temp_data,
#                            pca_object = pca_object,
#                            color_by = "batch") +
#   scale_fill_manual(values = batch_color)
# 
# plot
# 
# ggsave(plot,
#        filename = "pca_data_quality_plot.pdf",
#        width = 8,
#        height = 7)
# 
# #####GA
# qc_id <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "QC") %>%
#   pull(sample_id)
# 
# temp_data <-
#   urine_metabolomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% qc_id)
# 
# pca_object <-
#   temp_data %>%
#   scale_data() %>%
#   run_pca()
# 
# plot <-
#   massstat::pca_score_plot(
#     object = temp_data,
#     pca_object = pca_object,
#     color_by = "GA",
#     frame = FALSE
#   ) +
#   geom_point(size = 5, shape = 21) +
#   guides(fill = guide_colourbar(title = "GA (week)")) +
#   scale_fill_gradientn(colours = c(
#     alpha("#155F83FF", 1),
#     alpha("#155F83FF", 0.4),
#     alpha("#FFA319FF", 0.4),
#     alpha("#FFA319FF", 1)
#   ),
#   na.value = "#C71000FF")
# 
# plot
# 
# ggsave(plot,
#        filename = "pca_ga.pdf",
#        width = 8,
#        height = 7)
# 
# 
# 
# #
# #
# # ##pca plot age
# # plot <-
# #   ggplot(x, aes(PC1, PC2, fill = age)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(size = 5, shape = 21) +
# #   guides(fill = guide_colourbar(title = "Age")) +
# #   scale_fill_gradient(low = "blue", high = "red") +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     legend.position = c(1, 1),
# #     legend.justification = c(1, 1),
# #     legend.background = element_blank()
# #   ) +
# #   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
# #   #          colour = "#C71000FF", size = 5) +
# #   labs(
# #     x = paste("PC1 (", round(summary(pca_object)$importance[1, 1], 2), "%)", sep = ""),
# #     y = paste("PC2 (", round(summary(pca_object)$importance[1, 2], 2), "%)", sep = "")
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "pca_plot_age.pdf",
# #        width = 7,
# #        height = 7)
# #
# # ##pca plot ethinicty
# # color <-
# #   ggsci::pal_futurama()(7)[1:7]
# #
# # names(color) <- sort(unique(ethnicity))
# #
# # plot <-
# #   ggplot(x, aes(PC1, PC2, fill = ethnicity)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(size = 5, shape = 21) +
# #   scale_fill_manual(values = color) +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     legend.position = c(0, 1),
# #     legend.justification = c(0, 1),
# #     legend.background = element_blank()
# #   ) +
# #   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
# #   #          colour = "#C71000FF", size = 5) +
# #   labs(
# #     x = paste("PC1 (", round(summary(pca_object)$importance[1, 1], 2), "%)", sep = ""),
# #     y = paste("PC2 (", round(summary(pca_object)$importance[1, 2], 2), "%)", sep = "")
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "pca_plot_ethnicity.pdf",
# #        width = 7,
# #        height = 7)
# #
# # ##pca plot parity
# # plot <-
# #   ggplot(x, aes(PC1, PC2, fill = parity)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(size = 5, shape = 21) +
# #   scale_fill_gradient(low = alpha("red", 0.3), high = alpha("red", 1)) +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     legend.position = c(1, 1),
# #     legend.justification = c(1, 1),
# #     legend.background = element_blank()
# #   ) +
# #   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
# #   #          colour = "#C71000FF", size = 5) +
# #   labs(
# #     x = paste("PC1 (", round(summary(pca_object)$importance[1, 1], 2), "%)", sep = ""),
# #     y = paste("PC2 (", round(summary(pca_object)$importance[1, 2], 2), "%)", sep = "")
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "pca_plot_parity.pdf",
# #        width = 7,
# #        height = 7)
# #
# #
# # ##pca plot bmi
# # plot <-
# #   ggplot(x, aes(PC1, PC2, fill = bmi)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(size = 5, shape = 21) +
# #   scale_fill_gradient(low = alpha("blue", 0.3), high = "blue") +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     legend.position = c(1, 1),
# #     legend.justification = c(1, 1),
# #     legend.background = element_blank()
# #   ) +
# #   # annotate(geom = "point", x = x$PC1[x$GA==0], y = x$PC2[x$GA == 0],
# #   #          colour = "#C71000FF", size = 5) +
# #   labs(
# #     x = paste("PC1 (", round(summary(pca_object)$importance[1, 1], 2), "%)", sep = ""),
# #     y = paste("PC2 (", round(summary(pca_object)$importance[1, 2], 2), "%)", sep = "")
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "pca_plot_bmi.pdf",
# #        width = 7,
# #        height = 7)
# #
# #
# # ####for each participant
# # x <-
# #   data.frame(sample_id = colnames(subject_data),
# #              x,
# #              stringsAsFactors = FALSE)
# #
# # # x <-
# # #   x %>%
# # #   left_join(sample_info[,c("sample_id", "subject_id")],
# # #             by = c("sample_id"))
# #
# # x2 <-
# #   x %>%
# #   mutate(y = 0) %>%
# #   dplyr::filter(GA == 0)
# #
# #
# # plot <-
# #   x %>%
# #   dplyr::filter(GA != 0) %>%
# #   mutate(y = 0) %>%
# #   ggplot(aes(x = PC2, y = y, fill = GA)) +
# #   geom_point(shape = 21) +
# #   geom_point(
# #     mapping = aes(x = PC2, y = y),
# #     data = x2,
# #     fill = "#C71000FF",
# #     shape = 21
# #   ) +
# #   # annotate(geom = "point", x = x$PC1[x$GA==0], y = 0,
# #   #          colour = "#C71000FF") +
# #   ggrepel::geom_text_repel(aes(PC2, y = y, label = round(GA, 2)), size = 2.5) +
# #   scale_fill_gradientn(colours = c(
# #     alpha("#155F83FF", 1),
# #     alpha("#155F83FF", 0.4),
# #     alpha("#FFA319FF", 0.4),
# #     alpha("#FFA319FF", 1)
# #   )) +
# #   labs(y = "") +
# #   facet_wrap( ~ subject_id, ncol = 6) +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 8),
# #     axis.text = element_text(size = 8),
# #     legend.title = element_text(size = 8),
# #     legend.text = element_text(size = 8),
# #     legend.position = "top",
# #     legend.background = element_blank(),
# #     strip.background = element_rect(fill = "grey"),
# #     strip.text = element_text(color = "white"),
# #     axis.text.y = element_blank(),
# #     panel.grid.major.y = element_blank(),
# #     panel.grid.minor.y = element_blank(),
# #     axis.ticks.y = element_blank()
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "pca_for_each_person.pdf",
# #        width = 10,
# #        height = 8)
# #
# #
# # ##tsne analysis
# # tsne_object <- Rtsne::Rtsne(
# #   X = as.matrix(temp_data),
# #   dims = 2,
# #   perplexity = 30,
# #   verbose = TRUE
# # )
# #
# # Y <- tsne_object$Y
# # Y <-
# #   data.frame(Y, ga, subject_id, age, ethnicity, parity, bmi, stringsAsFactors = FALSE)
# #
# # ####t-sne for ga
# # plot <-
# #   ggplot(Y[Y$ga != 0, ], aes(X1, X2, fill = ga)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(size = 5, shape = 21) +
# #   guides(fill = guide_colourbar(title = "GA (week)")) +
# #   scale_fill_gradientn(colours = c(
# #     alpha("#155F83FF", 1),
# #     alpha("#155F83FF", 0.4),
# #     alpha("#FFA319FF", 0.4),
# #     alpha("#FFA319FF", 1)
# #   )) +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     legend.position = c(1, 1),
# #     legend.justification = c(1, 1),
# #     legend.background = element_blank()
# #   ) +
# #   annotate(
# #     geom = "point",
# #     x = Y$X1[Y$ga == 0],
# #     y = Y$X2[Y$ga == 0],
# #     fill = "#C71000FF",
# #     size = 5,
# #     shape = 21
# #   ) +
# #   labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2")
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "tsne_plot_ga.pdf",
# #        width = 7,
# #        height = 7)
# #
# #
# # ####t-sne for subject_id
# # color <-
# #   colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(Y$subject_id)))
# #
# # names(color) <- unique(Y$subject_id)
# #
# # plot <- Y %>%
# #   ggplot(aes(X1, X2)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(
# #     aes(fill = subject_id),
# #     size = 5,
# #     shape = 21,
# #     color = "black",
# #     alpha = 0.8
# #   ) +
# #   scale_fill_manual(values = color) +
# #   labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     strip.background = element_rect(fill = "#0099B47F"),
# #     strip.text = element_text(color = "white", size = 13)
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "tsne_plot_subject_id.pdf",
# #        width = 9,
# #        height = 7)
# #
# #
# #
# # ####t-sne for age
# # color <-
# #   colorRampPalette(colors = ggsci::pal_futurama(alpha = 0.5)(10))(length(unique(Y$subject_id)))
# #
# # names(color) <- unique(Y$subject_id)
# #
# # plot <- Y %>%
# #   ggplot(aes(X1, X2)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(
# #     aes(fill = age),
# #     size = 5,
# #     shape = 21,
# #     color = "black",
# #     alpha = 0.8
# #   ) +
# #   scale_fill_gradient(low = "blue", high = "red") +
# #   labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     strip.background = element_rect(fill = "#0099B47F"),
# #     strip.text = element_text(color = "white", size = 13)
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "tsne_plot_age.pdf",
# #        width = 9,
# #        height = 7)
# #
# #
# # ####t-sne for ethnicity
# # color <-
# #   ggsci::pal_futurama()(7)[1:7]
# #
# # names(color) <- sort(unique(ethnicity))
# #
# #
# # library(ggfortify)
# #
# # plot <- Y %>%
# #   ggplot(aes(X1, X2)) +
# #   geom_vline(xintercept = 0, linetype = 2) +
# #   geom_hline(yintercept = 0, linetype = 2) +
# #   geom_point(
# #     aes(fill = ethnicity),
# #     size = 5,
# #     shape = 21,
# #     color = "black",
# #     alpha = 0.8
# #   ) +
# #   scale_fill_manual(values = color) +
# #   labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
# #   theme_bw() +
# #   theme(
# #     legend.position = c(0, 1),
# #     legend.justification = c(0, 1),
# #     axis.title = element_text(size = 13),
# #     axis.text = element_text(size = 12),
# #     legend.title = element_text(size = 13),
# #     legend.text = element_text(size = 12),
# #     strip.background = element_rect(fill = "#0099B47F"),
# #     strip.text = element_text(color = "white", size = 13)
# #   )
# #
# # plot
# #
# # ggsave(plot,
# #        filename = "tsne_plot_ethnicity.pdf",
# #        width = 7,
# #        height = 7)
