setwd(r4projects::get_project_wd())
setwd("microbiome/")

library(tidyverse)

sample_info <- 
  readr::read_csv("metadata_for_Jasper.csv")

sample_info <- 
  sample_info %>% 
  dplyr::select(
    subject_id = pt_ID,
    delivery_ga = `Pregnancy period`,
    ga = Sample_GA,
    pp = Trimester
  ) %>% 
  dplyr::mutate(
    ga = as.numeric(ga),
    delivery_ga = as.numeric(delivery_ga)
  ) %>% 
  dplyr::mutate(
    class = case_when(
      pp == "Postpartum" ~ "After",
      TRUE ~ "Before"
    )
  ) %>% 
  dplyr::arrange(desc(delivery_ga)) %>% 
  dplyr::mutate(subject_id = factor(subject_id, levels = unique(subject_id)))
  
sample_info$ga[is.na(sample_info$ga)] <- 45

temp_data <- 
  sample_info
  
plot1 <- 
  temp_data %>% 
  ggplot(aes(x = ga, y = subject_id)) +
  # geom_vline(xintercept = 40, color = "black") +
  geom_point(aes(color = class), shape = 16, show.legend = FALSE) +
  scale_color_manual(values = c("Before" = "black", "After" = "#FB8072")) +
  geom_point(aes(x = delivery_ga, y = subject_id), shape = 17, color = "red") +
  # geom_rect(
  #   aes(
  #     xmin = 12,
  #     xmax = 26,
  #     ymin = 0,
  #     ymax = Inf
  #   ),
  #   data = data.frame(),
  #   fill = "#8DD3C7",
  #   alpha = 0.2,
  #   inherit.aes = FALSE) +
  labs(y = "Subject ID", x = "Gestational age (GA, weeks)") +
  scale_x_continuous(breaks = c(10,20,30,40,45), label = c(10,20,30,40,"PP")) +
  # scale_y_continuous(limits = c(-10, 50)) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)
  )

plot1


plot2 <- 
  temp_data %>% 
  ggplot(aes(x = ga)) +
  geom_histogram(binwidth = 1, 
                 fill = alpha(ggsci::pal_aaas()(10)[1], 0.8), 
                 color = ggsci::pal_aaas()(10)[1]) +
  geom_density(alpha = 1, aes(y = ..count..),
               fill = NA, 
               color = ggsci::pal_aaas()(10)[1]) +
  labs(y = "Sample number", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        # panel.grid.major.y = element_line(color = subject_line_color),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)
  )


plot2


plot3 <- 
  temp_data %>% 
  ggplot(aes(x = subject_id)) +
  geom_bar(width = 0.8, fill = alpha(ggsci::pal_aaas()(10)[1], 0.8)) +
  labs(y = "Sample number", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        # panel.grid.major.y = element_line(color = subject_line_color),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  coord_flip()


plot3

empty_space <- 
  temp_data %>% 
  ggplot(aes(x = subject_id)) +
  geom_bar(width = 0.8, fill = alpha(ggsci::pal_aaas()(10)[1], 0.8)) +
  labs(y = "", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # panel.grid.major.y = element_line(color = subject_line_color),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
  ) +
  coord_flip()

empty_space


library(patchwork)

plot <-
  {
    plot2 + empty_space + plot_layout(ncol = 2, widths = c(4, 1))
  } -
  {
    plot1 + plot3 + plot_layout(ncol = 2, widths = c(4, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 4))

plot

ggsave(plot, file = "sample_collection_light.pdf", width = 9, height = 7)





