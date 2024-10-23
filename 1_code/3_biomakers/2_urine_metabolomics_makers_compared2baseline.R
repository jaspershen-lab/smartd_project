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

dir.create(
  "3_data_analysis/3_biomakers/2_urine_metabolomics_makers_compared2baseline",
  recursive = TRUE
)
setwd("3_data_analysis/3_biomakers/2_urine_metabolomics_makers_compared2baseline")

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

sample_info$GA[is.na(sample_info$GA)] <- 45

###remove lm_fdr > 0.05
text_colour <-
  c(colorRampPalette(colors = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  ))(13), "red")

colnames(subject_data) == sample_info$sample_id

##ga range
library(plyr)
temp_data <-
  sample_info  %>%
  plyr::dlply(.(ga_range)) %>%
  lapply(function(x) {
    c(length(unique(x$subject_id)), length(unique(x$sample_id)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

library(ggradar)

colnames(temp_data) <- c("subject", "sample")

plot <-
  temp_data %>%
  tibble::rownames_to_column(var = "ga_range") %>%
  tidyr::pivot_longer(cols = -ga_range,
                      names_to = "class",
                      values_to = "number") %>%
  ggplot(aes(ga_range, number)) +
  geom_line(aes(group = class, color = class), size = 1) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_point(aes(group = class, fill = class),
             size = 5,
             shape = 21) +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm(alpha = 0.8) +
  scale_y_continuous(limits = c(5, 30)) +
  theme_bw() +
  labs(x = "", y = "#Subject/Sample") +
  theme(
    axis.text.x = element_text(size = 12, color = text_colour),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.border = element_blank()
  ) +
  coord_polar(start = -0.2)
plot
ggsave(plot,
       filename = "subject_sample_number_in_ga_range.pdf",
       width = 8,
       height = 7)


###combine different samples in one ga range together
library(plyr)

subject_data <-
  apply(subject_data, 1, function(x) {
    (x) / sd(x)
  })

subject_data2 <-
  subject_data %>%
  data.frame(., ga_range = sample_info$ga_range, stringsAsFactors = FALSE) %>%
  mutate(ga_range = factor(
    ga_range,
    levels = sample_info$ga_range %>% unique() %>% stringr::str_sort(numeric = TRUE)
  )) %>%
  plyr::dlply(.variables = .(ga_range))

# subject_data_mean <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, mean)
#   }) %>%
#   do.call(cbind, . )
#
# subject_data_sd <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, sd)
#   }) %>%
#   do.call(cbind, . )
#
# subject_data_sem <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, function(y){
#       sd(y)/(nrow(x) - 1)
#     })
#   }) %>%
#   do.call(cbind, . )
#
# save(subject_data_mean, file = "subject_data_mean")
# save(subject_data_sd, file = "subject_data_sd")
# save(subject_data_sem, file = "subject_data_sem")

load("subject_data_mean")
load("subject_data_sd")
load("subject_data_sem")

subject_data2 <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      dplyr::select(-ga_range)
  })

## find all the peaks in different time points
# fc_p_value <-
#   pbapply::pblapply(subject_data2[-1], function(x){
#     p_value <- lapply(1:ncol(x), function(idx){
#       t.test(x[,idx], subject_data2[[1]][,idx])$p.value
#     }) %>%
#       unlist() %>%
#       p.adjust(method = "fdr")
#
#     fc <- lapply(1:ncol(x), function(idx){
#       mean(x[,idx]) /mean(subject_data2[[1]][,idx])
#     }) %>%
#       unlist()
#
#     fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
#
#     data.frame(p_value, fc, stringsAsFactors = FALSE)
#   })
#
# save(fc_p_value, file = "fc_p_value")

load("fc_p_value")

names(fc_p_value)

fc_p_value$`(14,16]`

dir.create("marker_in_different_points")

# for(idx in 1:length(fc_p_value)){
#   cat(idx, " ")
#   plot <-
#     volcano_plot(fc = fc_p_value[[idx]][,2],
#                  p_value = fc_p_value[[idx]][,1],
#                  p.cutoff = 0.05, fc.cutoff = 1,
#                  theme = "light") +
#     scale_x_continuous(limits = c(-5, 5)) +
#     scale_y_continuous(limits = c(0, 15))
#
#   plot <-
#     plot +
#     labs(title = paste(names(fc_p_value)[idx], "(4,12]", sep = "/"))
#
#   plot
#
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_light.pdf", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
# }


##find markers for each time points
marker_each_point <-
  lapply(fc_p_value, function(x) {
    idx1 <- which(x$p_value < 0.05 & x$fc > 1)
    idx2 <- which(x$p_value < 0.05 & x$fc < 1)
    
    gene1 <-
      try(data.frame(
        gene_id = variable_info$variable_id[idx1],
        x[idx1, ],
        class = "increase",
        stringsAsFactors = FALSE
      ),
      silent = TRUE)
    
    if (class(gene1) == "try-error") {
      gene1 <- NULL
    }
    
    gene2 <-
      try(data.frame(
        gene_id = variable_info$variable_id[idx2],
        x[idx2, ],
        class = "decrease",
        stringsAsFactors = FALSE
      ),
      silent = TRUE)
    
    if (class(gene2) == "try-error") {
      gene2 <- NULL
    }
    
    rbind(gene1, gene2)
  })


marker_each_point[[10]]

names(marker_each_point)


##increase
temp_data_increase <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    peak_name1 <-
      x %>%
      dplyr::filter(class == "increase") %>%
      dplyr::pull(gene_id)
    
    if (length(peak_name1) == 0) {
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    
    lapply(marker_each_point, function(y) {
      if (is.null(y)) {
        return(NA)
      }
      
      peak_name2 <-
        y %>%
        dplyr::filter(class == "increase") %>%
        dplyr::pull(gene_id)
      
      if (length(peak_name2) == 0) {
        return(NA)
      }
      
      length(intersect(peak_name1, peak_name2)) / length(peak_name1)
      
    }) %>% unlist()
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data_increase <- temp_data_increase[-13, -13]

##increase
temp_data_decrease <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    peak_name1 <-
      x %>%
      dplyr::filter(class == "decrease") %>%
      dplyr::pull(gene_id)
    
    if (length(peak_name1) == 0) {
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    
    lapply(marker_each_point, function(y) {
      if (is.null(y)) {
        return(NA)
      }
      
      peak_name2 <-
        y %>%
        dplyr::filter(class == "decrease") %>%
        dplyr::pull(gene_id)
      
      if (length(peak_name2) == 0) {
        return(NA)
      }
      
      length(intersect(peak_name1, peak_name2)) / length(peak_name1)
      
    }) %>% unlist()
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    peak_name1_increase <-
      x %>%
      dplyr::filter(class == "increase") %>%
      dplyr::pull(gene_id)
    
    peak_name1_decrease <-
      x %>%
      dplyr::filter(class == "decrease") %>%
      dplyr::pull(gene_id)
    
    if (length(peak_name1_decrease) == 0 &
        length(peak_name1_increase) == 0) {
      x <- rep(NA, length(marker_each_point))
      names(x) <- names(marker_each_point)
      return(x)
    }
    
    lapply(marker_each_point, function(y) {
      if (is.null(y)) {
        return(NA)
      }
      
      peak_name2_decrease <-
        y %>%
        dplyr::filter(class == "decrease") %>%
        dplyr::pull(gene_id)
      
      peak_name2_increase <-
        y %>%
        dplyr::filter(class == "increase") %>%
        dplyr::pull(gene_id)
      
      if (length(peak_name2_decrease) == 0 &
          length(peak_name2_increase) == 0) {
        return(NA)
      }
      
      (length(intersect(
        peak_name1_decrease, peak_name2_decrease
      )) +
          length(intersect(
            peak_name1_increase, peak_name2_increase
          ))) / (length(peak_name1_decrease) + length(peak_name1_increase))
      
    }) %>% unlist()
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


temp_data_increase <- temp_data_increase[-13, -13]

temp_data_decrease <- temp_data_decrease[-13, -13]

temp_data <- temp_data[-13, -13]

# temp_data_increase[which(is.na(temp_data_increase), arr.ind = TRUE)] <- 0

library(corrplot)

corrplot(
  corr = as.matrix(temp_data_increase),
  method = "circle",
  type = "upper",
  is.corr = FALSE,
  tl.col = text_colour,
  diag = FALSE,
  na.label = "NA",
  addCoef.col = "grey",
  col =   colorRampPalette(colors = c(
    ggsci::pal_aaas()(10)[5], "white", ggsci::pal_aaas()(10)[6]
  ))(100)
)

corrplot(
  corr = as.matrix(temp_data_decrease),
  method = "circle",
  type = "upper",
  is.corr = FALSE,
  tl.col = text_colour,
  diag = FALSE,
  na.label = "NA",
  addCoef.col = "grey",
  col =   colorRampPalette(colors = c(
    ggsci::pal_aaas()(10)[6], "white", ggsci::pal_aaas()(10)[5]
  ))(100)
)


corrplot(
  corr = as.matrix(temp_data),
  method = "circle",
  type = "upper",
  is.corr = FALSE,
  tl.col = text_colour,
  diag = FALSE,
  na.label = "NA",
  addCoef.col = "grey"
)


# temp_data <-
#   lapply(1:length(marker_each_point), function(i){
#     x <- marker_each_point[[i]]
#     if(is.null(x)){
#       return(NA)
#     }
#     peak_name1_increase <-
#       x %>%
#       dplyr::filter(class == "increase") %>%
#       dplyr::pull(gene_id)
#
#     peak_name1_decrease <-
#       x %>%
#       dplyr::filter(class == "decrease") %>%
#       dplyr::pull(gene_id)
#
#     if(length(peak_name1_decrease) == 0 & length(peak_name1_increase) == 0){
#       return(NA)
#     }
#
#     peak_name2_decrease <-
#     lapply(marker_each_point[(i+1):(length(marker_each_point) - 1)], function(y){
#       if(is.null(y)){
#         return(NULL)
#       }
#         y %>%
#         dplyr::filter(class == "decrease") %>%
#         dplyr::pull(gene_id)
#     }) %>% unlist() %>%
#       unname()
#
#     peak_name2_increase <-
#       lapply(marker_each_point[(i+1):(length(marker_each_point) - 1)], function(y){
#         if(is.null(y)){
#           return(NULL)
#         }
#         y %>%
#           dplyr::filter(class == "increase") %>%
#           dplyr::pull(gene_id)
#       }) %>% unlist() %>%
#       unname()
#
#     decrease <-
#     lapply(peak_name1_decrease, function(z){
#       sum(z == peak_name2_decrease)
#     }) %>%
#       unlist()
#
#     increase <-
#     lapply(peak_name1_increase, function(z){
#       sum(z == peak_name2_increase)
#     }) %>%
#       unlist()
#
#     all <- c(decrease, increase)
#     all >= length((i+1):(length(marker_each_point) - 1)) - 1
#   })
#
# temp_data <- unlist(temp_data)
# temp_data <- temp_data[!is.na(temp_data)]
# sum(unlist(temp_data), na.rm = TRUE)/length(temp_data)


#####a sankey
marker_each_point %>%
  lapply(nrow) %>%
  unlist()

all_marker_name <-
  lapply(marker_each_point, function(x) {
    x$gene_id
  }) %>%
  unlist() %>%
  unique()

library(ggalluvial)

temp_data <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x <-
      data.frame(gene_id = all_marker_name, stringsAsFactors = FALSE) %>%
      left_join(x, by = "gene_id") %>%
      dplyr::select(gene_id, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
    
  })

temp_data <-
  purrr::map2(
    .x = temp_data,
    .y = names(temp_data),
    .f = function(x, y) {
      if (is.null(x)) {
        return(NULL)
      }
      data.frame(x, point = y, stringsAsFactors = FALSE)
    }
  )

temp_data <-
  do.call(rbind, temp_data)

temp_data %>%
  dplyr::group_by(point, class) %>%
  dplyr::summarise(n = n())

plot1 <-
  ggplot(
    temp_data,
    aes(
      x = point,
      y = freq,
      stratum = class,
      alluvium = gene_id,
      fill = class,
      label = class
    )
  ) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(
    values = c(
      "increase" = ggsci::pal_aaas()(10)[6],
      "decrease" = ggsci::pal_aaas()(10)[5],
      "no" = "grey"
    )
  ) +
  ggalluvial::geom_stratum(alpha = 1) +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 2
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot1

ggsave(
  plot1,
  file = file.path("marker_in_different_points", "gene_sankey_light.pdf"),
  width = 14,
  height = 7,
  bg = "transparent"
)

plot2 <-
  marker_each_point %>%
  lapply(function(x) {
    if (is.null(x)) {
      c(0, 0)
    } else{
      x$class %>% table() %>% as.numeric()
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  rownames_to_column(var = "point") %>%
  mutate(decrease = -V1, increase = V2) %>%
  mutate(point = factor(point, levels = point)) %>%
  ggplot(aes(point)) +
  geom_bar(aes(point, increase),
           stat = "identity",
           fill = ggsci::pal_aaas()(10)[6]) +
  geom_bar(aes(point, decrease),
           stat = "identity",
           fill = ggsci::pal_aaas()(10)[5]) +
  theme_bw() +
  labs(x = "", y = "Peak number") +
  theme(
    legend.position = "top",
    # panel.border = element_blank(),
    # panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    # axis.text.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot2

ggsave(
  plot2,
  file = file.path("marker_in_different_points", "barplot.pdf"),
  width = 14,
  height = 7,
  bg = "transparent"
)


##heatmap
all_marker_name
subject_data_mean
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c(ggsci::pal_aaas()(10)[5], "white", ggsci::pal_aaas()(10)[6]))
temp_data <-
  subject_data_mean[all_marker_name, ]

temp_data <-
  apply(temp_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

range(temp_data)

plot3 <-
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    border = TRUE,
    col = col_fun,
    row_km = 2
  )

plot3

library(ggplotify)

plot3 <- as.ggplot(plot3)
plot3
ggsave(
  plot3,
  file = file.path("marker_in_different_points", "heatmap.pdf"),
  width = 14,
  height = 7,
  bg = "transparent"
)


###upset plot show overlap of markers
marker_each_point_up <-
  purrr::map(
    .x = marker_each_point,
    # .y = names(marker_each_point),
    .f = function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      x <-
        x %>%
        dplyr::filter(class == "increase")
      
      if (nrow(x) == 0) {
        return(NULL)
      }
      x$gene_id
      # data.frame(name = x$gene_id, time = y, stringsAsFactors = FALSE)
    }
  )
# do.call(rbind, .)

remove_idx <- lapply(marker_each_point_up, is.null) %>%
  unlist() %>%
  which()

if (length(remove_idx) > 0) {
  marker_each_point_up <- marker_each_point_up[-remove_idx]
}

marker_each_point_down <-
  purrr::map(
    .x = marker_each_point,
    # .y = names(marker_each_point),
    .f = function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      x <-
        x %>%
        dplyr::filter(class == "decrease")
      
      if (nrow(x) == 0) {
        return(NULL)
      }
      x$gene_id
      # data.frame(name = x$gene_id, time = y, stringsAsFactors = FALSE)
    }
  )
# do.call(rbind, .)

remove_idx <- lapply(marker_each_point_down, is.null) %>%
  unlist() %>%
  which()

if (length(remove_idx) > 0) {
  marker_each_point_down <- marker_each_point_down[-remove_idx]
}

library(UpSetR)
library(ComplexHeatmap)

text_colour <- c(colorRampPalette(colors = c(
  alpha("#155F83FF", 1),
  alpha("#155F83FF", 0.4),
  alpha("#FFA319FF", 0.4),
  alpha("#FFA319FF", 1)
))(13), "red")

m = make_comb_mat(marker_each_point_up)

plot <-
  UpSet(
    m,
    set_order = names(marker_each_point_up),
    # comb_order = order(comb_size(m)),
    bg_col = c(ggsci::pal_aaas(alpha = 0.5)(10)[6], "white"),
    right_annotation = upset_right_annotation(
      m,
      # ylim = c(0, 30),
      gp = gpar(fill = text_colour[4:14]),
      annotation_name_side = "bottom",
      axis_param = list(side = "bottom")
    )
  )

plot <- ggplotify::as.ggplot(plot = plot)
plot
ggsave(plot,
       filename = "up_marker_upset_plot.pdf",
       width = 14,
       height = 10)

m = make_comb_mat(marker_each_point_down)

plot <-
  UpSet(
    m,
    set_order = names(marker_each_point_down),
    # comb_order = order(comb_size(m)),
    bg_col = c(ggsci::pal_aaas(alpha = 0.5)(10)[5], "white"),
    right_annotation = upset_right_annotation(
      m,
      # ylim = c(0, 30),
      gp = gpar(fill = text_colour[3:14]),
      annotation_name_side = "bottom",
      axis_param = list(side = "bottom")
    )
  )

plot <- ggplotify::as.ggplot(plot = plot)
plot

ggsave(plot,
       filename = "down_marker_upset_plot.pdf",
       width = 14,
       height = 10)
