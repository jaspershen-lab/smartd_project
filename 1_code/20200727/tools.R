plotLambdaVSdeviation <-
  function(object,
           xlab = "Log lambda",
           ylab = "Deviation ratio (%)") {
    data.frame(
      lambda = object$lambda,
      dev.ratio = object$dev.ratio,
      stringsAsFactors = FALSE
    ) %>%
      ggplot(aes(log(lambda), dev.ratio * 100)) +
      labs(x = xlab, y = ylab) +
      geom_point(size = 2) +
      geom_line() +
      theme_bw() +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 13))
  }


plotLambdaVScoefficients <-
  function(object,
           xlab = "Log lambda",
           ylab = "Coefficients") {
    beta <-
      object$beta %>%
      as.matrix() %>%
      t() %>%
      as_tibble() %>%
      mutate(lambda = object$lambda) %>%
      tidyr::gather(., key = "feature", value = "coef", -lambda)
    
    label_index <- seq(range(log(beta$lambda))[1],
                       range(log(beta$lambda))[2],
                       by = 1)
    label <-
      lapply(label_index, function(x) {
        c(log(object$lambda) - x) %>%
          abs() %>%
          which.min() %>%
          `[`(object$df + 1, .)
      }) %>%
      unlist()
    
    
    beta %>%
      ggplot(., aes(log(lambda), coef)) +
      geom_line(aes(colour = feature), show.legend = FALSE) +
      scale_x_continuous(
        position = "bottom",
        sec.axis = sec_axis(
          ~ .,
          name = "",
          breaks = label_index,
          labels = label
        )
      ) +
      scale_colour_manual(values = colorRampPalette(ggsci::pal_uchicago()(5))(600)) +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 13))
  }


plotLambdaVSerror <- function(object,
                              xlab = "Log lambda",
                              ylab = "Mean absolute error") {
  cvm <-
    data.frame(
      lambda = object$lambda,
      df = object$glmnet.fit$df,
      cvm = object$cvm,
      cvup = object$cvup,
      cvlo = object$cvlo,
      stringsAsFactors = FALSE
    )
  
  cvm %>%
    ggplot(., aes(log(lambda), cvm)) +
    geom_vline(xintercept = log(c(object$lambda.min,
                                  object$lambda.1se)),
               linetype = 2) +
    geom_errorbar(aes(ymin = cvlo, ymax = cvup), colour = "#155F83FF") +
    geom_point(size = 2, colour = "#FFA319FF") +
    scale_x_continuous(
      position = "bottom",
      sec.axis = sec_axis(
        trans = ~ .,
        breaks = log(cvm$lambda)[seq(1, 100, by = 7)],
        labels = cvm$df[seq(1, 100, by = 7)],
        name = ""
      )
    ) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13)
          # plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"))
          )
}


# plotSVMtunePerformance <-
#   function(object){
#     performance <-
#       object$performances
#
#     object$performances %>%
#     ggplot2::ggplot(aes(gamma, cost, colour = error)) +
#       geom_point()
#
#   }



setwd_project <- function() {
  currect_wd <-
    getwd()
  
  candidate_wd <-
    currect_wd %>%
    stringr::str_split("/") %>%
    unlist()
  
  if (length(candidate_wd) == 1) {
    candidate_wd <- currect_wd
  } else {
    candidate_wd <-
      lapply(2:length(candidate_wd), function(i) {
        paste(candidate_wd[1:i], collapse = "/")
      })
  }
  
  candidate_wd <-
    rev(candidate_wd)
  
  for (i in 1:length(candidate_wd)) {
    wd <- candidate_wd[[i]]
    file_name <-
      list.files(wd,
                 recursive = ifelse(wd == currect_wd, TRUE, FALSE),
                 full.names = TRUE)
    project_index <-
      grep(".Rproj", file_name)
    
    if (length(project_index) != 0) {
      project_wd <-
        file_name[project_index[1]] %>%
        stringr::str_split("/") %>%
        unlist() %>%
        head(-1) %>%
        paste(collapse = "/")
      cat(
        "The project name is:",
        file_name[project_index[1]] %>%
          stringr::str_split("/") %>%
          unlist() %>%
          tail(1),
        "\n"
      )
      cat("The project wd is:",
          project_wd, "\n")
      
      setwd(project_wd)
      break()
    }
  }
  
  if (length(project_index) == 0) {
    cat("There are no .Rproj in your file. No change for wd.\n")
  }
}


trans_ht <- function(x) {
  x <-
    stringr::str_replace(x, "cm", "")
  x_inch <- x[grep("'", x)]
  if (length(x_inch) > 0) {
    x_inch <-
      stringr::str_split(x_inch, "'") %>%
      lapply(function(y) {
        y[2] <- stringr::str_replace(y[2], '"', "")
        (as.numeric(y[1]) * 12 + as.numeric(y[2])) * 2.54
      }) %>%
      unlist()
    x[grep("'", x)] <-
      x_inch
  }
  x <- as.numeric(x)
  x
}


trans_wt <- function(x) {
  sapply(x, function(y) {
    if (is.na(y)) {
      return(NA)
    }
    if (stringr::str_detect(y, "kg")) {
      y <- stringr::str_replace(y, "kg", "") %>%
        as.numeric()
    } else {
      y <- as.numeric(y) * 0.4535922921969
    }
    y
  }) %>%
    unname()
}


cal_r2 <- function(predicted, true) {
  1 - sum((predicted - true) ^ 2) / sum((mean(true) - true) ^ 2)
}

sxt_rank <- function(x) {
  data.frame(
    rank = c(1:length(x)),
    order = order(x),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(order) %>%
    data.frame(x, ., stringsAsFactors = FALSE) %>%
    pull(rank)
}



volcano_plot <- function(fc,
                         p_value,
                         text,
                         p.cutoff = 0.05,
                         log2_fc = TRUE,
                         fc.cutoff = 2,
                         top = 5,
                         size_range = c(0.3, 5),
                         theme = c("light", "dark")) {
  theme <- match.arg(theme)
  if (log2_fc) {
    temp_data <- data.frame(
      fc = log(fc, 2),
      p_value = -log(p_value, 10),
      text,
      stringsAsFactors = FALSE
    )
  } else{
    temp_data <- data.frame(
      fc = fc,
      p_value = -log(p_value, 10),
      text,
      stringsAsFactors = FALSE
    )
  }
  
  temp_data <-
    temp_data %>%
    dplyr::mutate(class = case_when(
      p_value > -log(p.cutoff, 10) & fc > log(fc.cutoff, 2) ~ "Increase",
      p_value > -log(p.cutoff, 10) &
        fc < log(1 / fc.cutoff, 2) ~ "Decrease",
      TRUE ~ "No"
    ))
  

  top_text_up <-
    temp_data %>% 
      dplyr::arrange(desc(fc), desc(p_value)) %>% 
      head(top)
    
  top_text_down <-
    temp_data %>% 
      dplyr::arrange(desc(fc), p_value) %>% 
      tail(top)
  
  temp_data <-
  temp_data %>% 
    mutate(text = ifelse(text %in% c(top_text_up$text, top_text_down$text), text, ""))
  
  plot <-
    temp_data %>%
    ggplot(aes(fc, p_value)) +
    geom_hline(
      yintercept = -log(p.cutoff, 10),
      linetype = 2
    ) +
    geom_point(
      shape = 21,
      aes(fill = class,
          size = p_value),
      show.legend = FALSE,
      alpha = 0.5,
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_size_continuous(range = size_range) +
    scale_fill_manual(values = c("Increase" = ggsci::pal_aaas()(10)[2],
                                 "Decrease" = ggsci::pal_aaas()(10)[1],
                                 "No" = "grey")) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 13)
    ) +
    labs(x = "log2(Fold change)",
         y = "-log10(p value, FDR)") +
    ggrepel::geom_text_repel(
      data = temp_data,
      aes(label = text),
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey",
      segment.size = 0.5,
      segment.alpha = 0.5,
      size = 3
    )

  if (theme == "light") {
    plot <- plot + theme_classic()  +
      theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 13)
      )
  } else{
    plot <- plot + ggdark::dark_theme_classic() +
      theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 13)
      )
  }
}


###----------------------------------------------------------------------------
readPIUMet <-
  function(path = ".",
           variable_info,
           fc_p_table,
           text = TRUE,
           layout = "kk",
           size_range = c(3, 8),
           width_range = c(0.2, 0.5),
           marker_name = NULL) {
    annotation_result <-
      read.table(
        file.path(
          path,
          "peaks_putative_metabolites_w10.0_b2.0_mu0.0005_R1.txt"
        ),
        sep = "\t",
        header = TRUE
      )
    
    if (nrow(annotation_result) == 0) {
      cat("No result.\n")
      return(NULL)
    }
    
    annotation_result <-
      annotation_result %>%
      dplyr::mutate(mz =
                      stringr::str_replace(mz.peak, "m/z=", "") %>%
                      as.numeric() %>%
                      round(4)) %>%
      dplyr::mutate(mz2 = as.character(mz))
    
    if (!is.null(marker_name)) {
      marker <-
        variable_info %>%
        dplyr::filter(name %in% marker_name) %>%
        dplyr::mutate(polarity = case_when(
          stringr::str_detect(name, "POS") ~ "positive",
          stringr::str_detect(name, "NEG") ~ "negative"
        )) %>%
        dplyr::mutate(
          mz2 = case_when(
            polarity == "positive" ~ as.character(round(mz, 4) - 1),
            polarity == "negative" ~ as.character(round(mz, 4) + 1)
          )
        )
    } else{
      idx <- which(fc_p_table[, 1] < 0.05)
      marker <-
        variable_info[idx, ] %>%
        dplyr::mutate(polarity = case_when(
          stringr::str_detect(name, "POS") ~ "positive",
          stringr::str_detect(name, "NEG") ~ "negative"
        )) %>%
        dplyr::mutate(
          mz2 = case_when(
            polarity == "positive" ~ as.character(round(mz, 4) - 1),
            polarity == "negative" ~ as.character(round(mz, 4) + 1)
          )
        )
    }
    
    
    
    annotation_result <-
      annotation_result %>%
      dplyr::left_join(marker, by = "mz2") %>%
      dplyr::select(-c(mz.y, HMDB.ID.y)) %>%
      dplyr::rename(mz = mz.x, HMDB.ID = HMDB.ID.x)
    
    edge_attr <-
      read.table(
        file.path(path, "result_edge_frequency_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(edge = V1)
    
    edge_data <-
      read.table(
        file.path(path, "result_union_net_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(from = V1, to = V2) %>%
      dplyr::mutate(edge = paste(from, "(pp)", to, sep = " ")) %>%
      dplyr::left_join(edge_attr, by = "edge")
    
    node_data <-
      read.table(
        file.path(path, "result_node_frequency_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(node = V1,
                    node_class = V3,
                    HMDB_ID = V4) %>%
      dplyr::left_join(annotation_result[, c("name", "Metabolite.Name", "super.class")],
                       by = c("node" = "Metabolite.Name"))
    
    node <-
      node_data$node %>%
      stringr::str_replace("m/z=", "") %>%
      as.numeric() %>%
      round(4) %>%
      as.character()
    
    node <- marker$name[match(node, marker$mz2)]
    
    rename <-
      data.frame(name1 = node_data$node[!is.na(node)],
                 name2 = node[!is.na(node)])
    
    node_data$node <-
      sapply(node_data$node, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    edge_data$from <-
      sapply(edge_data$from, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    edge_data$to <-
      sapply(edge_data$to, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    
    edge_data$edge <-
      paste(edge_data$from, "(pp)", edge_data$to, sep = " ")
    
    node_data <-
      node_data %>%
      dplyr::select(-name) %>%
      dplyr::distinct()
    
    node_data$node_class[grep("Metabolite", node_data$node_class)] <-
      "Metabolite"
    
    graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(Degree = tidygraph::centrality_degree(mode = 'all'))
    
    fill <-
      c(
        "m/z Peak" = ggsci::pal_aaas()(10)[9],
        "Metabolite" = ggsci::pal_aaas()(10)[1],
        "Protein" = ggsci::pal_aaas(alpha = 0.1)(10)[2]
        # "Protein" = "tomato"
      )
    
    col <-
      c(
        "m/z Peak" = ggsci::pal_aaas()(10)[9],
        "Metabolite" = ggsci::pal_aaas()(10)[1],
        "Protein" = ggsci::pal_aaas()(10)[2]
      )
    
    shape = c("m/z Peak" = 24,
              "Metabolite" = 22,
              # "Metabolite_others" = 22,
              "Protein" = 21)
    
    require(ggraph)
    if (text) {
      plot <-
        ggraph(graph,
               layout = layout) +
        geom_edge_link(
          aes(edge_width = V3),
          alpha = 1,
          color = "black",
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(
            size = V2,
            fill = node_class,
            shape = node_class
          ),
          alpha = 1,
          show.legend = TRUE
        ) +
        scale_shape_manual(values = shape) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggraph::geom_node_text(aes(label = node,
                                   color = node_class),
                               repel = TRUE,
                               size = 3) +
        ggraph::scale_edge_width(range = width_range) +
        scale_size_continuous(range = size_range) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = col) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
    } else{
      plot <-
        ggraph(graph,
               layout = layout) +
        geom_edge_link(
          aes(edge_width = V3),
          alpha = 1,
          color = "black",
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(
            size = V2,
            fill = node_class,
            shape = node_class
          ),
          alpha = 1,
          show.legend = TRUE
        ) +
        scale_shape_manual(values = shape) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggraph::scale_edge_width(range = width_range) +
        scale_size_continuous(range = size_range) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = col) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
    }
    
    
    
    output_path <- file.path(path, "Result")
    dir.create(output_path)
    save(edge_data, file = file.path(output_path, "edge_data"))
    save(node_data, file = file.path(output_path, "node_data"))
    save(graph, file = file.path(output_path, "graph"))
    save(annotation_result, file = file.path(output_path, "annotation_result"))
    
    ggsave(
      plot,
      filename = file.path(output_path, "graph_plog.pdf"),
      width = 7,
      height = 7
    )
    
    plot
  }


plot_silhouette <- function(sil,
                            color = "red") {
  temp_data <-
    data.frame(
      cluster = sil[, 1],
      neighbor = sil[, 2],
      sil_width = sil[, 3],
      stringsAsFactors = FALSE
    )
  temp_data <-
    temp_data %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::arrange(desc(cluster), sil_width) %>%
    dplyr::mutate(index = 1:nrow(temp_data))
  
  plot <-
    temp_data %>%
    ggplot() +
    geom_bar(
      aes(
        x = sil_width,
        y = index,
        fill = cluster,
        color = cluster
      ),
      stat = "identity",
      show.legend = FALSE
    ) +
    # geom_segment(aes(x = 0, y = index, xend = sil_width, yend = index,
    #                  color = cluster)) +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    scale_x_continuous(expand = expansion(mult = c(0, .2))) +
    theme_classic() +
    labs(
      x = paste(
        "Silhousette width",
        "\nAverage silhousettle widht:",
        round(mean(temp_data$sil_width), 2)
      ),
      y = ""
    ) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank()
    )
  plot <-
    plot +
    ggplot2::annotate(
      geom = "text",
      x = 0,
      y = max(temp_data$index),
      label = paste("n =", nrow(temp_data)),
      color = "black",
      hjust = 0,
      vjust = -1,
      size = 4
    )
  
  cluster_num <- as.numeric(max(temp_data$cluster))
  title <-
    paste(cluster_num, "clusters Cj", "\n", "j: nj | avei<Cj Si")
  plot <-
    plot +
    ggplot2::annotate(
      geom = "text",
      x = max(temp_data$sil_width),
      y = max(temp_data$index),
      label = title,
      color = "black",
      hjust = 0,
      vjust = 0,
      size = 4
    )
  
  class <- temp_data$cluster %>%
    unique() %>%
    sort() %>%
    rev()
  
  cluster_num <-
    temp_data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(cluster)) %>%
    dplyr::pull(n)
  
  cluster_num <-
    sapply(1:length(cluster_num), function(x) {
      if (x == 1) {
        cluster_num[x] / 2
      } else{
        tail(cumsum(cluster_num[1:(x - 1)]), 1)  +  cluster_num[x] / 2
      }
    })
  
  
  for (i in 1:length(class)) {
    label <-
      paste(class[i],
            ":",
            sum(temp_data$cluster == class[i]),
            "|",
            round(mean(temp_data$sil_width[temp_data$cluster == class[i]]), 2))
    plot <-
      plot +
      ggplot2::annotate(
        geom = "text",
        x = max(temp_data$sil_width),
        y = cluster_num[i],
        label = label,
        color = "black",
        hjust = 0,
        vjust = 0,
        size = 4
      )
  }
  
  plot
}