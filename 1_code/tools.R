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
      theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)
      )
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
