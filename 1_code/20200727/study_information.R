#to avoind source
no_exist_function()

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(tidyverse)

# load data
load(
  "3_data_analysis/data_analysis20200108/urine_metabolome/sample_information/sample_information"
)

setwd("3_data_analysis/data_analysis20200108/study_information/")

##### the clinic information
sample_information$Patient_ID %>%
  unique() %>%
  length()

sample_information %>%
  dplyr::filter(is.na(GA)) %>%
  dplyr::select(Sample_ID)

library(ggExtra)

temp_data <-
  sample_information %>%
  mutate(GA2 = case_when(is.na(GA) ~ 45, !is.na(GA) ~ GA)) %>%
  mutate(class = case_when(is.na(GA) ~ "PP", !is.na(GA) ~ "Normal")) %>%
  mutate(Begin.Date = as.Date(EDD) - 280) %>%
  mutate(term.date = as.Date(DD) - Begin.Date) %>%
  mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  arrange(diff_day)


## which person is preterm
preterm21 <-
  temp_data %>%
  # mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::filter(diff_day >= 21) %>%
  dplyr::pull(Patient_ID) %>%
  unique()

preterm7 <-
  temp_data %>%
  # mutate(diff_day = as.Date(EDD) - as.Date(DD)) %>%
  dplyr::filter(diff_day >= 7 & diff_day < 21) %>%
  dplyr::pull(Patient_ID) %>%
  unique()

text_colour <- unique(temp_data$Patient_ID)
text_colour <- case_when(
  text_colour %in% preterm21 ~ "#58593FFF",
  text_colour %in% preterm7 ~ "#800000FF",
  TRUE ~ "#BEBEBE"
)

text_face <- unique(temp_data$Patient_ID)
text_face <- case_when(text_face %in% preterm21 ~ "plain",
                       text_face %in% preterm7 ~ "plain",
                       TRUE ~ "plain")


plot1 <-
  # temp_data %>%
  ggplot(data = temp_data, aes(x = GA2, y = factor(Patient_ID, level = unique(Patient_ID)))) +
  geom_point(
    aes(
      x = GA2,
      y = factor(Patient_ID, level = unique(Patient_ID)),
      colour = class
    ),
    data = temp_data,
    show.legend = FALSE,
    size = 2
  ) +
  scale_colour_manual(values = c("PP" = "#155F83FF", "Normal" = "#FFA319FF")) +
  labs(x = "Gestational age (GA, weeks)", y = "Subject ID") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(colour = "#8A9045FF", linetype = 1),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    axis.title = element_text(size = 13),
    axis.text.y = element_text(
      size = 12,
      colour = text_colour,
      face = text_face
    ),
    axis.text.x = element_text(size = 12)
  )

term_date <-
  temp_data %>%
  select(Patient_ID, term.date) %>%
  mutate(term.date = as.numeric(term.date / 7)) %>%
  distinct()

plot1 <-
  plot1 +
  ggplot2::annotate(
    geom = "point",
    shape = 17,
    colour = "black",
    size = 2,
    x = term_date$term.date,
    y = term_date$Patient_ID
  )

plot1

plot2 <-
  ggplot(temp_data, aes(x = GA2)) +
  geom_histogram(
    binwidth = 0.5,
    colour = "#8DD3C7",
    fill = "#8DD3C7"
  ) +
  labs(x = "GA (weeks)", y = "Sample number") +
  theme_classic() +
  scale_x_continuous(# limits = c(10, 46),
    name = NULL,
    labels = NULL,
    breaks = NULL) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  theme(
    # panel.border = element_blank(),
    # axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )

plot2

plot3 <-
  ggplot(temp_data, aes(x = factor(Patient_ID, levels = unique(Patient_ID)))) +
  geom_bar(width = 0.8, fill = "#8DD3C7") +
  labs(x = "GA (weeks)", y = "Sample number") +
  theme_bw() +
  scale_x_discrete(name = NULL,
                   label = NULL,
                   breaks = NULL) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  coord_flip() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )

plot3

space_plot <-
  ggplot(temp_data, aes(x = factor(Patient_ID, levels = unique(Patient_ID)))) +
  geom_bar(width = 0.8, fill = "#8A9045FF") +
  labs(x = "", y = "") +
  theme_classic() +
  scale_x_discrete(name = NULL,
                   label = NULL,
                   breaks = NULL) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  coord_flip() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

space_plot

library(patchwork)

plot <-
  {
    plot2 + space_plot + plot_layout(ncol = 2, widths = c(4, 1))
  } -
  {
    plot1 + plot3 + plot_layout(ncol = 2, widths = c(4, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 4))


plot


# ggsave(plot,
#        filename = "sample_colection_distribution.pdf",
#        width = 10, height = 7
# )

plot



#### a circos plot to show the information
info <-
  readxl::read_xlsx("../../patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")

info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste("SF", ., sep = "")


info <-
  info %>%
  dplyr::filter(ID %in% temp_data$Patient_ID)

info <-
  info[match(unique(temp_data$Patient_ID), info$ID), ]

info$ID == unique(temp_data$Patient_ID)

## we have
colnames(info)

###out put information for each person
range((as.numeric(info$Age)))
IQR((as.numeric(info$Age)))

birth_wr <- info$`Birth wt`
birth_wr[33] <- 3015 + 3170
birth_wr[34] <- 2315 + 2190
birth_wr[36] <- 1430 + 1630
birth_wr <- as.numeric(birth_wr)

a <- mean(birth_wr, na.rm = TRUE)
s <- sd(birth_wr, na.rm = TRUE)
n <- sum(!is.na(birth_wr))
xbar <- 6185
plot(density(birth_wr, na.rm = TRUE))
2 * (1 - pnorm(xbar, mean = a, sd = s / sqrt(n)))

range((as.numeric(birth_wr)), na.rm = TRUE)
IQR(as.numeric(birth_wr), na.rm = TRUE)
mean(as.numeric(birth_wr), na.rm = TRUE)

bmi <-
  trans_wt(info$Wt) / ((trans_ht(info$Ht)) / 100) ^ 2

a <- mean(bmi, na.rm = TRUE)
s <- sd(bmi, na.rm = TRUE)
n <- sum(!is.na(bmi))
xbar <- 57.23
plot(density(bmi, na.rm = TRUE))
2 * (1 - pnorm(xbar, mean = a, sd = s / sqrt(n)))

range((as.numeric(bmi)), na.rm = TRUE)
IQR(as.numeric(bmi), na.rm = TRUE)

parity <- info$Parity
parity <-
  parity %>%
  stringr::str_extract("G[0-9]{1}") %>%
  stringr::str_replace("G", "") %>%
  as.numeric()


a <- mean(parity, na.rm = TRUE)
s <- sd(parity, na.rm = TRUE)
n <- sum(!is.na(parity))
xbar <- max(parity, na.rm = TRUE)
plot(density(parity, na.rm = TRUE))
2 * (1 - pnorm(xbar, mean = a, sd = s / sqrt(n)))

range((as.numeric(parity)), na.rm = TRUE)
IQR(as.numeric(parity), na.rm = TRUE)
mean(parity)


#------------------------------------------------------------------------------
## metabolomics data
met_data <-
  readr::read_csv(
    "/Users/shenxt/projects/smartD/data_analysis20191125/annotation/RPLC/metabolite_table.csv"
  )

met_data <-
  met_data %>%
  dplyr::filter(Level == 1 | Level == 2) %>%
  select(dplyr::starts_with("X"), dplyr::starts_with("SFU"))

met_data <-
  met_data[, -unique(which(is.na(met_data), arr.ind = TRUE)[, 2])]

met_data <-
  apply(met_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  })


met_data <- t(met_data) %>% as.data.frame()


sample_info <-
  temp_data %>%
  dplyr::filter(Sample_ID %in% colnames(met_data))

met_data <-
  lapply(unique(sample_info$Patient_ID), function(id) {
    sample_id <- sample_info %>%
      dplyr::filter(Patient_ID == id) %>%
      pull(Sample_ID)
    met_data %>%
      select(sample_id) %>%
      apply(1, mean)
  })


met_data <-
  do.call(cbind, met_data)

colnames(met_data) <- unique(sample_info$Patient_ID)

#------------------------------------------------------------------------------
## due date
due_data <-
  (as.Date(info$DD) - (as.Date(info$EDD) - 280)) / 7
range(due_data)

#------------------------------------------------------------------------------
## age
age <-
  info$Age %>%
  as.numeric()
age
range(age)

#------------------------------------------------------------------------------
## ethinic
ethinic <-
  info$`Ethinic Group`

ethinic <-
  case_when(
    ethinic == "1" ~ "White",
    ethinic == "2" ~ "Black",
    ethinic == "3" ~ "Latina",
    ethinic == "4" ~ "Pacific Islander",
    ethinic == "5" ~ "Asian",
    ethinic == "4 (Asian)" ~ "Asian",
    ethinic == "Afr Am" ~ "Black",
    TRUE ~ ethinic
  )

#------------------------------------------------------------------------------
## BMI
bmi <-
  info$BMI

bmi <-
  trans_wt(info$Wt) / ((trans_ht(info$Ht)) / 100) ^ 2



bmi[is.na(bmi)] <- 0

#------------------------------------------------------------------------------
## parity
parity <- info$Parity
parity <-
  parity %>%
  stringr::str_extract("G[0-9]{1}") %>%
  stringr::str_replace("G", "") %>%
  as.numeric()


## sex and twins
info$Sex
sex <- info$Sex
sex <-
  case_when(
    is.na(sex) ~ "NA",
    sex == "F, F" ~ "F_F",
    sex == "M,M" ~ "M_M",
    sex == "M,M" ~ "M_M",
    sex == "M / F" ~ "M_F",
    TRUE ~ sex
  )
## IVF
ivf <- rep(NA, 36)
ivf[grep("ivf", stringr::str_to_lower(info$`Pregnancy dx`))] <- "YES"
ivf[grep("transfer", stringr::str_to_lower(info$Dating))] <- "YES"
ivf[is.na(ivf)] <- "NO"

## induction
induction <-
  info$Induction

induction <-
  case_when(
    is.na(induction) ~ "NA",
    induction == "Y" ~ "YES",
    induction == "Yes" ~ "YES",
    induction == "N" ~ "NO",
  )


## type o delivery
# spont<U+81EA><U+53D1>
# augment <U+589E><U+52A0><U+589E><U+5927>
# vacuum <U+771F><U+7A7A>
# forcep:<U+94B3><U+5B50>

# birth weight
birth_wt <-
  info$`Birth wt`
birth_wt[33] <- "3015;3170"

birth_wt <-
  sapply(birth_wt, function(x) {
    if (!is.na(x)) {
      x <- stringr::str_split(x, ";")[[1]] %>%
        stringr::str_trim() %>%
        as.numeric() %>%
        sum()
    } else {
      NA
    }
    x
  })

birth_wt <- unname(birth_wt) %>%
  as.numeric()





data.frame(class = "age",
           value = age,
           stringsAsFactors = FALSE) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot() +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "", y = "Age (years)") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10)
  )
##dark
data.frame(class = "age",
           value = age,
           stringsAsFactors = FALSE) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(fill = "transparent", color = "white") +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "Age", y = "Years") +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(color = "grey")
  )

# ggsave(filename = "age_dark.png", width = 2, height = 7, bg = "transparent")




data.frame(class = "bmi",
           value = bmi,
           stringsAsFactors = FALSE) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot() +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "BMI", y = "BMI") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10)
  )
##dark
data.frame(class = "bmi",
           value = bmi,
           stringsAsFactors = FALSE) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(fill = "transparent", color = "white") +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "BMI", y = "BMI") +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(color = "grey")
  )

# ggsave(filename = "bmi_dark.png", width = 2, height = 7, bg = "transparent")

data.frame(class = "birth_wt",
           value = birth_wt,
           stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "", y = "Birth weight (g)") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10)
  )


##dark
data.frame(class = "birth_wt",
           value = birth_wt,
           stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(fill = "transparent", color = "white") +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "Birth weight", y = "g") +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(color = "grey")
  )

# ggsave(filename = "birth_wt_dark.png", width = 2, height = 7, bg = "transparent")


data.frame(class = "parity",
           value = parity,
           stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "", y = "Birth weight (g)") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10)
  )


##dark
data.frame(class = "parity",
           value = parity,
           stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(fill = "transparent", color = "white") +
  geom_jitter(colour = alpha("orange", 0.5),
              size = 6,
              shape = 16) +
  labs(x = "Parity", y = "No.") +
  ggdark::dark_theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(color = "grey")
  )

# ggsave(filename = "parity_dark.png", width = 2, height = 7, bg = "transparent")



par(mar = c(2, 2, 2, 2))
par(xpd = TRUE)
library(circlize)
set.seed(999)

## 36 participants
n <- 36
df <- data.frame(
  factors = sample(letters[1], n, replace = TRUE),
  # factors = paste("a",1:36,sep = "_"),
  x = 1:36,
  y = 1
)

# info1 <- info
# info1 <- rbind(rep(NA, ncol(info1)), info1)
df <- cbind(df, info)

library(circlize)

circos.par(
  "track.height" = 0.06,
  start.degree = 90,
  clock.wise = FALSE,
  gap.after = 95,
  cell.padding = c(0, 0, 0, 0)
)
circos.initialize(factors = df$factors, x = df$x)

### label
circos.track(
  factors = df$factors,
  x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = NA,
  # bg.col = NA,
  track.height = 0.0001,
  panel.fun = function(x, y) {
    circos.points(
      x = x - 0.5,
      y = y,
      pch = 19,
      col = text_colour
    )
    circos.text(
      x = x - 0.5,
      y = y + uy(2, "mm"),
      niceFacing = TRUE,
      labels = df$ID,
      facing = "clockwise",
      adj = c(0, 0.5),
      col = text_colour,
      cex = 0.8
    )
  }
)

# metabolomics data
col_fun <- colorRamp2(c(-2, 0, 2), c("#008EA0FF", "white", "#C71000FF"))

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 256),
  bg.border = "#FFA319FF",
  bg.col = NA,
  track.height = 0.2,
  panel.fun = function(x, y) {
    # m2 = met_data[, order.dendrogram(dend)]
    m2 <- met_data
    col_mat <- col_fun(m2)
    nr <- nrow(m2)
    nc <- ncol(m2)
    for (i in 1:nr) {
      circos.rect(
        xleft = 1:nc - 1,
        ybottom = rep(nr - i, nc),
        xright = 1:nc,
        ytop = rep(nr - i + 1, nc),
        border = col_mat[i, ],
        col = col_mat[i, ]
      )
    }
  }
)




## due date
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(30, 42),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    circos.points(
      x = 1:36 - 0.5,
      y = due_data,
      cex = 1,
      pch = 19,
      col = "#8A9045FF"
      # fill = "#8A9045FF"
    )
    circos.lines(x = 1:36 - 0.5,
                 y = due_data,
                 col = "#8A9045FF")
  }
)



### age
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 40),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 0.8,
      # ybottom = rep(ylim[1], n_breaks - 1),
      ybottom = rep(0, 36),
      xright = breaks[-1] - 1.2,
      # ytop = rep(ylim[2], n_breaks - 1),
      ytop = age,
      col = "#FFA319FF"
      # border = NA
    )
  }
)


## ethinic
library(RColorBrewer)
# display.brewer.all()
# display.brewer.pal(n = 7, name = "Set3")
col <-
  brewer.pal(n = 7, name = "Set3")

temp_col <- ethinic
temp_col[ethinic == unique(ethinic)[1]] <- col[1]
temp_col[ethinic == unique(ethinic)[2]] <- col[2]
temp_col[ethinic == unique(ethinic)[3]] <- col[3]
temp_col[ethinic == unique(ethinic)[4]] <- col[4]
temp_col[ethinic == unique(ethinic)[5]] <- col[5]
temp_col[ethinic == unique(ethinic)[6]] <- col[6]
temp_col[ethinic == unique(ethinic)[7]] <- col[7]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 1,
      ybottom = rep(ylim[1], n_breaks - 1),
      xright = breaks[-1] - 1,
      ytop = rep(ylim[2], n_breaks - 1),
      col = temp_col
      # border = NA
    )
  }
)


## BMI
bmi[bmi == 0] <- NA
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 45),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 0.8,
      # ybottom = rep(ylim[1], n_breaks - 1),
      ybottom = rep(0, 36),
      xright = breaks[-1] - 1.2,
      # ytop = rep(ylim[2], n_breaks - 1),
      ytop = bmi,
      col = "#8A9045FF"
      # border = NA
    )
  }
)

## parity

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 10),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 0.8,
      # ybottom = rep(ylim[1], n_breaks - 1),
      ybottom = rep(0, 36),
      xright = breaks[-1] - 1.2,
      # ytop = rep(ylim[2], n_breaks - 1),
      ytop = parity,
      col = "#155F83FF"
      # border = NA
    )
  }
)


#### sex
col_fun <- colorRamp2(c(-2, 0, 2), c("green", "black", "red"))
sex[sex == "NA"] <- NA
temp_col <- sex
temp_col <-
  case_when(
    temp_col == "M" ~ "#FFA319FF",
    temp_col == "F" ~ "#8F3931FF",
    temp_col == "M_M" ~ "#C16622FF",
    temp_col == "F_F" ~ "#800000FF",
    temp_col == "M_F" ~ "#155F83FF",
    temp_col == "NA" ~ "#767676FF"
  )

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 1,
      ybottom = rep(ylim[1], n_breaks - 1),
      xright = breaks[-1] - 1,
      ytop = rep(ylim[2], n_breaks - 1),
      col = temp_col
      # border = NA
    )
  }
)


## ivf
temp_col <- ivf
temp_col <-
  case_when(temp_col == "YES" ~ "skyblue", temp_col == "NO" ~ "grey")

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 1,
      ybottom = rep(ylim[1], n_breaks - 1),
      xright = breaks[-1] - 1,
      ytop = rep(ylim[2], n_breaks - 1),
      col = temp_col
      # border = NA
    )
  }
)



## induction
temp_col <- induction
temp_col <-
  case_when(temp_col == "YES" ~ "red", temp_col == "NO" ~ "grey")

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 1,
      ybottom = rep(ylim[1], n_breaks - 1),
      xright = breaks[-1] - 1,
      ytop = rep(ylim[2], n_breaks - 1),
      col = temp_col
      # border = NA
    )
  }
)


## birth weight

range(birth_wt, na.rm = TRUE)
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 4400),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    breaks <- 1:37
    n_breaks <- length(breaks)
    circos.rect(
      xleft = breaks[-n_breaks] - 0.8,
      # ybottom = rep(ylim[1], n_breaks - 1),
      ybottom = rep(0, 36),
      xright = breaks[-1] - 1.2,
      # ytop = rep(ylim[2], n_breaks - 1),
      ytop = birth_wt,
      col = "#8A9045FF"
      # border = NA
    )
  }
)


circos.clear()



#### 20191030
setwd(r4projects::get_project_wd())
setwd("patient information/")
info <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
sample_information <- readr::read_csv("sample_information.csv")
info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste("SF", ., sep = "")


info <-
  info %>%
  dplyr::filter(ID %in% sample_information$Patient_ID)

## ID:participant ID
## Maternal DOB:participant ID
load(
  "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_dis"
)
load(
  "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_val"
)


dis_id <-
  sample_data_dis$Patient_ID %>%
  unique()

val_id <-
  sample_data_val$Patient_ID %>%
  unique()

dis_info <-
  info %>%
  dplyr::filter(ID %in% dis_id)

val_info <-
  info %>%
  dplyr::filter(ID %in% val_id)

## mother age
library(tidyverse)

dis_info <-
  dis_info %>%
  dplyr:::mutate(Age = as.numeric(Age))

dis_info %>%
  dplyr::select(Age) %>%
  mutate(Age = as.numeric(Age)) %>%
  summarise(mean = mean(Age), sd = sd(Age))

val_info <-
  val_info %>%
  dplyr:::mutate(Age = as.numeric(Age))

val_info %>%
  dplyr::select(Age) %>%
  mutate(Age = as.numeric(Age)) %>%
  summarise(mean = mean(Age), sd = sd(Age))


## BMI
dis_info %>%
  dplyr::select(Ht, Wt, BMI) %>%
  mutate(ht = trans_ht(Ht),
         wt = trans_wt(Wt),
         bmi = wt / (ht / 100) ^ 2) %>%
  summarise(mean.bmi = mean(bmi), sd.bmi = sd(bmi))


dis_info <-
  dis_info %>%
  dplyr::mutate(ht = trans_ht(Ht),
                wt = trans_wt(Wt),
                bmi = wt / (ht / 100) ^ 2)

val_info %>%
  dplyr::select(Ht, Wt, BMI) %>%
  mutate(ht = trans_ht(Ht),
         wt = trans_wt(Wt),
         bmi = wt / (ht / 100) ^ 2) %>%
  summarise(mean.bmi = mean(bmi), sd.bmi = sd(bmi))

val_info <-
  val_info %>%
  dplyr::mutate(ht = trans_ht(Ht),
                wt = trans_wt(Wt),
                bmi = wt / (ht / 100) ^ 2)


## GA
dis_info2 <-
  sample_data_dis %>%
  dplyr::filter(!is.na(day_0)) %>%
  select(-contains("_POS")) %>%
  select(-contains("_NEG")) %>%
  left_join(dis_info, ., by = c("ID" = "Patient_ID")) %>%
  mutate(
    day_0 = as.Date(day_0, "%m/%d/%Y"),
    edd = as.Date(EDD.x),
    begin = edd - 280
  ) %>%
  transmute(ga = as.Date(DD.x) - begin) %>%
  summarise(mean(as.numeric(ga)), sd(as.numeric(ga)))


dis_info2

val_info2 <-
  sample_data_val %>%
  dplyr::filter(!is.na(day_0)) %>%
  select(-contains("_POS")) %>%
  select(-contains("_NEG")) %>%
  left_join(val_info, ., by = c("ID" = "Patient_ID")) %>%
  mutate(
    day_0 = as.Date(day_0, "%m/%d/%Y"),
    edd = as.Date(EDD.x),
    begin = edd - 280
  ) %>%
  transmute(ga = as.Date(DD.x) - begin) %>%
  summarise(mean(as.numeric(ga)), sd(as.numeric(ga)))

val_info2

dis_info %>%
  select(`Birth wt`) %>%
  dplyr::filter(!is.na(`Birth wt`)) %>%
  pull(`Birth wt`) %>%
  stringr::str_split(";") %>%
  unlist() %>%
  as.numeric() %>%
  as_tibble() %>%
  summarise(mean = mean(value), sd = sd(value))


val_info %>%
  select(`Birth wt`) %>%
  dplyr::filter(!is.na(`Birth wt`)) %>%
  pull(`Birth wt`) %>%
  stringr::str_replace("/", ";") %>%
  stringr::str_split(";") %>%
  unlist() %>%
  as.numeric() %>%
  as_tibble() %>%
  summarise(mean = mean(value), sd = sd(value))



## gender of child
dis_info$Sex %>%
  as_tibble() %>%
  dplyr::filter(!is.na(value))


dis_info








####20200610
#### 20191030
setwd(r4projects::get_project_wd())
setwd("patient information/")
info <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
sample_information <- readr::read_csv("sample_information.csv")
info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste("SF", ., sep = "")


info <-
  info %>%
  dplyr::filter(ID %in% sample_information$Patient_ID)

## ID:participant ID
## Maternal DOB:participant ID
load(
  "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_dis"
)
load(
  "/Users/shenxt/projects/smartD/data_analysis20200108/prediction/metabolites/sample_data_val"
)


dis_id <-
  sample_data_dis$Patient_ID %>%
  unique()

val_id <-
  sample_data_val$Patient_ID %>%
  unique()

dis_info <-
  info %>%
  dplyr::filter(ID %in% dis_id)

val_info <-
  info %>%
  dplyr::filter(ID %in% val_id)

dis_info$class <- "Training"
val_info$class <- "Testing"

patient_info <- rbind(dis_info, val_info)

## mother age
library(tidyverse)

patient_info <-
  patient_info %>%
  dplyr:::mutate(Age = as.numeric(Age))

## BMI

patient_info <-
  patient_info %>%
  dplyr::mutate(ht = trans_ht(Ht),
                wt = trans_wt(Wt),
                bmi = wt / (ht / 100) ^ 2)

## GA

##ethnicity
patient_info <-
  patient_info %>%
  dplyr::mutate(
    ethinic =  case_when(
      ethinic == "1" ~ "White",
      ethinic == "2" ~ "Black",
      ethinic == "3" ~ "Latina",
      ethinic == "4" ~ "Pacific Islander",
      ethinic == "5" ~ "Asian",
      ethinic == "4 (Asian)" ~ "Asian",
      ethinic == "Afr Am" ~ "Black",
      TRUE ~ ethinic
    )
  )



patient_info <-
  patient_info %>%
  dplyr::mutate(
    parity = stringr::str_extract(Parity, "G[0-9]{1}") %>%
      stringr::str_replace("G", "") %>%
      as.numeric()
  )




xlsx::write.xlsx(patient_info, "patient_info.xlsx")
