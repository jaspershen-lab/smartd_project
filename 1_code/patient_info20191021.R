setwd(r4projects::get_project_wd())
setwd("patient information/")
all346_urine_info <-
  readr::read_csv("SmartD_all346urine.csv")

clinic_varibale_summary <-
  readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx", sheet = 1)

###This is batch data set
smartd_rplc_batch1 <-
  readr::read_csv(
    "E:/project/smartD/data_analysis20190828/prediction/identification_table/smartd_rplc_batch1.csv"
  )

smartd_rplc_batch2 <-
  readr::read_csv("E:/project/smartD/smartD_batch1_2/RPLC/POS_NEG/ms1_data_rplc_batch2.csv")

###we must get the samples with 'P' coresponding to which sample
all346_urine_info$sample_id_RPLC %>% stringr::str_sort(numeric = TRUE) %>%
  grep("SFU_B", ., value = TRUE)

##we have 198 sampels in total
colnames(smartd_rplc_batch2) %>% stringr::str_sort(numeric = TRUE) %>%
  grep("X", ., value = TRUE)

colnames(smartd_rplc_batch2) %>% stringr::str_sort(numeric = TRUE) %>%
  grep("P", ., value = TRUE)


###so it is very clear that X178 is P1, X179 is P2, and X180 is P3. X198 is P21. and so on
colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] <-
  colnames(smartd_rplc_batch2)[grep("P", colnames(smartd_rplc_batch2))] %>%
  stringr::str_replace("P", "") %>%
  as.numeric() %>%
  `+`(177) %>%
  paste("X", ., sep = "")

name_batch1 <-
  smartd_rplc_batch1 %>%
  select(-(name:Database)) %>%
  colnames()

name_batch2 <-
  smartd_rplc_batch2 %>%
  select(-(name:rt)) %>%
  select(-contains("blk")) %>%
  select(-contains("QC")) %>%
  colnames()

sample_id_RPLC <-
  all346_urine_info$sample_id_RPLC

sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] <-
  sample_id_RPLC[grep("SFU_B", sample_id_RPLC)] %>%
  gsub(pattern = "SFU_B",
       replacement = "",
       x = .) %>%
  paste("X", ., sep = "")

name_batch1 %in% sample_id_RPLC
name_batch2 %in% sample_id_RPLC


all346_urine_info$sample_id_RPLC <-
  sample_id_RPLC

all346_info <-
  all346_urine_info

all346_info <-
  all346_info %>%
  select(
    Patient_ID = PTID,
    Sample_ID = sample_id_RPLC,
    Visit,
    Date.Acquired,
    GA = Visit.GA
  ) %>%
  arrange(Patient_ID, Visit)

##some sample's GA are NA, should be removed
# all346_info %>%
#   filter(is.na(GA)) %>%
#   pull(Sample_ID)
#
# patient_info <-
#   readr::read_csv("patient_info.csv")
#
# patient_info$Sample_ID <-
#   patient_info$Sample_ID %>%
#   gsub(pattern = "SFU_B", replacement = "X", x = .)
#
# patient_info %>%
#   select(Sample_ID, `Visit GA`) %>%
#   arrange(Sample_ID)
#
# all346_info %>%
#   select(Sample_ID, GA) %>%
#   filter(stringr::str_starts(Sample_ID, "X")) %>%
#   arrange(Sample_ID)
#
#
# intersect(all346_info$Sample_ID, patient_info$Sample_ID)

###modify GA
GA <- all346_info$GA

GA <-
  (GA - trunc(GA)) * 10 / 7 + trunc(GA)

all346_info$GA <-
  GA

##patient_info
final_info <-
  readr::read_csv("FINALSMARTDiaphragm2_DATA_2019-07-01_1615.csv")

###change the Patient ID
all346_info$Patient_ID[-grep(pattern = "SF", all346_info$Patient_ID)] <-
  all346_info$Patient_ID[-grep(pattern = "SF", all346_info$Patient_ID)] %>%
  paste("SF", ., sep = "")

participant_id <-
  final_info$participant_id

participant_id <-
  gsub("sf", "SF", participant_id)

participant_id[-grep('SF', participant_id)] <-
  paste("SF", participant_id[-grep('SF', participant_id)], sep = "")


final_info$participant_id <-
  participant_id


setdiff(all346_info$Patient_ID, final_info$participant_id)
setdiff(final_info$participant_id, all346_info$Patient_ID)

final_info$redcap_event_name <-
  gsub("study_visit_|_arm_1", "", final_info$redcap_event_name) %>%
  as.numeric()


sample_info_191021 <-
  left_join(
    all346_info,
    final_info,
    by = c("Patient_ID" = "participant_id", "Visit" = "redcap_event_name")
  )

sample_info_191021 %>%
  filter(Sample_ID %in% name_batch1) %>%
  pull(GA)

sample_info_191021 %>%
  filter(Sample_ID %in% name_batch2) %>%
  pull(GA)

clinic_varibale_summary$ID <-
  stringr::str_to_upper(clinic_varibale_summary$ID)

clinic_varibale_summary$ID[-grep("SF", clinic_varibale_summary$ID)] <-
  paste("SF", clinic_varibale_summary$ID[-grep("SF", clinic_varibale_summary$ID)], sep = "")

sample_info_191021$Patient_ID %in% clinic_varibale_summary$ID

sample_info_191021 <-
  sample_info_191021 %>%
  left_join(clinic_varibale_summary, by = c("Patient_ID" = "ID"))



write.csv(sample_info_191021, "sample_info_191021.csv", row.names = FALSE)






# patient_info$Sample_ID
#
# info1 <-
#   sample_info_191021 %>%
#   filter(stringr::str_starts(Sample_ID, "X")) %>%
#   select(Patient_ID:GA) %>%
#   arrange(Sample_ID)
#
# info2 <-
#   patient_info %>%
#   filter(stringr::str_starts(Sample_ID, "SFU_B")) %>%
#   select(Patient_ID:`Visit GA`) %>%
#   arrange(Sample_ID) %>%
#   mutate(Sample_ID = stringr::str_replace(Sample_ID, "SFU_B", "X"))



#####the clinic information
sample_info_191021$Patient_ID %>%
  unique() %>%
  length()

sample_info_191021 %>%
  filter(is.na(GA)) %>%
  dplyr::select(Sample_ID)

library(ggExtra)


temp_data <-
  sample_info_191021 %>%
  mutate(GA2 = case_when(is.na(GA) ~ 45, !is.na(GA) ~ GA)) %>%
  mutate(class = case_when(is.na(GA) ~ "PP", !is.na(GA) ~ "Normal"))




##which person is preterm
preterm21 <-
  sample_info_191021 %>%
  mutate(diff_day = as.Date(sample_info_191021$EDD) - as.Date(sample_info_191021$DD)) %>%
  filter(diff_day >= 21) %>%
  pull(Patient_ID) %>%
  unique()


preterm7 <-
  sample_info_191021 %>%
  mutate(diff_day = as.Date(sample_info_191021$EDD) - as.Date(sample_info_191021$DD)) %>%
  filter(diff_day >= 7 & diff_day < 21) %>%
  pull(Patient_ID) %>%
  unique()



text_colour <- sort(unique(temp_data$Patient_ID))
text_colour <- case_when(
  text_colour %in% preterm21 ~ "#58593FFF",
  text_colour %in% preterm7 ~ "#800000FF",
  TRUE ~ "#BEBEBE"
)

text_face <- sort(unique(temp_data$Patient_ID))
text_face <- case_when(text_face %in% preterm21 ~ "bold",
                       text_face %in% preterm7 ~ "bold",
                       TRUE ~ "plain")


(
  plot1 <-
    temp_data %>%
    ggplot(aes(GA2, Patient_ID)) +
    geom_point(aes(colour = class), show.legend = FALSE, size = 2) +
    geom_rect(
      aes(
        xmin = 10,
        ymin = 0,
        xmax = 15,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_rect(
      aes(
        xmin = 20,
        ymin = 0,
        xmax = 25,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_rect(
      aes(
        xmin = 30,
        ymin = 0,
        xmax = 35,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_rect(
      aes(
        xmin = 40,
        ymin = 0,
        xmax = 45,
        ymax = Inf
      ),
      fill = "grey",
      alpha = 0.01,
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = match(preterm7, sort(unique(
        temp_data$Patient_ID
      ))),
      colour = "#800000FF",
      linetype = 1
    ) +
    geom_hline(
      yintercept = match(preterm21, sort(unique(
        temp_data$Patient_ID
      ))),
      colour = "#58593FFF",
      linetype = 1
    ) +
    scale_colour_manual(values = c(
      "PP" = "#155F83FF", "Normal" = "#FFA319FF"
      # "Preterm>7_days" = "#800000B2",
      # "Preterm>21_days" = "#350E20B2"
    )) +
    labs(x = "Gestational age (weeks)", y = "Subject ID") +
    scale_x_continuous(
      breaks = c(10, 20, 30, 40, 45),
      labels = c(10, 20, 30, 40, "PP"),
      limits = c(10, 46)
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(colour = "#8A9045FF", linetype = 2),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      axis.title = element_text(size = 15),
      axis.text.y = element_text(
        size = 11,
        colour = text_colour,
        face = text_face
      ),
      axis.text.x = element_text(size = 13)
    )
)

(
  plot2 <-
    ggplot(temp_data, aes(x = GA2)) +
    geom_histogram(
      binwidth = 0.5,
      colour = "white",
      fill = "#8A9045B2"
    ) +
    labs(x = "GA (weeks)", y = "Sample number") +
    theme_bw() +
    scale_x_continuous(
      limits = c(10, 46),
      name = NULL,
      labels = NULL,
      breaks = NULL
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13)
    )
)

(
  plot3 <-
    ggplot(temp_data, aes(x = Patient_ID)) +
    geom_bar(width = 0.8, fill = "#8A9045B2") +
    labs(x = "GA (weeks)", y = "Sample number") +
    theme_bw() +
    scale_x_discrete(
      name = NULL,
      label = NULL,
      breaks = NULL
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    coord_flip() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13)
    )
)


library(patchwork)

plot <-
  {
    plot2 + plot_spacer() + plot_layout(ncol = 2, widths = c(3, 1))
  } -
  {
    plot1 + plot3 + plot_layout(ncol = 2, widths = c(3, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 3))


plot


ggsave(plot,
       filename = "Sample_colection_distribution.pdf",
       width = 7,
       height = 7)

ggsave(plot,
       filename = "Sample_colection_distribution2.pdf",
       width = 10,
       height = 7)

plot


####a circos plot to show the information
info <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")

info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste('SF', ., sep = "")


info <-
  info %>%
  filter(ID %in% sample_info_191021$Patient_ID)

##we have
colnames(info)
#------------------------------------------------------------------------------
##metabolomics data
met_data <-
  readr::read_csv("E:/project/smartD/data_analysis20191125/annotation/RPLC/metabolite_table.csv")

met_data <-
  met_data %>%
  filter(Level == 1 | Level == 2) %>%
  select(dplyr::starts_with("X"), dplyr::starts_with("SFU"))

met_data <-
  met_data[, -unique(which(is.na(met_data), arr.ind = TRUE)[, 2])]

met_data <-
  apply(met_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  })


met_data <- t(met_data) %>% as.data.frame()


sample_info <-
  sample_info_191021 %>%
  filter(Sample_ID %in% colnames(met_data))

met_data <-
  lapply(unique(sample_info$Patient_ID), function(id) {
    sample_id <- sample_info %>%
      filter(Patient_ID == id) %>%
      pull(Sample_ID)
    met_data %>%
      select(sample_id) %>%
      apply(1, mean)
  })


met_data <-
  do.call(cbind, met_data)

colnames(met_data) <- unique(sample_info$Patient_ID)

#------------------------------------------------------------------------------
##due date
due_data <-
  (as.Date(info$DD) - (as.Date(info$EDD) - 280)) / 7
range(due_data)

#------------------------------------------------------------------------------
##age
age <-
  info$Age %>%
  as.numeric()
age
range(age)

#------------------------------------------------------------------------------
##ethinic
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
  )

#------------------------------------------------------------------------------
##BMI
bmi <-
  info$BMI
bmi[is.na(bmi)] <- 0

#------------------------------------------------------------------------------
##parity
parity <- info$Parity
parity <-
  parity %>%
  stringr::str_extract("G[0-9]{1}") %>%
  stringr::str_replace("G", "") %>%
  as.numeric()


##sex and twins
info$Sex
sex <- info$Sex
sex <-
  case_when(
    is.na(sex) ~ "NA",
    sex == "F, F" ~ 'F_F',
    sex == "M,M" ~ 'M_M',
    sex == "M,M" ~ 'M_M',
    sex == "M / F" ~ 'M_F',
    TRUE ~ sex
  )
##IVF
ivf <- rep(NA, 36)
ivf[grep("ivf", stringr::str_to_lower(info$`Pregnancy dx`))] <- "YES"
ivf[grep("transfer", stringr::str_to_lower(info$Dating))] <- "YES"
ivf[is.na(ivf)] <- "NO"

##induction
induction <-
  info$Induction

induction <-
  case_when(
    is.na(induction) ~ "NA",
    induction == "Y" ~ "YES",
    induction == "Yes" ~ "YES",
    induction == "N" ~ "NO",
  )


##type o delivery
#spont自发
#augment 增加增大
#vacuum 真空
#forcep:钳子

#birth weight
birth_wt <-
  info$`Birth wt`
birth_wt[19] <- "3015;3170"

birth_wt <-
  sapply(birth_wt, function(x) {
    if (!is.na(x)) {
      x <- stringr::str_split(x, ";")[[1]] %>%
        as.numeric() %>%
        mean()
    } else{
      x
    }
    x
    
  })

birth_wt <- unname(birth_wt) %>%
  as.numeric()


par(mar = c(2, 2, 2, 2))
par(xpd = TRUE)
library(circlize)
set.seed(999)

##36 participants
n <- 36
df <- data.frame(
  factors = sample(letters[1], n, replace = TRUE),
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

###label
circos.track(
  factors = df$factors,
  x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = NA,
  bg.col = NA,
  track.height = 0.0001,
  panel.fun = function(x, y) {
    circos.points(
      x = x - 0.5,
      y = y,
      pch = 19,
      col = "black"
    )
    circos.text(
      x = x - 0.5,
      y = y + uy(2, "mm"),
      niceFacing = TRUE,
      labels = df$ID,
      facing = "clockwise",
      adj = c(0, 0.5),
      col = "black",
      cex = 0.8
    )
  }
)

#metabolomics data
col_fun = colorRamp2(c(-2, 0, 2), c("#155F837F", "white", "#8000007F"))

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
    m2 = met_data
    col_mat = col_fun(m2)
    nr = nrow(m2)
    nc = ncol(m2)
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


##due date
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(30, 42),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.points(
      x = 1:36 - 0.5,
      y =  due_data,
      cex = 1,
      pch = 19,
      col = "#8A9045FF"
      # fill = "#8A9045FF"
    )
    circos.lines(x = 1:36 - 0.5,
                 y =  due_data,
                 col = "#8A9045FF")
  }
)



###age
circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 40),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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


##ethinic
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
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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


##BMI
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
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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

##parity

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 10),
  # bg.border = NA,
  bg.col = NA,
  track.height = 0.06,
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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


####sex
col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))
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
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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


##ivf
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
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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



##induction
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
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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


##birth weight

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
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = 1:37
    n_breaks = length(breaks)
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



####20191030
setwd(r4projects::get_project_wd())
setwd("patient information/")
info <- readxl::read_xlsx("SmartD_ClinicalVariables_PartiallySummarized.xlsx")
sample_info_191021 <- readr::read_csv("sample_info_191021.csv")
info$ID <-
  info$ID %>%
  stringr::str_replace("sf", "") %>%
  paste('SF', ., sep = "")


info <-
  info %>%
  dplyr::filter(ID %in% sample_info_191021$Patient_ID)
##ID:participant ID
##Maternal DOB:participant ID
load(
  "E:/project/smartD/data_analysis20191015/prediction/identification_table/sample_data_dis"
)
load(
  "E:/project/smartD/data_analysis20191015/prediction/identification_table/sample_data_val"
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

##mother age
library(tidyverse)

dis_info %>%
  dplyr::select(Age) %>%
  mutate(Age = as.numeric(Age)) %>%
  summarise(mean = mean(Age), sd = sd(Age))

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


val_info %>%
  dplyr::select(Ht, Wt, BMI) %>%
  mutate(ht = trans_ht(Ht),
         wt = trans_wt(Wt),
         bmi = wt / (ht / 100) ^ 2) %>%
  summarise(mean.bmi = mean(bmi), sd.bmi = sd(bmi))



##GA
dis_info2 <-
  sample_data_dis %>%
  filter(!is.na(day_0)) %>%
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
  filter(!is.na(day_0)) %>%
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
  filter(!is.na(`Birth wt`)) %>%
  pull(`Birth wt`) %>%
  stringr::str_split(";") %>%
  unlist() %>%
  as.numeric() %>%
  as_tibble() %>%
  summarise(mean = mean(value), sd = sd(value))


val_info %>%
  select(`Birth wt`) %>%
  filter(!is.na(`Birth wt`)) %>%
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
  filter(!is.na(value))
