###20191125 use all the identified metabolites (Level 1 and Level 2),
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
setwd("data_analysis20191125/metabolite_table/correlation_network/")
rm(list = ls())
metabolite_table <- readr::read_csv("metabolite_table.csv")

###we only use the level 1 and level 2 identifications
#1-3 means identification confidence level, 4 means MS2 spectra match si not good enough, and 5 means this compound
#may be not in human body
metabolite_table <- 
  metabolite_table %>% 
  filter(Level == 1 | Level == 2)

dim(metabolite_table)

##remove QC_DL and Blank Samples
metabolite_table <- 
  metabolite_table %>% 
  select(., -contains("BLK")) %>% 
  select(., -contains("QC")) %>% 
  select(., -contains("blk"))

####remove the duplicated metabolites according to score
metabolite_table <- 
  metabolite_table %>% 
  group_by(., Compound.name) %>% 
  filter(., SS == max(SS)) %>% 
  ungroup()

dim(metabolite_table)


metabolite_tags <- 
  metabolite_table %>% 
  select(name:Database)

####give new information for metabolites
kegg.id <- 
  lapply(metabolite_tags$Compound.name, function(x){
    chemicalTranslate::transID(query = x, from = "Chemical name", to = "KEGG", top = 1)
  }) %>% 
  do.call(rbind, .)

hmdb.id <- 
  lapply(metabolite_tags$Compound.name, function(x){
    chemicalTranslate::transID(query = x, from = "Chemical name",
                               to = "Human Metabolome Database", top = 1)
  }) %>% 
  do.call(rbind, .)

ID_trans <-
  data.frame(
    KEGG = kegg.id$KEGG,
    HMDB = hmdb.id$`Human Metabolome Database`,
    stringsAsFactors = FALSE
  )

hmdb_kegg <- 
  metabolite_tags %>% 
  select(KEGG.ID, HMDB.ID) %>% 
  data.frame(ID_trans, stringsAsFactors = FALSE)

hmdb_kegg <-
  apply(hmdb_kegg, 1, function(x){
    hmdb <- x[c(2,4)] 
    hmdb <- hmdb[!is.na(hmdb)]
    kegg <- x[c(1,3)]
    kegg <- kegg[!is.na(kegg)]
    
    hmdb <- ifelse(length(hmdb) == 0, NA, hmdb[1])
    kegg <- ifelse(length(kegg) == 0, NA, kegg[1])
    
    c(hmdb, kegg)
    
  }) %>% 
  t() %>% 
  as_tibble()


colnames(hmdb_kegg) <- 
  c("HMDB.ID", "KEGG.ID")

hmdb_kegg$HMDB.ID <- 
  hmdb_kegg$HMDB.ID %>% 
  sapply(function(x) {
    if(is.na(x)) {
      return(NA) 
    }
    if(nchar(x) == 9){
      pre <- stringr::str_extract(x, "[A-Za-z]{1,}")
      end <- stringr::str_extract(x, "[0-9]{1,}")
      end <- paste("00", end, sep = "")
      x <- paste(pre, end, sep = "")
    }
    x 
  }
    ) %>% 
  unname()


cbind(metabolite_tags$HMDB.ID, hmdb_kegg[,1])

metabolite_tags[,c("HMDB.ID", 'KEGG.ID')] <-
  hmdb_kegg


###should use classyfire to calculate the class of all the metabolites
inchikey1 <- 
  metabolite_tags$HMDB.ID %>% 
  lapply(function(x){
    chemicalTranslate::transID(query = x, from = "Human Metabolome Database", to = "InChIKey", top = 1)
  }) %>% 
  do.call(rbind, .)


inchikey2 <- 
  metabolite_tags$KEGG.ID %>% 
  lapply(function(x){
    chemicalTranslate::transID(query = x, from = "KEGG", to = "InChIKey", top = 1)
  }) %>% 
  do.call(rbind, .)


inchikey3 <- 
  metabolite_tags$name %>% 
  lapply(function(x){
    chemicalTranslate::transID(query = x, from = "Chemical name", to = "InChIKey", top = 1)
  }) %>% 
  do.call(rbind, .)

inchikey <- cbind(inchikey1 = inchikey1$InChIKey,
                  inchikey2 = inchikey2$InChIKey,
                  inchikey3 = inchikey3$InChIKey) %>% 
  as_tibble() %>% 
  apply(1, function(x){
    x <- x[!is.na(x)]
    ifelse(length(x) == 0, NA, x[1])
  })

metabolite_class <- vector(mode = "list", length = length(inchikey))



metabolite_class <- 
  lapply(inchikey, function(x){
    Sys.sleep(time = 5)
    result <- metflow2::get_metclass(inchikey = x, sleep = 5)
  })


idx <- 
  lapply(metabolite_class, class) %>% unlist() %>% 
  `==`("logical") %>% 
  which()

inchikey[idx]


super_class <- lapply(metabolite_class, function(x){
  if(is.na(x)) return(NA)
  x@classification_info %>% 
    dplyr::filter(name == "Superclass") %>% 
    pull(value)
}) %>% 
  unlist()


class <- lapply(metabolite_class, function(x){
  if(is.na(x)) return(NA)
  x@classification_info %>% 
    dplyr::filter(name == "Class") %>% 
    pull(value)
}) %>% 
  unlist()

sub_class <- lapply(metabolite_class, function(x){
  if(is.na(x)) return(NA)
  x@classification_info %>% 
    dplyr::filter(name == "Subclass") %>% 
    pull(value)
}) %>% 
  unlist()

sub_class[which(sub_class == "Not available")] <- NA


metabolite_tags <- data.frame(metabolite_tags,
                              super_class,
                              class, stringsAsFactors = FALSE)


metabolite_table <- 
  metabolite_table %>% 
  select(-c(name:Database))

metabolite_table <- metabolite_table %>%
  select(contains("SFU"))

##sample sampels may have NAs
remove_idx <- 
which(is.na(metabolite_table), arr.ind = TRUE)[,2] %>% 
  unique()

colnames(metabolite_table)[remove_idx]

metabolite_table[,remove_idx]

metabolite_table <- metabolite_table[,-remove_idx]

dim(metabolite_table)

###get mean value for each person
sample_info_191021 <- 
  readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")


sample_info <- 
  sample_info_191021 %>% 
  filter(stringr::str_detect(Sample_ID, "SFU"))


metabolite_table <-
  metabolite_table %>% 
  select(one_of(sample_info$Sample_ID))


sample_info <-
  sample_info %>% 
  filter(Sample_ID %in% colnames(metabolite_table))

metabolite_table <-
lapply(unique(sample_info$Patient_ID), function(x){
  sample_id <- 
    sample_info %>% 
    filter(Patient_ID == x) %>% 
    pull(Sample_ID)
  
  metabolite_table[,sample_id] %>% 
    apply(1, mean)
  
}) %>% 
  do.call(cbind, .)

colnames(metabolite_table) <-
  unique(sample_info$Patient_ID)


##add clinical information
##due date
info <- 
  readxl::read_xlsx("E:/project/smartD/patient information/SmartD_ClinicalVariables_PartiallySummarized.xlsx")
info <- 
  info %>% 
  mutate(ID = stringr::str_replace(ID, "sf", "")) %>% 
  mutate(ID = paste("SF", ID, sep = ""))

info <- 
  info %>% 
  filter(ID %in% colnames(metabolite_table))


info$ID == colnames(metabolite_table)

due_date <- 
  (as.Date(info$DD) - (as.Date(info$EDD) - 280))/7 %>% 
  as.numeric()

##age
age <- 
  info$Age %>% 
  as.numeric()
age

##ethinic
ethinic <- 
  info$`Ethinic Group`

ethinic <- 
  case_when(ethinic == "1" ~ "White",
            ethinic == "Caucasian" ~ "White",
            ethinic == "2" ~ "Black",
            ethinic == "3" ~ "Latina",
            ethinic == "4" ~ "Pacific Islander",
            ethinic == "5" ~ "Asian",
            ethinic == "4 (Asian)" ~ "Asian",
            ethinic == "Asian" ~ "Asian",
            ethinic == "Afr Am" ~ "Black",
            TRUE ~ ethinic
  )

##BMI
ht <- trans_ht(info$Ht)/100
wt <- trans_wt(info$Wt)
bmi <- wt/(ht^2)

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
  case_when(is.na(sex) ~ "NA", 
            sex == "F, F" ~ 'F_F',
            sex == "M,M" ~ 'M_M',
            sex == "M,M" ~ 'M_M',
            sex == "M / F" ~ 'M_F',
            TRUE ~ sex
  )
##IVF
ivf <- rep(NA, nrow(info))
ivf[grep("ivf", stringr::str_to_lower(info$`Pregnancy dx`))] <- "YES"
ivf[grep("transfer", stringr::str_to_lower(info$Dating))] <- "YES"
ivf[is.na(ivf)] <- "NO"

##induction
induction <- 
  info$Induction

induction <- 
  case_when(is.na(induction) ~ "NA",
            induction == "Y" ~ "YES",
            induction == "Yes" ~ "YES",
            induction == "N" ~ "NO",
  )


birth_wt <- 
  info$`Birth wt`
birth_wt[15] <- "3015;3170"
birth_wt <- 
  sapply(birth_wt, function(x){
    if(!is.na(x)){
      x <- stringr::str_split(x, ";")[[1]] %>% 
        as.numeric() %>% 
        mean()  
    }else{
      x
    }
    x
  })

birth_wt <- unname(birth_wt) %>% 
  as.numeric()

clinical_data <- rbind(
  due_date,
  age,
  ethinic,
  bmi,
  parity,
  sex,ivf,
  induction,
  birth_wt
)

colnames(clinical_data) <- 
  colnames(metabolite_table)

##due_date, age, bmi, birth_wt
clinical_data <- 
  clinical_data[c("due_date", "age", "bmi", "birth_wt"),] %>% 
  apply(2, as.numeric)

rownames(clinical_data) <- 
  c("due_date", "age", "bmi", "birth_wt")


###calculate correlation matrix
# ###cross sectional network
# setwd("cross_sectional_network/")
# 
# ##patient information
# sfu1_148 <- 
#   readr::read_csv("E:/project/smartD/patient information/SFU1-148_GA.csv")
# 
# patient_info <- 
#   readr::read_csv("E:/project/smartD/data_analysis20190828/patient_info/patient_info.csv")
# 
# sfu1_148 <- 
#   sfu1_148 %>% 
#   mutate(subject_id = as.character(subject_id), visit = as.character(visit)) %>% 
#   arrange(subject_id)
# 
# patient_info <- 
#   patient_info %>% 
#   mutate(Patient_ID = as.character(Patient_ID), Visit = as.character(Visit)) %>% 
#   arrange(Patient_ID)

####data for correlation network analysis
cross_metabolite_table <- metabolite_table

####log 10 and scale
cross_metabolite_table <- 
  log(cross_metabolite_table, 10) %>% 
  as.data.frame(cross_metabolite_table)

cross_metabolite_table <-
  apply(cross_metabolite_table, 1, function(x){
    (x - mean(x))/sd(x)
  })

cross_metabolite_table <- 
  cross_metabolite_table %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "Sample_name")

colnames(cross_metabolite_table)[-1] <- 
  metabolite_tags$Compound.name

clinical_data <- 
t(clinical_data) %>% 
  as.data.frame()

rownames(clinical_data) == cross_metabolite_table$Sample_name

cross_table <-
  cbind(clinical_data, cross_metabolite_table) %>% 
  select(Sample_name, everything())


# ##add the subejct ID to each sample
# ###patient and sampel information
# sample_info_191021 <- 
#   readr::read_csv("E:/project/smartD/patient information/sample_info_191021.csv")
# 
# ##the postpurm samples may have no patient inforamtion
# metabolite_table_cross <- 
#   metabolite_table_cross %>% 
#   inner_join(sample_info_191021[,c(1:5)], by = c("Sample_name" = "Sample_ID"))
# 
# metabolite_table_cross <- 
#   metabolite_table_cross %>% 
#   select(Sample_name, Patient_ID:GA, everything())
# 
# write.csv(metabolite_table_cross, "metabolite_table_cross.csv", 
#           row.names = FALSE)

###calculate correlation
library(plyr)
# metabolite_table_cross2 <-
#   metabolite_table_cross %>% 
#   select(-c(Sample_name))
  # plyr::dlply(., .(Patient_ID)) %>% 
  # lapply(., function(x){
  #   x <- 
  #     x %>% 
  #     select(., -c(Patient_ID)) %>% 
  #     apply(., 2, mean)
  # }) %>% 
  # do.call(rbind, .) %>% 
  # as_tibble()
  
  
cross_cor <- 
  cross_table %>% 
  select(-Sample_name) %>% 
  as.matrix() %>% 
  cor(., method = "spearman", use = "na.or.complete")

rownames(cross_cor) <- 
  colnames(cross_cor)

##set lower tri as o
cross_cor[lower.tri(cross_cor)] <- 2

cross_cor <- 
  as_tibble(cross_cor)

rownames(cross_cor) <- 
  colnames(cross_cor)

cross_cor <- 
  cross_cor %>% 
  rownames_to_column(., var = "name1") %>% 
  gather(., key = "name2", value = "Correlation", -name1) %>% 
  distinct() %>% 
  filter(., Correlation != 1 & Correlation != 2) %>% 
  arrange(., desc(abs(Correlation))) %>% 
  filter(abs(Correlation) > 0.5)

cross_cor_p <-
  pbapply::pbapply(cross_cor, 1, function(x) {
    peak1 <- as.character(x[1])
    peak2 <- as.character(x[2])
    int1 <-
      pull(cross_table, var = peak1)
    int2 <-
      pull(cross_table, var = peak2)
    p <-  cor.test(int1, int2, alternative = "two.sided",
                   method = "spearman",  use = "na.or.complete")$p.value
    p
  })

sum(cross_cor_p < 0.05)
###p adjustment
cross_cor_p2 <- 
  p.adjust(cross_cor_p, method = "BH")

cross_cor <- 
  cross_cor %>% 
  mutate(P.adjust = cross_cor_p2) %>% 
  filter(P.adjust < 0.05 & abs(Correlation) > 0.6)

dim(cross_cor)


###############################################################################
###############################################################################
###############################################################################
cross_cor <- 
  cross_cor %>% 
  dplyr::rename(from = name1, to = name2)

######circos figure
# library(circlize)
# cross_correlation_data2 <- 
#   cross_correlation_data[1:100,]
# 
# ggplot(data = cross_correlation_data2, aes(x = Peak.name1, y = Peak.name2)) +
#   geom_tile(aes(fill = Correlation)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   theme_bw()
# 
# Metabolite <- 
#   unique(c(cross_correlation_data2$Compound.name1, cross_correlation_data2$Compound.name2))
# Metabolite
# cross_temp_data <- data.frame(factor = "Metabolomics", 
#                         x = 1:length(Metabolite),
#                         y = 1,
#                         Metabolite = Metabolite)
# 
# 
# ##use circlize
# 
# library(circlize)
# par(mar = c(2,2,2,2))
# circos.par(start.degree = 75, clock.wise = TRUE, gap.after = 6.3)
# circos.initialize(factors = cross_temp_data$factor, x = cross_temp_data$x)
# circos.track(factors = cross_temp_data$factor, 
#              x = cross_temp_data$x, 
#              y = cross_temp_data$y - 1,
# 
#              ylim = c(0, 1),
#              bg.border = NA,
#              bg.col = NA,
#              track.height = 0.4,
#              # track.margin = c(0,0,0,0),
#              panel.fun = function(x, y){
#                circos.points(x = x, y = y, pch = 19, col = "black")
#                circos.text(x = x, y = y + uy(2, "mm"), 
#                            niceFacing = TRUE, 
#                            labels = cross_temp_data$Metabolite,
#                            facing = "clockwise", 
#                            adj = c(0, 0.5),
#                            col = "black", cex = 0.8)
#   
# })
# 
# col_fun = colorRamp2(c(-1, 0, 1), c("#3C5488FF", "white", "#E64B35FF"))
# 
# for(i in 1:nrow(cross_correlation_data2)){
#   cat(i, " ")
#   point1 <- 
#     match(cross_correlation_data[i,]$Compound.name1, cross_temp_data$Metabolite)
#   point2 <- 
#     match(cross_correlation_data[i,]$Compound.name2, cross_temp_data$Metabolite)
#   correlation <- cross_correlation_data[i,]$Correlation
#     col <- col_fun(x = correlation)
#   circos.link(sector.index1 = "Metabolomics", point1 = point1, 
#               sector.index2 = "Metabolomics", point2 = point2, 
#               h = 1, col = col, lwd = 2)
# }
# 
# circos.clear()


###----------------------------------------------------------------------------
#top 100 correlation network
##use ggraph
library(ggraph)
library(tidygraph)

cross_cor$P.adjust[which(cross_cor$P.adjust == 0)] <- 
min(cross_cor$P.adjust[cross_cor$P.adjust != 0])


range(-log(cross_cor$P.adjust, 10))

cross_cor2 <- 
  cross_cor[1:100,]

metabolite <- 
  unique(c(cross_cor2$from, cross_cor2$to))

node_attr <- 
  data.frame(name = metabolite, stringsAsFactors = FALSE) %>% 
  left_join(metabolite_tags, by = c("name" = "Compound.name"))

node_attr$node_class <- rep(NA, nrow(node_attr))

node_attr$node_class <- 
  case_when(is.na(node_attr$Total.score) ~ "Clinical information",
            TRUE ~ 'Metabolite')


node_attr$node_class[node_attr$node_class == "Metabolite"] <- 
  node_attr$super_class[node_attr$node_class == "Metabolite"]

node_attr$node_class[is.na(node_attr$node_class)] <- "Other"

node_attr$Total.score[is.na(node_attr$Total.score)] <- 
  mean(node_attr$Total.score, na.rm = TRUE)

##add metaboliet class information to the node

node_attr$sort <- 
  case_when(
    node_attr$node_class == "Clinical information" ~ 1,
    node_attr$node_class == "Benzenoids" ~ 2,
    node_attr$node_class == "Lipids and lipid-like molecules" ~ 3,
    node_attr$node_class == "Organic acids and derivatives" ~ 4,
    node_attr$node_class == "Organic nitrogen compounds" ~ 5,
    node_attr$node_class == "Organic oxygen compounds" ~ 6,
    node_attr$node_class == "Organoheterocyclic compounds" ~ 7,
    node_attr$node_class == "Phenylpropanoids and polyketides" ~ 8,
    node_attr$node_class == "Other" ~ 9
    )

node_attr <- node_attr %>% 
  arrange(sort)

node_attr$node_class <- factor(node_attr$node_class, 
                               levels = unique(node_attr$node_class))

cross_graph <-
  tidygraph::tbl_graph(nodes = node_attr, 
                       edges = cross_cor2, 
                       directed = FALSE)


##angle for label
metabolite <- igraph::V(cross_graph)$name
id = 1:length(metabolite)
angle <- 360 * (id - 0.5)/length(metabolite)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

set_graph_style(family = "Arial")
extrafont::loadfonts()
plot <- 
ggraph(cross_graph, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(edge_colour = Correlation, 
                    edge_width = -log(P.adjust, 10))) +
  scale_edge_colour_gradient2(low = "#155F83FF", 
                              mid = "white",
                              high = "#800000FF") +
  scale_edge_width_continuous(range = c(0.2,1.2)) +
  geom_node_point(aes(size = Total.score, colour = node_class)) +
  scale_colour_manual(values = c("Clinical information" = "#FF6F00FF",
                                 "Benzenoids" = "#C71000FF", 
                                 "Lipids and lipid-like molecules" = "#008EA0FF",
                                 "Organic acids and derivatives" = "#8A4198FF",
                                 "Organic nitrogen compounds" = "#5A9599FF",
                                 "Organic oxygen compounds" = "#FF6348FF",
                                 "Organoheterocyclic compounds" = "#84D7E1FF",
                                 "Phenylpropanoids and polyketides" = "#FF95A8FF",
                                 "Other" = "#3D3B25FF")) +
  # ggsci::scale_color_futurama(alpha = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous(range = c(2,7),
                        guide = guide_legend(override.aes = list(colour = "black"))) +
  geom_node_text(aes(x = x * 1.05,
                     y = y * 1.05,
                     label = name, 
                     colour = node_class), 
                 angle = angle, 
                 hjust = hjust,
                 # colour = "black",
                 size = 3.5) +
  # ggdark::dark_theme_void() +
  theme_void() +
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
  
plot

# plot + ggdark::dark_theme_void()
ggsave(plot, filename = "top100_correlation_network.pdf")

######network community analysis
##https://github.com/IOR-Bioinformatics/PCSF/
##modularity定义是指网络中连接社区结构内部顶点的边所占的比例与另外一个随机网络中丽娜姐社区结构内部顶部
#的边所占比例的期望值相减得到的差值.
library(igraph)
###construct igraph project of correlation data
node_name <- unique(c(cross_cor$from, 
                      cross_cor$to))

node_attr <- 
  tibble(name = node_name) %>% 
  left_join(metabolite_tags, by = c("name" = "Compound.name"))

node_attr$node_class <- rep(NA, nrow(node_attr))

node_attr$node_class <- 
  case_when(is.na(node_attr$Total.score) ~ "Clinical information",
            TRUE ~ 'Metabolite')

node_attr$node_class[node_attr$node_class == "Metabolite"] <- 
  node_attr$super_class[node_attr$node_class == "Metabolite"]

node_attr$node_class[is.na(node_attr$node_class)] <- "Other"

node_attr$Total.score[is.na(node_attr$Total.score)] <- 
  mean(node_attr$Total.score, na.rm = TRUE)

##add metaboliet class information to the node

node_attr$sort <- 
  case_when(
    node_attr$node_class == "Clinical information" ~ 1,
    node_attr$node_class == "Alkaloids and derivatives" ~ 2,
    node_attr$node_class == "Benzenoids" ~ 3,
    node_attr$node_class == "Lipids and lipid-like molecules" ~ 4,
    node_attr$node_class == "Nucleosides, nucleotides, and analogues" ~ 5,
    node_attr$node_class == "Organic acids and derivatives" ~ 6,
    node_attr$node_class == "Organic nitrogen compounds" ~ 7,
    node_attr$node_class == "Organic oxygen compounds" ~ 8,
    node_attr$node_class == "Organoheterocyclic compounds" ~ 9,
    node_attr$node_class == "Organosulfur compounds" ~ 10,
    node_attr$node_class == "Phenylpropanoids and polyketides" ~ 11,
    node_attr$node_class == "Other" ~ 12
  )

node_attr <- node_attr %>% 
  arrange(sort)

node_attr$node_class <- factor(node_attr$node_class, 
                               levels = unique(node_attr$node_class))

node_attr$Total.score[is.na(node_attr$Total.score)] <- 
  mean(node_attr$Total.score, na.rm = TRUE)

cross_graph <- 
  graph_from_data_frame(d = cross_cor, 
                        directed = FALSE, 
                        vertices = node_attr
                        )

# coords <- layout.davidson.harel(cross_correlation_data_igraph)
library(ggraph)

plot <- 
ggraph(cross_graph, 
       layout = "auto", circular = TRUE) +
  geom_edge_link2(aes(edge_colour = Correlation, 
                      edge_width = -log(P.adjust, 10))) +
  scale_edge_colour_gradient2(low = "#155F83FF", 
                              mid = "white", 
                              high = "#800000FF") +
  scale_edge_width_continuous(range = c(0.8,1.5)) +
  geom_node_point(aes(size = Total.score,
                      colour = node_class)) +
  scale_colour_manual(values = c("Clinical information" = "#FF6F00FF",
                                 "Alkaloids and derivatives" = "#ADE2D0FF",
                                 "Benzenoids" = "#C71000FF", 
                                 "Lipids and lipid-like molecules" = "#008EA0FF",
                                 "Nucleosides, nucleotides, and analogues" = "#1A5354FF",
                                 "Organic acids and derivatives" = "#8A4198FF",
                                 "Organic nitrogen compounds" = "#5A9599FF",
                                 "Organic oxygen compounds" = "#FF6348FF",
                                 "Organoheterocyclic compounds" = "#84D7E1FF",
                                 "Organosulfur compounds" = "#3F4041FF",
                                 "Phenylpropanoids and polyketides" = "#FF95A8FF",
                                 "Other" = "#3D3B25FF")) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous(range = c(2,7),
                        guide = guide_legend(override.aes = list(colour = "black"))) +
  # geom_node_text(aes(x = x * 1.05,
  #                    y = y * 1.05,
  #                    label = name), 
  #                angle = angle, 
  #                hjust = hjust,
  #                colour = "black",
  #                size = 3.5) +
  theme_void() 
 # expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
plot

ggsave(plot, filename = "whole_correlation_network.pdf")

# lay <- layout.auto(graph = cross_correlation_data_igraph)
# plot(cross_correlation_data_igraph, 
#      vertex.shape = "circle",
#      vertex.label = NA,
#      vertex.color = "#FFA319FF",
#      vertex.size = vertex_attr(cross_correlation_data_igraph, "SS") * 10,
#      # edge.color = "#767676FF",
#      edge.color = ifelse(edge_attr(cross_correlation_data_igraph, "Correlation") > 0,
#                          "#80000099", "#155F8399"),
#      edge.width = abs(edge_attr(cross_correlation_data_igraph, "Correlation")) * 5,
#      layout = lay
#      )

##legend
# size_legend <-  
# quantile(x = vertex_attr(cross_correlation_data_igraph, "SS") * 10 %>% 
#              sort(), 
#            probs = c(0, 0.5, 1)) %>% 
#   round()
#   
# legend("topleft", 
#        legend = size_legend, 
#        bty = "n", 
#        pt.cex = size_legend/3,
#        pch = 19, 
#        # col = "black",
#        col = "#FFA319FF")

cross_subnetworks <- 
  cluster_edge_betweenness(graph = cross_graph, 
                           weights = abs(edge_attr(cross_graph,
                                                   "Correlation")))

par(mar = c(5,5,4,2))

plot <- 
ggplot(
  data.frame(index = 1:length(cross_subnetworks$modularity),
             modu = cross_subnetworks$modularity, stringsAsFactors = FALSE),
  aes(index, modu) 
) +
  geom_vline(xintercept = which.max(cross_subnetworks$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line(colour = "black") +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

plot <-
  plot + 
  ggplot2::annotate(geom = "point", 
                    x = which.max(cross_subnetworks$modularity),
                    y = max(cross_subnetworks$modularity), 
                    size = 3, 
                    colour = "#FFA319FF") +
    annotate(geom = "text", 
             x = which.max(cross_subnetworks$modularity),
             y = max(cross_subnetworks$modularity), 
             label = paste("(",  which.max(cross_subnetworks$modularity),
                           ",", 
                           max(cross_subnetworks$modularity) %>% round(3),
                           ")"),
             size = 5,
             colour = "#FFA319FF"
             )
  
plot  
ggsave(plot, filename = "modularity.pdf", width = 7, height = 7)

membership(cross_subnetworks)
modularity(cross_subnetworks)
length(cross_subnetworks)
sizes(cross_subnetworks)
algorithm(cross_subnetworks)
merges(cross_subnetworks)
crossing(communities = cross_subnetworks, graph = cross_correlation_data_igraph)
code_len(cross_subnetworks)
# is_hierarchical(cross_subnetworks)
# plot(as.dendrogram(cross_subnetworks))
# plot(as.hclust(cross_subnetworks))

which(sort(table(membership(communities = cross_subnetworks)),
           decreasing = TRUE) > 2)
which(sort(table(membership(communities = cross_subnetworks)), 
           decreasing = TRUE) > 5)

library(ggraph)
#> Loading required package: ggplot2
library(tidygraph)

temp_id <-
  sizes(cross_subnetworks) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 5) %>% 
  pull(., subnetwork_ID)

temp_subnetwork <-
lapply(temp_id, function(x){
  subgraph(graph = cross_graph, 
           v = names(membership(cross_subnetworks))[membership(cross_subnetworks) == x])
})

###node number
node_number <- 
lapply(temp_subnetwork, function(x){
vertex_attr(graph = x, name = "name") %>% 
    length()
}) %>% 
  unlist() 
  # sum()

###edge number
lapply(temp_subnetwork, function(x){
as_edgelist(x) %>% nrow()
}) %>% 
  unlist() %>% 
  mean()


sizes(cross_subnetworks) %>% 
  as_tibble() %>% 
  dplyr::rename(., subnetwork_ID = `Community sizes`, size = n) %>% 
  arrange(., desc(size)) %>% 
  filter(., size > 5)


###there are 12 cross_subnetworkss that with node bigger than 5.
##the first subnetwork
cross_subnetwork <- 
  subgraph(graph = cross_graph, 
           v = names(membership(cross_subnetworks))[membership(cross_subnetworks) == 9])

cross_subnetwork

metabolite <- igraph::V(cross_subnetwork)$name
id = 1:length(metabolite)
angle <- 360 * (id - 0.5)/length(metabolite)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)


plot <- 
  ggraph(cross_subnetwork, layout = "auto", circular = FALSE) +
  geom_edge_link(aes(edge_colour = Correlation, 
                    edge_width = -log(P.adjust, 10))) +
  scale_edge_colour_gradient2(low = "#155F83FF", 
                              mid = "white",
                              high = "#800000FF") +
  scale_edge_width_continuous(range = c(0.2,1.2)) +
  geom_node_point(aes(size = Total.score, colour = node_class)) +
  scale_colour_manual(values = c("Clinical information" = "#FF6F00FF",
                                 "Alkaloids and derivatives" = "#ADE2D0FF",
                                 "Benzenoids" = "#C71000FF", 
                                 "Lipids and lipid-like molecules" = "#008EA0FF",
                                 "Nucleosides, nucleotides, and analogues" = "#1A5354FF",
                                 "Organic acids and derivatives" = "#8A4198FF",
                                 "Organic nitrogen compounds" = "#5A9599FF",
                                 "Organic oxygen compounds" = "#FF6348FF",
                                 "Organoheterocyclic compounds" = "#84D7E1FF",
                                 "Organosulfur compounds" = "#3F4041FF",
                                 "Phenylpropanoids and polyketides" = "#FF95A8FF",
                                 "Other" = "#3D3B25FF")) +
  # ggsci::scale_color_futurama(alpha = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous(range = c(2,7),
                        guide = guide_legend(override.aes = list(colour = "black"))) +
  geom_node_text(aes(x = x * 1.05,
                     y = y * 1.05,
                     label = name, 
                     colour = node_class), 
                 # angle = angle, 
                 # hjust = hjust,
                 size = 3.5,
                 repel = TRUE) +
  # ggdark::dark_theme_void() +
  theme_void() +
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
# expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
plot


save(plot, file = "cross_subnetwork9.pdf")

#pathway enrichment for each cross_subnetwork
cross_subnetwork_node_info <-
  vertex_attr(cross_subnetwork) %>% 
  bind_cols()

###calculate the betweness and degree and all other
cross_betweenness <- betweenness(graph = cross_subnetwork)
cross_degree <- igraph::degree(graph = cross_subnetwork)
cross_closeness <- closeness(graph = cross_subnetwork)

cross_importance <- 
  tibble(cross_betweenness, cross_degree, cross_closeness) %>% 
  mutate_all(., function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  mutate(.,mean = apply(., 1 ,mean)) %>% 
  mutate(name = vertex_attr(graph = cross_subnetwork, name = "name")) %>% 
  select(name, everything()) %>% 
  arrange(., desc(mean))

metabolite_tags$super_class[is.na(metabolite_tags$super_class)] <- 
  "Other"


cols <- c("Clinical information" = "#FF6F00FF",
          "Alkaloids and derivatives" = "#ADE2D0FF",
          "Benzenoids" = "#C71000FF", 
          "Lipids and lipid-like molecules" = "#008EA0FF",
          "Nucleosides, nucleotides, and analogues" = "#1A5354FF",
          "Organic acids and derivatives" = "#8A4198FF",
          "Organic nitrogen compounds" = "#5A9599FF",
          "Organic oxygen compounds" = "#FF6348FF",
          "Organoheterocyclic compounds" = "#84D7E1FF",
          "Organosulfur compounds" = "#3F4041FF",
          "Phenylpropanoids and polyketides" = "#FF95A8FF",
          "Other" = "#3D3B25FF")

temp_data <- 
cross_importance[c(1:30),] %>% 
  arrange(., mean) %>% 
  left_join(metabolite_tags, by = c("name" = "Compound.name")) %>% 
  mutate(name = factor(name, levels = name))

temp_data %>% 
  ggplot(.,aes(x = name, y = mean)) +
  geom_point(size = 2, aes(colour = super_class)) +
  scale_colour_manual(values = c("Clinical information" = "#FF6F00FF",
                                 "Alkaloids and derivatives" = "#ADE2D0FF",
                                 "Benzenoids" = "#C71000FF", 
                                 "Lipids and lipid-like molecules" = "#008EA0FF",
                                 "Nucleosides, nucleotides, and analogues" = "#1A5354FF",
                                 "Organic acids and derivatives" = "#8A4198FF",
                                 "Organic nitrogen compounds" = "#5A9599FF",
                                 "Organic oxygen compounds" = "#FF6348FF",
                                 "Organoheterocyclic compounds" = "#84D7E1FF",
                                 "Organosulfur compounds" = "#3F4041FF",
                                 "Phenylpropanoids and polyketides" = "#FF95A8FF",
                                 "Other" = "#3D3B25FF")) +
  labs(x = "", 
       y = "Mean (Betweenness, Degree, Closeness)") +
  geom_segment(aes(x = name, y = 0, 
                   xend = name, 
                   yend = mean,
                   colour = super_class)) +
  # scale_colour_manual(values = c("Yes" = "#FFA319FF", "No" = "#767676FF")) +
  theme_bw() +
  coord_flip() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 10, colour = cols[temp_data$super_class]), 
        legend.position = c(1, 0),
        legend.justification = c(1,0),
        axis.text.x = element_text(), legend.background = element_blank())

###betweenness, degree and closness
# test <- 
cross_importance[c(1:30),] %>% 
  select(., name:mean) %>% 
  arrange(., mean) %>% 
  gather(., key = "Class", value = "Value", -name) %>% 
  mutate(name = factor(name, levels = name[1:30])) %>% 
  filter(., Class != "mean") %>% 
  ggplot(.,aes(x = name, 
               y = Value)) +
  geom_point(size = 2, 
             aes(colour = Class)
             # colour = "#155F83FF"
             ) +
  geom_segment(aes(x = name, y = 0, 
                   xend = name, 
                   yend = Value,
                   colour = Class)
               # colour = "#155F83FF"
               ) +
  scale_color_manual(values = c(
    "cross_betweenness" = "#FFA319FF",
    "cross_closeness" = "#8A9045FF",
    "cross_degree" = "#C16622FF"
  )) +
  facet_grid(. ~ Class, 
             labeller = as_labeller(
               c("cross_betweenness" = "Betweenness",
                 "cross_closeness" = "Closeness",
                 "cross_degree" = "Degree"
                 )
             )
             ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 10), legend.position = "none", 
        axis.text.x = element_text(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 13)) +
  coord_flip()

###find the nodes with top 3 imporance
cross_importance[1,]$name

# cross_temp_subgraph <- 
# neighbors(graph = cross_subnetwork11, 
#           v = cross_importance11[1,]$peak.name) %>% 
#   names() %>% 
#   `c`(., cross_importance11[1,]$peak.name) %>% 
#   unique() %>% 
#   subgraph(cross_subnetwork11, v = .)
# 
# 
# cross_temp_subgraph <- 
# set_vertex_attr(graph = cross_temp_subgraph, name = "label", 
#                 value = ifelse(vertex_attr(cross_temp_subgraph, "name") == cross_importance11[1,]$peak.name, "Yes", "No"))
# 
# ##get top 20 correlations
# E(cross_temp_subgraph)
# cross_temp_subgraph <- 
# as_long_data_frame(cross_temp_subgraph) %>% 
#   select(., from_name, to_name, Correlation) %>%
#   filter(., from_name == cross_importance11[1,]$peak.name | to_name == cross_importance11[1,]$peak.name) %>% 
#   mutate(Correlation = abs(Correlation)) %>% 
#   arrange(., desc(Correlation)) %>% 
#   top_n(., 20)
#   
#   
# cross_temp_subgraph <- 
#   subgraph(graph = cross_subnetwork11, 
#            v = unique(c(cross_temp_subgraph$from_name, cross_temp_subgraph$to_name)))
# 
# 
# cross_temp_subgraph <- 
#   set_vertex_attr(graph = cross_temp_subgraph, name = "label", 
#                   value = ifelse(vertex_attr(cross_temp_subgraph, "name") == cross_importance11[1,]$peak.name, "Yes", "No"))
# 
# which(vertex_attr(cross_temp_subgraph, "name") == cross_importance11[1,]$peak.name)
# 
# lay <- layout_as_tree(cross_temp_subgraph, root = 18)
# 
#   plot(cross_temp_subgraph,
#        vertex.color = ifelse(vertex_attr(cross_temp_subgraph, "label")=="Yes", "#FFA31999", "#8A904599"),
#        vertex.label = vertex_attr(cross_temp_subgraph, "Compound.name"),
#        vertex.label.color = "black",
#        vertex.size = vertex_attr(cross_temp_subgraph, "SS") * 20,
#        vertex.label.dist = 1,
#        
#        # edge.color = "#767676FF",
#        edge.color = ifelse(edge_attr(cross_temp_subgraph, "Correlation") > 0,
#                            "#80000099", "#155F8399"),
#        edge.width = abs(edge_attr(cross_temp_subgraph, "Correlation")) * 5,
#        layout = -lay[, 2:1]
#        )
  
 
  
  
  
  
  
  
  
  
  
  






  








