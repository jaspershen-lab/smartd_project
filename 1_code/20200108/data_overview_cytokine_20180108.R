#-------------------------------------------------------------------------------
##RPLC pos
#-------------------------------------------------------------------------------
setwd(r4projects::get_project_wd())
rm(list=ls())
setwd("data_analysis20200108/data_preparation_for_analysis/cytokine/")
load("cytokine_pheno")
load("cytokine_table")
load("cytokine_tags")

setwd(r4projects::get_project_wd())
setwd("patient information/")
sample_info <- readr::read_csv("sample_info_191021.csv")
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_overview/cytokine/")
library(tidyverse)

###log
subject_data <- 
  cytokine_table %>% 
  select(one_of(sample_info$Sample_ID))

sample_info <-
  sample_info %>% 
  filter(Sample_ID %in% colnames(subject_data))
  
colnames(subject_data) == sample_info$Sample_ID

subject_data <- 
  log(subject_data)

subject_data <- 
  t(
    apply(subject_data, 1, function(x){
      (x - mean(x))/sd(x)
    })
  )


subject_data <-
  subject_data %>% 
  as_tibble()

subject_data2 <- t(subject_data)
subject_data2 <- as_tibble(subject_data2)
rownames(subject_data2)

subject_data2 <- 
  subject_data2 %>% 
  mutate(GA = sample_info$GA)



##PCA analysis
temp_data <- 
  subject_data2 %>% 
  select(-c(GA))

pca_object <- 
  prcomp(x = 
           temp_data)

library(ggfortify)

x <- pca_object$x

x <- x[,1:2]

x <- data.frame(x, 
                GA = sample_info$GA, 
                stringsAsFactors = FALSE)

ggplot(x[x$GA != 0,], aes(PC1, PC2, colour = GA)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 5) +
  # scale_colour_gradient(low = alpha("#155F83FF", 0.1),
  #                       high = alpha("#FFA319FF", 1)) +
  guides(colour = guide_colourbar(title = "GA (week)")) +
  scale_colour_gradientn(colours = c(
    alpha("#155F83FF", 1),
    alpha("#155F83FF", 0.4),
    alpha("#FFA319FF", 0.4),
    alpha("#FFA319FF", 1)
  )) +
  # scale_colour_brewer() +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        # legend.position = "top",
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15))+
annotate(geom = "point", x = x$PC1[x$GA == 0], y = x$PC1[x$GA == 0], colout = "black")




