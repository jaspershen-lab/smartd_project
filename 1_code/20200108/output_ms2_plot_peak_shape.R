sxtTools::setwd_project()
setwd("data_analysis20200108/")
rm(list = ls())
load("data_preparation_for_analysis/metabolite_table")
load("data_preparation_for_analysis/metabolite_tags")
load("data_preparation_for_analysis/peak_table")


sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/GA_prediction/")
marker1 <- readr::read_csv("marker_rf_final.csv")

sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/time_to_due_prediction/remove_cs/")
marker2 <- readr::read_csv("marker_rf_final.csv")


marker <- rbind(marker1, marker2) %>% 
  distinct(name, .keep_all = TRUE)

marker <-
  marker %>% 
  dplyr::mutate(Model = case_when(
    (name %in% marker1$name) & (name %in% marker2$name) ~ "GA&Sampling time to delivery",
    (name %in% marker1$name) & !(name %in% marker2$name) ~ "GA",
    !(name %in% marker1$name) & (name %in% marker2$name) ~ "Sampling time to delivery"
  ))


# marker <- 
#   marker %>% 
#   dplyr::left_join(metabolite_tags, by = "name")

xlsx::write.xlsx(x = marker, file = "marker_paper.xlsx")

sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/all_marker_ms2_shape/peak_shape/")
xlsx::write.xlsx2(x = marker %>% select(name, mz, rt) %>% as.data.frame(), 
            file = "marker.xlsx", row.names = FALSE)


peak_data_pos <-
  metflow2::extractPeaks(
    path = ".",
    ppm = 15,
    threads = 4,
    rt.tolerance = 20,
    is.table = "marker.xlsx"
  )



plot <- vector(mode = "list", length = nrow(marker))
for(i in 1:nrow(marker)){
  cat(i, " ")
  plot[[i]] <- metflow2::showPeak(object = peak_data_pos, 
                                  peak.index = i, 
                                  alpha = 0.5, 
                                  interactive = FALSE)
  
  plot[[i]] <- 
    plot[[i]] +
    labs(title = marker$Compound.name.x[i]) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 7),
          axis.ticks = element_blank()) 
}


library(cowplot)
plot <- 
plot_grid(plot[[1]], plot[[2]],
          plot[[3]], plot[[4]],
          plot[[5]], plot[[6]],
          plot[[7]], plot[[8]],
          plot[[9]], plot[[10]],
          plot[[11]], plot[[12]],
          plot[[13]], plot[[14]],
          plot[[15]], plot[[16]],
          plot[[17]], plot[[18]],
          plot[[19]], plot[[20]],
          plot[[21]], plot[[22]],
          plot[[23]], plot[[24]],
          labels = letters[1:24], 
          ncol = 6)

plot

ggsave(plot, filename = "peak_shape.pdf", width = 7, height = 7)

ggsave(plot, filename = "peak_shape.png", width = 7, height = 7)


###out put MS2 plot
sxtTools::setwd_project()
setwd("data_analysis20200108/prediction/metabolites/RF/all_marker_ms2_shape/MS2_plot/")
load("result.pRPLC.nce25")
load("result.pRPLC.nce50")
load("HMDB.metabolite.data")
load("hmdbDatabase0.0.1")
load("massbankDatabase0.0.1")
load("metlinDatabase0.0.1")
load("monaDatabase0.0.1")
load("msDatabase_hilic0.0.1")
load("msDatabase_rplc0.0.1")
load("nistDatabase0.0.1")
load("orbitrapDatabase0.0.1")

i <- 1
temp_name <- 
marker$name[i]
temp_compund <- marker$Compound.name.x[i]
temp_compund
temp_database <- marker$Database.x[i]
temp_idx <- which(names(result.pRPLC.nce25) == temp_database)

plot <- 
metID::ms2plot(object = result.pRPLC.nce50[[temp_idx]], 
               database = get(temp_database), 
               which.peak = temp_name)

plot <- plot +
  labs(title = paste(temp_name, temp_database, temp_compund, sep = "/")) +
  theme(plot.title = element_text(size = 15))

plot
ggsave(filename = paste(temp_name, "pdf", sep = "."), width = 8, height = 6)


##dark
plot <- 
  metID::ms2plot(object = result.pRPLC.nce50[[temp_idx]], 
                 database = get(temp_database), 
                 which.peak = temp_name, col.exp = "white")

plot <- 
  plot +
  labs(title = paste(temp_name, temp_database, temp_compund, sep = "/")) +
  ggdark::dark_theme_bw() +
  theme(plot.title = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "grey")
        )
plot

plot <- plotly::ggplotly(plot)

ggsave(filename = paste(temp_name, "png", sep = "."), 
       width = 8, height = 6, bg = "transparent")

export(p, file='image.png')

















plot <- vector(mode = "list", length = 24)


i <- 24
temp_name <- 
  marker$name[i]
temp_compund <- marker$Compound.name.x[i]
temp_compund
temp_database <- marker$Database.x[i]
temp_idx <- which(names(result.pRPLC.nce25) == temp_database)

plot[[i]] <- 
  metID::ms2plot(object = result.pRPLC.nce50[[temp_idx]], 
                 database = get(temp_database), 
                 which.peak = temp_name)

# plot[[i]] <- plot[[i]] +
#   labs(title = paste(temp_name, temp_database, temp_compund, sep = "/")) +
#   theme(plot.title = element_text(size = 15))

plot[[i]]


save(plot, file = "plot")

for(i in 1:24){
  plot[[i]] <-
    plot[[i]] +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 8))
}



plot_new <- 
  plot_grid(plot[[1]], plot[[2]],
            plot[[3]], plot[[4]],
            plot[[5]], plot[[6]],
            plot[[7]], plot[[8]],
            plot[[9]], plot[[10]],
            plot[[11]], plot[[12]],
            labels = letters[1:12], 
            ncol = 3)

plot_new

ggsave(plot_new, filename = "MS2_plot1.pdf", width = 7, height = 7)


plot_new <- 
  plot_grid(plot[[13]], plot[[14]],
            plot[[15]], plot[[16]],
            plot[[17]], plot[[18]],
            plot[[19]], plot[[20]],
            plot[[21]], plot[[22]],
            plot[[23]], plot[[24]],
            labels = letters[13:24], 
            ncol = 3)

plot_new

ggsave(plot_new, filename = "MS2_plot2.pdf", width = 7, height = 7)


  
plot_new <- 
  plot_grid(plot[[1]], plot[[2]],
            plot[[3]], plot[[4]],
            plot[[5]], plot[[6]],
            plot[[7]], plot[[8]],
            plot[[9]], plot[[10]],
            plot[[11]], plot[[12]],
            plot[[13]], plot[[14]],
            plot[[15]], plot[[16]],
            plot[[17]], plot[[18]],
            plot[[19]], plot[[20]],
            plot[[21]], plot[[22]],
            plot[[23]], plot[[24]],
            labels = letters[1:24], 
            ncol = 6)

plot_new

ggsave(plot_new, filename = "MS2_plot.pdf", width = 7, height = 7)

ggsave(plot_new, filename = "MS2_plot.pdf", width = 7, height = 7)