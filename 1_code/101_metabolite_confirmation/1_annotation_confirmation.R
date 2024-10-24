# setwd(r4projects::get_project_wd())
# setwd("data_analysis20200108/prediction/metabolites/RF/GA_prediction")
# 
# marker_rf <- readr::read_csv("marker_rf.csv")
# 
# setwd("annotation_confirmation/")
# load("hmdbDatabase0.0.1")
# load("massbankDatabase0.0.1")
# load("metlinDatabase0.0.1")
# load("monaDatabase0.0.1")
# load("msDatabase_rplc0.0.1")
# load("nistDatabase0.0.1")
# load("orbitrapDatabase0.0.1")
# load("result.pRPLC.nce50")
# load("result.pRPLC.nce25")
# load("result.nRPLC.nce50")
# load("result.nRPLC.nce25")
# 
# 
# ###POS
# database_name <-
#   marker_rf %>%
#   filter(stringr::str_detect(name, "POS")) %>%
#   pull(Database) %>%
#   unique()
# 
# for (i in 1:length(database_name)) {
#   cat(i, " ")
#   temp_database <- database_name[i]
#   idx <- lapply(result.pRPLC.nce25, function(x)
#     x@database) %>%
#     unlist() %>%
#     unname() %>%
#     `==`(temp_database) %>%
#     which()
#   temp_name <- marker_rf %>%
#     filter(Database == temp_database) %>%
#     filter(stringr::str_detect(name, "POS")) %>%
#     pull(name)
#   
#   if (length(temp_name) == 0) {
#     next()
#   } else{
#     metID::ms2plot(
#       object = result.pRPLC.nce25[[idx]],
#       which.peak = temp_name,
#       database = get(temp_database),
#       path = file.path(
#         "annotation_confirmation",
#         paste("POS", temp_database, "NCE25", sep = "_")
#       ),
#       show.plot = FALSE
#     )
#     
#     metID::ms2plot(
#       object = result.pRPLC.nce50[[idx]],
#       which.peak = temp_name,
#       database = get(temp_database),
#       path =  file.path(
#         "annotation_confirmation",
#         paste("POS", temp_database, "NCE50", sep = "_")
#       ),
#       ,
#       show.plot = FALSE
#     )
#   }
#   
# }
# 
# 
# 
# ###NEG
# database_name <-
#   marker_rf %>%
#   filter(stringr::str_detect(name, "NEG")) %>%
#   pull(Database) %>%
#   unique()
# 
# for (i in 1:length(database_name)) {
#   cat(i, " ")
#   temp_database <- database_name[i]
#   idx <- lapply(result.pRPLC.nce25, function(x)
#     x@database) %>%
#     unlist() %>%
#     unname() %>%
#     `==`(temp_database) %>%
#     which()
#   temp_name <- marker_rf %>%
#     filter(Database == temp_database) %>%
#     filter(stringr::str_detect(name, "NEG")) %>%
#     pull(name)
#   
#   if (length(temp_name) == 0) {
#     next()
#   } else{
#     metID::ms2plot(
#       object = result.nRPLC.nce25[[idx]],
#       which.peak = temp_name,
#       database = get(temp_database),
#       path = paste("NEG", temp_database, "NCE25", sep = "_")
#     )
#     
#     metID::ms2plot(
#       object = result.nRPLC.nce50[[idx]],
#       which.peak = temp_name,
#       database = get(temp_database),
#       path = paste("NEG", temp_database, "NCE50", sep = "_")
#     )
#   }
#   
# }
