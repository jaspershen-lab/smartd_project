library(KEGGgraph)
library(KEGGREST)
library(KEGGlincs)
library(tidyverse)


FoxO_KGML <- get_KGML("map00564")
slot(FoxO_KGML, "pathwayInfo")

image_link <- slot(slot(FoxO_KGML, "pathwayInfo"), "image")
download.file(url = image_link, destfile = basename(image_link), mode = "wb")



##get the ID of all hsa pathways
KEGGREST::listDatabases()
path_ID <- 
KEGGREST::keggList(organism = "hsa", database = "pathway") %>% 
  names() %>% 
  unique() %>% 
  stringr::str_replace_all(., "path:", "")

dir.create(path = "kegg_pathway")
setwd("kegg_pathway/")


node_info <- vector(mode = "list", length = length(path_ID))
edge_info <- vector(mode = "list", length = length(path_ID))

for(i in i:length(path_ID)){
  cat(i, " ")
  temp_kgml <- 
    get_KGML(pathwayid = path_ID[i])
  if(is.na(temp_kgml)){
    next()
  }
  # slot(object = temp_kgml, name = "pathwayInfo")
  #Download a static pathway image (png file) to working directory:
  temp_image_link <- slot(slot(temp_kgml, "pathwayInfo"), "image")
  download.file(temp_image_link, basename(temp_image_link), mode = "wb")
  temp_kegg_mappings <- expand_KEGG_mappings(KGML_file = temp_kgml)
  temp_kegg_edges <- expand_KEGG_edges(temp_kgml, temp_kegg_mappings)
  #Modify existing data sets; specify as nodes and edges
  temp_node_mapping_info <- node_mapping_info(temp_kegg_mappings)
  temp_edge_mapping_info <- edge_mapping_info(temp_kegg_edges)
  #Create an igraph object
  # kegg_graph <- get_graph_object(temp_node_mapping_info, temp_edge_mapping_info)
  
  edge_info[[i]] <- 
    temp_kegg_edges
  node_info[[i]] <- temp_node_mapping_info

  temp_kegg_nodes <- temp_node_mapping_info

}


save(node_info, file = "node_info")
save(edge_info, file = "edge_info")



path_name <- KEGGREST::keggList(organism = "hsa", database = "pathway") %>% 
  unname()

path_name_id <- 
  paste(path_ID, path_name, sep = ";")


for(i in 1:length(path_name_id)){
  cat(i, " ")
  if(!is.null(node_info[[i]])){
    temp <- node_info[[i]]
    temp[,3] <- lapply(temp[,3], function(x) paste(x, collapse = ";")) %>% unlist()
    temp[,4] <- lapply(temp[,4], function(x) paste(x, collapse = ";")) %>% unlist()
    temp[,13] <- lapply(temp[,13], function(x) paste(x, collapse = ";")) %>% unlist()
    temp[,14] <- lapply(temp[,14], function(x) paste(x, collapse = ";")) %>% unlist()
    readr::write_csv(temp[,-c(6,7,8,9,10,11,12,13,17,18)], paste(path_ID[i], "_node_info.csv", sep = ""))
  }
  
  if(!is.null(edge_info[[i]])){
    temp <- edge_info[[i]]
    # temp[,3] <- lapply(temp[,3], function(x) paste(x, collapse = ";")) %>% unlist()
    # temp[,4] <- lapply(temp[,4], function(x) paste(x, collapse = ";")) %>% unlist()
    # temp[,13] <- lapply(temp[,13], function(x) paste(x, collapse = ";")) %>% unlist()
    # temp[,14] <- lapply(temp[,14], function(x) paste(x, collapse = ";")) %>% unlist()
    readr::write_csv(temp, paste(path_ID[i], "_edge_info.csv", sep = ""))
  }
  
  
}


for(i in 1:length(path_ID)){
  cat(i, " ")
  file <- grep(path_ID[i], dir(), value = TRUE)
  dir.create(path_ID[i])
  file.copy(from = file, to = path_ID[i]) 
}

unlink(grep("csv", dir(), value = TRUE))

library(igraph)


as_edgelist(GO)



vertex_attr(GO)



library(tidyverse)
library(ggforce)
library(nycflights13)

p <- airports %>%
  filter(lon < 0, tzone != "\\N") %>%
  ggplot(aes(lon, lat, color = tzone)) + 
  geom_point(show.legend = FALSE)  




