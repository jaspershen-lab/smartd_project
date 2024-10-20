sxtTools::setwd_project()
setwd("data_analysis20200108/urine_metabolome/kegg_pathway/")

###pathway analysis for module
load("kegg_hsa_pathway_database")
path_id <- lapply(kegg_hsa_pathway_database, function(x){
  x$ENTRY
}) %>% 
  unlist() %>% 
  unname()

path_name <- lapply(kegg_hsa_pathway_database, function(x){
  x$NAME
}) %>% 
  unlist() %>% 
  unname()

hsa_pathway <- lapply(kegg_hsa_pathway_database, function(x){
  unique(names(x$COMPOUND))
})

names(hsa_pathway) <- paste(path_name, path_id, sep = ";")

##remove cancer pathway
path_class <- 
  lapply(kegg_hsa_pathway_database, function(x){
    x <- x$CLASS
    if(is.null(x)){
      return(NA)
    }
    return(x)
  })  %>% 
  unlist()

remove_idx1 <- which(is.na(path_class))
remove_idx2 <- 
  stringr::str_detect(string = path_class, pattern = "Human Diseases;") %>% 
  which()

remove_idx <- unique(c(remove_idx1, remove_idx2))

hsa_pathway <- hsa_pathway[-remove_idx]

#remove null
remove_idx <- 
  lapply(hsa_pathway, is.null) %>% 
  unlist() %>% 
  which()

if(length(remove_idx) > 0){
  hsa_pathway <- hsa_pathway[-remove_idx]
}

save(hsa_pathway, file = "hsa_pathway")



###disease pathway
sxtTools::setwd_project()
setwd("data_analysis20200108/urine_metabolome/kegg_pathway/")

###pathway analysis for module
load("kegg_hsa_pathway_database")
path_id <- lapply(kegg_hsa_pathway_database, function(x){
  x$ENTRY
}) %>% 
  unlist() %>% 
  unname()

path_name <- lapply(kegg_hsa_pathway_database, function(x){
  x$NAME
}) %>% 
  unlist() %>% 
  unname()

hsa_disease_pathway <- lapply(kegg_hsa_pathway_database, function(x){
  unique(names(x$COMPOUND))
})

names(hsa_disease_pathway) <- paste(path_name, path_id, sep = ";")

##remain cancer pathway
path_class <- 
  lapply(kegg_hsa_pathway_database, function(x){
    x <- x$CLASS
    if(is.null(x)){
      return(NA)
    }
    return(x)
  })  %>% 
  unlist()

remove_idx1 <- which(is.na(path_class))
remove_idx2 <- 
  stringr::str_detect(string = path_class, pattern = "Human Diseases;") %>% 
  `!`() %>% 
  which()

remove_idx <- unique(c(remove_idx1, remove_idx2))

hsa_disease_pathway <- hsa_disease_pathway[-remove_idx]

#remove null
remove_idx <- 
  lapply(hsa_disease_pathway, is.null) %>% 
  unlist() %>% 
  which()

if(length(remove_idx) > 0){
  hsa_disease_pathway <- hsa_disease_pathway[-remove_idx]
}

save(hsa_disease_pathway, file = "hsa_disease_pathway")
