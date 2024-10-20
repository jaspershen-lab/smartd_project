#get the identification results
#####for metabolite identification result
##load data
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_cleaning/RPLC/POS/")
load("rplc_pos_6")

setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_cleaning/RPLC/NEG/")
load("rplc_neg_6")

##get the identification results
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/metabolite_identification/RPLC/POS/")
identification_table_pos <- readr::read_csv("identification.table.new.csv")

setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/metabolite_identification/RPLC/NEG/")
identification_table_neg <- readr::read_csv("identification.table.new.csv")

##first set the work directory to project folder
setwd(r4projects::get_project_wd())
setwd("data_analysis20200108/data_preparation_for_analysis/")

identification_table_pos <- 
  identification_table_pos %>% 
  select(name:rt, MS2.spectrum.name:Database)


identification_table_neg <- 
  identification_table_neg %>% 
  select(name:rt, MS2.spectrum.name:Database)

metabolite_table <- 
  rbind(identification_table_pos, identification_table_neg)

###we only use the level 1 and level 2 identifications
#1-3 means identification confidence level, 4 means MS2 spectra match si not good enough, and 5 means this compound
#may be not in human body
metabolite_table <- 
  metabolite_table %>% 
  filter(Level == 1 | Level == 2)

dim(metabolite_table)

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
  pbapply::pblapply(metabolite_tags$Compound.name, function(x){
    metflow2::transID(query = x, from = "Chemical name", to = "KEGG", top = 1)
  }) %>% 
  do.call(rbind, .)

hmdb.id <- 
  pbapply::pblapply(metabolite_tags$Compound.name, function(x){
    metflow2::transID(query = x, from = "Chemical name",
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
  pbapply::pblapply(function(x){
    metflow2::transID(query = x, from = "Human Metabolome Database", to = "InChIKey", top = 1)
  }) %>% 
  do.call(rbind, .)

inchikey2 <- 
  metabolite_tags$KEGG.ID %>% 
  pbapply::pblapply(function(x){
    metflow2::transID(query = x, from = "KEGG", to = "InChIKey", top = 1)
  }) %>% 
  do.call(rbind, .)

inchikey3 <- 
  metabolite_tags$name %>% 
  pbapply::pblapply(function(x){
    metflow2::transID(query = x, from = "Chemical name", to = "InChIKey", top = 1)
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
  pbapply::pblapply(inchikey, function(x){
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
                              class, 
                              sub_class,
                              stringsAsFactors = FALSE)

###get the metabolite table
peak_table_pos <- rplc_pos_6@ms1.data[[1]]
peak_table_neg <- rplc_neg_6@ms1.data[[1]]

##remove Blank and QC_DL
peak_table_pos <- 
  peak_table_pos %>% 
  select(-dplyr::starts_with('BLK')) %>% 
  select(-dplyr::starts_with('blk'))


peak_table_neg <- 
  peak_table_neg %>% 
  select(-dplyr::starts_with('BLK')) %>% 
  select(-dplyr::starts_with('blk'))


setdiff(colnames(peak_table_pos), colnames(peak_table_neg))

setdiff(colnames(peak_table_neg), colnames(peak_table_pos))

intersect_name <- 
  intersect(colnames(peak_table_pos), colnames(peak_table_neg))

peak_table_pos <- 
  peak_table_pos %>% 
  select(one_of(intersect_name))

peak_table_neg <- 
  peak_table_neg %>% 
  select(one_of(intersect_name))


colnames(peak_table_pos) == colnames(peak_table_neg)

peak_table <- 
  rbind(peak_table_pos, peak_table_neg)


readr::write_csv(peak_table, "peak_table.csv")
save(peak_table, file = "peak_table")


metabolite_table <- 
  peak_table %>% 
  filter(name %in% metabolite_tags$name) %>% 
  arrange(name)

metabolite_tags <- 
  metabolite_tags %>% 
  filter(name %in% metabolite_table$name) %>% 
  arrange(name)

metabolite_table$name == metabolite_tags$name

readr::write_csv(metabolite_table, "metabolite_table.csv")
save(metabolite_table, file = "metabolite_table")

readr::write_csv(metabolite_tags, "metabolite_tags.csv")
save(metabolite_tags, file = "metabolite_tags")
# metabolite_table <- metabolite_table %>%
#   select(contains("SFU"))

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

