setwd("/home/gleison/Documents/metastasisevolution")


library(tidyr)
library(dplyr)
library(httr)
library(vroom)
library(here)
library(UpSetR)

################################## CREATING DATABASE ###############################################
# Function to search genes by GOids filtering by exoerimental evidences
query_GO <- function(GO_id) {
  url <- paste0("https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class%3D%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22", GO_id,"%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&fq=type:%22protein%22&fq=evidence_subset_closure_label:%22experimental%20evidence%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=*%3A*")
  response <- tryCatch({
    vroom(url, delim = "\t", col_names = FALSE)
  }, error = function(e) {
    return(data.frame())
  })
  
  if (ncol(response) == 0) {
    return(character(0))
  }
  
  genes <- response[[1]] 
  return(unique(as.vector(genes)))
}

# Create a list with GOids and their respective genes 
get_genes_by_GOid <- function(go_ids) {
  genes_per_term <- list()
  for (go_id in go_ids) {
    genes_per_term[[go_id]] <- query_GO(go_id)
  }
  return(genes_per_term)
}

# Create gene data frame from list
list_to_df <- function(genesPerTerm, ids) {
  result <- data.frame(external_gene_name = character(), Signature = character(), stringsAsFactors = FALSE)
  
  for (go_id in names(genesPerTerm)) {
    genes <- genesPerTerm[[go_id]]
    signature <- ids$Description[ids$GO_id == go_id]
    
    if (length(signature) > 0 && length(genes) > 0) {
      temp <- data.frame(external_gene_name = genes, Signature = rep(signature, length(genes)), stringsAsFactors = FALSE)
      result <- rbind(result, temp)
    }
  }
  
  return(result)
}



cell_adhesion <- vroom::vroom("_dev/cell_adhesion.csv")
cell_adhesion <- separate(cell_adhesion, Signature, into = c("GO_id", "Description"), sep = 10, extra = "merge")

# Get genes per GOid
genesPerTerm <- get_genes_by_GOid(cell_adhesion$GO_id)

# Create gene data frame
result <- data.frame(external_gene_name = character(), Signature = character(), stringsAsFactors = FALSE)

result <- list_to_df(genesPerTerm, cell_adhesion) %>%
  unique()

result |> 
  vroom::vroom_write(file = here("_dev/cell_adhesion_genes.csv"), delim = ",")

head(result)

################################## ORTHOLOGY ROOTING ###############################################

library(GeneBridge)
library(geneplast.data)
library(readr)
library(dplyr)
library(purrr)
library(biomaRt)
library(magrittr)
library(KEGGREST)
library(ape)
library(tidyverse)
library(data.table)
library(stringi)
library(AnnotationHub)
library(sourcetools)
library(here)

# Function to download needed files
download_if_missing <- function(url, filename = basename(url)) {
  filename <- here::here("assets/", filename)
  if (!file.exists(filename)) {
    download.file(url, filename)
  }
}

# get IDs from STRING DB
get_string_ids <- function(genes_hgnc, species_id = "9606") {
  
  req <- RCurl::postForm(
    "https://string-db.org/api/tsv/get_string_ids",
    identifiers = paste(genes_hgnc, collapse = "%0D"),  
    echo_query = "1",
    species = species_id, 
    .opts = list(ssl.verifypeer = FALSE)
  )
  
  map_ids <- read.table(text = req, sep = "\t", header = TRUE, quote = "") %>%
    dplyr::select(-queryIndex) %>%
    unique()
  
  map_ids$stringId <- substring(map_ids$stringId, 6, 1000)
  
  return(map_ids)
}

# Get STRING interactions
get_network_interaction <- function(map_ids, stringId, species_id = "9606") {
  
  identifiers <- map_ids %>% pull(stringId) %>% na.omit %>% paste0(collapse="%0d") 
  
  req2 <- RCurl::postForm(
    "https://string-db.org/api/tsv/network",
    identifiers = identifiers, 
    required_core = "0", 
    species = species_id,
    .opts = list(ssl.verifypeer = FALSE)
  )
  
  int_network <- read.table(text = req2, sep = "\t", header = TRUE)
  
  int_network <- unique(int_network)
  
  return(int_network)
}

## Recomputing scores
combine_scores <- function(dat, evidences = "all", confLevel = 0.4) {
  if(evidences[1] == "all"){
    edat<-dat[,-c(1,2,ncol(dat))]
  } else {
    if(!all(evidences%in%colnames(dat))){
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat<-dat[,evidences]
  }
  if (any(edat > 1)) {
    edat <- edat/1000
  }
  edat<-1-edat
  sc<- apply(X = edat, MARGIN = 1, FUN = function(x) 1-prod(x))
  dat <- cbind(dat[,c(1,2)],combined_score = sc)
  idx <- dat$combined_score >= confLevel
  dat <-dat[idx,]
  return(dat)
}

# Metastasis Genes data frame
df <- vroom::vroom("_dev/cell_adhesion_genes.csv")

# NCBI Eukaryotes 
load('assets/string_eukaryotes.rda')

## Table with Orthologous Groups and their proteins
#download_if_missing("https://stringdb-static.org/download/COG.mappings.v11.0.txt.gz")

cogs <- fread(
  "assets/COG.mappings.v11.0.txt.gz",
  header           = F,
  stringsAsFactors = F,
  skip             = 1,
  sep              = "\t",
  col.names        = c("taxid.string_id","cog_id"),
  select           = c(1,4),
  quote            = ""
)

# Query Phylotree and OG data
ah <- AnnotationHub()
meta <- query(ah, "geneplast")
load(meta[["AH83116"]])
cogdata <- dplyr::select(cogdata, "protein_id", "ssp_id", "og_id" = "cog_id")

# Query by STRING API
map_ids <- get_string_ids(df$external_gene_name) 

# Spliting first column into taxid and string_id
separated_ids <- cogs %$% stri_split_fixed(taxid.string_id, pattern = ".", n = 2, simplify = T)

cogs[["taxid"]]     <- separated_ids[, 1]
cogs[["string_id"]] <- separated_ids[, 2]

# Freeing up some memory
rm(separated_ids)

# keeping only eukaryotes
cogs %<>% dplyr::select(-taxid.string_id) %>% 
  filter(taxid %in% string_eukaryotes[["taxid"]])
gc()

cogdata2 <- cogs %>% 
  dplyr::select(protein_id = string_id, ssp_id = taxid, og_id = cog_id)
cogdata2 <- as.data.frame(cogdata2)

ogdata <- unique(rbind(cogdata, cogdata2))

# Subsetting cogs of interest - METASTASIS GENES
gene_cogs <- cogs %>%
  filter(string_id %in% map_ids[["stringId"]]) %>%
  dplyr::select(-taxid) %>%
  group_by(string_id) %>%
  summarise(n = n(), cog_id = paste(cog_id, collapse = " / "))

## Proteins with multiple COGs
gene_cogs %>% filter(n > 1)

# Resolving main proteins
gene_cogs_resolved <- tribble(
  ~string_id, ~cog_id,
  "ENSP00000239462",	"KOG1225",
  "ENSP00000264808",	"KOG1721",	
  "ENSP00000265562",	"KOG2220",	
  "ENSP00000323856",	"COG5069",		
  "ENSP00000359085",  "KOG3512",		
  "ENSP00000361467",  "KOG3528",		
  "ENSP00000363317",	"KOG3528",		
  "ENSP00000383042",  "KOG0685",		
  "ENSP00000418112",	"COG5599",		
  "ENSP00000422533",	"KOG3545"
)

# Removing unresolved cases and adding manual assignments
gene_cogs %<>%
  filter(n == 1) %>%
  dplyr:: select(-n) %>%
  bind_rows(gene_cogs_resolved)

# Exporting for package use
gene_cogs |> 
  vroom::vroom_write(file = here("results/orthology_data/gene_cogs.csv"), delim = ",")

map_ids |> 
  vroom::vroom_write(file = here("results/orthology_data/map_ids.csv"), delim = ",")

cogs |> 
  vroom::vroom_write(file = here("results/orthology_data/cogs.csv"), delim = ",")