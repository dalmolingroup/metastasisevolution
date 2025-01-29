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
download_if_missing("https://stringdb-static.org/download/COG.mappings.v11.0.txt.gz")

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
  "ENSP00000239462", "KOG1225",     
  "ENSP00000265131", "KOG1225",     
  "ENSP00000265562", "KOG2220", 
  "ENSP00000359085", "KOG3512", 
  "ENSP00000361467", "KOG3528",    
  "ENSP00000418112", "COG5599"    
)

# Removing unresolved cases and adding manual assignments
gene_cogs %<>%
  filter(n == 1) %>%
  dplyr:: select(-n) %>%
  bind_rows(gene_cogs_resolved)

# Exporting for package use
#gene_cogs |> 
#  vroom::vroom_write(file = here("results/orthology_data/gene_cogs.csv"), delim = ",")

#map_ids |> 
#  vroom::vroom_write(file = here("results/orthology_data/map_ids.csv"), delim = ",")

#cogs |> 
#  vroom::vroom_write(file = here("results/orthology_data/cogs.csv"), delim = ",")

# Get proteins interaction
string_edgelist <- get_network_interaction(map_ids)

# Recomputing scores
string_edgelist <- combine_scores(string_edgelist, 
                                  evidences = c("ascore", "escore", "dscore"), 
                                  confLevel = 0.4)

colnames(string_edgelist) <- c("stringId_A", "stringId_B", "combined_score")

# Remove o species id
string_edgelist$stringId_A <- substring(string_edgelist$stringId_A, 6, 1000)
string_edgelist$stringId_B <- substring(string_edgelist$stringId_B, 6, 1000)

# How many edgelist proteins are absent in gene_ids? (should return 0)
setdiff(
  string_edgelist %$% c(stringId_A, stringId_B),
  map_ids %>% pull(stringId)
) 

# Exporting for package use
# string_edgelist |> 
#   vroom::vroom_write(file = here("results/orthology_data/string_edgelist.csv"), delim = ",")

## Run GeneBridge
ogr <- newBridge(ogdata=ogdata, phyloTree=phyloTree, ogids = gene_cogs$cog_id, refsp="9606")

ogr <- runBridge(ogr, penalty = 2, threshold = 0.5, verbose = TRUE)

ogr <- runPermutation(ogr, nPermutations=1000, verbose=FALSE)

res <- getBridge(ogr, what="results")

# res |> 
#   vroom::vroom_write(file = here("results/orthology_data/genebridge_result.csv"), delim = ",")

#save(ogr, file = "../results/orthology_data/genebridge_ogr.RData")

## Naming the rooted clades and getting the final results table
CLADE_NAMES <- "https://raw.githubusercontent.com/dalmolingroup/neurotransmissionevolution/ctenophora_before_porifera/analysis/geneplast_clade_names.tsv"

lca_names <- read_table(CLADE_NAMES)

groot_df <- res %>%
  tibble::rownames_to_column("cog_id") %>%
  dplyr::select(cog_id, root = Root) %>%
  inner_join(lca_names) %>%
  inner_join(gene_cogs) 

# groot_df |> 
#   vroom::vroom_write(file = here("results/orthology_data/groot_df.csv"), delim = ",")
# 
# ## Create
 nodelist <- data.frame(node = unique(c(string_edgelist$stringId_A, string_edgelist$stringId_B)))
# 
 merged_paths <- merge(nodelist, groot_df, by.x = "node", by.y = "string_id")
# 
# merged_paths |> 
#   vroom::vroom_write(file = here("results/orthology_data/merged_paths.csv"), delim = ",")

net <- get_network_interaction(merged_paths, "node")
net <- combine_scores(net, evidences = c("ascore", "escore", "dscore"), confLevel = 0.4)

net <-  net %>%
  separate(stringId_A,
           into = c("ncbi_taxon_id", "stringId_A"),
           sep = "\\.") %>%
  separate(stringId_B,
           into = c("ncbi_taxon_id", "stringId_B"),
           sep = "\\.") 

network_filtered <- net %>%
  dplyr::select(stringId_A, stringId_B) |>
  distinct()

pivotada <- df %>% 
  dplyr::select(external_gene_name, Signature) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::pivot_wider(
    id_cols = external_gene_name,
    names_from = Signature,
    values_from = n,
    values_fn = list(n = length),
    values_fill = list(n = 0),
  )

source_statements <-
  colnames(pivotada)[2:length(pivotada)]

nodelist <-
  data.frame(node = unique(c(network_filtered$stringId_A, network_filtered$stringId_B))) %>%
  left_join(merged_paths, by = c("node" = "node")) %>%
  left_join(map_ids, by = c("node" = "stringId")) %>%
  left_join(pivotada, by = c("queryItem" = "external_gene_name"))

# Network Metrics
connected_nodes <- rle(sort(c(network_filtered[,1], network_filtered[,2])))
connected_nodes <- data.frame(count=connected_nodes$lengths, node=connected_nodes$values)
connected_nodes <- left_join(nodelist, connected_nodes, by = c("node" = "node"))

# nodelist |> 
#   vroom::vroom_write(file = here("results/orthology_data/nodelist.csv "), delim = ",")
# 
# connected_nodes |> 
#   vroom::vroom_write(file = here("results/orthology_data/connected_nodes.csv"), delim = ",")
# 
# network_filtered |> 
#   vroom::vroom_write(file = here("results/orthology_data/network_filtered.csv"), delim = ",")

################################ PLOTING ROOTS #######################################################

library(ggplot2)
library(ggraph)
library(dplyr)
library(tidyr)
library(igraph)
library(purrr)
library(vroom)
library(paletteer)
library(easylayout)
library(UpSetR)
library(tinter)

color_mappings <- c(
  "cell adhesion involved in sprouting angiogenesis"   = "#06141FFF"
  ,"cell adhesion mediated by integrin"                = "#742C14FF"
  ,"cell adhesion mediator activity"                   = "#3D4F7DFF"
  ,"cell-substrate adhesion"                           = "#E48C2AFF"
  ,"negative regulation of cell adhesion"              ="#72874EFF"
  ,"positive regulation of cell adhesion"              ="#046E8FFF"
  ,"protein complex involved in cell adhesion"         = "red"
  ,"regulation of cell adhesion"                       = "green"
)

subset_graph_by_root <-
  function(geneplast_result, root_number, graph) {
    filtered <- geneplast_result %>%
      filter(root >= root_number) %>%
      pull(node)
    
    induced_subgraph(graph, which(V(graph)$name %in% filtered))
  }

adjust_color_by_root <- function(geneplast_result, root_number, graph) {
  filtered <- geneplast_result %>%
    filter(root == root_number) %>%
    pull(node)
  
  V(graph)$color <- ifelse(V(graph)$name %in% filtered, "black", "gray")
  return(graph)
}

# Configure graph collors by genes incrementation
subset_and_adjust_color_by_root <- function(geneplast_result, root_number, graph) {
  subgraph <- subset_graph_by_root(geneplast_result, root_number, graph)
  adjusted_graph <- adjust_color_by_root(geneplast_result, root_number, subgraph)
  return(adjusted_graph)
}

plot_network <- function(graph, title, nodelist, xlims, ylims, legend = "none") {
  
  # Generate color map
  source_statements <-
    colnames(nodelist)[10:length(nodelist)]
  
  color_mappings <- c(
    "cell adhesion involved in sprouting angiogenesis"   = "#06141FFF"
    ,"cell adhesion mediated by integrin"                = "#742C14FF"
    ,"cell adhesion mediator activity"                   = "#3D4F7DFF"
    ,"cell-substrate adhesion"                           = "#E48C2AFF"
    ,"negative regulation of cell adhesion"              ="#72874EFF"
    ,"positive regulation of cell adhesion"              ="#046E8FFF"
    ,"protein complex involved in cell adhesion"         = "red"
    ,"regulation of cell adhesion"                       = "green"
  )
  
  vertices <- igraph::as_data_frame(graph, "vertices")
  
  ggraph:: ggraph(graph,
                  "manual",
                  x = V(graph)$x,
                  y = V(graph)$y) +
    ggraph::geom_edge_link0(edge_width = 0.2, color = "#90909020") +
    ggraph::geom_node_point(ggplot2::aes(color = I(V(graph)$color)), size = 0.5) +
    scatterpie::geom_scatterpie(
      aes(x=x, y=y, r=18),
      cols = source_statements,
      data = vertices[rownames(vertices) %in% V(graph)$name[V(graph)$color == "black"],],
      colour = NA,
      pie_scale = 1
    ) +
    geom_node_text(aes(label = ifelse(V(graph)$color == "black", V(graph)$queryItem, NA)), 
                   nudge_x = 1, nudge_y = 1, size = 0.5, colour = "#BFBEBF") +
    ggplot2::scale_fill_manual(values = color_mappings, drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(limits = xlims) +
    ggplot2::scale_y_continuous(limits = ylims) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend,
      legend.key.size = ggplot2::unit(0.5, 'cm'),
      legend.key.height = ggplot2::unit(0.5, 'cm'),
      legend.key.width = ggplot2::unit(0.5, 'cm'),
      legend.title = ggplot2::element_text(size=6),
      legend.text = ggplot2::element_text(size=6),
      panel.border = ggplot2::element_rect(
        colour = "#161616",
        fill = NA,
        linewidth = 1
      ),
      plot.title = ggplot2::element_text(size = 8, face = "bold")
    ) +
    ggplot2::guides(
      color = "none",
      fill = "none"
    ) +
    ggplot2::labs(fill = "Source:", title = title)
}

source_statements <- colnames(nodelist)[10:length(nodelist)]

upset(dplyr::select(as.data.frame(nodelist), 
                   "queryItem", 
                   "cell adhesion involved in sprouting angiogenesis"
                   ,"cell adhesion mediated by integrin"
                   ,"cell adhesion mediator activity"
                   ,"cell-substrate adhesion"
                   ,"negative regulation of cell adhesion"
                   ,"positive regulation of cell adhesion"
                   ,"protein complex involved in cell adhesion"
                   ,"regulation of cell adhesion"),
     nsets = 50, nintersects = NA,
     #sets.bar.color = c("#06141FFF", "#72874EFF", "#3D4F7DFF",
     #                   "#742C14FF","#046E8FFF", "#E48C2AFF"), 
     mainbar.y.label = "Biological Process \nIntersections",
     sets.x.label = "Set Size")

graph <-
  graph_from_data_frame(network_filtered, directed = FALSE, vertices = nodelist)

#layout <- easylayout::easylayout(graph)
#V(graph)$x <- layout[, 1]
#V(graph)$y <- layout[, 2]

#save(graph, file="_dev/graph_celladhesion.Rdata")

ggraph(graph, "manual", x = V(graph)$x, y = V(graph)$y) +
  geom_edge_link0(color = "#90909020") +
  geom_node_point(aes(color = -root), size = 2) +
  theme_void() +
  theme(legend.position = "left")



geneplast_roots <- merged_paths |>
  arrange(root)

buffer <- c(-50, 50)
xlims <- ceiling(range(V(graph)$x)) + buffer
ylims <- ceiling(range(V(graph)$y)) + buffer

roots <- unique(geneplast_roots$root) %>%
  set_names(unique(geneplast_roots$clade_name))

# Subset graphs by LCAs
subsets <-
  map(roots, ~ subset_and_adjust_color_by_root(geneplast_roots, .x, graph))

# Plot titles
titles <- names(roots)

plots <-
  map2(
    subsets,
    titles,
    plot_network,
    nodelist = nodelist,
    xlims = xlims,
    ylims = ylims,
    legend = "right"
  ) %>%
  discard(is.null)

design <-
  '123'

#net_all_roots <-
patchwork::wrap_plots(
  plots$Metamonada, plots$Choanoflagellata, plots$Actinopterygii,
  design = design,
  nrow = 6,
  ncol = 6
)














calculate_cumulative_genes <- function(nodelist) {
  
  # Obter todas as categorias possíveis de clade_name
  all_clades <- node_annotation %>%
    arrange(desc(root)) %>%
    dplyr:: select(clade_name) %>%
    unique()
  
  # Definir as colunas de interesse
  process_columns <- c("queryItem", "root", "clade_name", 
                       "cell adhesion involved in sprouting angiogenesis"
                       ,"cell adhesion mediated by integrin"
                       ,"cell adhesion mediator activity"
                       ,"cell-substrate adhesion"
                       ,"negative regulation of cell adhesion"
                       ,"positive regulation of cell adhesion"
                       ,"protein complex involved in cell adhesion"
                       ,"regulation of cell adhesion")
  
  # Calcular o cumulativo agrupando por clade_name
  cumulative_genes <- nodelist %>%
    arrange(desc(root)) %>%
    dplyr::select(all_of(process_columns)) %>%
    group_by(clade_name, root) %>%
    summarise(count_genes = n(), .groups = "drop") %>%
    arrange(desc(root)) %>%
    mutate(cumulative_sum = cumsum(count_genes)) %>%
    right_join(all_clades, by = "clade_name") %>%
    fill(cumulative_sum, .direction = "down")
  
  return(cumulative_genes)
}

calculate_cumulative_bp <- function(nodelist) {
  
  # Obter todas as categorias possíveis de clade_name
  all_clades <- node_annotation %>%
    arrange(desc(root)) %>%
    dplyr:: select(clade_name) %>%
    unique()
  
  # Definir as colunas de interesse
  process_columns <- c("queryItem", "root", "clade_name", 
                       "cell adhesion involved in sprouting angiogenesis"
                       ,"cell adhesion mediated by integrin"
                       ,"cell adhesion mediator activity"
                       ,"cell-substrate adhesion"
                       ,"negative regulation of cell adhesion"
                       ,"positive regulation of cell adhesion"
                       ,"protein complex involved in cell adhesion"
                       ,"regulation of cell adhesion")
  
  # Calcular a soma cumulativa para cada processo biológico
  cumulative_bp <- nodelist %>%
    dplyr::select(all_of(process_columns)) %>%
    distinct(root, queryItem, .keep_all = TRUE) %>%
    mutate(across(all_of(process_columns[-c(1:3)]), ~ as.numeric(.))) %>%
    group_by(root, clade_name) %>%
    summarise(across(all_of(process_columns[-c(1:3)]), 
                     ~ sum(. , na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(desc(root)) %>%
    mutate(across(all_of(process_columns[-c(1:3)]), ~ cumsum(.))) %>%
    right_join(all_clades, by = "clade_name") %>%
    fill(everything(), .direction = "down")
  
  return(cumulative_bp)
}

node_annotation <- nodelist %>%
  inner_join(gene_cogs, by = c("node" = "string_id", "cog_id")) %>%
  inner_join(df, by = c("queryItem" = "external_gene_name")) %>%
  distinct(queryItem, cog_id, Signature, root, clade_name)

cumulative_genes <- calculate_cumulative_genes(nodelist) 
cumulative_bp <- calculate_cumulative_bp(nodelist)

cumulative_data <- left_join(cumulative_genes, cumulative_bp)


long_data <- cumulative_data %>%
  pivot_longer(cols = 5:12, 
               names_to = "Process", 
               values_to = "Value")

#a <-
ggplot() +
  # Gráfico de barras para cumulative_sum
  geom_bar(data = cumulative_data, 
           aes(x = factor(clade_name, levels = clade_name), y = cumulative_sum), 
           stat = "identity", fill = "darkgray", colour = NA) +
  geom_text(data = cumulative_data, 
            aes(x = factor(clade_name, levels = clade_name), y = cumulative_sum, label = cumulative_sum), 
            vjust = -0.5, size = 3, color = "darkgray") +
  scale_color_manual(values = color_mappings) +
  
  labs(x = "Clade Name", y = "Cumulative Sum", 
       title = "Cumulative Sum and Biological Processes",
       fill = "Cumulative Sum",
       color = "Biological Processes") +
  
  theme_main +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#b <- 
ggplot() +
  # Gráfico de barras para cumulative_sum
  geom_bar(data = cumulative_data, 
           aes(x = factor(-root), y = cumulative_sum), 
           stat = "identity", fill = "darkgray", colour = NA) +
  geom_text(data = cumulative_data, 
            aes(x = factor(-root), y = cumulative_sum, label = cumulative_sum), 
            vjust = -0.5, linewidth = 3, color = "darkgray") +
  
  # Gráfico de linhas para os processos biológicos
  geom_line(data = long_data, 
            aes(x = factor(-root), y = Value, color = Process, group = Process), 
            linewidth = 1) +
  
  # Usar a paleta de cores definida
  scale_color_manual(values = color_mappings) +
  
  labs(x = "Clade Name", y = "Cumulative Sum", 
       title = "Cumulative Sum and Biological Processes",
       fill = "Cumulative Sum",
       color = "Biological Processes") +
  
  theme_main +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
