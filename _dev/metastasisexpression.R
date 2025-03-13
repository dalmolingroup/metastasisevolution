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

library("ComplexHeatmap")
library(dplyr)

metastasis_genes <- vroom::vroom("gene_publication.csv", delim = ",")
gene_counts <- metastasis_genes %>%
  count(gene) %>%
  filter(n > 1)

i <- get_string_ids(unique(metastasis_genes$gene))

nodelist <- vroom::vroom("../results/orthology_data/nodelist.csv")

sum(i$stringId %in% nodelist$node)

list_intersect <- nodelist[nodelist$preferredName %in% i$preferredName,]

inter <- list(paper = i$preferredName,
              mydata = nodelist$preferredName)

lt <- make_comb_mat(inter, mode = "intersect")
lt

UpSet(lt, comb_order = order(comb_size(lt)), top_annotation = upset_top_annotation(lt, add_numbers = TRUE),
      right_annotation = upset_right_annotation(lt, add_numbers = TRUE))


cumulative_genes <- calculate_cumulative_genes(list_intersect) 
cumulative_bp <- calculate_cumulative_bp(list_intersect)

cumulative_data <- left_join(cumulative_genes, cumulative_bp)


long_data <- cumulative_data %>%
  pivot_longer(cols = 5:10, 
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
  scale_color_manual(values = annotation_colors) +
  
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
            vjust = -0.5, size = 3, color = "darkgray") +
  
  # Gráfico de linhas para os processos biológicos
  geom_line(data = long_data, 
            aes(x = factor(-root), y = Value, color = Process, group = Process), 
            size = 1) +
  
  # Usar a paleta de cores definida
  scale_color_manual(values = annotation_colors) +
  
  labs(x = "Clade Name", y = "Cumulative Sum", 
       title = "Cumulative Sum and Biological Processes",
       fill = "Cumulative Sum",
       color = "Biological Processes") +
  
  theme_main +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

      