# Instale pacotes necessários se ainda não os tiver
# install.packages("igraph")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("enrichplot")
# BiocManager::install("biomaRt")

library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)

# Suponha que 'subsets' seja sua lista de objetos igraph
# subsets <- list(
#   # Exemplo de dados
#   Dermoptera = igraph::graph_from_literal(ENSP00000005260-ENSP00000379350, ENSP00000005260-ENSP00000366283),
#   Afrotheria = igraph::graph_from_literal(ENSP00000231449-ENSP00000252999, ENSP00000239462-ENSP00000258341)
# )

# Função para extrair proteínas de um objeto igraph
extract_proteins <- function(graph) {
  V(graph)$name
}

# Função para traduzir ENSP para ENSG
translate_ensp_to_ensg <- function(proteins) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensp_to_ensg <- getBM(filters = "ensembl_peptide_id", 
                        attributes = c("ensembl_peptide_id", "ensembl_gene_id"), 
                        values = proteins, 
                        mart = mart)
  ensp_to_ensg <- setNames(ensp_to_ensg$ensembl_gene_id, ensp_to_ensg$ensembl_peptide_id)
  ensg_ids <- ensp_to_ensg[proteins]
  ensg_ids[!is.na(ensg_ids)]
}

# Lista para armazenar os resultados do enriquecimento
enrichment_results <- list()

# Loop para realizar o enriquecimento funcional em cada subconjunto
for (name in names(subsets)) {
  proteins <- extract_proteins(subsets[[name]])
  genes <- translate_ensp_to_ensg(proteins)
  
  # Realize o enriquecimento funcional (GO BP)
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05)
  
  enrichment_results[[name]] <- ego
}

# Crie um dataframe combinado para evolução temática
combined_results <- do.call(rbind, lapply(names(enrichment_results), function(name) {
  result <- enrichment_results[[name]]
  if (is.null(result)) return(NULL)
  result_df <- as.data.frame(result)
  result_df$Clade <- name
  result_df
}))

combined_results  <- inner_join(combined_results, lca_names, by = c("Clade" = "clade_name")) 




# Ordene lca_names de acordo com a sequência desejada
ordered_clades <- lca_names$clade_name
#ordered_clades <- c(ordered_clades[length(ordered_clades)], ordered_clades[-length(ordered_clades)]) # Move "Metamonada" para o primeiro lugar

# Transforme Clade em um fator com a ordem especificada
combined_results$Clade <- factor(combined_results$Clade, levels = ordered_clades)

# Seleciona as colunas necessárias
df <- combined_results[, c("Clade", "Description", "qvalue")]

# Verifica a estrutura dos dados
head(df)

# Converte os dados para o formato longo usando reshape2
df_melt <- melt(df, id.vars = c("Clade", "Description"), measure.vars = "qvalue")

# Crie o heatmap usando ggplot2
ggplot(df_melt, aes(x = Description, y = Clade, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = "Heatmap de Enriquecimento Funcional KEGG",
       x = "Clado",
       y = "Termo KEGG",
       fill = "q-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 0.8))




























# Função para traduzir ENSEMBL Peptide IDs para Entrez Gene IDs
translate_ensp_to_entrez <- function(proteins) {
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensp_to_entrez <- getBM(filters = "ensembl_peptide_id", 
                          attributes = c("ensembl_peptide_id", "entrezgene_id"), 
                          values = proteins, 
                          mart = mart)
  ensp_to_entrez <- setNames(ensp_to_entrez$entrezgene_id, ensp_to_entrez$ensembl_peptide_id)
  ensg_ids <- ensp_to_entrez[proteins]
  ensg_ids[!is.na(ensg_ids)]
}

## Enriquecimento funcional KEGG
library(clusterProfiler)

kegg_enrichment_results <- list()

for (name in names(subsets)) {
  proteins <- extract_proteins(subsets[[name]])
  genes <- translate_ensp_to_entrez(proteins)
  
  # Realize o enriquecimento funcional (KEGG)
  kegg_enrich <- enrichKEGG(
    gene = genes,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH"
  )
  
  kegg_enrichment_results[[name]] <- kegg_enrich
}

# Crie um dataframe combinado para evolução temática
combined_results_kegg <- do.call(rbind, lapply(names(kegg_enrichment_results), function(name) {
  result <- kegg_enrichment_results[[name]]
  if (is.null(result)) return(NULL)
  result_df <- as.data.frame(result)
  result_df$Clade <- name
  result_df
}))

# Supondo que combined_results e lca_names já existam e tenham as colunas corretas
combined_results_kegg <- inner_join(combined_results_kegg, lca_names, by = c("Clade" = "clade_name"))
combined_results_kegg2 <- subset(combined_results_kegg, category != "Human Diseases" )                                                                        

combined_results_kegg2 <- a


library(ggplot2)
library(reshape2)
library(dplyr)

# Supondo que combined_results_kegg e lca_names sejam os seus dataframes

# Ordene lca_names de acordo com a sequência desejada
ordered_clades <- lca_names$clade_name
#ordered_clades <- c(ordered_clades[length(ordered_clades)], ordered_clades[-length(ordered_clades)]) # Move "Metamonada" para o primeiro lugar

# Transforme Clade em um fator com a ordem especificada
combined_results_kegg2$Clade <- factor(combined_results_kegg2$Clade, levels = ordered_clades)

# Seleciona as colunas necessárias
df <- combined_results_kegg2[, c("Clade", "Description", "qvalue")]

# Verifica a estrutura dos dados
head(df)

# Converte os dados para o formato longo usando reshape2
df_melt <- melt(df, id.vars = c("Clade", "Description"), measure.vars = "qvalue")

# Crie o heatmap usando ggplot2
ggplot(df_melt, aes(x = Description, y = Clade, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = "Heatmap de Enriquecimento Funcional KEGG",
       x = "Clado",
       y = "Termo KEGG",
       fill = "q-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
        




combined_results_kegg <- inner_join(combined_results_kegg, lca_names, by = c("Clade" = "clade_name"))
combined_results_kegg2 <- subset(combined_results_kegg, category != "Human Diseases")











combined_results_GO <- subset(combined_results, ID %in% metastasis_ids$GO_id)

library(reshape2)
library(dplyr)
library(ggplot2)

# Ordene lca_names de acordo com a sequência desejada
ordered_clades <- lca_names$clade_name

# Transforme Clade em um fator com a ordem especificada
combined_results_GO$Clade <- factor(combined_results_GO$Clade, levels = ordered_clades)

# Seleciona as colunas necessárias
df <- combined_results_GO[, c("Clade", "Description", "qvalue")]

# Contar o número de clados em que cada descrição aparece
description_counts <- df %>%
  group_by(Description) %>%
  summarize(CladeCount = n_distinct(Clade)) %>%
  arrange(CladeCount, Description)

# Reordenar as descrições com base no número de clados
df$Description <- factor(df$Description, levels = description_counts$Description)

# Converte os dados para o formato longo usando reshape2
df_melt <- melt(df, id.vars = c("Clade", "Description"), measure.vars = "qvalue")

# Crie o heatmap usando ggplot2
ggplot(df_melt, aes(x = Description, y = Clade, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = "Heatmap de Enriquecimento Funcional KEGG",
       x = "Termo KEGG",
       y = "Clado",
       fill = "q-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        
        
        
# Instale e carregue as bibliotecas necessárias
install.packages("reshape2")
install.packages("tidyr")
library(reshape2)
library(dplyr)
library(tidyr)

# Supondo que combined_results_kegg e lca_names sejam os seus dataframes

# Ordene lca_names de acordo com a sequência desejada
ordered_clades <- lca_names$clade_name
ordered_clades <- c(ordered_clades[length(ordered_clades)], ordered_clades[-length(ordered_clades)]) # Move "Metamonada" para o primeiro lugar

# Filtre apenas os clados presentes em combined_results_kegg
valid_clades <- intersect(ordered_clades, unique(combined_results_kegg$Clade))

# Transforme Clade em um fator com a ordem especificada
combined_results_kegg$Clade <- factor(combined_results_kegg$Clade, levels = valid_clades)

# Crie uma matriz onde as linhas são os termos KEGG, as colunas são os clados e os valores são os q-values
heatmap_data <- combined_results_kegg %>%
  dplyr::select(Clade, Description, qvalue) %>%
  spread(Clade, qvalue)

# Verifique se todas as colunas estão na ordem correta e complete com NA onde necessário
heatmap_data <- heatmap_data %>%
dplyr::  select(Description, valid_clades)

# Converta a coluna Description para row names e remova-a do dataframe
rownames(heatmap_data) <- heatmap_data$Description
heatmap_data <- heatmap_data %>% dplyr::select(-Description)

# Converta para matriz
heatmap_matrix <- as.matrix(heatmap_data)

# Verifique as dimensões da matriz e os nomes das colunas
dim(heatmap_matrix)
colnames(heatmap_matrix)

# Crie o heatmap
heatmap(heatmap_matrix, Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "red"))(100),
        scale = "none", margins = c(10, 10), xlab = "Clado", ylab = "Termo KEGG",
        main = "Heatmap de Enriquecimento Funcional KEGG")








library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)

# Cores fornecidas
annotation_colors <- c(
  "cell adhesion"                               = "#06141FFF", 
  "extracellular matrix organization"           = "#742C14FF", 
  "epithelial to mesenchymal transition"        = "#3D4F7DFF",
  "mesenchymal to epithelial transition"        = "#CD4F38FF",
  "regulation of metallopeptidase activity"     = "#E48C2AFF",
  "cell junction organization"                  = "#72874EFF",
  "cellular extravasation"                      = "#046E8FFF"
)

# Preparar os dados
combined_results_GO <- subset(combined_results, ID %in% metastasis_ids$GO_id)

# Ordene lca_names de acordo com a sequência desejada
ordered_clades <- lca_names$clade_name

# Transforme Clade em um fator com a ordem especificada
combined_results_GO$Clade <- factor(combined_results_GO$Clade, levels = ordered_clades)

# Seleciona as colunas necessárias
df <- combined_results_GO[, c("Clade", "Description", "qvalue", "ID")]


# Adicione as assinaturas da tabela metastasis_ids
df <- left_join(df, metastasis_ids, by = c("ID" = "GO_id"))
df <- dplyr:: select(df, "Clade", "ID", Description = "Description.x", "qvalue", "Hallmark", "Signature")

# Contar o número de clados em que cada descrição aparece
description_counts <- df %>%
  group_by(Description) %>%
  summarize(CladeCount = n_distinct(Clade)) %>%
  arrange(CladeCount, Description)

# Reordenar as descrições com base no número de clados
df$Description <- factor(df$Description, levels = description_counts$Description)

# Converte os dados para o formato longo usando reshape2
df_melt <- melt(df, id.vars = c("Clade", "Description", "Signature"), measure.vars = "qvalue")

# Crie o heatmap usando ggplot2
#ggplot(df_melt, aes(x = Description, y = Clade, fill = value)) +
# geom_tile() +
#  scale_fill_gradient(low = "pink", high = "red") +
#  labs(title = "Heatmap de Enriquecimento Funcional KEGG",
#       x = "Termo KEGG",
#       y = "Clado",
#       fill = "q-value") +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_melt <- df_melt %>%
  mutate(color = case_when(
    Signature == "cell adhesion" ~ "#06141FFF",
    Signature == "extracellular matrix organization" ~ "#742C14FF",
    Signature == "epithelial to mesenchymal transition" ~ "#3D4F7DFF",
    Signature == "mesenchymal to epithelial transition" ~ "#CD4F38FF",
    Signature == "regulation of metallopeptidase activity" ~ "#E48C2AFF",
    Signature == "cell junction organization" ~ "#72874EFF",
    Signature == "cellular extravasation" ~ "#046E8FFF",
    TRUE ~ "#FFFFFF" # Default color if no match
  ))

# Cria o heatmap com as cores baseadas na assinatura
ggplot(df_melt, aes(x = Description, y = Clade, fill = color, alpha = value)) +
  geom_tile() +
  scale_fill_identity() + # Usamos scale_fill_identity para usar as cores diretamente
  scale_alpha_continuous(range = c(1,0.01)) + # Ajuste o intervalo de alpha conforme necessário
  labs(title = "Heatmap de Enriquecimento Funcional KEGG",
       x = "Termo KEGG",
       y = "Clado",
       fill = "Signature",
       alpha = "q-value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
























library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)

# Função para clarear as cores
lighten <- function(color, factor=0.6) {
  col <- col2rgb(color)
  col <- col + (255 - col) * factor
  rgb(t(col), maxColorValue=255)
}

# Cores fornecidas
annotation_colors <- c(
  "cell adhesion"                               = "#06141FFF", 
  "extracellular matrix organization"           = "#742C14FF", 
  "epithelial to mesenchymal transition"        = "#3D4F7DFF",
  "mesenchymal to epithelial transition"        = "#CD4F38FF",
  "regulation of metallopeptidase activity"     = "#E48C2AFF",
  "cell junction organization"                  = "#72874EFF",
  "cellular extravasation"                      = "#046E8FFF"
)

# Clarear as cores fornecidas
light_colors <- sapply(annotation_colors, lighten)
color_scale <- c(light_colors, annotation_colors)

# Preparar os dados
combined_results_GO <- subset(combined_results, ID %in% metastasis_ids$GO_id)

# Ordene lca_names de acordo com a sequência desejada
ordered_clades <- lca_names$clade_name

# Transforme Clade em um fator com a ordem especificada
combined_results_GO$Clade <- factor(combined_results_GO$Clade, levels = ordered_clades)

# Seleciona as colunas necessárias
df <- combined_results_GO[, c("Clade", "Description", "qvalue", "ID")]

# Adicione as assinaturas da tabela metastasis_ids
df <- left_join(df, metastasis_ids, by = c("ID" = "GO_id"))
df <- dplyr::select(df, "Clade", "ID", Description = "Description.x", "qvalue", "Hallmark", "Signature")

# Contar o número de clados em que cada descrição aparece
description_counts <- df %>%
  group_by(Signature, Description) %>%
  summarize(CladeCount = n_distinct(Clade)) %>%
  arrange(Signature, CladeCount, Description)

# Reordenar as descrições com base na assinatura e clado
df$Description <- factor(df$Description, levels = description_counts$Description[description_counts$Signature %in% unique(df$Signature)])

# Converte os dados para o formato longo usando reshape2
df_melt <- melt(df, id.vars = c("Clade", "Description", "Signature"), measure.vars = "qvalue")

# Mapear as cores das assinaturas para os valores de qvalue
df_melt <- df_melt %>%
  mutate(color = case_when(
    Signature == "cell adhesion" ~ "#06141FFF",
    Signature == "extracellular matrix organization" ~ "#742C14FF",
    Signature == "epithelial to mesenchymal transition" ~ "#3D4F7DFF",
    Signature == "mesenchymal to epithelial transition" ~ "#CD4F38FF",
    Signature == "regulation of metallopeptidase activity" ~ "#E48C2AFF",
    Signature == "cell junction organization" ~ "#72874EFF",
    Signature == "cellular extravasation" ~ "#046E8FFF",
    TRUE ~ "#FFFFFF" # Default color if no match
  ))

# Crie o heatmap com as cores baseadas na assinatura e escala de alpha para qvalue
ggplot(df_melt, aes(x = Description, y = Clade, fill = color, alpha = value)) +
  geom_tile() +
  scale_fill_identity(name = "Signature") +  # Usamos scale_fill_identity para usar as cores diretamente
  scale_alpha_continuous(range = c(1, 0.1)) +  # Ajuste o intervalo de alpha conforme necessário
  labs(title = "Heatmap de Enriquecimento Funcional KEGG",
       x = "Assinatura",
       y = "Termo KEGG",
       fill = "Signature",
       alpha = "q-value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))
