library(dplyr)
library(RedeR)
library(igraph)

setwd("/home/gleison/Documents/metastasisevolution/")

load("results/plots/graph")
nodelist <- vroom::vroom("results/orthology_data/nodelist.csv") %>%
  mutate(node_color = case_when(
    clade_name == "Metamonada" ~ "#B3E6C6FF",
    clade_name == "Choanoflagellata" ~ "#FFCCCCFF",
    clade_name == "Actinopterygii" ~ "#99CCFFFF",
    TRUE ~ "#999999FF"
  )) %>%
  mutate(node_size = case_when(
    clade_name %in% c("Metamonada","Choanoflagellata","Actinopterygii") ~ 50,
    TRUE ~ 35
  ))

startRedeR()

resetRedeR()

g3 <- subg(g = graph, dat = nodelist[nodelist$root %in% 20:37, ], refcol = 1, connected = F, maincomp = F, transdat = T)
g3  <- att.setv(g=g3, from="queryItem", to="nodeLabel")
#g3 <- att.setv(g=g3, from = "node_color", to="nodeColor")
#V(g3)$nodeColor <- ifelse(V(g3)$clade_name  %in% c("Metamonada","Choanoflagellata", "Actinopterygii"), "black", "gray")
V(g3)$nodeColor <- V(g3)$node_color
V(g3)$nodeSize <- V(g3)$node_size

g2 <- subg(g = graph, dat = nodelist[nodelist$root %in% 30:37, ], refcol = 1, connected = F, maincomp = F, transdat = T)
g2  <- att.setv(g=g2, from="queryItem", to="nodeLabel")
#g2 <- att.setv(g=g2, from = "node_color", to="nodeColor")
#V(g2)$nodeColor <- ifelse(V(g2)$clade_name  %in% c("Metamonada","Choanoflagellata", "Actinopterygii"), "black", "gray")
V(g2)$nodeColor <- V(g2)$node_color
V(g2)$nodeSize <- V(g2)$node_size

g1 <- subg(g = graph, dat = nodelist[nodelist$root %in% 37, ], refcol = 1, connected = F, maincomp = F, transdat = T)
g1  <- att.setv(g=g1, from="queryItem", to="nodeLabel")
#g1 <- att.setv(g=g1, from = "node_color", to="nodeColor")
#V(g1)$nodeColor <- ifelse(V(g1)$clade_name  %in% c("Metamonada","Choanoflagellata", "Actinopterygii"), "black", "gray")
V(g1)$nodeColor <- V(g1)$node_color
V(g1)$nodeSize <- V(g1)$node_size

#g1$nestAlias <- "Human-Actinopterygii"
#g2$nestAlias <- "LCA"
# g3$nestAlias <- "LCA"

#N1 <- addGraphToRedeR( g1, gcoord=c(10,25), gscale=20, isNested=TRUE, theme='tm1', zoom=30)
#N2 <- addGraphToRedeR( g2, gcoord=c(20,70), gscale=50, isNested=TRUE, theme='tm1', zoom=30)

N3 <- addGraphToRedeR( g3, gcoord=c(0,0), gscale=80, isNested=FALSE, theme='tm1', zoom=30)
N3 <- nestNodes(nodes = V(g3)$name, parent = N3, theme = 'tm1', gatt = list(nestLabel = "Human-Actinopterygii",
                                                                            nestShape = "Circle",
                                                                            nestSize = 3000,
                                                                            nestLabelSize = 80,
                                                                            nestLabelColor = "#99CCFFFF",
                                                                            nestLabelCoords = c(x=-10, y=-2000),
                                                                            nestLineType = "SOLID",
                                                                            nestLineWidth = 15,
                                                                            nestLineColor = "#99CCFFFF"))
  

#N4 <- nestNodes( nodes=V(g1)$name, parent=N2, theme='tm1', status = "transparent")
N5 <- nestNodes(nodes=V(g2)$name, gcoord = c(40, 57), parent=N3, theme='tm1', gatt = list(nestLabel = "Human-Choanoflagellata",
                                                                       nestShape = "Circle",
                                                                       nestSize = 2300,
                                                                       nestLabelSize = 80,
                                                                       nestLabelColor = "#FFCCCCFF",
                                                                       nestLabelCoords = c(x=-9, y=-1780),
                                                                       nestLineType = "SOLID",
                                                                       nestLineWidth = 15,
                                                                       nestLineColor = "#FFCCCCFF"))


nestNodes( nodes=V(g1)$name, gcoord = c(41, 62), parent=N5, theme='tm1', gatt = list(nestLabel = "Human-Metamonada",
                                                               nestShape = "Circle",
                                                               nestSize = 1600,
                                                               nestLabelSize = 80,
                                                               nestLabelColor = "#B3E6C6FF",
                                                               nestLabelCoords = c(x=-7, y=-1500),
                                                               nestLineType = "SOLID",
                                                               nestLineWidth = 15,
                                                               nestLineColor = "#B3E6C6FF"))

mergeOutEdges( nlevels=2)

# RedeR force-directed parameters
 relaxRedeR(
p1 = 200,
p2 = 500,
p3 = 200,
p4 = 350,
p5 = 50,
p6 = 75,
p7 = 100,
p8 = 300,
p9 = 1500
)
 
a <- getGraphFromRedeR()
b <- graph_from_edgelist(a, vertices = nodelist, directed = F)



#########################################################################################################3

library(dplyr)
library(tidyr)

nodelist <- vroom::vroom("results/orthology_data/nodelist.csv")

BP_by_clade <- function(nodelist) {
  # Selecionar colunas de processos biológicos
  process_cols <- c("cellular extravasation", "cell junction organization", 
                    "cell adhesion", "extracellular matrix organization", 
                    "epithelial to mesenchymal transition", 
                    "regulation of metallopeptidase activity")
  
  # Transformar em formato longo
  long_df <- nodelist %>%
    pivot_longer(cols = all_of(process_cols), 
                 names_to = "process", 
                 values_to = "participation") %>%
    filter(participation == 1)  # Considerar apenas genes que participam do processo
  
  # Contagem de genes por processo e clado
  summary_df <- long_df %>%
    group_by(clade_name, process) %>%
    summarise(genes_in_process = n(), .groups = 'drop')
  
  # Calcular o total de genes por processo (todos os clados)
  total_genes_per_process <- long_df %>%
    group_by(process) %>%
    summarise(total_genes_in_process = n(), .groups = 'drop')
  
  # Juntar as informações e calcular a proporção
  final_df <- summary_df %>%
    left_join(total_genes_per_process, by = "process") %>%
    mutate(proportion = genes_in_process / total_genes_in_process)
  
  return(final_df)
}

df_proportion <- BP_by_clade(nodelist)


plot_donut_chart <- function(data, clado_selecionado) {
    # Filtrar o clado de interesse
    clado_data <- data %>%
    filter(clade_name == clado_selecionado)
   
    # Criar o gráfico do tipo donut
    ggplot(clado_data, aes(x = 2, y = proportion, fill = process)) +
       geom_col(width = 0.1, color = "white") +
       coord_polar(theta = "y") +
       xlim(0.5, 2.5) +
       theme_void() +
       theme(legend.position = "right") +
       labs(title = paste("Distribuição dos Processos Biológicos em", clado_selecionado),
            fill = "Processo Biológico") +
       scale_fill_brewer(palette = "Set3")
}

plot_stacked_bar_chart <- function(data) {
  # Filtrar os clados de interesse
  clados_interesse <- c("Metamonada", "Actinopterygii", "Choanoflagellata")
  
  filtered_data <- data %>%
    filter(clade_name %in% clados_interesse)
  
  # Criar o gráfico de barras empilhadas
  ggplot(filtered_data, aes(x = clade_name, y = proportion, fill = process)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    theme_minimal() +
    labs(
      title = "Distribuição dos Processos Biológicos por Clado",
      x = "Clado",
      y = "Proporção de Genes",
      fill = "Processo Biológico"
    ) +
    scale_fill_brewer(palette = "Paired") +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
}

# Exemplo de uso:
# resultado <- contar_processos_por_clado(nodelist)
# plot_stacked_bar_chart(resultado)

 # Exemplo de uso:
# resultado <- contar_processos_por_clado(nodelist)
 plot_donut_chart(df_proportion, "Metamonada")
 plot_donut_chart(df_proportion, "Choanoflagellata")
 plot_donut_chart(df_proportion, "Actinopterygii")
