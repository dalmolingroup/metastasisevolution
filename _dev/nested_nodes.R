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
# #
V(g3)$nodeColor <- V(g3)$node_color
V(g3)$nodeSize  <- V(g3)$node_size
# #
# V(g3)$cog_id    <- V(g3)$cog_id
# V(g3)$root      <- V(g3)$root
# V(g3)$clade     <- V(g3)$clade_name
# V(g3)$queryItem <- V(g3)$queryItem
# V(g3)$taxonId   <- V(g3)$ncbiTaxonId
# V(g3)$taxonName <- V(g3)$taxonName
# V(g3)$preferredName <- V(g3)$preferredName
# V(g3)$annotation    <- V(g3)$annotation
# #
# V(g3)$cell_extravasation               <- V(g3)$`cellular extravasation`
# V(g3)$cell_junction_organization       <- V(g3)$`cell junction organization`
# V(g3)$cell_adhesion                    <- V(g3)$`cell adhesion`
# V(g3)$extracellular_matrix_organization <- V(g3)$`extracellular matrix organization`
# V(g3)$epithelial_mesenchymal_transition <- V(g3)$`epithelial to mesenchymal transition`
# V(g3)$regulation_metallopeptidase      <- V(g3)$`regulation of metallopeptidase activity`



g2 <- subg(g = graph, dat = nodelist[nodelist$root %in% 30:37, ], refcol = 1, connected = F, maincomp = F, transdat = T)
g2  <- att.setv(g=g2, from="queryItem", to="nodeLabel")
# #
V(g2)$nodeColor <- V(g2)$node_color
V(g2)$nodeSize  <- V(g2)$node_size
# #
# V(g2)$cog_id    <- V(g2)$cog_id
# V(g2)$root      <- V(g2)$root
# V(g2)$clade     <- V(g2)$clade_name
# V(g2)$queryItem <- V(g2)$queryItem
# V(g2)$taxonId   <- V(g2)$ncbiTaxonId
# V(g2)$taxonName <- V(g2)$taxonName
# V(g2)$preferredName <- V(g2)$preferredName
# V(g2)$annotation    <- V(g2)$annotation
# #
# V(g2)$cell_extravasation               <- V(g2)$`cellular extravasation`
# V(g2)$cell_junction_organization       <- V(g2)$`cell junction organization`
# V(g2)$cell_adhesion                    <- V(g2)$`cell adhesion`
# V(g2)$extracellular_matrix_organization <- V(g2)$`extracellular matrix organization`
# V(g2)$epithelial_mesenchymal_transition <- V(g2)$`epithelial to mesenchymal transition`
# V(g2)$regulation_metallopeptidase      <- V(g2)$`regulation of metallopeptidase activity`

g1 <- subg(g = graph, dat = nodelist[nodelist$root %in% 37, ], refcol = 1, connected = F, maincomp = F, transdat = T)
g1  <- att.setv(g=g1, from="queryItem", to="nodeLabel")
#
V(g1)$nodeColor <- V(g1)$node_color
V(g1)$nodeSize  <- V(g1)$node_size
# #
# V(g1)$cog_id    <- V(g1)$cog_id
# V(g1)$root      <- V(g1)$root
# V(g1)$clade     <- V(g1)$clade_name
# V(g1)$queryItem <- V(g1)$queryItem
# V(g1)$taxonId   <- V(g1)$ncbiTaxonId
# V(g1)$taxonName <- V(g1)$taxonName
# V(g1)$preferredName <- V(g1)$preferredName
# V(g1)$annotation    <- V(g1)$annotation
# #
# V(g1)$cell_extravasation               <- V(g1)$`cellular extravasation`
# V(g1)$cell_junction_organization       <- V(g1)$`cell junction organization`
# V(g1)$cell_adhesion                    <- V(g1)$`cell adhesion`
# V(g1)$extracellular_matrix_organization <- V(g1)$`extracellular matrix organization`
# V(g1)$epithelial_mesenchymal_transition <- V(g1)$`epithelial to mesenchymal transition`
# V(g1)$regulation_metallopeptidase      <- V(g1)$`regulation of metallopeptidase activity`

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
 

 a <- getGraphFromRedeR(type = "node")
### ATENÇÃO: NECESSÁRIO ARRANJAR OS ATRIBUTOS
 # Ordenar os vértices dos dois grafos pelo nodeLabel
 # a <- a %>% induced_subgraph(order(V(a)$nodeLabel))
 # g3 <- g3 %>% induced_subgraph(order(V(g3)$nodeLabel))
 index_match <- match(V(a)$nodeLabel, V(g3)$nodeLabel)
 
 # Setting graph attributes
 V(a)$cog_id        <- V(g3)$cog_id[index_match]
 V(a)$root          <- V(g3)$root[index_match]
 V(a)$clade         <- V(g3)$clade_name[index_match]
 V(a)$queryItem     <- V(g3)$queryItem[index_match]
 V(a)$taxonId       <- V(g3)$ncbiTaxonId[index_match]
 V(a)$taxonName     <- V(g3)$taxonName[index_match]
 V(a)$preferredName <- V(g3)$preferredName[index_match]
 V(a)$annotation    <- V(g3)$annotation[index_match]
 #
 V(a)$`cellular extravasation`                  <- V(g3)$`cellular extravasation`[index_match]
 V(a)$`cell junction organization`              <- V(g3)$`cell junction organization`[index_match]
 V(a)$`cell adhesion`                           <- V(g3)$`cell adhesion`[index_match]
 V(a)$`extracellular matrix organization`       <- V(g3)$`extracellular matrix organization`[index_match]
 V(a)$`epithelial to mesenchymal transition`    <- V(g3)$`epithelial to mesenchymal transition`[index_match]
 V(a)$`regulation of metallopeptidase activity` <- V(g3)$`regulation of metallopeptidase activity`[index_match]
 
 ## Plots
ggraph(a, "manual", x = V(a)$x, y = V(a)$y) +
  geom_edge_link0(color = "#90909020") +  
  geom_node_point(aes(color = -root), size = 2) +  
  theme_void() + 
  theme(legend.position = "left") +
  guides(color = "none")

#########################################################################################################3
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

nodelist <- vroom::vroom("results/orthology_data/nodelist.csv")

BP_by_clade <- function(data) {
  data %>%
    pivot_longer(
      cols = c(`cellular extravasation`, `cell junction organization`, 
               `cell adhesion`, `extracellular matrix organization`, 
               `epithelial to mesenchymal transition`, 
               `regulation of metallopeptidase activity`),
      names_to = "process",
      values_to = "presence"
    ) %>%
    filter(presence == 1) %>%
    group_by(clade_name, process, root) %>% 
    summarise(count = n(), .groups = "drop") %>%
    mutate(total_genes = sum(count), .by = process) %>%
    mutate(proportion = count / total_genes)
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
  clados_interesse <- c("Metamonada", "Choanoflagellata", "Actinopterygii")
  
  filtered_data <- data %>%
    filter(clade_name %in% clados_interesse) %>%
    mutate(clade_name = factor(clade_name, levels = clados_interesse)) 
  
  # Criar o gráfico de barras empilhadas
  ggplot(filtered_data, aes(x = process, y = proportion, fill = clade_name)) +
    geom_bar(position = "fill", stat = "identity", color = NA, width = 0.7) +
    theme_minimal() +
    labs(
      title = "Percent BP by Clades",
      x="",
      y = "Genes Proportion",
      fill = "Biological Processes"
    ) +
    scale_fill_manual(values = c("Metamonada" = "#B3E6C6", 
                                 "Choanoflagellata" = "#FFCCCC", 
                                 "Actinopterygii" = "#99CCFF")) +
    scale_y_continuous(labels = scales::percent) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
}

# Exemplo de uso:
 resultado <- BP_by_clade(nodelist)
 plot_stacked_bar_chart(resultado)

 # Exemplo de uso:
# resultado <- contar_processos_por_clado(nodelist)
 plot_donut_chart(df_proportion, "Metamonada")
 plot_donut_chart(df_proportion, "Choanoflagellata")
 plot_donut_chart(df_proportion, "Actinopterygii")

 
 
 
 
 
 

 
 
 
 plot_stacked_bar_chart <- function(data) {
   # Definir os grupos de interesse
   clados_interesse <- c("Metamonada", "Choanoflagellata", "Actinopterygii")
   
   # Criar as categorias "Outros" com base nos intervalos de root
   filtered_data <- data %>%
     mutate(clade_group = case_when(
       clade_name == "Metamonada" ~ "Metamonada",
       root >= 31 & root <= 36 ~ "Intermediate M-C",
       clade_name == "Choanoflagellata" ~ "Choanoflagellata",
       root >= 21 & root <= 29 ~ "Intermediate C-A",
       clade_name == "Actinopterygii" ~ "Actinopterygii",
       TRUE ~ NA_character_
     )) %>%
     filter(!is.na(clade_group)) %>%
     mutate(clade_group = factor(clade_group, 
                                 levels = c("Metamonada", "Intermediate M-C", 
                                            "Choanoflagellata", "Intermediate C-A", 
                                            "Actinopterygii")))
   
   # Ajustar a largura das barras: mais finas para "Outros"
   filtered_data <- filtered_data %>%
     mutate(bar_width = ifelse(grepl("Intermediate", clade_group), 0.4, 0.8)) %>%
     mutate(alpha = ifelse(grepl("Intermediate", clade_group), 0.55, 1))
   
   # Criar o gráfico de barras empilhadas com larguras ajustadas
   ggplot(filtered_data, aes(x = clade_group, y = proportion, fill = process, 
                             width = bar_width)) +
     geom_bar(position = "fill", stat = "identity", color = NA, alpha = filtered_data$alpha) +
     theme_minimal() +
     labs(
       title = "Biological Process Distribution",
       x="",
       y = "Genes Proportion",
       fill = "Biological Processes"
     ) +
     scale_fill_manual(values = c(
       "cell adhesion"                               = "#06141FFF"
       ,"extracellular matrix organization"          = "#742C14FF"
       ,"epithelial to mesenchymal transition"       = "#3D4F7DFF"
       ,"regulation of metallopeptidase activity"    = "#E48C2AFF"
       ,"cell junction organization"                  ="#72874EFF"
       ,"cellular extravasation"                      ="#046E8FFF"
     )) + 
     scale_y_continuous(labels = scales::percent) +
     theme(
       axis.text = element_text(size = 12),
       axis.title = element_text(size = 14),
       axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
       legend.position = "right",
       plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
       panel.grid = element_blank()
     ) +
     guides(alpha = "none")
 }

 # Exemplo de uso:
  resultado <- BP_by_clade(nodelist)
  plot_stacked_bar_chart(resultado)
  
  
  plot_piedonut_chart <- function(data) {
    # Definir a ordem dos clados (do centro para fora)
    clados_ordem <- c("Metamonada", "Choanoflagellata", "Actinopterygii")
    
    # Definir paleta de cores fixa para os processos biológicos
    process_colors <- c(
      "cell adhesion"                               = "#06141FFF",
      "extracellular matrix organization"          = "#742C14FF",
      "epithelial to mesenchymal transition"       = "#3D4F7DFF",
      "regulation of metallopeptidase activity"    = "#E48C2AFF",
      "cell junction organization"                 = "#72874EFF",
      "cellular extravasation"                     = "#046E8FFF"
    )
    
    # Garantir que a ordem dos processos seja consistente
    data <- data %>%
      filter(clade_name %in% clados_ordem) %>%
      mutate(
        ring = case_when(
          clade_name == "Metamonada" ~ 1,
          clade_name == "Choanoflagellata" ~ 2,
          clade_name == "Actinopterygii" ~ 3
        ),
        process = factor(process, levels = names(process_colors))  # Ordem fixa dos processos
      )
    
    # Criar o gráfico piedonut
    ggplot(data, aes(x = ring, y = proportion, fill = process)) +
      geom_col(width = 1, color = "white") +         # Barras preenchendo o anel
      coord_polar(theta = "y") +                     # Transformação para gráfico circular
      xlim(0.5, 3.5) +                               # Espaço radial para os 3 anéis
      theme_void() +                                 # Remover fundos e eixos
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      labs(
        title = "Distribuição dos Processos Biológicos",
        fill = "Processo Biológico"
      ) +
      scale_fill_manual(values = process_colors)     # Paleta de cores fixa
  }
  
  # Exemplo de chamada da função
   plot_piedonut_chart(resultado)
  
   
   ######################################################################################
   library(phyper)
   nodelist <- vroom::vroom("results/orthology_data/nodelist.csv") %>% 
     filter(., root >= 20)
   
   
   # N -> Número de genes total = Universo
   # M -> Número de genes do set do Gene Ontoloy
   # n -> Número de genes da sua lista de interesse
   # k -> Número de genes na interseção entre M e n
   N <- length(nodelist$queryItem)
   M <- sum(nodelist$`cell adhesion`)
   n <- length(nodelist$clade_name == "Choanoflagellata")
   k <- sum(filter(nodelist, clade_name == "Choanoflagellata")$`cell adhesion`)
   
   a <- 1 - phyper(k - 1, M, N - M, n)
   