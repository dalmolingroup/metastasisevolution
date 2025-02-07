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
