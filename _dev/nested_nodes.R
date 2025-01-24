load("results/plots/graph")
nodelist <- vroom::vroom("results/orthology_data/nodelist.csv")



library(RedeR)
library(igraph)

#rdp <- RedPort()
#calld(rdp)
startRedeR()

resetRedeR()

g1 <- subg(g = graph, dat = nodelist[nodelist$root %in% 37, ], refcol = 1)

g1  <- att.setv(g=g1, from="queryItem", to="nodeLabel")
V(g1)$nodeColor <- ifelse(V(g1)$clade_name  %in% c("Metamonada","Choanoflagellata", "Actinopterygii"), "black", "gray")


g2 <- subg(g = graph, dat = nodelist[nodelist$root %in% 30:37, ], refcol = 1)

g2  <- att.setv(g=g2, from="queryItem", to="nodeLabel")
V(g2)$nodeColor <- ifelse(V(g2)$clade_name  %in% c("Metamonada","Choanoflagellata", "Actinopterygii"), "black", "gray")

g3 <- subg(g = graph, dat = nodelist[nodelist$root %in% 20:37, ], refcol = 1)

g3  <- att.setv(g=g3, from="queryItem", to="nodeLabel")
V(g3)$nodeColor <- ifelse(V(g3)$clade_name  %in% c("Metamonada","Choanoflagellata", "Actinopterygii"), "black", "gray")

g1$nestAlias <- "LCA"
g2$nestAlias <- "LCA"
g3$nestAlias <- "LCA"

N1 <- addGraphToRedeR( g1, gcoord=c(10,25), gscale=20, isNested=TRUE, theme='tm1', zoom=30)
N2 <- addGraphToRedeR( g2, gcoord=c(20,70), gscale=50, isNested=TRUE, theme='tm1', zoom=30)
N3 <- addGraphToRedeR( g3, gcoord=c(70,55), gscale=80, isNested=TRUE, theme='tm1', zoom=30)

N4 <- nestNodes( nodes=V(g1)$name, parent=N2, theme='tm2')
N5 <- nestNodes( nodes=V(g2)$name, parent=N3, theme='tm2')
nestNodes( nodes=V(g1)$name, parent=N5, theme='tm3')

mergeOutEdges( nlevels=2)


relax(rdp, p1=100, p2=100, p3=5, p4=150, p5=5, p8=10, p9=20)
selectNodes(rdp,"RET")



library("RedeR")
library("igraph")

resetRedeR()

#############################################################3
g1 <- subg(g = graph, dat = nodelist[nodelist$root %in% 37, ], refcol = 1)
g1  <- att.setv(g=g1, from="queryItem", to="nodeLabel")

g2 <- subg(g = graph, dat = nodelist[nodelist$root %in% 31:37, ], refcol = 1)
g2  <- att.setv(g=g2, from="queryItem", to="nodeLabel")

g3 <- subg(g = graph, dat = nodelist[nodelist$root %in% 30:37, ], refcol = 1)
g3  <- att.setv(g=g3, from="queryItem", to="nodeLabel")

g4 <- subg(g = graph, dat = nodelist[nodelist$root %in% 21:37, ], refcol = 1)
g4  <- att.setv(g=g4, from="queryItem", to="nodeLabel")

g5 <- subg(g = graph, dat = nodelist[nodelist$root %in% 20:37, ], refcol = 1)
g5  <- att.setv(g=g5, from="queryItem", to="nodeLabel")

N1 <- addGraphToRedeR( g1, gcoord=c(10,25), gscale=20, isNested=TRUE, theme='tm1', zoom=30)
N2 <- addGraphToRedeR( g2, gcoord=c(20,70), gscale=50, isNested=TRUE, theme='tm1', zoom=30)
N3 <- addGraphToRedeR( g3, gcoord=c(70,55), gscale=80, isNested=TRUE, theme='tm1', zoom=30)
N4 <- addGraphToRedeR( g4, gcoord=c(70,55), gscale=80, isNested=TRUE, theme='tm1', zoom=30)
N5 <- addGraphToRedeR( g5, gcoord=c(70,55), gscale=80, isNested=TRUE, theme='tm1', zoom=30)

N4 <- nestNodes( nodes=V(g1)$name, parent=N2, theme='tm2')
N5 <- nestNodes( nodes=V(g2)$name, parent=N3, theme='tm2')
N6 <- nestNodes( nodes=V(g3)$name, parent=N4, theme='tm2')


N9 <- nestNodes( nodes=V(g1)$name, parent=N5, theme='tm1')
N10 <- nestNodes( nodes=V(g2)$name, parent=N9, theme='tm1')
N11 <- nestNodes( nodes=V(g3)$name, parent=N10 , theme='tm1')
nestNodes( nodes=V(g4)$name, parent=N11 , theme='tm1')
