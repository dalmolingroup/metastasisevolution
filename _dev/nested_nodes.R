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


g2 <- subg(g = graph, dat = nodelist[nodelist$root %in% 30:37, ], refcol = 1)

g2  <- att.setv(g=g2, from="queryItem", to="nodeLabel")


g3 <- subg(g = graph, dat = nodelist[nodelist$root %in% 20:37, ], refcol = 1)

g2  <- att.setv(g=g2, from="queryItem", to="nodeLabel")

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

g1 <- subg(g=graph, dat=nodelist[nodelist$root %in% 32], refcol=1)
g1 <- att.setv(g=g1, from="queryItem" to="nodeLabel")
V(g1)$nodeColor <- ifelse(clade[V(g1), ] == 'Metamonada', "red", "gray")
#V(g1)$nodeColor <- ifelse(V(graph)$name %in% filtered, "black", "gray")