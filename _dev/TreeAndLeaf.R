color_palette <- c(
   "cell adhesion"                              = "#06141F"
  ,"extracellular matrix organization"          = "#742C14"
  ,"epithelial to mesenchymal transition"       = "#3D4F7D"
  ,"regulation of metallopeptidase activity"    = "#E48C2A"
  ,"cell junction organization"                 = "#72874E"
  ,"cellular extravasation"                     = "#046E8F"
)


metastasis_ids <- metastasis_ids %>%
  mutate(color = color_palette[Signature])

mSim <- mgoSim(metastasis_ids$GO_id, metastasis_ids$GO_id, semData = semData, measure = "Wang", combine = NULL)

# Setting nodes sizes
size <- 100
aux <- sort(unique(size))
names(aux) <- as.character(1:length(aux))
sizeRanks <- as.factor(names(aux[match(size, aux)]))
sizeIntervals <- 1
sizeMultiplier <- 7
sizeBase <- 100
metastasis_ids$size <- (sizeBase + (as.numeric(sizeRanks) * sizeMultiplier))

# Terms Clustering
hc <- hclust(dist(mSim), "average")
tal <- treeAndLeaf(hc)

tal <- att.mapv(g = tal, dat = metastasis_ids, refcol = 3)
pal <- brewer.pal(9, "OrRd")
tal <- att.setv(g = tal, from = "Description", to = "nodeLabel")
tal <- att.setv(g = tal, from = "Signature", to = "nodeColor")
tal <- att.setv(
  g = tal,
  from = "color",
  to = "nodeColor",
  pal = 1,
  cols = c(
     "#046E8F"
    ,"#06141F"
    ,"#3D4F7D"
    ,"#72874E"
    ,"#742C14"
    ,"#E48C2A"
  )#,
#  nquant = 5
)
tal <- att.setv(
  g = tal,
  from = "size",
  to = "nodeSize",
  xlim = c(
    (sizeBase + sizeMultiplier),
    (sizeBase + (sizeMultiplier * sizeIntervals)),
    sizeMultiplier
  )#,
 # nquant = sizeIntervals
)
tal <- att.addv(tal, "nodeFontSize", value = 15, index = V(tal)$isLeaf)
tal <- att.adde(tal, "edgeWidth", value = 3)
tal <- att.adde(tal, "edgeColor", value = "gray80")
tal <- att.addv(tal, "nodeLineColor", value = "gray80")

V(tal)$nodeColor

resetRedeR()
addGraphToRedeR(tal)

startRedeR()

"#06141F"
