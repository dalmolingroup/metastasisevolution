# =============================================================================
# Metastasis Suppressor Genes (MSG) GO Enrichment Analysis with TreeAndLeaf
# =============================================================================
# Author: Gleison M. Azevedo
# Description: Performs GO enrichment of metastasis suppressor genes across
#              3 evolutionary groups, clusters by semantic similarity, and
#              visualizes using TreeAndLeaf/RedeR
# =============================================================================

# -----------------------------------------------------------------------------
# Load Libraries
# -----------------------------------------------------------------------------
library(TreeAndLeaf)
library(RedeR)
library(igraph)
library(tidyverse)
library(here)
library(vroom)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(classInt)
library(UpSetR)

# -----------------------------------------------------------------------------
# Configuration: Define Evolutionary Group Colors
# -----------------------------------------------------------------------------
group_colors <- c(
  "Unicellular" = "#E69F00", # Orange
  "Multicellular non-vertebrate" = "#0072B2", # Blue
  "Vertebrates" = "#009E73" # Green
)

# -----------------------------------------------------------------------------
# Read Data
# -----------------------------------------------------------------------------
# Nodelist with rooting information and biological process annotations
nodelist <- vroom::vroom(here("results/orthology_data/nodelist.csv"))

# Metastasis suppressor genes reference list
msg_reference <- vroom::vroom(
  here("assets/metastasis_supressors_nature.csv"),
  delim = "\t"
)

# -----------------------------------------------------------------------------
# Filter Metastasis Suppressor Genes
# -----------------------------------------------------------------------------
# Filter nodelist to keep only metastasis suppressor genes
msg_nodelist <- nodelist %>%
  filter(`metastasis suppressor` == 1) %>%
  dplyr::select(cog_id, root, clade_name, string_id, queryItem)

cat("Total metastasis suppressor genes in nodelist:", nrow(msg_nodelist), "\n")

# -----------------------------------------------------------------------------
# Assign Evolutionary Groups Based on Root Values
# -----------------------------------------------------------------------------
# Group 1: Unicellular - before Choanoflagellata (exclusive), root > 30
# Group 2: Multicellular non-vertebrate - Choanoflagellata to Tunicata, 21 <= root <= 30
# Group 3: Vertebrates - from Actinopterygii (inclusive), root <= 20

msg_nodelist <- msg_nodelist %>%
  mutate(
    evolutionary_group = case_when(
      root > 30 ~ "Unicellular",
      root >= 21 & root <= 30 ~ "Multicellular non-vertebrate",
      root <= 20 ~ "Vertebrates"
    )
  ) %>%
  mutate(
    group_color = group_colors[evolutionary_group]
  )

# Summary of genes per group
cat("\n--- Gene Distribution by Evolutionary Group ---\n")
print(table(msg_nodelist$evolutionary_group))

# -----------------------------------------------------------------------------
# GO Enrichment Analysis per Evolutionary Group
# -----------------------------------------------------------------------------
cat("\n--- Performing GO Enrichment Analysis ---\n")

# Function to perform GO enrichment for a group
perform_go_enrichment <- function(genes, group_name) {
  cat(sprintf("  Processing %s (%d genes)...\n", group_name, length(genes)))

  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (nrow(gene_ids) == 0) {
    cat(sprintf("    No genes mapped for %s\n", group_name))
    return(NULL)
  }

  # Perform GO enrichment
  ego <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat(sprintf("    No enriched terms for %s\n", group_name))
    return(NULL)
  }

  result <- as.data.frame(ego) %>%
    mutate(evolutionary_group = group_name)

  cat(sprintf("    Found %d enriched GO terms\n", nrow(result)))
  return(result)
}

# Get unique genes per group
genes_unicellular <- msg_nodelist %>%
  filter(evolutionary_group == "Unicellular") %>%
  pull(queryItem) %>%
  unique()

genes_multicellular <- msg_nodelist %>%
  filter(evolutionary_group == "Multicellular non-vertebrate") %>%
  pull(queryItem) %>%
  unique()

genes_vertebrates <- msg_nodelist %>%
  filter(evolutionary_group == "Vertebrates") %>%
  pull(queryItem) %>%
  unique()

# Perform enrichment for each group
go_unicellular <- perform_go_enrichment(genes_unicellular, "Unicellular")
go_multicellular <- perform_go_enrichment(
  genes_multicellular,
  "Multicellular non-vertebrate"
)
go_vertebrates <- perform_go_enrichment(genes_vertebrates, "Vertebrates")

# Combine all results
go_combined <- bind_rows(
  go_unicellular,
  go_multicellular,
  go_vertebrates
) %>%
  filter(!is.na(ID))

if (nrow(go_combined) == 0) {
  stop("No enriched GO terms found in any group.")
}

cat(sprintf(
  "\nTotal enriched GO terms across all groups: %d\n",
  nrow(go_combined)
))

# Save enrichment results
if (!is.null(go_unicellular) && nrow(go_unicellular) > 0) {
  vroom::vroom_write(
    go_unicellular,
    here("results/paper_data/MSG_GO_enrichment_unicellular.csv"),
    delim = ","
  )
}

if (!is.null(go_multicellular) && nrow(go_multicellular) > 0) {
  vroom::vroom_write(
    go_multicellular,
    here("results/paper_data/MSG_GO_enrichment_multicellular.csv"),
    delim = ","
  )
}

if (!is.null(go_vertebrates) && nrow(go_vertebrates) > 0) {
  vroom::vroom_write(
    go_vertebrates,
    here("results/paper_data/MSG_GO_enrichment_vertebrates.csv"),
    delim = ","
  )
}

vroom::vroom_write(
  go_combined,
  here("results/paper_data/MSG_GO_enrichment_combined.csv"),
  delim = ","
)

cat("\nEnrichment results saved to results/paper_data/\n")

# -----------------------------------------------------------------------------
# UpSet Plot - GO Term Overlap Between Groups
# -----------------------------------------------------------------------------
cat("\n--- Creating UpSet Plot for GO Term Overlap ---\n")

# Create list of GO terms per group
go_terms_list <- list(
  "Unicellular" = if (!is.null(go_unicellular) && nrow(go_unicellular) > 0) {
    go_unicellular$ID
  } else {
    character(0)
  },
  "Multicellular" = if (
    !is.null(go_multicellular) && nrow(go_multicellular) > 0
  ) {
    go_multicellular$ID
  } else {
    character(0)
  },
  "Vertebrates" = if (!is.null(go_vertebrates) && nrow(go_vertebrates) > 0) {
    go_vertebrates$ID
  } else {
    character(0)
  }
)

# Remove empty groups
go_terms_list <- go_terms_list[sapply(go_terms_list, length) > 0]

if (length(go_terms_list) >= 2) {
  # Create UpSet plot
  upset_plot <- upset(
    fromList(go_terms_list),
    order.by = "freq",
    decreasing = TRUE,
    mb.ratio = c(0.6, 0.4),
    number.angles = 0,
    text.scale = 1.5,
    point.size = 3.5,
    line.size = 1.5,
    mainbar.y.label = "GO Terms Intersection",
    sets.x.label = "GO Terms per Group",
    sets.bar.color = c(
      group_colors["Unicellular"],
      group_colors["Multicellular non-vertebrate"],
      group_colors["Vertebrates"]
    )[names(go_terms_list)]
  )

  # Save UpSet plot
  pdf(here("results/plots/MSG_GO_upset_plot.pdf"), width = 10, height = 7)
  print(upset_plot)
  dev.off()

  png(
    here("results/plots/MSG_GO_upset_plot.png"),
    width = 1000,
    height = 700,
    res = 100
  )
  print(upset_plot)
  dev.off()

  cat("UpSet plot saved to results/plots/MSG_GO_upset_plot.pdf and .png\n")

  # Print overlap statistics
  cat("\n--- GO Term Overlap Statistics ---\n")

  # Find terms in multiple groups
  all_terms <- unique(go_combined$ID)
  terms_in_groups <- sapply(all_terms, function(term) {
    groups <- go_combined %>%
      filter(ID == term) %>%
      pull(evolutionary_group) %>%
      unique()
    paste(groups, collapse = " & ")
  })

  overlap_summary <- table(terms_in_groups)
  cat("Distribution of GO terms across groups:\n")
  print(overlap_summary)

  # List terms that appear in multiple groups
  shared_terms <- go_combined %>%
    group_by(ID, Description) %>%
    summarise(
      groups = paste(unique(evolutionary_group), collapse = " & "),
      n_groups = n_distinct(evolutionary_group),
      .groups = "drop"
    ) %>%
    filter(n_groups > 1) %>%
    arrange(desc(n_groups))

  if (nrow(shared_terms) > 0) {
    cat(sprintf(
      "\nGO terms appearing in multiple groups: %d\n",
      nrow(shared_terms)
    ))
    cat("\nShared GO terms:\n")
    print(as.data.frame(shared_terms), row.names = FALSE)

    # Save shared terms
    vroom::vroom_write(
      shared_terms,
      here("results/paper_data/MSG_GO_shared_terms.csv"),
      delim = ","
    )
    cat("\nShared terms saved to: results/paper_data/MSG_GO_shared_terms.csv\n")
  } else {
    cat("\nNo GO terms appear in multiple groups.\n")
  }
} else {
  cat("Not enough groups with enriched terms for UpSet plot.\n")
}

# -----------------------------------------------------------------------------
# Semantic Similarity Clustering
# -----------------------------------------------------------------------------
cat("\n--- Calculating GO Semantic Similarity ---\n")

# Get unique GO term IDs
go_terms <- unique(go_combined$ID)
cat(sprintf(
  "Calculating similarity for %d unique GO terms...\n",
  length(go_terms)
))

# Create semantic data object
semData <- godata(OrgDb = org.Hs.eg.db, ont = "BP")

# Calculate pairwise semantic similarity
mSim <- mgoSim(
  go_terms,
  go_terms,
  semData = semData,
  measure = "Wang",
  combine = NULL
)

# Handle NA values by replacing with 0
mSim[is.na(mSim)] <- 0

# Create distance matrix (convert similarity to distance)
dist_matrix <- as.dist(1 - mSim)

# Hierarchical clustering
hc <- hclust(dist_matrix, method = "average")

cat("Semantic similarity clustering complete.\n")

# -----------------------------------------------------------------------------
# Prepare Data for Visualization
# -----------------------------------------------------------------------------
# Calculate ratio of genes in GO term
go_combined <- go_combined %>%
  mutate(
    path_length = as.integer(sapply(strsplit(GeneRatio, "/"), "[", 2)),
    ratio = Count / path_length
  )

# For each GO term, determine the primary group (if term appears in multiple groups)
go_term_info <- go_combined %>%
  group_by(ID) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(
    ID,
    Description,
    Count,
    ratio,
    pvalue,
    p.adjust,
    evolutionary_group
  )

# Organize node sizes
size <- go_term_info$Count
aux <- sort(unique(size))
names(aux) <- as.character(1:length(aux))
sizeRanks <- as.factor(names(aux[match(size, aux)]))
sizeIntervals <- min(7, length(unique(size)))
sizeMultiplier <- 15
sizeBase <- 50
go_term_info$size <- (sizeBase + (as.numeric(sizeRanks) * sizeMultiplier))

# Add group color
go_term_info$group_color <- group_colors[go_term_info$evolutionary_group]

# -----------------------------------------------------------------------------
# Create TreeAndLeaf Object
# -----------------------------------------------------------------------------
cat("\n--- Creating TreeAndLeaf Visualization ---\n")

tal <- treeAndLeaf(hc)

# Map attributes to the graph
tal <- att.mapv(g = tal, dat = as.data.frame(go_term_info), refcol = 1)

# Color palette for ratio
pal <- brewer.pal(9, "OrRd")

# Set node labels (GO term descriptions)
tal <- att.setv(g = tal, from = "Description", to = "nodeAlias")

# Set node colors based on ratio
tal <- att.setv(
  g = tal,
  from = "ratio",
  to = "nodeColor",
  cols = pal,
  nquant = 5
)

# Set node sizes based on gene count
tal <- att.setv(
  g = tal,
  from = "size",
  to = "nodeSize",
  xlim = c(
    (sizeBase + sizeMultiplier),
    (sizeBase + (sizeMultiplier * sizeIntervals)),
    sizeMultiplier
  ),
  nquant = sizeIntervals
)

# Set font size for labels
tal <- att.addv(tal, "nodeFontSize", value = 15, index = V(tal)$isLeaf)

# Set edge attributes
tal <- att.adde(tal, "edgeWidth", value = 3)
tal <- att.adde(tal, "edgeColor", value = "gray80")

# -----------------------------------------------------------------------------
# Set nodeLineColor Based on Evolutionary Group
# -----------------------------------------------------------------------------
# Get leaf nodes
leaf_indices <- which(V(tal)$isLeaf)

# Set nodeLineColor for each leaf based on its group
for (i in leaf_indices) {
  node_name <- V(tal)$name[i]

  # Find the group for this GO term
  node_group <- go_term_info %>%
    filter(ID == node_name) %>%
    pull(evolutionary_group)

  if (length(node_group) > 0 && !is.na(node_group)) {
    V(tal)$nodeLineColor[i] <- group_colors[node_group]
  } else {
    V(tal)$nodeLineColor[i] <- "gray50"
  }
}

# Set line width for the border
tal <- att.addv(tal, "nodeLineWidth", value = 4, index = V(tal)$isLeaf)

# For non-leaf nodes, set gray color
non_leaf_indices <- which(!V(tal)$isLeaf)
V(tal)$nodeLineColor[non_leaf_indices] <- "gray80"
V(tal)$nodeLineWidth[non_leaf_indices] <- 1

# -----------------------------------------------------------------------------
# Display with RedeR
# -----------------------------------------------------------------------------
cat("\n--- Starting RedeR Visualization ---\n")
cat("Group Colors (nodeLineColor):\n")
for (group in names(group_colors)) {
  cat(sprintf("  %s: %s\n", group, group_colors[group]))
}

# Start RedeR
rdp <- RedPort()
calld(rdp)

# Reset and add graph
resetd(rdp)
addGraph(rdp, tal, layout = NULL)

# Add legends
addLegend.color(
  obj = rdp,
  tal,
  title = "Gene Ratio",
  position = "topright"
)

addLegend.size(
  obj = rdp,
  tal,
  title = "Gene Count",
  position = "bottomright"
)

# Relax the network layout
relax(rdp, p1 = 21, p2 = 70, p3 = 255, p4 = 100, p5 = 50)

cat("\n--- Visualization Complete ---\n")
cat("The TreeAndLeaf visualization is displayed in RedeR.\n")
cat("\nLegend:\n")
cat("  - Node color (fill): Gene ratio (proportion of genes in GO term)\n")
cat("  - Node size: Number of genes\n")
cat("  - Node border (nodeLineColor): Evolutionary group\n")
cat("    - Orange (#E69F00): Unicellular\n")
cat("    - Blue (#0072B2): Multicellular non-vertebrate\n")
cat("    - Green (#009E73): Vertebrates\n")

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n--- Analysis Summary ---\n")
cat(sprintf("Total MSG genes analyzed: %d\n", nrow(msg_nodelist)))
cat(sprintf(
  "  Unicellular: %d genes\n",
  sum(msg_nodelist$evolutionary_group == "Unicellular")
))
cat(sprintf(
  "  Multicellular non-vertebrate: %d genes\n",
  sum(msg_nodelist$evolutionary_group == "Multicellular non-vertebrate")
))
cat(sprintf(
  "  Vertebrates: %d genes\n",
  sum(msg_nodelist$evolutionary_group == "Vertebrates")
))

cat(sprintf("\nTotal enriched GO terms: %d\n", nrow(go_combined)))
cat("GO terms per group:\n")
print(table(go_term_info$evolutionary_group))

# Export MSG nodelist with group assignment
msg_nodelist %>%
  vroom::vroom_write(
    file = here("results/orthology_data/MSG_nodelist_grouped.csv"),
    delim = ","
  )

cat(
  "\nMSG nodelist saved to: results/orthology_data/MSG_nodelist_grouped.csv\n"
)
cat("GO enrichment results saved to: results/paper_data/\n")

# =============================================================================
# KEGG Pathway Enrichment Analysis with TreeAndLeaf
# =============================================================================
cat("\n")
cat("===========================================================================\n")
cat("KEGG Pathway Enrichment Analysis\n")
cat("===========================================================================\n")

# -----------------------------------------------------------------------------
# KEGG Enrichment per Evolutionary Group
# -----------------------------------------------------------------------------
cat("\n--- Performing KEGG Enrichment Analysis ---\n")

# Function to perform KEGG enrichment for a group
perform_kegg_enrichment <- function(genes, group_name) {
  cat(sprintf("  Processing %s (%d genes)...\n", group_name, length(genes)))

  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (nrow(gene_ids) == 0) {
    cat(sprintf("    No genes mapped for %s\n", group_name))
    return(NULL)
  }

  # Perform KEGG enrichment
  ekegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    cat(sprintf("    No enriched KEGG pathways for %s\n", group_name))
    return(NULL)
  }

  # Convert to readable gene symbols
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

  result <- as.data.frame(ekegg) %>%
    mutate(evolutionary_group = group_name)

  cat(sprintf("    Found %d enriched KEGG pathways\n", nrow(result)))
  return(result)
}

# Perform KEGG enrichment for each group
kegg_unicellular <- perform_kegg_enrichment(genes_unicellular, "Unicellular")
kegg_multicellular <- perform_kegg_enrichment(
  genes_multicellular,
  "Multicellular non-vertebrate"
)
kegg_vertebrates <- perform_kegg_enrichment(genes_vertebrates, "Vertebrates")

# Combine all KEGG results
kegg_combined <- bind_rows(
  kegg_unicellular,
  kegg_multicellular,
  kegg_vertebrates
) %>%
  filter(!is.na(ID))

if (nrow(kegg_combined) > 0) {
  cat(sprintf(
    "\nTotal enriched KEGG pathways across all groups: %d\n",
    nrow(kegg_combined)
  ))

  # Save KEGG enrichment results
  if (!is.null(kegg_unicellular) && nrow(kegg_unicellular) > 0) {
    vroom::vroom_write(
      kegg_unicellular,
      here("results/paper_data/MSG_KEGG_enrichment_unicellular.csv"),
      delim = ","
    )
  }

  if (!is.null(kegg_multicellular) && nrow(kegg_multicellular) > 0) {
    vroom::vroom_write(
      kegg_multicellular,
      here("results/paper_data/MSG_KEGG_enrichment_multicellular.csv"),
      delim = ","
    )
  }

  if (!is.null(kegg_vertebrates) && nrow(kegg_vertebrates) > 0) {
    vroom::vroom_write(
      kegg_vertebrates,
      here("results/paper_data/MSG_KEGG_enrichment_vertebrates.csv"),
      delim = ","
    )
  }

  vroom::vroom_write(
    kegg_combined,
    here("results/paper_data/MSG_KEGG_enrichment_combined.csv"),
    delim = ","
  )

  cat("KEGG enrichment results saved to results/paper_data/\n")

  # ---------------------------------------------------------------------------
  # KEGG UpSet Plot
  # ---------------------------------------------------------------------------
  cat("\n--- Creating KEGG UpSet Plot ---\n")

  kegg_terms_list <- list(
    "Unicellular" = if (!is.null(kegg_unicellular) && nrow(kegg_unicellular) > 0)
      kegg_unicellular$ID else character(0),
    "Multicellular" = if (!is.null(kegg_multicellular) && nrow(kegg_multicellular) > 0)
      kegg_multicellular$ID else character(0),
    "Vertebrates" = if (!is.null(kegg_vertebrates) && nrow(kegg_vertebrates) > 0)
      kegg_vertebrates$ID else character(0)
  )

  kegg_terms_list <- kegg_terms_list[sapply(kegg_terms_list, length) > 0]

  if (length(kegg_terms_list) >= 2) {
    pdf(here("results/plots/MSG_KEGG_upset_plot.pdf"), width = 10, height = 7)
    print(upset(
      fromList(kegg_terms_list),
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = 1.5,
      point.size = 3.5,
      line.size = 1.5,
      mainbar.y.label = "KEGG Pathways Intersection",
      sets.x.label = "Pathways per Group"
    ))
    dev.off()

    png(
      here("results/plots/MSG_KEGG_upset_plot.png"),
      width = 1000, height = 700, res = 100
    )
    print(upset(
      fromList(kegg_terms_list),
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = 1.5,
      point.size = 3.5,
      line.size = 1.5,
      mainbar.y.label = "KEGG Pathways Intersection",
      sets.x.label = "Pathways per Group"
    ))
    dev.off()

    cat("KEGG UpSet plot saved to results/plots/\n")
  }

  # ---------------------------------------------------------------------------
  # KEGG Clustering and TreeAndLeaf
  # ---------------------------------------------------------------------------
  cat("\n--- Creating KEGG TreeAndLeaf Visualization ---\n")

  # Get unique KEGG pathway IDs
  kegg_ids <- unique(kegg_combined$ID)
  cat(sprintf("Creating clustering for %d unique KEGG pathways...\n", length(kegg_ids)))

  # Calculate distance based on gene overlap (Jaccard similarity)
  # Create a matrix of gene lists per pathway
  kegg_gene_matrix <- matrix(0, nrow = length(kegg_ids), ncol = length(kegg_ids))
  rownames(kegg_gene_matrix) <- kegg_ids
  colnames(kegg_gene_matrix) <- kegg_ids

  for (i in seq_along(kegg_ids)) {
    for (j in seq_along(kegg_ids)) {
      if (i <= j) {
        genes_i <- unlist(strsplit(
          kegg_combined$geneID[kegg_combined$ID == kegg_ids[i]][1], "/"
        ))
        genes_j <- unlist(strsplit(
          kegg_combined$geneID[kegg_combined$ID == kegg_ids[j]][1], "/"
        ))

        intersection <- length(intersect(genes_i, genes_j))
        union <- length(unique(c(genes_i, genes_j)))

        if (union > 0) {
          jaccard <- intersection / union
        } else {
          jaccard <- 0
        }

        kegg_gene_matrix[i, j] <- jaccard
        kegg_gene_matrix[j, i] <- jaccard
      }
    }
  }

  # Convert similarity to distance
  kegg_dist <- as.dist(1 - kegg_gene_matrix)

  # Handle zero variation case
  if (all(kegg_dist == 0)) {
    # Add small random noise
    kegg_dist <- kegg_dist + runif(length(kegg_dist), 0, 0.001)
  }

  # Hierarchical clustering
  hc_kegg <- hclust(kegg_dist, method = "average")

  # Create TreeAndLeaf object
  tal_kegg <- treeAndLeaf(hc_kegg)

  # Prepare KEGG data for visualization
  kegg_combined <- kegg_combined %>%
    mutate(
      path_length = as.integer(sapply(strsplit(GeneRatio, "/"), "[", 2)),
      ratio = Count / path_length
    )

  kegg_term_info <- kegg_combined %>%
    group_by(ID) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(ID, Description, Count, ratio, pvalue, p.adjust, evolutionary_group)

  # Organize node sizes
  kegg_size <- kegg_term_info$Count
  kegg_aux <- sort(unique(kegg_size))
  names(kegg_aux) <- as.character(1:length(kegg_aux))
  kegg_sizeRanks <- as.factor(names(kegg_aux[match(kegg_size, kegg_aux)]))
  kegg_sizeIntervals <- min(7, length(unique(kegg_size)))
  kegg_term_info$size <- (sizeBase + (as.numeric(kegg_sizeRanks) * sizeMultiplier))

  kegg_term_info$group_color <- group_colors[kegg_term_info$evolutionary_group]

  # Map attributes to KEGG graph
  tal_kegg <- att.mapv(g = tal_kegg, dat = as.data.frame(kegg_term_info), refcol = 1)

  # Set node labels
  tal_kegg <- att.setv(g = tal_kegg, from = "Description", to = "nodeAlias")

  # Set node colors based on ratio
  pal_kegg <- brewer.pal(9, "YlOrRd")
  tal_kegg <- att.setv(
    g = tal_kegg,
    from = "ratio",
    to = "nodeColor",
    cols = pal_kegg,
    nquant = 5
  )

  # Set node sizes
  tal_kegg <- att.setv(
    g = tal_kegg,
    from = "size",
    to = "nodeSize",
    xlim = c(
      (sizeBase + sizeMultiplier),
      (sizeBase + (sizeMultiplier * kegg_sizeIntervals)),
      sizeMultiplier
    ),
    nquant = kegg_sizeIntervals
  )

  # Set font size
  tal_kegg <- att.addv(tal_kegg, "nodeFontSize", value = 15, index = V(tal_kegg)$isLeaf)

  # Set edge attributes
  tal_kegg <- att.adde(tal_kegg, "edgeWidth", value = 3)
  tal_kegg <- att.adde(tal_kegg, "edgeColor", value = "gray80")

  # Set nodeLineColor based on evolutionary group
  kegg_leaf_indices <- which(V(tal_kegg)$isLeaf)

  for (i in kegg_leaf_indices) {
    node_name <- V(tal_kegg)$name[i]
    node_group <- kegg_term_info %>%
      filter(ID == node_name) %>%
      pull(evolutionary_group)

    if (length(node_group) > 0 && !is.na(node_group)) {
      V(tal_kegg)$nodeLineColor[i] <- group_colors[node_group]
    } else {
      V(tal_kegg)$nodeLineColor[i] <- "gray50"
    }
  }

  tal_kegg <- att.addv(tal_kegg, "nodeLineWidth", value = 4, index = V(tal_kegg)$isLeaf)

  kegg_non_leaf_indices <- which(!V(tal_kegg)$isLeaf)
  V(tal_kegg)$nodeLineColor[kegg_non_leaf_indices] <- "gray80"
  V(tal_kegg)$nodeLineWidth[kegg_non_leaf_indices] <- 1

  # Display KEGG TreeAndLeaf in RedeR
  cat("\n--- Starting KEGG RedeR Visualization ---\n")

  # Add to existing RedeR session
  Sys.sleep(2)  # Wait for previous visualization
  addGraph(rdp, tal_kegg, layout = NULL)

  addLegend.color(
    obj = rdp,
    tal_kegg,
    title = "KEGG Gene Ratio",
    position = "topleft"
  )

  relax(rdp, p1 = 21, p2 = 70, p3 = 255, p4 = 100, p5 = 50)

  cat("KEGG TreeAndLeaf visualization added to RedeR.\n")

  # KEGG Summary
  cat("\n--- KEGG Analysis Summary ---\n")
  cat(sprintf("Total enriched KEGG pathways: %d\n", nrow(kegg_combined)))
  cat("KEGG pathways per group:\n")
  print(table(kegg_term_info$evolutionary_group))

} else {
  cat("\nNo KEGG pathways were enriched in any group.\n")
}

cat("\n=== Analysis Complete ===\n")
