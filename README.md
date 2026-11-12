# Reconstruction of the Evolutionary Landscape of Biological Processes Involved in the Early Stages of the Metastatic Cascade

This repository hosts the code, data and results supporting the article:

> Azevedo GM, Farias E, Cavalcante JVF, Ferraz RS, William B, Coelho DM and Dalmolin RJS (2026) Reconstruction of the evolutionary landscape of biological processes involved in the early stages of the metastatic cascade. *Genetics and Molecular Biology* 49(Suppl 2):e20250197. https://doi.org/10.1590/1678-4685-GMB-2025-0197

## Abstract

Metastasis is not a *de novo* functional module innovation but rather the pathological redeployment of deeply conserved biological programs. Here, we reconstruct the evolutionary landscape of biological functions involved in the early stages of the metastatic cascade, including cell adhesion, extracellular matrix organization (ECM), regulation of metallopeptidase activity, cell junction organization, epithelial–mesenchymal transition (EMT), and cellular extravasation alongside the physiological constraints that suppress them (Metastasis Suppressor Genes). Using phyletic pattern reconstruction of 787 orthologs across 473 Clusters of Orthologous Groups with the Bridge algorithm, we identified a divergence between the evolutionary timelines of metastasis-enabling programs and metastasis suppressors. Our results indicate that the molecular systems associated with the structural capacity for tumor dissemination are evolutionary ancient. ECM organization traces to the Human–Discoba last common ancestor, EMT to Ctenophora, cell adhesion to Sauropsida, and cellular extravasation to Actinopterygii. In contrast, mechanisms responsible for suppressing genomically unstable cells during continuous tissue renewal emerged later, peaking in the vertebrate lineage. These findings support the Serial Atavism Model and suggest that metastasis arises from the progressive erosion of recently evolved regulatory constraints, allowing the reactivation of ancestral cellular programs that predate complex multicellularity.

## Repository Structure

```
├── analysis        # Quarto documents with the main analyses (01–06)
├── assets          # Input data required to run the analyses
└── results
    ├── enrichment         # MSG functional enrichment (GO BP and KEGG)
    ├── exclusive_genes    # Rooting of single-process (exclusive) genes
    ├── Figures            # Article figures
    ├── metastasis_genes   # Gene sets retrieved from Gene Ontology
    ├── orthology_data     # STRING/COG mapping and GeneBridge rooting results
    ├── paper_data         # Tables and data used in the article
    └── rooting_analysis   # Rooting distribution and hypergeometric tests
```

### Analyses

| Script | Description |
|---|---|
| `01_metastais_GO_query.qmd` | Gene Ontology query and construction of the metastasis-associated gene set |
| `02_orthology_rooting_analysis.qmd` | STRING/COG mapping and evolutionary rooting with GeneBridge |
| `03_Plotting_roots.qmd` | Rooting distribution, statistical tests and visualization |
| `04_exclusive_genes_analysis.qmd` | Rooting of genes exclusive to a single biological process |
| `05_MSG_enrichment.qmd` | Functional enrichment of metastasis suppressor genes per evolutionary group |
| `06_MSG_enrichment_cumulative.qmd` | Cumulative functional enrichment of metastasis suppressor genes |
