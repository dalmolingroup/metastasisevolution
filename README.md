# Reconstruction of the Evolutionary Landscape of Biological Processes Involved in the Early Stages of the Metastatic Cascade

This repository hosts the analysis and resources supporting the study of the evolutionary origins of biological processes underlying the early stages of the metastatic cascade.

Metastasis is a multistep process through which tumor cells disseminate from the primary site and colonize distant tissues. While the molecular mechanisms of metastasis are well studied, the evolutionary history of the genes and pathways enabling these processes remains underexplored. This project addresses that gap using phylogenetic reconstruction.

## Abstract

This study investigates the evolutionary trajectory of orthologous genes involved in six fundamental biological processes associated with the initial metastatic cascade: cell adhesion, extracellular matrix (ECM) organization, epithelial-mesenchymal transition (EMT), cell junction organization, metalloproteinase regulation, and cellular extravasation.

Using Gene Ontology (GO)-based selection, 668 protein-coding genes were mapped across 476 eukaryotic species and rooted with the GeneBridge algorithm. The results reveal a layered evolutionary assembly of metastasis-related modules:

-   **Human-Discoba LCA**: Early enrichment of ECM organization genes linked to aggregative multicellularity.\
-   **Human-Porifera LCA**: Emergence of junctional machinery for epithelial integrity.
-   **Human-Choanoflagellata**: Genes involved in cell-cell and extracellular matrix adhesion were rooted.\
-   **Human-Actinopterygii**: Genes tied to cell extravasation and the major histocompatibility complex were predominantly rooted.

Findings suggest that metastasis does not represent a novel cancer-specific innovation, but rather the pathological reconfiguration of deeply conserved biological programs, progressively assembled during eukaryotic and metazoan evolution.

## Repository Structure
├── _dev                     
├── analysis                                 
├── results                  
│   ├── Figures              
│   ├── metastasis_genes     
│   ├── orthology_data       
│   └── plots                
└── assets                

- **_dev**  
  Contains scripts for analyses currently under development.  

- **analysis**  
  Includes Quarto documents with the main analyses.  

- **assets**  
  Provides required data files to run the scripts.  

- **results**  
  Stores the outputs of the analyses. Subdirectories include:  
  - **Figures**: Visualization outputs generated during the project.  
  - **metastasis_genes**: Results focusing on genes implicated in metastasis-related pathways.  
  - **orthology_data**: Results of orthologous genes in the evolutionary analyses.  
  - **plots**: Generated plots for exploratory or final reporting purposes.  

## Notes  

This repository is a research compendium ensuring reproducibility and extensibility of the analyses. The project is still under development, and results may evolve as the research progresses.  
