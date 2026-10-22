# Plasticity Associated with Adoption of Social Roles in Clown Anemonefish

Code and data associated with study exploring how juvenile clown anemonefish (*Amphiprion percula*) undergo coordinated changes in growth, feeding, and agonistic behaviors, and gene expression as they adopt dominant, subordinate, or solitary social roles.

**Citation:** Vizer LF, Alvarado D, Bove C, Herrera M, Hughes A, Thompson K, Bogdanowicz SM, Laudet V, Davies SW, Buston PM (2026). Plasticity Associated with Adoption of Social Roles in Clown Anemonefish. *eLife* 15:RP109871. https://doi.org/10.7554/eLife.109871.3  

If you use the code or data in this repository, please cite both the paper above and the archived release: https://doi.org/10.5281/zenodo.22817747

[![DOI](https://zenodo.org/badge/1000510336.svg)](https://doi.org/10.5281/zenodo.22817747)

----
## Repository Structure

- `Code/` – R scripts for all analyses and figure generation  
- `Data/` – Processed data files used in analyses  
- `Manuscript_Figures/` – Final figures of main text and supplementary materials
- `Genotyping/` – Genotyping files and related metadata  

## Raw Data

Raw FASTQ files are publicly available via the [NCBI SRA BioProject PRJNA1280023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1280023).

## Key Analysis Scripts

- **Growth, behavior, and food intake analysis:**  
  `Code/percula_phenotypic_analysis_20250428.Rmd`

- **Gene expression analysis:**  
  `Code/percula_phenotypic_plasticity_GE.Rmd`

- **Pathway enrichment analysis:**  
  `Code/pathway_analysis20250229.R`

Additional supporting scripts are included in the `Code/` directory and are referenced within the main R Markdown files.
