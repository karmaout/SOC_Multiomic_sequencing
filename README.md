# SOC_Multiomic_sequencing

# Abstract

Fear, while crucial for survival, is a component of a myriad of psychiatric illnesses in its extreme. Persistent fear memories can form through processes such as second-order conditioning (SOC), during which a second-order conditioned stimulus (CS2) acquires significance by associating with a first-order conditioned stimulus (CS1).

The neural circuitry underlying SOC, particularly the roles of sensory cortices, remains poorly understood. We explored the mechanisms of olfactory SOC in rats, focusing on the basolateral amygdala (BLA) and posterior piriform cortex (pPC). Our results demonstrate that NMDAR-dependent plasticity in both regions is essential for SOC. 

The BLA mediates the CS2-CS1 association, while the pPC, receiving inputs from the locus coeruleus and BLA, is critical for memory encoding and retrieval. Transcriptomic analysis of Fos+ ensembles in the two structures shows distinct gene activation in excitatory neurons, highlighting key pathways involved in learning and memory.

These findings highlight the pPC’s role in integrating sensory and affective information, essential for long-term threat memory consolidation.
This repository provides scripts used for our study. 

The inventory contains the codes for analyzing the Multiomi-seq data of the study, the intermediate data files are generated and deposited in Zenodo inventory: https://zenodo.org/uploads/14014437

| **Name**                                | **Description**                                                                                          |
|----------------------------------------- |--------------------------------------------------------------------------------------------------------- |
| `ATAC_motif_discovery.R`                | R script for discovering motifs in snATAC-seq data.                                                      |
| `ATAC_processing.R`                     | R script for preprocessing and analyzing snATAC-seq data.                                                |
| `BLA_QC&clustering.R`                   | R script for performing quality control and clustering analysis on BLA cells.                            |
| `BLA_neurons_subsetting_analysis.R`     | R script for subsetting BLA neurons and performing targeted analysis on specific groups or conditions.   |
| `Gene_ontology_analysis.R`              | R script for conducting Gene Ontology (GO) analysis to identify biological processes associated with DEGs. |
| `LOESS_combination.py`                  | Python script for combining and running LOESS regression analyses on integrated data.                    |
| `Motif_LOESS_regression.py`             | Python script for performing LOESS regression analysis on motif data.                                    |
| `PC_QC&clustering.R`                    | R script for quality control and clustering of PC cells.                                                 |
| `PC_neurons_subsetting_analysis.R`      | R script for subsetting PC neurons and performing targeted analysis on specific groups or conditions.    |
