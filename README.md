# ctDNA Analysis - Williams et al. Nature
[![DOI](https://zenodo.org/badge/1013093862.svg)](https://doi.org/10.5281/zenodo.15798215)

This repository contains code for generating results from Williams et al. Nature: **Tracking clonal evolution during treatment
in ovarian cancer using cell-free DNA** ([Paper Link Placeholder])

## Data Requirements

All analyses are performed in R. Before running the analysis scripts, you need to download the required data. Once you've downloaded the data, change the config file to point to the root of the synapse directory and change the path to the tables file.

### Synapse Data
Download the project data from Synapse:
```bash
# Install synapseclient if not already installed
pip install synapseclient

# Login to Synapse (you'll need to create an account and get access to the project)
synapse login

# Download the project data
synapse get syn66399325 -r
```

### Supplementary Tables
Download the supplementary tables from the published paper (available from the journal website or supplementary materials).

## Docker Environment

A Docker image with all required R dependencies is available:
```bash
docker pull marcjwilliams1/signals:latest
```

## Main figures

**Figure01_OV-004.Rmd** - scWGS from patient OV-004, showing clone specific copy number differences, SVs and SNVs in chromosomes 8 and 17.

**Figure02_detectingclonalSVs.Rmd** - Analyzes error rates and detection performance of structural variants (SVs) and SNVs in ctDNA, comparing duplex vs uncorrected sequencing methods. Generates correlation plots between TP53 and clonal SV tumor fractions.

**Figure03_subclonalsvs.Rmd** - Creates integrated visualizations of phylogenetic trees, copy number profiles, and subclonal SV frequencies across multiple patient samples. Produces comprehensive plots showing clone-specific genomic alterations and their detection in cfDNA.

**Figure04_mullerplots.Rmd** - Generates longitudinal Muller plots showing clone frequency dynamics over time, integrated with clinical data including treatments, CA-125 levels, and recurrence events. Also analyzes clone-specific genetic alterations like ERBB2 amplifications.

**Figure05_scrna.Rmd** - Analyzes single-cell RNA sequencing data to compare pathway activities and gene expression between tumor clones. Creates UMAP visualizations and violin plots showing differential expression of key genes and pathways across clones.

**Figure06_wrightfisher.Rmd** - Implements Wright-Fisher population genetics simulations to model clone frequency changes over time and test for neutral evolution. Compares observed clone dynamics against neutral expectations using CA-125 levels as a proxy for tumor burden.