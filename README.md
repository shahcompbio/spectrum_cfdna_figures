# ctDNA Analysis - Williams et al. Nature

This repository contains code for generating results from Williams et al. Nature: **[Paper Title Placeholder]** ([Paper Link Placeholder])

## Data Requirements

All analyses are performed in R. Before running the analysis scripts, you need to download the required data:

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

**Figure02_detectingclonalSVs.Rmd** - Analyzes error rates and detection performance of structural variants (SVs) and SNVs in ctDNA, comparing duplex vs uncorrected sequencing methods. Generates correlation plots between TP53 and clonal SV tumor fractions.

**Figure03_subclonalsvs.Rmd** - Creates integrated visualizations of phylogenetic trees, copy number profiles, and subclonal SV frequencies across multiple patient samples. Produces comprehensive plots showing clone-specific genomic alterations and their detection in cfDNA.

**Figure04_mullerplots.Rmd** - Generates longitudinal Muller plots showing clone frequency dynamics over time, integrated with clinical data including treatments, CA-125 levels, and recurrence events. Also analyzes clone-specific genetic alterations like ERBB2 amplifications.

**Figure05_scrna.Rmd** - Analyzes single-cell RNA sequencing data to compare pathway activities and gene expression between tumor clones. Creates UMAP visualizations and violin plots showing differential expression of key genes and pathways across clones.

**Figure06_wrightfisher.Rmd** - Implements Wright-Fisher population genetics simulations to model clone frequency changes over time and test for neutral evolution. Compares observed clone dynamics against neutral expectations using CA-125 levels as a proxy for tumor burden.