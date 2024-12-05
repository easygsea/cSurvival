# cSurvival: a framework with an integrated survival database to study gene set, gene, locus, and drug dependencies and interactions in cancers

This repository contains essential scripts for:

* cSurvival's algorithms, pipelines, and visualizations (built on R and R Shiny)
* Data extraction and processing steps to establish cSurvival's database (shell and R scripts)

## Quick start

Step 1. Install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/)

Step 2. cSurvival cannot start without its underlying database. Download example [TCGA-LUAD data sets (.zip)](https://tau.cmmt.ubc.ca/cSurvival/project_data.zip), decompress the files, and put the entire TCGA-LUAD folder into the designated place under www/project_data

![folders](https://tau.cmmt.ubc.ca/cSurvival/Presentation1.jpg)

Step 3. To start cSurvival, open its global.R file in RStudio, click **Run App** on top right. If first time use, it may take a while to install required packages.

# Citation
Cheng, X., Liu, Y., Wang, J., Chen, Y., Robertson, A.G., Zhang, X., Jones, S.J. and Taubert, S., 2022. cSurvival: a web resource for biomarker interactions in cancer outcomes and in cell lines. Briefings in Bioinformatics, 23(3), p.bbac090.
https://doi.org/10.1093/bib/bbac090
