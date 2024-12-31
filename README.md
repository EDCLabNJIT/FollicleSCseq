# Single cell RNA-seq reveals that granulosa cells are a target of phthalate toxicity in the ovary

## Abstract

abstract here

[link to publication](ADD%20LINK%20HERE)

## Citation

citation here

## Contents

#### wrapper.R

A pipeline to run the scripts for all analysis done.

#### scRNAhelper

Acts a package that runs each step of the analysis through functions to be called by `wrapper.R`

#### morrismarkers

Marker genes compiled from [Morris et al.](https://doi.org/10.7554/elife.77239)

#### DAVIDresults

output from using DAVID for gene enrichment, as well as lists of pathways for the figure used. Only for plotting

## Setup

-   [Maurine atlas Seurat file](https://singlecell.broadinstitute.org/single_cell/study/SCP1914/a-single-cell-atlas-of-the-cycling-murine-ovary) from [Morris et al.](https://doi.org/10.7554/elife.77239)

    -   Named `ovary_0.rds`

-   10x Genomics Cellranger output. Accessible [here](ADD%20LINK%20HERE)

    -   Named `raw_feature_bc_matrix.h5`
