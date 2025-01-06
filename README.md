# Single cell RNA-seq reveals that granulosa cells are a target of phthalate toxicity in the ovary

## Abstract

Phthalates are known endocrine disrupting chemicals and ovarian toxicants that are used widely in consumer products. Phthalates have been shown to exert ovarian toxicity on multiple endpoints, altering transcription of genes responsible for normal ovarian function. However, the molecular mechanisms by which phthalates act on the ovary are not well understood. In this study, we hypothesized that phthalates specifically target granulosa cells within the ovarian follicle. To test our hypothesis, we cultured whole mouse antral follicles for 96 hours in the presence of vehicle or 10 µg/mL of a phthalate metabolite mixture. At the end of the culture period, follicles were dissociated into single cell suspensions and subjected to single cell RNA sequencing. We used markers from published studies to identify cell type clusters, the largest of which were granulosa and theca/stroma cells. We further identified sub-populations of granulosa, theca, and stromal cells and analyzed differentially expressed genes between the phthalate treatment and control. Granulosa cells, specifically mural granulosa cells, had the most differentially expressed genes. Pathway analysis of differentially expressed genes from the overall granulosa cell cluster revealed disruption of cell cycle and mitosis, whereas pathway analysis of the mural granulosa cell subcluster identified terms related to translation, ribosome, and endoplasmic reticulum. Our findings suggest that phthalates have both broad impacts on cell types and specific impacts on cellular subtypes, emphasizing the complexity of phthalate toxicity and highlighting how bulk sequencing can mask effects on vulnerable cell types.

[link to publication](https://doi.org/10.1093/toxsci/kfaf001)

## Citation



Erik Mattson, Genoa R Warner, Single cell RNA-seq reveals that granulosa cells are a target of phthalate toxicity in the ovary, Toxicological Sciences, 2025, kfaf001, https://doi.org/10.1093/toxsci/kfaf001


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
