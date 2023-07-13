# sarcoma_differentiation_landscape

This GitHub repository contains all relevant analysis code for the publication:

Mapping the Single-cell Differentiation Landscape of Osteosarcoma 
Danh Dinh Truong*, Corey Weistuch*, Kevin Murgas*, Joseph Deasy, Antonios G. Mikos, Allen Tannenbaum+, and Joseph Ludwig+
Under submission.

All code was developed in joint effort by Danh Truong, Kevin Murgas, and Corey Weistuch. Any questions about code implementation are welcome.

A summary of the code used in the following analysis steps:
1. Preparation of RNA-sequencing data (including quality control and pre-processing):
    - analysis_OSPDX.Rmd - prepare single-cell RNA-seq data of 3 osteosarcoma patient-derived xenografts (citation: Truong et al. BMC Cancer 2023)
    - analysis_MTL.Rmd - prepare single-cell RNA-seq data of osteogenesis and adipogenesis differentiation experiment (experimentally generated for this study)
    - prep_GSE160625_chondro.R - prepare single-cell RNA-seq data from a chondrogenesis differentiation experiment (citation: Wu et al. Nat Commun 2021)
    - prep_GSE152048_OS11.R - prepare single-cell RNA-seq data from 11 adult osteosarcoma tumors (citation: Zhou et al. Nat Commun 2020)
    - prep_targetOS.R - prepare bulk RNA-seq data from the TARGET OS pediatric sarcoma study
2. Gene expression data SCTransform and quantile normalization:
    - mtl_chondro_SCT_qnorm.R - apply SCTransform and merge MTL and chondrogenesis datasets, then quantile normalize merged data
    - ospdx_SCT_qnorm.R - apply SCTransform then quantile normalize OSPDX data
    - os11_SCT_qnorm.R - apply SCTransform then quantile normalize OS11 data
3. Archetype analysis via Non-negative Matrix Factorization (NMF):
    - *NOTE*: Expression data needs to be normalized before computing archetypes. For this study, we apply Seurat's SCTransform (based on negative binomial regression) and then quantile normalization.
    - Archetype_Analysis.m - performs all archetype computations presented in the above manuscript
    - NMF.m - computes the archetype scores [H] and weights [W] from input data [V] and a specified number of archetypes [k]
    - NMF_New_Weights.m - computes new archetype scores [H] from input weights [W], data [V], and the appropriate number of archetypes [k]
    - Archetype_Computer_MTL.m - computes the archetype scores [H] associated with our MTL-trained archetypes from new data [Data]
4. Figure generation:
    - figure2.R
    - figure3_4.R
    - figure5_targetos_survival.R
