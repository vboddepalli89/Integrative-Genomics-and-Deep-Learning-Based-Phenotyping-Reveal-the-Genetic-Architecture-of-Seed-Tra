# Integrative Genomics and Deep Learning-Based Phenotyping of Mungbean Seed Traits

This repository contains the statistical and genomic analysis pipelines for the manuscript: **"Integrative Genomics and Deep Learning-Based Phenotyping Reveal the Genetic Architecture of Seed Traits in Mungbean (Vigna radiata L.)"** submitted to *The Plant Genome*.

This project integrates high-throughput image-based phenotyping with Genome-Wide Association Studies (GWAS) and Genomic Prediction (GP) to dissect the genetic architecture of seed size, shape, and weight in the Iowa Mungbean Diversity (IMD) panel.

## 🔗 Deep Learning Phenotyping Pipeline (SAM)
The raw seed images were processed and traits were extracted using a zero-shot Segment Anything Model (SAM) pipeline. 
* **The complete image processing pipeline, custom Python scripts, and configuration files are hosted separately in the Baskar Group repository here:** [MungSeedAI - GitHub](https://github.com/baskargroup/MungSeedAI/tree/master).
* --

## 🧬 Genomic & Statistical Analysis Scripts

This repository contains the R scripts used for all downstream phenotypic and genomic analyses. They are numbered sequentially to reflect the workflow:

* **`01_Descriptive_Statistics.R`**: Performs data quality control, outlier removal, Box-Cox transformations, and extracts Best Linear Unbiased Estimators (BLUEs) using mixed linear models (`lme4`). It also calculates broad-sense heritability ($H^2$) and phenotypic variance components.
* **`02_GWAS_MultiModel_GAPIT.Rmd`**: Conducts multi-model GWAS using the `GAPIT3` package in R. Evaluates three frequentist models: BLINK, MLM, and FarmCPU to control for false discoveries.
* **`03_GWAS_SVEN.R`**: Conducts Bayesian variable selection GWAS using the Selection of Variables with Embedded Screening (SVEN) model via the `bravo` package to identify true marker-trait associations.
* **`04_Post_GWAS_calculations.R`**: Calculates the phenotypic variance explained (PVE), allele substitution effects, and standard errors for the significant SNPs identified across the GWAS models. 
* **`05_Genomic_prediction_GBLUP_analysis.R`**: An integrated pipeline using `rrBLUP` to evaluate Genomic Prediction Accuracies (GPA). It performs 100 cycles of 10-fold cross-validation for both standard **gBLUP** and marker-assisted **gBLUP+SNPs** models, and computes the top-10 Selection Coincidence Index (CI).

## 💻 Software & Dependencies

The scripts were written and executed using **R (v4.0+)**. The following primary R packages are required to run the pipelines:
* **Phenotypic Analysis:** `lme4`, `lmerTest`, `emmeans`, `car`, `MASS`, `tidyr`, `dplyr`
* **GWAS:** `GAPIT3`, `bravo`, `Matrix`
* **Genomic Prediction:** `rrBLUP`
* **Visualization:** `ggplot2`, `patchwork`, `qqman`

## 📊 Data Availability
* **Genotypic Data:** The raw Genotyping-by-Sequencing (GBS) data aligned to the Weilv-9 reference genome is available at the NCBI Sequence Read Archive (SRA) under BioProject PRJNAXXXXXX *(Update with actual link prior to publication)*.
* **Phenotypic Data:** Processed phenotypic BLUEs and extracted traits are provided in the supplementary files of the manuscript. 

## ✉️ Contact
For questions regarding the genomic analysis, please contact **Venkata Naresh Boddepalli** (naresh@iastate.edu). For questions regarding the SAM image extraction pipeline, please refer to the [MungSeedAI repository](https://github.com/baskargroup/MungSeedAI/tree/master).
