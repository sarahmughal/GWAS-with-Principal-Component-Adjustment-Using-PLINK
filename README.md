# GWAS-with-Principal-Component-Adjustment-Using-PLINK

This project implements an end-to-end genome-wide association study (GWAS) pipeline using PLINK and R, including quality control, minor allele frequency assessment, Hardy–Weinberg equilibrium testing, LD pruning, and principal component analysis (PCA). The analysis compares naive and ancestry-adjusted models to demonstrate how population structure can confound association results and how PCA-based correction improves inference.

---

- Set up a reproducible project structure for GWAS analysis using PLINK and R
- Perform quality control, including minor allele frequency (MAF) estimation and Hardy–Weinberg equilibrium (HWE) testing
- Apply linkage disequilibrium (LD) pruning to remove correlated SNPs
- Conduct principal component analysis (PCA) to assess population structure
- Visualize genetic clustering and evaluate evidence of ancestry stratification
- Run a naive genome-wide association study (GWAS) using linear regression
- Adjust GWAS models for ancestry using principal components as covariates
- Compare naive and PC-adjusted results using Manhattan and Q-Q plots to assess inflation and confounding

---

## Data Sources

The dataset consists of simulated genotype data in PLINK format (.bed/.bim/.fam) and a continuous phenotype. The data were generated to illustrate key GWAS concepts, including population structure and ancestry confounding.

---

## Project Structure

```text
project/
├── software/      # plink binary goes here
├── data/          # .bed, .bim, .fam
├── output/        # PCA results, GWAS logs
└── scripts/       # scripts
