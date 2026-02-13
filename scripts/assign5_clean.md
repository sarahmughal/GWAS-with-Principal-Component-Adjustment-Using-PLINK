---
title: "Molecular & Genetic Epidemiology - Assignment 5"
author: "Sarah Mughal"
date: last-modified
format:
  docx:
    toc: true
    number-sections: true
    highlight-style: github
---
<!--
cd scripts
quarto render assign5_clean.md
-->

Due by noon on the day before next class.

(I'm a little late with getting this out this time, so please just let me know if you need a little more time)

# Instructions

- **Both tracks**: Choose your preferred language (python, R, or stata) and run a GWAS using the simulated data.
- Please answer anything marked QUESTION, and fill-in code marked TODO.
- I was trying to simulate something special with the data, which made it a little hard for me to get the correct allele frequency distribution. Sorry.


# Suggested data structure

A suggested data structure for this assignment (that will work with the code), is as follows:

```
project/
├── software/      # plink binary goes here
├── data/          # .bed, .bim, .fam
├── output/        # PCA results, GWAS logs
└── scripts/       # scripts
```

# python env & libraries

I don't think this one is very dependent on this choice.

```
pyenv local 3.12.10
```

# Simulating the data

<!--
```
python assign5_sim_data_v06.py
```
-->

<!--
```
# quick check on the data
library(data.table)
bim <- fread("data/homework.bim")
table(bim[[1]])
```
-->


# Option 1: Analyzing it in plink + R

```
# -- load R libraries --
library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)


# -- setup file structure --
# You should be in the scripts folder (needed for markdown)
try(dir.create("../output"))
try(dir.create("../software"))
# scripts should have this file, and data should have the data


# -- download plink 1.9 --
if (!file.exists("plink")) {
  message("[*] Downloading PLINK 1.9...")
  curwd <- getwd()
  dir.create("../software")
  setwd("../software")

  # For Mac (change to 'linux-x86_64' or 'win64' as needed)
  download.file("https://s3.amazonaws.com/plink1-assets/plink_mac_20231211.zip", "plink.zip")
  unzip("plink.zip")

  # Probably need to do on mac/linux
  system("chmod u+x plink") # make executable
  system("ls -lh")

  setwd(curwd)
}


# -- QC --
# allele frequencies - run in plink
system("../software/plink --bfile ../data/homework --freq --out ../output/stats_freq")

# - and then load into R
freq <- fread("../output/stats_freq.frq")
head(freq)

# TODO: Create a histogram of MAFs

# QUESTION: What is the range of the frequency you see in these plots, and why?

# Delete the dataset to save memory, sprinkle throughout if need be
rm(freq)

# Hardy-Weinberg Equilibrium
system("../software/plink --bfile ../data/homework --hardy --out ../output/stats_hwe")

hwe <- fread("../output/stats_hwe.hwe", header=TRUE)

# TODO: Create a histogram of the P-values.


# QUESTION: What is the median value?
# QUESTION: Any guesses why this might be happening?

# NOTE: You would usually filtering for Minor Allele Frequency (MAF), and we'll do that here. You'd do a different cutoff depending on what you are doing

# -- PCA --
# - We must 'prune' for LD before PCA, otherwise LD blocks will skew the clusters.
# - Others will argue to remove HLA and other regions of long LD
#   - Our simulated data is not quite that accurate
# - Note that MAF is set to 5% here
system("../software/plink --bfile ../data/homework --maf 0.05 --indep-pairwise 50 5 0.2 --out ../output/pruning")
# - Then run PCA (this takes a little longer, but we have a pretty small dataset)
system("../software/plink --bfile ../data/homework --extract ../output/pruning.prune.in --pca 10 --out ../output/pca_results")


# - load in the PCA eigenvectors
# Note: PLINK .eigenvec files usually have FID and IID in the first two columns
pcs <- read.table("../output/pca_results.eigenvec", header=FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

# TODO: Plot at least PC1 vs. PC2 here


# Look at the phenotype we have
phe <- fread("../data/homework.fam")
head(phe)
phe <- phe[, c(1,2,6)]
colnames(phe) <- c("FID", "IID", "Pheno")
# TODO histogram of phenotype# TODO: summarize the phenotype

fwrite(phe, file="../output/homework.phe", quote=FALSE, sep="\t")

# -- The GWAS Comparison --

# GWAS A: Naive (No covariates)
system("../software/plink --bfile ../data/homework --maf 0.005 --linear --pheno ../output/homework.phe --mpheno 1  --out ../output/gwas_naive")

gwas_res <- fread("../output/gwas_naive.assoc.linear", header=TRUE)
table(gwas_res$CHR)

manhattan(gwas_res,
          chr="CHR", bp="BP", p="P", snp="SNP",
          main="GWAS Results (raw)",
          col=c("orange", "blue"), # Alternating colors
          suggestiveline = FALSE,
          genomewideline = -log10(5e-8), # Red line
          cex = 0.6) # Size of points

qq(gwas_res$P, main="Q-Q: Naive")

# GWAS B: GWAS adjusted for ancestry (this takes longer)
fwrite(pcs, "../output/covariates.txt", row.names=FALSE, quote=FALSE, sep="\t")
system("../software/plink --bfile ../data/homework --linear --pheno ../output/homework.phe --mpheno 1 --covar ../output/covariates.txt --hide-covar --out ../output/gwas_corrected")

gwas_res <- fread("../output/gwas_corrected.assoc.linear")
table(gwas_res$CHR)

## TODO Q-Q PLOT, MANHATTEN
```

# Analyzing it in python + R

```
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


# -- setup file structure --
# you should be in the scripts folder
directories = ["../output", "../software"]
for d in directories:
    os.makedirs(d, exist_ok=True)


# -- download PLINK 1.9 --
plink_path = "../software/plink"
if not os.path.exists(plink_path):
    print("[*] Downloading PLINK 1.9...")
    # Change URL to 'linux-x86_64' or 'win64' if necessary
    url = "https://s3.amazonaws.com/plink1-assets/plink_mac_20231211.zip"
    subprocess.run(["curl", "-L", url, "-o", "../software/plink.zip"], check=True)
    subprocess.run(["unzip", "-o", "../software/plink.zip", "-d", "../software"], check=True)
    subprocess.run(["chmod", "u+x", plink_path], check=True)


# -- QC --
# allele frequencies - run in plink
subprocess.run([plink_path, "--bfile", "../data/homework", "--freq", "--out", "../output/stats_freq"], check=True)


freq = pd.read_csv("../output/stats_freq.frq", sep=r'\s+')

# TODO: Create a histogram of MAFs

# QUESTION: What is the range of the frequency you see in these plots, and why?


# Hardy-Weinberg Equilibrium
subprocess.run([plink_path, "--bfile", "../data/homework", "--hardy", "--out", "../output/stats_hwe"], check=True)

# TODO: Create a histogram of the HWE p-values
hwe = pd.read_csv("../output/stats_hwe.hwe", sep=r'\s+')

# QUESTION: What is the median value?
# QUESTION: Any guesses why this might be happening?

# NOTE: You would usually filtering for Minor Allele Frequency (MAF), and we'll do that here. You'd do a different cutoff depending on what you are doing

# -- PCA --
# - We must 'prune' for LD before PCA, otherwise LD blocks will skew the clusters.
# - Others will argue to remove HLA and other regions of long LD
#   - Our simulated data is not quite that accurate
# - Note that MAF is set to 5% here
subprocess.run([plink_path, "--bfile", "../data/homework", "--maf", "0.05",
                "--indep-pairwise", "50", "5", "0.2", "--out", "../output/pruning"], check=True)

# - Then run PCA (this takes a little longer, but we have a pretty small dataset)
subprocess.run([plink_path, "--bfile", "../data/homework", "--extract", "../output/pruning.prune.in",
                "--pca", "10", "--out", "../output/pca_results"], check=True)

# Load PCA Results
pcs = pd.read_csv("../output/pca_results.eigenvec", sep=r'\s+', header=None)
pcs.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, 11)]

# Plot PC1 vs PC2
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pcs, x='PC1', y='PC2', alpha=0.4, color='darkslategrey', s=20)
plt.title("PCA: Genetic Population Structure\nThree distinct clusters revealed")
plt.show()

# -- Phenotype prep --
phe = pd.read_csv("../data/homework.fam", sep=r'\s+', header=None)
phe = phe[[0, 1, 5]]
phe.columns = ["FID", "IID", "Pheno"]

# TODO histogram of phenotype
# TODO: summarize the phenotype

phe.to_csv("../output/homework.phe", sep="\t", index=False)


# -- The GWAS Comparison --

# there's no standard package for manhattan plots, so we
#  provide code to create them
# honestly, the R code result looked pretty hideous,
#  I would always just do this anyway
def plot_gwas(df, title):
    # Manhattan Plot Prep
    df['-log10P'] = -np.log10(df['P'])
    df = df.sort_values(['CHR', 'BP'])
    df['pos'] = range(len(df)) # Simple index for X-axis

    # Manhattan
    colors = ['orange', 'blue']
    for i, chr_num in enumerate(df['CHR'].unique()):
        subset = df[df['CHR'] == chr_num]
        plt.scatter(subset['pos'], subset['-log10P'], c=colors[i % 2], s=5)
    plt.axhline(-np.log10(5e-8), color='red', linestyle='--')
    plt.show()

def plot_qq(p_values, title="Q-Q Plot"):
    # 1. Clean data: Remove NaNs and ensure p > 0
    p_values = p_values[~np.isnan(p_values) & (p_values > 0)]
    n = len(p_values)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(np.sort(p_values))
    plt.figure(figsize=(6, 6))
    plt.scatter(expected, observed, s=12, color='blue', alpha=0.5, edgecolors='none')
    #max_val = max(max(expected), max(observed))
    plt.plot([0, max(expected)], [0, max(expected)], color='red', linestyle='--')
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.tight_layout()
    plt.show()


# GWAS A: Naive
subprocess.run([plink_path, "--bfile", "../data/homework", "--maf", "0.005",
                "--linear", "--pheno", "../output/homework.phe", "--out", "../output/gwas_naive"], check=True)

gwas = pd.read_csv("../output/gwas_naive.assoc.linear", sep=r'\s+')
plot_gwas(gwas, "GWAS Results (Raw)")
plot_qq(gwas['P'], "GWAS Results (Raw)")

# GWAS B: Corrected
pcs.to_csv("../output/covariates.txt", sep="\t", index=False)
subprocess.run([plink_path, "--bfile", "../data/homework", "--maf", "0.005",
                "--linear", "--pheno", "../output/homework.phe",
                "--covar", "../output/covariates.txt", "--hide-covar", "--out", "../output/gwas_corrected"], check=True)

gwas = pd.read_csv("../output/gwas_corrected.assoc.linear", sep=r'\s+')
## TODO Q-Q PLOT, MANHATTEN
```


# Analyzing it in STATA + PLINK

- You must manually download plink, and extract it to the software subfolder: https://www.cog-genomics.org/plink/

```
* -- Setup file structure --
* Assuming you are in the 'scripts' folder
capture mkdir "../output"
capture mkdir "../software"

* -- download plink 1.9 --
* Note: Stata is not great at unzipping files across OS types.
* This assumes YOU HAVE DOWNLOADED PLINK to the ../software folder manually
* if the 'shell' command below fails.

* -- QC --

* The allele frequencies - run in plink (might need to do plink.exe in windows)
* Running inside the stata window can look a little wierd at times
shell ../software/plink --bfile ../data/homework --freq --out ../output/stats_freq

* Sorry, I simulated these not quite right. It should be a bit more of a gentle triangle sloping distribution...

* - and then load into Stata
import delimited "../output/stats_freq.frq", delimiters(" ", collapse) varnames(1) clear
list in 1/5

* TODO: Create a histogram of MAFs


* QUESTION: What is the range of the frequency you see in these plots, and why?

* Hardy-Weinberg Equilibrium
shell ../software/plink --bfile ../data/homework --hardy --out ../output/stats_hwe

import delimited "../output/stats_hwe.hwe", delimiters(" ", collapse) varnames(1) clear

* TODO: Create a histogram of the P-values.

* QUESTION: What is the median value?


* QUESTION: Any guesses why this might be happening?

* NOTE: You would usually filtering for Minor Allele Frequency (MAF), we'll do that here.


* -- PCA --
* - Note that MAF is set to 5% here
  shell ../software/plink --bfile ../data/homework --maf 0.05 --indep-pairwise 50 5 0.2 --out ../output/pruning
shell ../software/plink --bfile ../data/homework --extract ../output/pruning.prune.in --pca 10 --out ../output/pca_results

* Load the PCA results (Eigenvectors)
import delimited "../output/pca_results.eigenvec", delimiters(" ", collapse) varnames(1) clear
* Rename columns (PLINK outputs V1=FID, V2=IID, V3=PC1...)
  rename v1 fid
  rename v2 iid
  forval i = 1/10 {
      local j = `i' + 2
      rename v`j' pc`i'
  }

* TODO: Plot at least PC1 vs. PC2 here


* Look at the phenotype we have
* In Stata, we'll save the PCs to merge later or use them as a covariate file
export delimited using "../output/covariates.txt", delimiter(" ") replace

import delimited "../data/homework.fam", delimiters(whitespace, collapse) clear
list in 1/5
keep v1 v2 v6
rename v1 fid
rename v2 iid
rename v6 pheno
list in 1/5

* TODO histogram of phenotype

* TODO: summarize the phenotype

export delimited fid iid pheno using "../output/homework.phe", delimiter(tab) replace


* -- The GWAS Comparison --

* GWAS A: Naive (No covariates)
* In stata, you CANNOT see the progress bar. This takes some time to run. Alternatively run from a command prompt.
shell ../software/plink --bfile ../data/homework --maf 0.005 --linear --pheno ../output/homework.phe --mpheno 1 --out ../output/gwas_naive

import delimited "../output/gwas_naive.assoc.linear", delimiters(whitespace, collapse) varnames(1) clear
* Generate -log10(p) for plotting
gen logp = -log10(p)

* TODO: Manhattan Plot (Naive)
* Stata Manhattan plots are complex, this is a simplified version
sort CHR BP
* simple index across all
gen simplepos = _n
twoway (scatter logp simplepos, msize(tiny)), legend(off)

* TODO: Q-Q Plot
sort p
gen expected_p = _n / (_N + 1)
gen exp_logp = -log10(expected_p)
twoway (scatter logp exp_logp, msize(tiny) mc(navy)) (line exp_logp exp_logp, lcolor(red)), legend(off) aspect(1)


* GWAS B: GWAS adjusted for ancestry
* We use the exported PCs from earlier
shell ../software/plink --bfile ../data/homework --linear --pheno ../output/homework.phe --covar ../output/pca_results.eigenvec --hide-covar --out ../output/gwas_corrected

import delimited "../output/gwas_corrected.assoc.linear", delimiters(whitespace, collapse) varnames(1) clear
gen logp = -log10(p)

* TODO: Manhattan Plot (Corrected)

* TODO: Q-Q Plot
```
