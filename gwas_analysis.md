1. Load R Libraries
```{r}
library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)
```

2. Setup file structure
```{r}
try(dir.create("../output"))
try(dir.create("../software"))
```

3. Download PLINK 1.9
```{r}
if (!file.exists("../software/plink")) {
  message("[*] Downloading PLINK 1.9...")
  curwd <- getwd()
  dir.create("../software", showWarnings = FALSE)
  setwd("../software")

  download.file(
    "https://s3.amazonaws.com/plink1-assets/plink_mac_20231211.zip",
    "plink.zip"
  )
  unzip("plink.zip")

  system("chmod u+x plink")
  system("ls -lh")

  setwd(curwd)
}
```

4. Allele Frequencies
```{r}
system("../software/plink --bfile ../data/homework --freq --out ../output/stats_freq")
```

5. Load into R
```{r}
freq <- fread("../output/stats_freq.frq")
head(freq)
```

TODO: Create a histogram of MAFs
```{r}
ggplot(freq, aes(x = MAF)) +
  geom_histogram(bins = 50) +
  labs(title = "MAF Distribution", x = "Minor Allele Frequency", y = "Count")
```

QUESTION: What is the range of the frequency you see in these plots, and why?

ANSWER: The minor allele frequencies range approximately from 0.05 to 0.35, with most variants centered around ~0.18–0.20. Unlike real genetic data, which is typically skewed toward rare variants, this distribution appears bell-shaped. This pattern is likely due to the simulated nature of the dataset, where allele frequencies were generated from a specific distribution rather than reflecting natural population genetics processes.

6. Delete the dataset to save memory (optional)
```{r}
rm(freq)
```

7. Hardy-Weinberg Equilibrium
```{r}
system("../software/plink --bfile ../data/homework --hardy --out ../output/stats_hwe")
hwe <- fread("../output/stats_hwe.hwe")
```

# TODO: Create a histogram of the P-Values
```{r}
ggplot(hwe, aes(x = P)) +
  geom_histogram(bins = 50) +
  labs(title = "HWE P-value Distribution")

median(hwe$P, na.rm = TRUE)
```

QUESTION: What is the median P-value?
ANSWER: The median Hardy-Weinberg P-value is approximately 5.017e-13

QUESTION: Any guesses why this might be happening?
ANSWER: The extremely small median P-value suggests widespread deviation from Hardy–Weinberg equilibrium across SNPs. In real data, this might indicate genotyping errors, population stratification, inbreeding, or selection. However, since this dataset is simulated, it is likely that the genotype data were generated in a way that does not strictly follow Hardy–Weinberg assumptions, leading to systematic deviation.

8. Linkage Disequilibrium Pruning
Because nearby SNPs on the same chromosome tend to be inherited together (linkage disequilibrium), we prune correlated variants prior to PCA to ensure that the principal components reflect genome-wide ancestry structure rather than local chromosomal linkage effects.
```{r}
system("../software/plink --bfile ../data/homework --maf 0.05 --indep-pairwise 50 5 0.2 --out ../output/pruning")
```

9. PCA on Pruned SNPs
```{r}
system("../software/plink --bfile ../data/homework --extract ../output/pruning.prune.in --pca 10 --out ../output/pca_results")
```

10. Load PCA Eigenvectors
```{r}
pcs <- read.table("../output/pca_results.eigenvec", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))
```

TODO: Plot at least PC1 vs. PC2
The PCA plot shows three clearly separated clusters along PC1 and PC2, indicating strong population structure within the simulated dataset. This suggests the presence of multiple genetically distinct subgroups. Such structure can confound association analyses if not properly adjusted for, which justifies including principal components as covariates in the GWAS.
```{r}
library(ggplot2)

ggplot(pcs, aes(PC1, PC2)) +
  geom_point(alpha = 0.6, size = 1) +
  labs(title = "PCA: PC1 vs PC2") +
  theme_classic()

# How much variance does each PC explain
eigenvals <- scan("../output/pca_results.eigenval")
eigenvals / sum(eigenvals)
```

11. Load & Prepare Phenotype
```{r}
phe <- fread("../data/homework.fam", header = FALSE)
head(phe)

# Keep FID, IID, phenotype column (column 6)
phe <- phe[, .(FID = V1, IID = V2, Pheno = V6)]
```

TODO: Histogram of Phenotype
The phenotype appears approximately normally distributed, supporting the use of linear regression for the GWAS analysis.
```{r}
ggplot(phe, aes(x = Pheno)) +
  geom_histogram(bins = 50) +
  labs(title = "Phenotype Distribution",
       x = "Phenotype",
       y = "Count") +
  theme_minimal()
```

TODO: Summarize the Phenotype
The phenotype is continuous, with a mean of approximately 3.4 and a range from -51.5 to 74.3. The distribution appears roughly symmetric, supporting the use of linear regression in the GWAS analysis.
```{r}
summary(phe$Pheno)
```

12. Write Phenotype File for PLINK
```{r}
fwrite(phe, file="../output/homework.phe", quote=FALSE, sep="\t")
```

13. GWAS A: Naive (No Covariates)
```{r}
system("../software/plink --bfile ../data/homework --maf 0.005 --linear --pheno ../output/homework.phe --mpheno 1 --out ../output/gwas_naive")
```

14. Load Naive Results + Sanity Check
```{r}
gwas_res <- fread("../output/gwas_naive.assoc.linear", header = TRUE)
table(gwas_res$CHR)
```

15. Manhattan + QQ Plot for Naive
```{r}
manhattan(gwas_res,
          chr="CHR", bp="BP", p="P", snp="SNP",
          main="GWAS Results (raw)",
          col=c("orange", "blue"),
          suggestiveline = FALSE,
          genomewideline = -log10(5e-8),
          cex = 0.6)

qq(gwas_res$P, main="Q-Q: Naive")
```

16. Write Covariates for PLINK
```{r}
covar <- pcs[, c("FID", "IID", paste0("PC", 1:5))]
fwrite(covar, "../output/covariates.txt", row.names = FALSE, quote = FALSE, sep = "\t")
```

17. GWAS B: GWAS Adjusted for Ancestry (Covariates)
```{r}
system("../software/plink --bfile ../data/homework --linear --pheno ../output/homework.phe --mpheno 1 --covar ../output/covariates.txt --hide-covar --out ../output/gwas_corrected")
```

18. Load Adjusted Results + Sanity Check
```{r}
gwas_res <- fread("../output/gwas_corrected.assoc.linear")
table(gwas_res$CHR)
```

TODO: Manhattan + QQ Plot for Adjusted
The naive GWAS showed strong inflation, meaning many SNPs appeared significant even when they likely were not. The Q–Q plot clearly deviated from what we would expect under the null hypothesis, suggesting widespread false positives caused by population structure. After adjusting for ancestry using principal components, the Q–Q plot moved much closer to the expected line, indicating that much of the inflation was due to confounding rather than true genetic effects. Although some strong signals remained, the overall number of false positives was reduced, demonstrating the importance of correcting for population structure in GWAS.
```{r}
manhattan(gwas_res,
          chr="CHR", bp="BP", p="P", snp="SNP",
          main="GWAS Results (PC-adjusted)",
          col=c("orange", "blue"),
          suggestiveline = FALSE,
          genomewideline = -log10(5e-8),
          cex = 0.6)

qq(gwas_res$P, main="Q-Q: Corrected")
```
