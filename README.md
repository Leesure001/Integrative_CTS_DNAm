# Integrative Cell-Type-Specific Differential Methylation

This pipeline demonstrates how to perform integrative cell-type-specific differential methylation analysis. The main steps are:
* (a) cell type proportions estimation.
* (b) application of cell-type-specific differential methylation methods. In this pipeline, we include [CellDMC](https://github.com/sjczheng/EpiDISH), [TCA](https://github.com/cozygene/TCA/tree/master), [TOAST](https://github.com/ziyili20/TOAST/tree/master), [CeDAR](https://github.com/ziyili20/TOAST/blob/master/R/cedar.R), [HIRE](https://github.com/XiangyuLuo/HIREewas).
* (c) extracting results from the output of these methods to be combined based on average p-value or min p-value approach described in our paper.

## Install and load required packages
Load required functions.
```r
source("Integrative_CTS_DNAm.R")
```

Install `Cran` packages.
```r
cran_packages <- c("TCA", "devtools", "BiocManager", "assertthat", "matrixStats", "ggVennDiagram", "ggplot2", "ggpubr", "GGally", "gridExtra", "reshape2", "grid")
install_cran_pkgs(pkgs = cran_packages)
```

Install `Bioconductor` packages.
```r
bioc_packages <- c("EpiDISH", "TOAST", "ComplexHeatmap")
install_bioconductor_pkgs(pkgs = bioc_packages)
```

Install `HIREewas`.
```r
devtools::install_github("XiangyuLuo/HIREewas")
```

Run test and load all packages.
```r
check_installed(pkgs = c(bioc_packages, cran_packages, "HIREewas"))
lapply(c(bioc_packages, cran_packages, "HIREewas"), require, character.only = TRUE)
```

## Estimating cell type proportions
Here we demonstrate how to estimate the cell type proportions using [HEpiDISH](https://github.com/sjczheng/EpiDISH). Users can also use other methods to estimate the cell type proportions.
```r
load('sample_450k.RData')
```
This sample data was simulated based on `450K` array benchmarking dataset. It has 10,000 CpGs for 100 samples.

```r
library(EpiDISH)
frac_epidish <- hepidish(beta.m = sample_450k, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m[,c(1, 3, 5)], h.CT.idx = 3, method = 'RPC')
frac_epidish <- frac_epidish[, c("Epi", "Fib", "CD4T", "Mono", "B")]
```

## Load phenotype of interest matrix
```r
load('sample_450k_phe.RData')
```
This simulated phenotype matrix is for 100 samples, with 50 cases and 50 controls.

## Run cell type specific differential methylation analysis methods
### Run CellDMC
Run `CellDMC`.
```r
celldmc_output <- CellDMC(sample_450k, sample_450k_phe, frac_epidish)
```

Check `CellDMC` output.
```r
celldmc_dmct <- celldmc_output$dmct
celldmc_coe <- celldmc_output$coe
```

### Run TCA without refitting cell type porportions
Run `TCA`.
```r
sample_450k_phe_tca <- t(sample_450k_phe)
sample_450k_phe_tca[, 1] <- factor(sample_450k_phe_tca[, 1])
tca_output <- tca(X = sample_450k, W = frac_epidish, C1 = sample_450k_phe_tca)
```

Check `TCA` output.
```r
dmct_tca_estimate <- tca_output$gammas_hat
dmct_tca_pvalue <- tca_output$gammas_hat_pvals
```

### Run TOAST
Run `TOAST`.
```r
design <- data.frame(phe = as.factor(sample_450k_phe))
Design_out <- makeDesign(design, frac_epidish)
fitted_model <- fitModel(Design_out, sample_450k)
toast_output <- csTest(fitted_model, coef = "phe", cell_type = NULL)
```

Check `TOSAT` output.
```r
head(toast_output[[1]])
```

### Run CeDAR
Run `CeDAR`.
```r
design <- data.frame(phe = as.factor(sample_450k_phe))
cedar_output <- cedar(Y_raw = sample_450k, prop = frac_epidish, design.1 = design,
                   factor.to.test = 'phe',cutoff.tree = c('pval',0.01),
                   cutoff.prior.prob = c('pval',0.01))
```

Check `CeDAR` output.
```r
head(cedar_output$tree_res$full$pp)
```

### Run HIRE (initialized using inputting cell type proportions)
Run `HIRE`.
```r
hire_output <- HIRE_V2(sample_450k, sample_450k_phe, num_celltype = ncol(frac_epidish), props_true = t(frac_epidish))
```

Check `HIRE` output.
```r
head(hire_output$beta_t)
head(hire_output$pvalues)
```

## Extract results from each method's output
The original methods' output have different formats, we need to processe them so that they are in the same format to compare and integrate the results from different methods.

### CellDMC
```r
celldmc_coe <- celldmc_output$coe
```
Each **coe** list consists of matrices named after cell types. For each cell type, it will return a matrix, row names are CpGs, column names are `Estimate`, `p` and `adjP`. The default adjusted method is `FDR`. 

```r
dmct.l_celldmc <- get_dmct(coe = celldmc_coe, adjPThresh = 0.05) # Adjust adjPThresh for other FDR control.
dmct_celldmc <- dmct.l_celldmc$dmct
dmct_count_celldmc <- dmct.l_celldmc$dmct_count
```
**dmct.l**: DMCT list, each list consists of 2 matrices, **dmct** and **dmct_count**.
* **dmct** shows the indicators for each cell type:
  * 0 means this CpG is not differentially methylated in that cell type.
  *	1 means this CpG is hypermethylated in that cell type.
  *	-1 means this CpG is hypomethylated in that cell type.
  *	`DMC` column means the result across all cell types, if one CpG is differentially methylated in at least one cell type, then its value will be 1, otherwise, 0.
* The **dmct_count** matrix shows the number of differentially methylated CpGs in each cell type, `dmct_total` is the sum of `dmct_hyper` and `dmct_hypo`.

```r
all(dmct_celldmc == celldmc_output$dmct) # This should be TRUE under FDR 0.05.
celldmc_result <- extract_cpgs(dmct.l_celldmc)
```
Final result will list all significant CpGs, including:
* **dmct** and **dmct_count**, they are the same as DMCT list above.
* A vector **dmc** consists of all differentially methylated CpGs across all cell types.
* Lists for each cell type, each list includes 3 vectors:
  * `_all` shows differentially methylated CpGs in this cell type.
  * `_hyper` shows hypermethylated CpGs in this cell type.
  * `_hypo` shows hypomethylated CpGs in this cell type.

### TCA
```r
tca_coe <- generate_tca_coe(tca_output)
dmct.l_tca <- get_dmct(coe = tca_coe, adjPThresh = 0.05)
dmct_tca <- dmct.l_tca$dmct
dmct_count_tca <- dmct.l_tca$dmct_count
tca_result <- extract_cpgs(dmct.l_tca)
```

### TOAST
```r
toast_coe <- generate_toast_coe(toast_output)
dmct.l_toast <- get_dmct(coe = toast_coe, adjPThresh = 0.05)
dmct_toast <- dmct.l_toast$dmct
dmct_count_toast <- dmct.l_toast$dmct_count
toast_result <- extract_cpgs(dmct.l_toast)
```

### CeDAR
```r
cedar_coe <- generate_cedar_coe(cedar_output)
```
For `CeDAR` **coe** list, we apply a `Bayesian FDR` approach to identify cell-type-specific differentially methylated CpGs, as outlined by [Newton et al](https://pubmed.ncbi.nlm.nih.gov/15054023/). For those CpGs with smaller p values, if the sum of their p values is smaller than FDR control threshold, then we set them as significant CpGs. Therefore, there are no adjusted p values for this method. In the `CeDAR` **coe** list, for the column `adjP`, we set the values for those significant CpGs as 0.01 (this value should be less than FDR control threshold, the default FDR control threshold is 0.05, if you need change it, please modify R function **generate_cedar_coe**), 1 for those insignificant CpGs. For the column `p`, we save `1 – pp` in this column.

```r
dmct.l_cedar <- get_dmct(coe = cedar_coe, adjPThresh = 0.05)
dmct_cedar <- dmct.l_cedar$dmct
dmct_count_cedar <- dmct.l_cedar$dmct_count
cedar_result <- extract_cpgs(dmct.l_cedar)
```

### HIRE
```r
hire_coe <- generate_hire_v2_coe(hire_output)
dmct.l_hire <- get_dmct(coe = hire_coe, adjPThresh = 0.05)
dmct_hire <- dmct.l_hire$dmct
dmct_count_hire <- dmct.l_hire$dmct_count
hire_result <- extract_cpgs(dmct.l_hire)
```
**Note:** `HIRE` requires a large amount of memory for large sample sizes with full EPIC array (e.g., sample size over 200). Users should ensure they have sufficient memory when running `HIRE` in such cases.

## Integrative analysis with average p-value method
### Combine 5 methods: CellDMC, TCA, TOAST, CeDAR, HIRE
```r
ave_pv_5_coe <- generate_overall_coe_v1(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 5, adjPMethod = "fdr")
```
For the **coe** lists of combined methods, the column `Estimate` is copied the coefficients from `CellDMC` model, it wouldn’t be used in the analysis.

```r
dmct.l_ave_pv_5 <- get_dmct(coe = ave_pv_5_coe, adjPThresh = 0.05)
dmct_ave_pv_5 <- dmct.l_ave_pv_5$dmct
dmct_count_ave_pv_5 <- dmct.l_ave_pv_5$dmct_count
ave_pv_5 <- extract_cpgs(dmct.l_ave_pv_5)
```

### Combine 4 methods: CellDMC, TCA, TOAST, CeDAR
```r
ave_pv_4_coe <- generate_overall_coe_v1(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 4, adjPMethod = "fdr")
dmct.l_ave_pv_4 <- get_dmct(coe = ave_pv_4_coe, adjPThresh = 0.05)
dmct_ave_pv_4 <- dmct.l_ave_pv_4$dmct
dmct_count_ave_pv_4 <- dmct.l_ave_pv_4$dmct_count
ave_pv_4 <- extract_cpgs(dmct.l_ave_pv_4)
```

### Combine 3 methods: CellDMC, TCA, TOAST
```r
ave_pv_31_coe <- generate_overall_coe_v1(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 31, adjPMethod = "fdr")
dmct.l_ave_pv_31 <- get_dmct(coe = ave_pv_31_coe, adjPThresh = 0.05)
dmct_ave_pv_31 <- dmct.l_ave_pv_31$dmct
dmct_count_ave_pv_31 <- dmct.l_ave_pv_31$dmct_count
ave_pv_31 <- extract_cpgs(dmct.l_ave_pv_31)
```

### Combine 3 methods: CellDMC, TCA, CeDAR
```r
ave_pv_32_coe <- generate_overall_coe_v1(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 32, adjPMethod = "fdr")
dmct.l_ave_pv_32 <- get_dmct(coe = ave_pv_32_coe, adjPThresh = 0.05)
dmct_ave_pv_32 <- dmct.l_ave_pv_32$dmct
dmct_count_ave_pv_32 <- dmct.l_ave_pv_32$dmct_count
ave_pv_32 <- extract_cpgs(dmct.l_ave_pv_32)
```

### Combine 2 methods: CellDMC, TCA
```r
ave_pv_2_coe <- generate_overall_coe_v1(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 2, adjPMethod = "fdr")
dmct.l_ave_pv_2 <- get_dmct(coe = ave_pv_2_coe, adjPThresh = 0.05)
dmct_ave_pv_2 <- dmct.l_ave_pv_2$dmct
dmct_count_ave_pv_2 <- dmct.l_ave_pv_2$dmct_count
ave_pv_2 <- extract_cpgs(dmct.l_ave_pv_2)
```

## Integrative analysis with minimum p-value method
### Combine 5 methods: CellDMC, TCA, TOAST, CeDAR, HIRE
```r
min_pv_5_coe <- generate_overall_coe_v2(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 5, adjPMethod = "fdr")
dmct.l_min_pv_5 <- get_dmct(coe = min_pv_5_coe, adjPThresh = 0.05)
dmct_min_pv_5 <- dmct.l_min_pv_5$dmct
dmct_count_min_pv_5 <- dmct.l_min_pv_5$dmct_count
min_pv_5 <- extract_cpgs(dmct.l_min_pv_5)
```

### Combine 4 methods: CellDMC, TCA, TOAST, CeDAR
```r
min_pv_4_coe <- generate_overall_coe_v2(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 4, adjPMethod = "fdr")
dmct.l_min_pv_4 <- get_dmct(coe = min_pv_4_coe, adjPThresh = 0.05)
dmct_min_pv_4 <- dmct.l_min_pv_4$dmct
dmct_count_min_pv_4 <- dmct.l_min_pv_4$dmct_count
min_pv_4 <- extract_cpgs(dmct.l_min_pv_4)
```

### Combine 3 methods: CellDMC, TCA, TOAST
```r
min_pv_31_coe <- generate_overall_coe_v2(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 31, adjPMethod = "fdr")
dmct.l_min_pv_31 <- get_dmct(coe = min_pv_31_coe, adjPThresh = 0.05)
dmct_min_pv_31 <- dmct.l_min_pv_31$dmct
dmct_count_min_pv_31 <- dmct.l_min_pv_31$dmct_count
min_pv_31 <- extract_cpgs(dmct.l_min_pv_31)
```

### Combine 3 methods: CellDMC, TCA, CeDAR
```r
min_pv_32_coe <- generate_overall_coe_v2(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 32, adjPMethod = "fdr")
dmct.l_min_pv_32 <- get_dmct(coe = min_pv_32_coe, adjPThresh = 0.05)
dmct_min_pv_32 <- dmct.l_min_pv_32$dmct
dmct_count_min_pv_32 <- dmct.l_min_pv_32$dmct_count
min_pv_32 <- extract_cpgs(dmct.l_min_pv_32)
```

### Combine 2 methods: CellDMC, TCA
```r
min_pv_2_coe <- generate_overall_coe_v2(celldmc_coe, tca_coe, toast_coe, cedar_coe, hire_coe, select = 2, adjPMethod = "fdr")
dmct.l_min_pv_2 <- get_dmct(coe = min_pv_2_coe, adjPThresh = 0.05)
dmct_min_pv_2 <- dmct.l_min_pv_2$dmct
dmct_count_min_pv_2 <- dmct.l_min_pv_2$dmct_count
min_pv_2 <- extract_cpgs(dmct.l_min_pv_2)
```

## Generate VennDiagram to compare the results
```r
library(ggVennDiagram)
library(ggplot2)

my_list <- list(CellDMC = celldmc_result$epi$epi_all, TCA = tca_result$epi$epi_all, HIRE = hire_result$epi$epi_all, TOAST = toast_result$epi$epi_all, CeDAR = cedar_result$epi$epi_all)
ggvenn_plot1 <- ggVennDiagram(
  my_list,
  label = "count"
) + 
scale_fill_gradient(low = "white", high = "dodgerblue") +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
ggtitle("Epi")

my_list <- list(CellDMC = celldmc_result$fib$fib_all, TCA = tca_result$fib$fib_all, HIRE = hire_result$fib$fib_all, TOAST = toast_result$fib$fib_all, CeDAR = cedar_result$fib$fib_all)
ggvenn_plot2 <- ggVennDiagram(
  my_list,
  label = "count"
) + 
scale_fill_gradient(low = "white", high = "dodgerblue") +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
ggtitle("Fib")

my_list <- list(CellDMC = celldmc_result$cd4t$cd4t_all, TCA = tca_result$cd4t$cd4t_all, HIRE = hire_result$cd4t$cd4t_all, TOAST = toast_result$cd4t$cd4t_all, CeDAR = cedar_result$cd4t$cd4t_all)
ggvenn_plot3 <- ggVennDiagram(
  my_list,
  label = "count"
) + 
scale_fill_gradient(low = "white", high = "dodgerblue") +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
ggtitle("CD4T")

my_list <- list(CellDMC = celldmc_result$mono$mono_all, TCA = tca_result$mono$mono_all, HIRE = hire_result$mono$mono_all, TOAST = toast_result$mono$mono_all, CeDAR = cedar_result$mono$mono_all)
ggvenn_plot4 <- ggVennDiagram(
  my_list,
  label = "count"
) + 
scale_fill_gradient(low = "white", high = "dodgerblue") +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
ggtitle("Mono")

my_list <- list(CellDMC = celldmc_result$b$b_all, TCA = tca_result$b$b_all, HIRE = hire_result$b$b_all, TOAST = toast_result$b$b_all, CeDAR = cedar_result$b$b_all)
ggvenn_plot5 <- ggVennDiagram(
  my_list,
  label = "count"
) + 
scale_fill_gradient(low = "white", high = "dodgerblue") +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
ggtitle("B")

library(ggpubr)
p1 <- ggarrange(ggvenn_plot1, ggvenn_plot2, ggvenn_plot3, ggvenn_plot4, ggvenn_plot5, ncol = 3, nrow = 2)
p1
```

## Generate UpSet plot to compare the results
```r
library(GGally)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(ggVennDiagram)
library(ComplexHeatmap)
library(grid)

fdr1 <- 0.05
my_list <- getList(celldmc_coe, tca_coe, hire_coe, toast_coe, cedar_coe, FDR = fdr1)
plot.all <- list()
for(i in seq_along(my_list)){
    m1 <- make_comb_mat(my_list[[i]])
    plot.all[[i]] <- grid.grabExpr(draw(UpSet(m1, column_title = names(my_list)[i], bg_col = "#F0F0FF", bg_pt_col = "#CCCCFF", set_order = 1:length(my_list[[i]]))))
}

pdf(paste0('sample_450k_UpSet.pdf'), width = 12, height = 8)
print(ggarrange(plot.all[[1]], plot.all[[2]], plot.all[[3]], plot.all[[4]], plot.all[[5]], ncol = 3, nrow = 2, labels = LETTERS[1:5]))
dev.off()
```

## References
If you use this pipeline, please cite:

S. Li and P. Kuan (2025). A systematic evaluation of cell-type-specific differential methylation analysis in bulk tissue (manuscript under revision) and the corresponding cell type specific algorithms used.
