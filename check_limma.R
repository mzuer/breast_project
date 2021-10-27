# df_limma from phenopath code
# compare with my way of setting the design
# -> yield same results, looks ok
set.seed(123)
require(phenopath)
sim <- simulate_phenopath()
suppressPackageStartupMessages(library(SummarizedExperiment))
exprs_mat <- t(sim$y)
pdata <- data.frame(x = sim$x)
sce <- SummarizedExperiment(assays = list(exprs = exprs_mat), 
                            colData = pdata)
cnts <- abs(assay(sce, "exprs"))
labs <- c(rep(1, 50), rep(0, 50))
design <- model.matrix(~labs, pdata)
dge <- DGEList(cnts)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
results <- decideTests(fit)
df_limma <- data.frame(coef = fit$coefficients[,2], 
                       pval = fit$p.value[,2])
df_limma <- df_limma[order(df_limma$pval, decreasing = FALSE),]

design2 <- model.matrix(~ labs)
dge <- DGEList(cnts)
dge <- calcNormFactors(dge)
v <- voom(dge, design2, plot=FALSE)
fit <- lmFit(v, design2)
eb_fit <- eBayes(fit)
DE_topTable <- topTable(eb_fit, coef=ncol(v$design), number=Inf, sort.by="p") ## if not 0+ in design -> coef=2


stopifnot(rownames(DE_topTable) == rownames(df_limma))
stopifnot(DE_topTable$P.Value == df_limma$pval)


# > head(DE_topTable)
# logFC AveExpr     t  P.Value adj.P.Val    B
# 25 -1.74    14.6 -19.0 1.35e-36  1.23e-35 72.7
# 27 -1.74    14.6 -19.0 1.44e-36  1.23e-35 72.6
# 21 -1.74    14.6 -19.0 1.45e-36  1.23e-35 72.6
# 26 -1.74    14.6 -19.0 1.49e-36  1.23e-35 72.6
# 30 -1.74    14.6 -19.0 1.54e-36  1.23e-35 72.6
# 29  1.90    14.6  18.7 4.66e-36  2.00e-35 71.5
# > head(df_limma)
# coef     pval
# 25 -1.74 1.35e-36
# 27 -1.74 1.44e-36
# 21 -1.74 1.45e-36
# 26 -1.74 1.49e-36
# 30 -1.74 1.54e-36
# 29  1.90 4.66e-36
# > 

```{r limma}
# retained_fnames <- rownames(sc_tumour[rowSums(counts(sc_tumour)) > 20, ])
dge <- DGEList(counts(sce))
dge <- calcNormFactors(dge)
design <- model.matrix(~ x, pData(sce))
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
results <- decideTests(fit)
# vennDiagram(results)
```

And merge the results with phenotime:
  
  ```{r merge-with-us}
int_df_limma <- dplyr::rename(df_beta, feature = gene)
qvals <- p.adjust(fit$p.value[,2], method = 'BH')
df_limma <- data_frame(coef = fit$coefficients[,2], 
                       pval = fit$p.value[,2],
                       qval = qvals,
                       mu = pcavi$m_mu,
                       feature = featureNames(sce)) %>% 
  left_join(int_df_limma, by = "feature") %>% 
  dplyr::mutate(log10qval = -log10(qval),
                limma_sig = qval < 0.05) 

