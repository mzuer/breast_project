
# Rscript viper_cond1cond2_networks_proteo_Johansson.R


cond1 <- "HER2"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"

cond1 <- "LumBHer2"
cond2 <- "LumA"
annotCol <- "PAM50.subtype_merged"

cond1 <- "LumB"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"


# how many rows to put in the summary tables
nSum <- 25

plotType <- "png"
myWidth <- 400
myHeight <- 400
myWidthGG <- myHeightGG <- 7

pvalThresh_plot <- 0.05
fcThresh_plot <- 0.1

library(TCGAbiolinks)
library(viper)
library('org.Hs.eg.db')
library(gplots)
library(ggplot2)
library(foreach)
library(igraph)
source("breast_utils.R")
library(CEMiTool)

outFolder <- file.path("VIPER_COND1COND2_NETWORKS_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

computeNullModel <- F
runMSviper <- F
runViper <- F

library(aracne.networks)
data("regulonbrca")

gmt_file <- file.path("c5.go.bp.v7.4.symbols.gmt")

### not enough samples to run ARACne ??!

####################################
### retrieve proteo data
####################################

proteo_dt <- read.delim("data/johansson_data_relative_ratios_to_pool.csv", sep=",")
stopifnot(!duplicated(proteo_dt$gene_symbol))
rownames(proteo_dt) <- proteo_dt$gene_symbol
proteo_dt$gene_symbol <- proteo_dt$gene_symbol <- proteo_dt$ensembl_id <- NULL 

####################################
### first retrieve PAM50 annotation data
####################################


annot_dt <- read.delim("data/johansson_tumor_annot.csv", sep=",")

annot_dt$PAM50.subtype_merged <- annot_dt$PAM50.subtype
annot_dt$PAM50.subtype_merged[annot_dt$PAM50.subtype_merged == "LumB" |
                                annot_dt$PAM50.subtype_merged == "HER2"] <- "LumBHer2"


samp_annot <- setNames(as.character(annot_dt[,paste0(annotCol)]), as.character(annot_dt$Tumor.ID))
samp_annot_all <- setNames(as.character(annot_dt[,paste0("PAM50.subtype")]), as.character(annot_dt$Tumor.ID))
stopifnot(!duplicated(annot_dt$Tumor.ID))

stopifnot(names(samp_annot) %in% colnames(proteo_dt))

####################################
### load ARACne data
####################################

# regulfile <- file.path("brca_tcga_rnaseq851_signalomeregulon.rda")
# regul_ <- load(regulfile)
# # this loads expression data: dset + regulon: regul
# regul

# expr_cols <- colnames(dset)
# expr_cols_short <- substr(expr_cols, 1, 12)
# stopifnot(!duplicated(expr_cols_short)) # not true

# samp_cols<- substr(expr_cols, 14, 15)
# norm_samps <- paste0(11:30)
# all_entrez <- rownames(dset)

# regulon object, 2020 regulators, 18135 targets, 509011 interactions
## the one loaded from aracne.networks has 
# Object of class regulon with 6054 regulators, 19359 targets and 331919 interactions
# 
brca_regul <- regulonbrca
brca_regul_tfs <- names(brca_regul)
brca_regul_targets <- unlist(lapply(brca_regul, function(x) names(x[["tfmode"]])),use.names = FALSE)
all_entrez <- unique(c(brca_regul_tfs, brca_regul_targets))

# for compatibility with previous code
dset <- data.frame(x=1:length(all_entrez),y="", stringsAsFactors = FALSE)
rownames(dset) <- all_entrez

# 
# all_symbs1 <- mapIds(org.Hs.eg.db, all_entrez, 'SYMBOL', 'ENTREZID')
# all_symbs2 <- getSYMBOL(all_entrez, data='org.Hs.eg')
# return same results
# all(na.omit(all_symbs1)==na.omit(all_symbs2))
# TRUE
all_symbs <- mapIds(org.Hs.eg.db, all_entrez, 'SYMBOL', 'ENTREZID')
# 106 duplicated
all_symbs_nona <- na.omit(all_symbs)
stopifnot(!duplicated(all_symbs_nona))

new_rownames <- setNames(all_symbs[rownames(dset)], rownames(dset))
# replace the NA with the entrezID
new_rownames[is.na(new_rownames)] <- names(new_rownames)[is.na(new_rownames)]
stopifnot(!duplicated(new_rownames))
stopifnot(names(new_rownames) == rownames(dset))
rownames(dset) <- new_rownames

stopifnot(names(brca_regul) %in% names(new_rownames))
brca_regul_gs <- brca_regul
names(brca_regul_gs) <- new_rownames[names(brca_regul)]
stopifnot(!is.na(names(brca_regul_gs)))
brca_regul <- lapply(brca_regul_gs, function(x) {
  stopifnot(names(x[["tfmode"]]) %in% names(new_rownames))
  names(x[["tfmode"]]) <- new_rownames[names(x[["tfmode"]])]
  stopifnot(!is.na(names(x[["tfmode"]])))
  x
})

outFile <- file.path(outFolder, "brca_regul.Rdata")
save(brca_regul, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

length(brca_regul)
# 6054
brca_regul_filt <- brca_regul[names(brca_regul) %in% rownames(proteo_dt)]
length(brca_regul_filt)
# 2978
brca_regul_filt2 <- lapply(brca_regul_filt, function(x){
  xtfm <- x[["tfmode"]]
  tok <- which(names(xtfm) %in% rownames(proteo_dt))
  xtfm2 <- xtfm[tok]
  xlkd <- x[["likelihood"]]
  xlkd2 <- xlkd[tok]
  list(
    tfmode=xtfm2, likelihood=xlkd2
  )
} )
# 
length(brca_regul_filt2)
# 2978
brca_regul_filt3 <- Filter(function(x) length(x[["tfmode"]]) > 0, brca_regul_filt2) 
length(brca_regul_filt3)
# 2976

length(unique(unlist(lapply(brca_regul, function(x)names(x[["tfmode"]])))))
#19359
length(unique(unlist(lapply(brca_regul_filt, function(x)names(x[["tfmode"]])))))
# 19188
length(unique(unlist(lapply(brca_regul_filt2, function(x)names(x[["tfmode"]])))))
# 8967
length(unique(unlist(lapply(brca_regul_filt3, function(x)names(x[["tfmode"]])))))
# 8967

####################################
### prepare expression data
####################################


exprDT <- proteo_dt

stopifnot(any(colnames(exprDT) %in% names(samp_annot)))
dim(exprDT)

cond1_samps <- names(samp_annot)[samp_annot==cond1]
stopifnot(length(cond1_samps) > 0)
length(cond1_samps)
# 562
cond1_samps <- cond1_samps[cond1_samps %in% colnames(exprDT)]  
stopifnot(length(cond1_samps) > 0)
length(cond1_samps)
# 362

cond2_samps <- names(samp_annot)[samp_annot==cond2]
stopifnot(length(cond2_samps) > 0)
length(cond2_samps)
# 209
cond2_samps <- cond2_samps[cond2_samps %in% colnames(exprDT)]  
stopifnot(length(cond2_samps) > 0)
length(cond2_samps)
# 167

x <- as.matrix(exprDT)
stopifnot(dim(x) == dim(exprDT))
exprDT <- x
cond1_dt <- exprDT[, colnames(exprDT) %in% cond1_samps]
stopifnot(dim(cond1_dt) > 0)

cond2_dt <- exprDT[, colnames(exprDT) %in% cond2_samps]
stopifnot(dim(cond2_dt) > 0)

cond12_dt <- cbind(cond1_dt, cond2_dt)

### check how many go in the same direction in prot data
all_tfms_dt <- do.call(rbind, lapply(1:length(brca_regul_filt3), function(x) { 
  tfm <- brca_regul_filt3[[x]][["tfmode"]]
  data.frame(g1=names(brca_regul_filt3)[x],
             g2=names(tfm),
             tfm = as.numeric(tfm),
             stringsAsFactors = FALSE)
}))
stopifnot(all_tfms_dt$g1 %in% rownames(cond12_dt))
stopifnot(all_tfms_dt$g2 %in% rownames(cond12_dt))

all_tfms_dt$prot_corr <-foreach(i = 1:nrow(all_tfms_dt), .combine='c') %dopar% {
  g1 <- all_tfms_dt$g1[i]
  g2 <- all_tfms_dt$g2[i]
  if(g1 %in% rownames(cond12_dt) & g2 %in% rownames(cond12_dt)) {
    cor(as.numeric(cond12_dt[paste0(g1),]), as.numeric(cond12_dt[paste0(g2),]))
  } else {
    NA
  }
}
stopifnot(!is.na(all_tfms_dt))

outFile <- file.path(outFolder, paste0("cmp_tfmodeARACNe_proteoCorr.", plotType))
do.call(plotType, list(height=myHeight, width=myWidth))
plot(x=all_tfms_dt$tfm,
     y=all_tfms_dt$prot_corr,
     xlab="ARACNe TFm",
     ylab="proteo corr.")
mtext(side=3, text=paste0("proteo corr: ", cond1, "+", cond2, " data"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

### try with brca_regul3
brca_regul <- brca_regul_filt3

outFile <- file.path(outFolder, "mrs_brca_regul.Rdata")
save(brca_regul, file=outFile)
cat(paste0("... written: ", outFile, "\n"))
# stop("---ok\n")
####################################
### Master Regulator Analysis performed by msVIPER
####################################

#### Generating the gene expression signatures

# viper:: rowTtest that efficiently performs Student’s t-test for each row of a dataset. 
# can take an ‘ExpressionSet’ object as argument 
# can also take two matrixes as arguments, the first one containing the ‘test’ samples and the second the‘reference’ samples
# by default a 2-tail test.

# here we compare lumA as test and lumB as reference

cond12_signature <- rowTtest(x=cond1_dt, y=cond2_dt)
# output: List of Student-t-statistic (statistic) and p-values (p.value)
# statistics : # genes x 1 matrix
# pvalue : # genes x 1 matrix


# While we could define the Gene Expression Signature (GES) by using the t-statistic, to be consistent
# with the z-score based null model for msVIPER (see section 6.2), we will estimate z-score values for the GES:
  
cond12z_signature <- (qnorm(cond12_signature$p.value/2, lower.tail = FALSE) *
                     sign(cond12_signature$statistic))[, 1]

#### Generating the null model

# uniform distribution of the targets on the GES is not a good prior for msVIPER. Given the high degree
# of co-regulation in transcriptional networks, the assumption of statistical independence of gene expression
# s unrealistic an can potentially lead to p-value underestimates. To account for the correlation structure
# between genes, we define a null model for msVIPER by using a set of signatures obtained after permuting
# the samples at random. The function ttestNull performs such process by shuffling the samples among the
# ‘test’ and ‘reference’ sets, according to the re-sampling mode and number of permutations indicated by the
# parameters repos and per, respectively 
# As output, the ttestNull function produces a numerical matrix of z-scores, with genes/probes in rows and
# permutation iterations in columns, than can be used as null model for the msVIPER analysis

# repos	= Logical, whether the permutations should be performed with reposition
# ! WARNING: this takes some time (5-10min)
if(computeNullModel) {
  cat(paste0("... start compute null model msVIPER..."))
  cond12_nullmodel <- ttestNull(x=cond1_dt, y=cond2_dt, per = 1000,
                               repos = TRUE, verbose = FALSE)
  cat(paste0("... done !\n"))
  outFile <- file.path(outFolder, "cond12_nullmodel.Rdata")
  save(cond12_nullmodel, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "cond12_nullmodel.Rdata")
  cond12_nullmodel <- get(load(outFile))
}

# we will use the an apropriate cell context-specific regulatory network provided for brca

#### msViper analysis
# we have everything for the master regulatory analysis
outFile <- file.path(outFolder, "cond12_mrs.Rdata")
if(runMSviper) {
  cat(paste0("... start compute msVIPER..."))
  cond12_mrs <- msviper(cond12z_signature, brca_regul, cond12_nullmodel, verbose = FALSE)
  cat(paste0("... done !\n"))
  save(cond12_mrs, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  cond12_mrs <- get(load(outFile))  
}


#### result summary and vizulization
mrs_summary <- summary(cond12_mrs)
# summary of the results
mrs_summary <- mrs_summary[order(mrs_summary$FDR, mrs_summary$p.value, decreasing=FALSE),]
# mrs_summary$Regulon_symb <- all_symbs[paste0(mrs_summary$Regulon)]
  
# add the targets
stopifnot(mrs_summary$Regulon %in% names(brca_regul))
mrs_tgts <- sapply(mrs_summary$Regulon, function(x) paste0(names(brca_regul[[x]][["tfmode"]]), collapse=","))
names(mrs_tgts) <- NULL
outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_mrs_summary.txt"))
write.table(mrs_summary, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

mrs_summary$targets <- mrs_tgts
outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_mrs_summary_withTgts.txt"))
write.table(mrs_summary, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

# results vizualization
outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_mrs_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(cond12_mrs, cex = .7)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

mrs_summary_all <- summary(cond12_mrs,length(cond12_mrs$regulon))
mrs_summary_all <- mrs_summary_all[order(mrs_summary_all$FDR, mrs_summary_all$p.value, decreasing=FALSE),]

outFile <- file.path(outFolder, "mrs_summary_all.Rdata")
save(mrs_summary_all, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

de_topTable_dt <- get(load("STRINGDB_COND1_COND2_PROTEO_JOHANSSON/test_LumB_vs_ref_LumA/DE_topTable.Rdata"))

#### 
nToPlot <- 10
for(i in 1:nToPlot) {
  i_reg <- mrs_summary_all$Regulon[i]
  stopifnot(i_reg %in% names(brca_regul))
  moi_genes <- names(brca_regul[[paste0(i_reg)]][["tfmode"]])
  moi_genes <- c(i_reg, moi_genes)
  stopifnot(moi_genes %in% de_topTable_dt$gene)
  if(length(moi_genes) >1) {
    stopifnot(moi_genes %in% rownames(cond12_dt))
    # p <- plot_mymodule(module_genes=c(moi_genes), prot_dt=cond12_dt,
    #                    cond1_s=cond1_samps, cond2_s=cond2_samps,
    #                    meanExprThresh=1, absLog2fcThresh=0.5 )
    # use logFC and threshold from the DE data
    p <- plot_mymodule_v2(module_genes=c(moi_genes), 
                          prot_dt=cond12_dt, de_dt=de_topTable_dt,
                       cond1_s=cond1_samps, cond2_s=cond2_samps,
                       pvalThresh=pvalThresh_plot, absLogFCThresh=fcThresh_plot ) + 
    ggtitle(paste0(i_reg, " reg."), subtitle=paste0(length(moi_genes),
                                                    " tgts (before pval<=", pvalThresh_plot, " & absLogFC>=", fcThresh_plot, ")"))
  }
  if(!is.na(p)) {
    outFile <- file.path(outFolder, paste0(i_reg, "_network_fc_corr.", plotType))
    ggsave(p, filename=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
  }
}

#### leading edge analysis
# msVIPER infers the relative activity of a regulatory gene based on the enrichment of its most closely-
#   regulated targets on a given GES, but does not identify which are the target genes enriched in the GES.
# a method called leading-edge analysis to identify the genes driving the
# enrichment of a gene-set on a GES based on Gene Set Enrichment Analysis (GSEA)

ledge_cond12_mrs <- ledge(cond12_mrs)
summary(ledge_cond12_mrs)

ledge_summary <- summary(ledge_cond12_mrs, nSum)
# ledge_summary <- ledge_summary[order(ledge_summary$FDR, ledge_summary$p.value, decreasing=FALSE),]
# ledge_summary$Regulon_symb <- all_symbs[paste0(ledge_summary$Regulon)]

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_ledgeSummary.txt"))
write.table(ledge_summary, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


########################################## ORA

go_fdr_thresh <- 0.1
go_mrs_dt <- mrs_summary_all[mrs_summary_all$FDR <= go_fdr_thresh,]

cat(paste0("... doing GO for FDR <= ", go_fdr_thresh, " modules :\n"))
cat(paste0("... ", nrow(go_mrs_dt), "/", nrow(mrs_summary_all), "\n"))

gmt_in <- read_gmt(gmt_file)
message("Using all genes in GMT file as universe.")
allgenes <- unique(gmt_in[, "gene"])


### mz added: do enrichment for the ms
# getMethod("mod_ora", "CEMiTool")
# CEMiTool:::ora)

cat(paste0("... start computing ORA\n"))

mods <- lapply(go_mrs_dt$Regulon, function(x) names(brca_regul[[x]][["tfmode"]]))
names(mods) <- go_mrs_dt$Regulon

res_list <- lapply(names(mods), CEMiTool:::ora, gmt_in, allgenes, mods)

if (all(lapply(res_list, nrow) == 0)) {
  warning("Enrichment is NULL. Either your gmt file is inadequate or your modules really aren't enriched for any of the pathways in the gmt file.")
  # return(cem)
  ora_res <- NULL
} else {
  names(res_list) <- names(mods)
  ora_res <- lapply(names(res_list), function(x) {
    if (nrow(res_list[[x]]) > 0) {
      as.data.frame(cbind(x, res_list[[x]]))
    }
  })
  ora_res_dt <- do.call(rbind, ora_res)
  names(ora_res_dt)[names(ora_res_dt) == "x"] <- "Module"
  rownames(ora_res_dt) <- NULL
}
go_signifThresh <- 0.05
ora_res_dt_signif <- ora_res_dt[ora_res_dt$p.adjust <= go_signifThresh,]

outFile <- file.path(outFolder, "ora_res_dt_signif.Rdata")
save(ora_res_dt_signif, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("... keeping GO with adj. pval <= ", go_signifThresh, "\n"))
cat(paste0("... ", nrow(ora_res_dt_signif), "/", nrow(ora_res_dt), "\n"))

ora_res_dt_maxSignif <- do.call(rbind, by(ora_res_dt,ora_res_dt$Module, function(x) x[which.min(x$p.adjust),]))
grep("prolif", ora_res_dt_signif$Description)

go_signif_countdt <- data.frame(table(ora_res_dt_signif$ID))
go_signif_countdt <- go_signif_countdt[order(go_signif_countdt$Freq, decreasing = TRUE),]
colnames(go_signif_countdt) <- c("Description", "Freq")

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_goFreq_dt.txt"))
write.table(go_signif_countdt, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


plot_ora_single(head(x, n = n), pv_cut = pv_cut, 
                graph_color = mod_cols[unique(x$Module)], title = unique(x$Module), 
                ...)
  
> CEMiTool:::plot_ora_single
function (es, ordr_by = "p.adjust", max_length = 50, pv_cut = 0.05, 
          graph_color = "#4169E1", title = "Over Representation Analysis") 
{
  comsub <- function(x) {
    d_x <- strsplit(x[c(1, length(x))], "")
    der_com <- match(FALSE, do.call("==", d_x)) - 1
    return(substr(x, 1, der_com + 1))
  }
  es[, "GeneSet"] <- es[, "ID"]
  ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length)
  ovf_data <- es[ovf_rows, "GeneSet"]
  test <- strtrim(ovf_data, max_length)
  dupes <- duplicated(test) | duplicated(test, fromLast = TRUE)
  if (sum(dupes) > 0) {
    test[dupes] <- ovf_data[dupes]
    test[dupes] <- comsub(test[dupes])
    max_length <- max(nchar(test))
  }
  es[ovf_rows, "GeneSet"] <- paste0(strtrim(test, max_length), 
                                    "...")
  es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)
  lvls <- es[order(es[, ordr_by], decreasing = TRUE), "GeneSet"]
  es[, "GeneSet"] <- factor(es[, "GeneSet"], levels = lvls)
  es[, "alpha"] <- 1
  es[es[, ordr_by] > pv_cut, "alpha"] <- 0
  es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
  my_squish <- function(...) {
    return(scales::squish(..., only.finite = FALSE))
  }
  y_axis <- paste("-log10(", ordr_by, ")")
  pl <- ggplot(es, aes_string(x = "GeneSet", y = y_axis, alpha = "alpha", 
                              fill = y_axis)) + geom_bar(stat = "identity") + theme(axis.text = element_text(size = 8), 
                                                                                    legend.title = element_blank()) + coord_flip() + scale_alpha(range = c(0.4, 
                                                                                                                                                           1), guide = "none") + labs(y = "-log10(adjusted p-value)", 
                                                                                                                                                                                      title = title, x = "") + geom_hline(yintercept = -log10(pv_cut), 
                                                                                                                                                                                                                          colour = "grey", linetype = "longdash") + scale_fill_gradient(low = "gray", 
                                                                                                                                                                                                                                                                                        high = graph_color, limits = c(2, 5), oob = my_squish)
  res <- list(pl = pl, numsig = sum(es[, ordr_by] < pv_cut, 
                                    na.rm = TRUE))
  return(res)
}
<bytecode: 0x55a3f7d4f458>
  <enviro




stop("---ok\n")

####################################
### Beyond msVIPER
####################################

#### shadow analysis

# A regulator may appear to be significantly activated based on its regulon’s analysis, simply because several
# of its targets may also be regulated by a bona fide activated TF (shadow effect)
# msVIPER and VIPER (section 8) address this issue by penalizig the contribution
# of the pleotropically regulated targets to the enrichment score. 
# a post-hoc shadow analysis can still be applied to the msVIPER results with the function shadow
# this function takes a class ‘msviper’ object, and performs a shadow analysis on a selected number of top MRs indicated by the
# argument regulators, which can be used to indicate either the enrichment p-value cutoff, the number of top
# MRs, or the names of the MRs to consider in the analysis

cond12_mrshadow <- shadow(cond12_mrs, regulators = 25, verbose = FALSE)
summary(cond12_mrshadow)

shadow_summary <- summary(cond12_mrshadow, nSum)
# shadow_summary <- shadow_summary[order(shadow_summary$FDR, shadow_summary$p.value, decreasing=FALSE),]
# shadow_summary$Regulon_symb <- all_symbs[paste0(shadow_summary$Regulon)]

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_shadowSummary.txt"))
fileConn<-file(outFile)
writeLines(shadow_summary[["Shadow.pairs"]], fileConn)
close(fileConn)
cat(paste0("... written: ", outFile, "\n"))

#### synergy analysis

# to predict synergistic interactions between regulators we first compute the enrichment of co-regulons, defined
# as the intersection between regulons. We expect that a combination of regulators will synergistically regulate
# a gene expression signature if their co-regulon show a significantly higher enrichment on the signature than
# the union of the corresponding regulons

# msviperCombinatorial = Co-regulon analysis; computes the enrichment of all
# co-regulons, generated from a selected number of MRs (indicated by the regulators parameter), on the GES
# as an example, we compute the enrichment of the co-regulons for the top 25 regulators,
cond12_mrscomb <- msviperCombinatorial(cond12_mrs, regulators = 25, verbose = FALSE)


# synergy analysis = comparison between the enrichment of the co-regulon versus the union of the corresponding regulons
# implemented by the function msviperSynergy
# input: msviperCombinatorial and the number of permutations used to compute the p-values (default 1000)

cond12_synmrs <- msviperSynergy(cond12_mrscomb, verbose = FALSE)
summary(cond12_synmrs)

synergy_summary <- summary(cond12_synmrs, nSum)

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_synergySummary.txt"))
write.table(synergy_summary, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_synergy_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(cond12_synmrs, 25, cex = .7)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



####################################
### VIPER
####################################

#### viper analysis

# VIPER is the extension of msVIPER to **single sample-based analysis**. 
# It effectively transforms a gene expression matrix to a regulatory protein activity matrix. 
# The simplest implementation of VIPER is based on single-sample gene expression signatures obtained by
# scaling the probes or genes – subtracting the mean and dividing by the standard devition of each row. 
# A gene expression matrix or ‘ExpressionSet’ object and appropriate regulatory network are the minimum set of parameters
# required to perform a VIPER analysis with the viper function.

# eset should be gene x samples matrix
cond12_vpres <- viper(cbind(cond1_dt, cond2_dt), brca_regul, verbose = FALSE)

dim(cond12_vpres)
# # regulatory genes x sample

# The differential activity of regulatory proteins between groups of samples, can be obtained
# by any hypothesis testing statistical method, like for example the Student’s t-test:
  

tmp <- rowTtest(x=cond12_vpres[, colnames(cond1_dt)], y=cond12_vpres[,colnames(cond2_dt)])

da_dt <- data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2),
          "p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:nSum], ]

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_viper_diffActivity_dt.txt"))
write.table(da_dt, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))



#### null model for VIPER


# VIPER computes the normalized enrichment score (NES) analytically, based on the assumption that in the
# null situation, the target genes are uniformly distributed on the gene expression signature. Because the
# extensive co-regulation of gene expression taking place in the cell, this assumption never holds true, and
# this is the reason why a null model based on sample permutations is used in msVIPER to estimate NES.
# The same approach can also be used for VIPER, given that a set of samples is used as reference for the
# analysis. We can generate a set of GESs based on a set of reference samples, and the corresponding null
# model based on sample permutations, with the function viperSignature.

# the number of permutations for the null model can be defined by the per argument, whose default value is 1,000
# two matrixes as arguments, the first one containing the expression data for all the ‘test’ samples, 
# and the second corresponding to the ‘reference’ samples. 
cond12_vpsig <- viperSignature(cond1_dt, cond2_dt, verbose = FALSE)

if(runViper){
  cat(paste0("... start compute VIPER..."))
  cond12_vpres <- viper(cond12_vpsig, brca_regul, verbose = FALSE)  
  cat(paste0("... done !"))
  outFile <- file.path(outFolder, "cond12_vpres.Rdata")
  save(cond12_vpres, file = outFile)
} else{
  outFile <- file.path(outFolder, "cond12_vpres.Rdata")
  cond12_vpres <- get(load(outFile))
}

#### sample clustering
# because VIPER expresses activity for all the regulatory proteins in the same scale – normalized enrich-
#   ment score –, euclidean distance is an appropriate measure of similarity between samples
# we can perform an unsupervised hierarchical cluster analysis of the samples in a similar way we would do
# it in the case of gene expression data
cond12_vpres_dd <- dist(t(cond12_vpres), method = "euclidean")
# heatmap(as.matrix(cond12_vpres_dd), Rowv = as.dendrogram(hclust(cond12_vpres_dd, method = "average")), symm = T)

all_labs <- labels(cond12_vpres_dd)
stopifnot(all_labs %in% names(samp_annot_all))
all_labs_annot <- samp_annot_all[paste0(all_labs)]
all_labs_cols <- as.numeric(as.factor(all_labs_annot))
all_labs_cols <- ifelse(all_labs_cols == 1, "black", "forestgreen")
stopifnot(!is.na(all_labs_cols))

# all_labs_cols[1:10] ="forestgreen"

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_vprSampleSimilarityHeatmap_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.8, width=myWidth*1.8))
heatmap.2(t(as.matrix(cond12_vpres_dd)), 
          # Rowv = as.dendrogram(hclust(cond12_vpres_dd, method = "average")),
          dendrogram="col",
          trace="none",
          Colv="Rowv",
          key=FALSE,
          ColSideColors = all_labs_cols,
          colCol = all_labs_cols, 
          colRow = all_labs_cols, 
          symm = T)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))






### => this shows only lumA samples because I did the comparison with lumA as test over lumB as reference

#########3 similarity between columns of a gene expression or VIPER-predicted activity matrix
# We have developed, and included in the viper package, a function to compute the similarity between the
# columns of a gene expression or VIPER-predicted activity matrix. It follows the same concept as the two-tail
# Gene Set Enrichment Analysis (GSEA)[7], but it is based on the aREA algorithm[12]. The viperSimilarity
# function takes an expression or activity matrix as input, and generates a matrix of similarity scores between
# sample pairs, in the form of a ‘similarityDistance’ class object.

cond12_dd <- viperSimilarity(cond12_vpres)

# We can use the generic function scale to ‘scale’ the similary matrix in the rage [-1; 1], and the resulting
# matrix will be analogous to a correlation matrix. In this case, identical signatures will produce a similarity
# score equal to 1, while perfectly reversed signatures will produce similarity scores equal to -1. Orthogonal
# signatures will be characterized by similarity scores close to zero. As for other matrixes of similarity, the
# ‘signatureDistance’ class object can be transformed into a ‘distance’ class object with the method as.dist,
# which in turn can be used to perform, for example, cluster analysis of the samples (Fig. 5).

all_labs <- colnames(cond12_dd)
stopifnot(all_labs %in% names(samp_annot_all))
all_labs_annot <- samp_annot_all[paste0(all_labs)]
all_labs_cols <- as.numeric(as.factor(all_labs_annot))
all_labs_cols <- ifelse(all_labs_cols == 1, "black", "forestgreen")
stopifnot(!is.na(all_labs_cols))

# all_labs_cols[1:10] ="forestgreen"

outFile <- file.path(outFolder, paste0(cond1, "_vs_", cond2, "_vprViperSimilarityHeatmap_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.8, width=myWidth*1.8))
heatmap.2(t(as.matrix(cond12_dd)), 
          # Rowv = as.dendrogram(hclust(cond12_vpres_dd, method = "average")),
          dendrogram="col",
          trace="none",
          Colv="Rowv",
          key=FALSE,
          ColSideColors = all_labs_cols,
          colCol = all_labs_cols, 
          colRow = all_labs_cols, 
          symm = T)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



stop("---ok \n")


################################# TRASH

expr_dt <- dset
expr_samps <- colnames(expr_dt)
expr_samps_short <- substr(expr_samps, 1, 12)
sum(annot_dt$cgc_case_id %in% expr_samps_short)



annot_dt <- get(load("../ovarian_project/tcga_data/DOWNLOAD_TCGA_BRCA_LUMALUMB_RECOUNT2/tcga_sampleAnnot.Rdata"))

cond1_samps <- annot_dt$cgc_case_id[annot_dt$PAM50 == "lumA"]
stopifnot(length(cond1_samps) > 0)

cond2_samps <- annot_dt$cgc_case_id[annot_dt$PAM50 == "lumB"]
stopifnot(length(cond2_samps) > 0)


dld_dt <- get(load("../ovarian_project/tcga_data/DOWNLOAD_TCGA_BRCA_LUMALUMB_RECOUNT2/all_counts_onlyPF_0.6.Rdata"))



We have shown that analysis of TF targets inferred by the ARACNe algorithm[1, 13], using the Master
Regulator Inference algorithm (MARINA)[4], is effective in identifying drivers of specific cellular phenotypes
which could be experimentally validated[4, 6]. While VIPER exploits the same principle as MARINA, it
implements a dedicated algorithm specially formulated to estimate regulon activity, which takes into account
the regulator mode of action, the regulator-target gene interaction confidence and the pleiotropic nature of
each target gene regulation. In addition, while especially straightforward for TFs, VIPER effectively extends
to signal transduction proteins. For this, we extended the concept of regulon to include the transcriptional
targets that are most directly affected by the protein’s activity, based on maximization of information transfer
over all alternative paths

VIPER is provided in this package in two flavors: a multiple sample version (msVIPER) designed for
gene expression signatures based in multiple samples or expression profiles, and the single sample version
(VIPER), which estimates relative protein activity on a sample-by-sample basis, thus allowing transformation
of a typical gene expression matrix (i.e. multiple mRNA profiled across multiple samples) into a protein
activity matrix, representing the relative activity of each protein in each sample.



Sequencing data and activity inference RNA-Seq raw gene counts were downloaded from the TCGA firehose web site
(http://gdac.broadinstitute.org, 2016-01-28 release), transformed to Reads Per Kilobase of transcript, per Million mapped reads (RPKM), 
using the average transcript length for each gene and log2 transformed. Transcriptome-wide expression signatures were computed 
by two non-parametric transformations. First, each column (tumor sample) was rank transformed and scaled between 0 and 1. 
Then each row (gene) was rank transformed and scaled between 0 and 1. Finally, regulatory protein activity was measured by 
the VIPER algorithm (Alvarez et al., 2016), using tis- sue-matched ARACNE regulons (Giorgi et al., 2016; Lachmann et al., 2016)
(see Figure S1B). Systematic experimental validation has confirmed that VIPER can accurately measure differential activity 
for > 80% of transcrip-
  tional regulator proteins, when R40% of the genes in a regulon represent bona fide targets of the protein (Alvarez et al., 2016). 
In addition, multiple studies have experimentally validated that > 70% of ARACNe-inferred targets represent bona fide, physical 
tran- scriptional targets—e.g., by Chromatin Immunoprecipitation (ChIP) and RNAi-mediated silencing, followed by gene expression
profiling (Alvarez et al., 2016; Basso et al., 2005; Carro et al., 2010; Lefebvre et al., 2010)—thus fulfilling the VIPER 
requirements for accurate protein measurement. The results of the VIPER analysis are reported as a Normalized Enrichment Scores (NES) values of a protein targets in differentially expressed genes with respect to the centroid of TCGA, as assessed by aREA (see below). This has been shown to accurately characterize differential protein activity. Positive NES values (shown as a red gradient) indicate increased protein activity while negative NES values (shown as a blue gradient) indicate decreased protein activity.
Genomic

I would definitely use the BRCA network from aracne.networks package (bioconductor). The network should be tissue-matched to the gene expression signatures you want to analyze, but does not need to be generated from the same dataset.

If you plan to run ARACNe, I’d run it using all genes annotated as TFs, neve a subset of some “interesting” TFs. The DPI step in ARACNe requires an unbiased and genome-wide representation of regulators in the network.

ARACNe infers direct regulatory relationships between transcriptional regulator proteins and target genes.
Uses gene expression profile data and a predefined list of gene regulators (e.g. TFs)
Reconstructs context-specific transcriptional networks in multiple tissue types.

3 steps
1) MI threshold estimation
2) bootstrapping/MI network reconstruction
3) building consensus network

VIPER outputs a list of genes whose activity levels are highly correlated with the expression of their regulons in the input network.

VIPER analysis entails the construction of networks between TF and target genes using ARACNE (Algorithm for the Reconstruction of 
Accurate Cellular Networks) [33], and then for each TF, a regulon is defined comprising the complete set of predicted downstream 
target genes for that TF. The regulons are then scored for differential enrichment between experimental conditions or on a 
single-sample basis using gene set analysis methods, to infer the biological activity of each TF. A limitation of VIPER is that hundreds
(or thousands) of samples are required for accurate network construction using ARACNE, however, in many applications, it is possible to 
use large data sets available in the public domain derived from the same tissue type to construct a tissue-specific network suitable for 
VIPER.

require(viper)
# download 4.10.21 from http://dx.doi.org/10.6084/m9.figshare.695962
regulfile <- file.path("brca_tcga_rnaseq851_signalomeregulon.rda")
regul_ <- load(regulfile)

expr_dt <- dset
expr_samps <- colnames(expr_dt)
expr_samps_short <- substr(expr_samps, 1, 12)
sum(annot_dt$cgc_case_id %in% expr_samps_short)
# 477/685 available

vpres <- viper(dset, regul, verbose = FALSE) 
# viper output: A matrix of inferred activity for each regulator gene in the network across all samples



The differential activity of regulatory proteins between groups of samples, for example between germinal
center B-cell and Na ̈ıve B-cells, can be obtained by any hypothesis testing statistical method, like for example
the Student’s t-test:
  > tmp <- rowTtest(vpres, "description", c("CB", "CC"), "N")
> data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2),
             + "p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:10], ]
Gene t p.value
TOP2A TOP2A 20.77 2.35e-11
WHSC1 WHSC1 17.44 2.13e-10
MYBL2 MYBL2 16.85 3.25e-10
ZNF274 ZNF274 -16.59 3.96e-10
BCL6 BCL6 16.46 4.36e-10
ZNF23 ZNF23 -16.43 4.47e-10
ZNF32 ZNF32 -16.36 4.70e-10
MORC2 MORC2 16.27 5.03e-10
ZNF101 ZNF101 -16.07 5.87e-10
MYB MYB 15.96 6.40e-10

> signature <- rowTtest(dset, "description", c("CB", "CC"), "N")
It can also take two matrixes as arguments, the first one containing the ‘test’ samples and the second the
‘reference’ samples.
While we could define the Gene Expression Signature (GES) by using the t-statistic, to be consistent
with the z-score based null model for msVIPER (see section 6.2), we will estimate z-score values for the
GES:
  > signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) *
                    + sign(signature$statistic))[, 1]
To account for the correlation structure
between genes, we define a null model for msVIPER by using a set of signatures obtained after permuting
the samples at random. The function ttestNull performs such process by shuffling the samples among the
‘test’ and ‘reference’ sets, according to the re-sampling mode and number of permutations indicated by the
parameters repos and per, respectively.
> nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000,
                         + repos = TRUE, verbose = FALSE)
As output, the ttestNull function produces a numerical matrix of z-scores, with genes/probes in rows and
permutation iterations in columns, than can be used as null model for the msVIPER analysis

he last element required by msVIPER that we are still missing is an apropriate cell context-specific reg-
  ulatory network. We have included a B-cell regulatory network in the bcellViper package, and additional
networks described in [12] for human B-cell, glioma and breast carcinoma can be obtained from figshare
(Table 1
  
  The msVIPER analysis is performed by the msVIPER function. It requires a GES, regulon object and
  null model as arguments, and produces an object of class ‘msVIPER’, containing the GES, regulon and
  estimated enrichment, including the Normalized Enrichment Score (NES) and p-value, as output.
  > mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
  
  The reults can be summarized by the generic function summary, which takes the msviper object and
  either the number of top regulators to report or a specific set of regulators to list. The default for this
  parameter is the top 10 master regulators (MRs).
  > summary(mrs)
  
  A graphics representation of the results (msVIPER plot) can be obtained by the generic function plot
  (shown in Fig. 1). It takes the msviper object and either, the number of top differentially active regulators,
  or the names of the regulators to include in the plot as arguments. The default behavior is to plot the top
  10 most differentially active MRs.
  > plot(mrs, cex = .7)
  
  6.3.1 Leading-edge analysis
  msVIPER infers the relative activity of a regulatory gene based on the enrichment of its most closely-
    regulated targets on a given GES, but does not identify which are the target genes enriched in the GES.
  Subramanian et al. [14] proposed a method called leading-edge analysis to identify the genes driving the
  enrichment of a gene-set on a GES based on Gene Set Enrichment Analysis (GSEA). We implemented the
  leading-edge analysis in the ledge function of the viper package. The function only has a ‘msviper’ class
  object as argument and generates an updated ‘msviper’ object that now includes a ‘ledge’ slot.
  > mrs <- ledge(mrs)
  > summary(mrs)


names(regul) <- sapply(strsplit(names(regul), split = ' - '), head, 1)
mrs <- msviper(ges = mySignature, regulon = viper_regulon, minsize = 4, ges.filter = F, verbose = TRUE)
# Save the results as a dataframe
TF_activities <- data.frame(Regulon = names(mrs$es$nes),
                            Size = mrs$es$size[ names(mrs$es$nes) ], 
                            NES = mrs$es$nes, 
                            p.value = mrs$es$p.value, 
                            FDR = p.adjust(mrs$es$p.value, method = 'fdr'))

Which TFs have the most significant change in activity?
  
  Which genes contributed most to these TF activity changes?
  
  mrs <- ledge(mrs)
# See the results
summary(mrs)

ARACNE shows promise in identifying direct transcriptional interactions in mammalian cellular networks

ARACNE reconstructs the network exactly (asymptotically) if the effect of loops in the network topology is negligible, 
and we show that the algorithm works well in practice, even in the presence of numerous loops and complex topologies
ability to infer validated transcriptional targets of
