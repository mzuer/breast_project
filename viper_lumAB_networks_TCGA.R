
# Rscript viper_lumAB_networks_TCGA.R

library(TCGAbiolinks)
library(viper)
library('org.Hs.eg.db')

outFolder <- file.path("VIPER_LUMAB_NETWORKS_TCGA")
dir.create(outFolder, recursive=TRUE)

computeNullModel <- FALSE
runMSviper <- TRUE

library(aracne.networks)
data("regulonbrca")

####################################
### first retrieve PAM50 annotation data
####################################
cancer <- "BRCA"
PlatformCancer <- "IlluminaHiSeq_RNASeqV2"
dataType <- "rsem.genes.results"
pathCancer <- "TCGAData/miRNA"
# get subtype information
dataSubt <- TCGAquery_subtype(tumor = cancer)

stopifnot(!duplicated(dataSubt$patient))
samp_annot <- setNames(dataSubt$BRCA_Subtype_PAM50, dataSubt$patient)


####################################
### load ARACne data
####################################

regulfile <- file.path("brca_tcga_rnaseq851_signalomeregulon.rda")
regul_ <- load(regulfile)
# this loads expression data: dset + regulon: regul
regul
# regulon object, 2020 regulators, 18135 targets, 509011 interactions
brca_regul <- regul

expr_cols <- colnames(dset)
expr_cols_short <- substr(expr_cols, 1, 12)
# stopifnot(!duplicated(expr_cols_short)) # not true

samp_cols<- substr(expr_cols, 14, 15)
norm_samps <- paste0(11:30)


all_entrez <- rownames(dset)

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

####################################
### prepare expression data
####################################

expr_dt <- dset[, !samp_cols %in% norm_samps]
stopifnot(!is.na(expr_dt))

# remove samples with dup vials -> keep most var
tmp_Dt <- as.data.frame(t(expr_dt))
tmp_Dt$samp_short <- substr(rownames(tmp_Dt), 1, 12)

dt1 <- tmp_Dt[!tmp_Dt$samp_short %in% tmp_Dt$samp_short[duplicated(tmp_Dt$samp_short)],]
dup_dt <- tmp_Dt[tmp_Dt$samp_short %in% tmp_Dt$samp_short[duplicated(tmp_Dt$samp_short)],]

x=dup_dt[dup_dt$samp_short==dup_dt$samp_short[1],]

dt2 <- do.call(rbind, by(dup_dt, dup_dt$samp_short, function(x) {
  xtmp <- x
  xtmp$samp_short <- NULL
  outtmp <- xtmp[which.max(apply(xtmp, 1, var)),, drop=FALSE]
  stopifnot(nrow(outtmp) == 1)
  outtmp <- as.data.frame(outtmp)
  outtmp$samp_short <- unique(x$samp_short)
  outtmp
}))
stopifnot(colnames(dt1) == colnames(dt2))
expr_dt_nodup <- rbind(dt1,dt2)
stopifnot(!duplicated(expr_dt_nodup$samp_short))
rownames(expr_dt_nodup) <- expr_dt_nodup$samp_short
expr_dt_nodup$samp_short <- NULL
exprDT <- t(expr_dt_nodup)

stopifnot(any(colnames(exprDT) %in% names(samp_annot)))
dim(exprDT)

lumA_samps <- names(samp_annot)[samp_annot=="LumA"]
length(lumA_samps)
# 562
lumA_samps <- lumA_samps[lumA_samps %in% colnames(exprDT)]  
length(lumA_samps)
# 362

lumB_samps <- names(samp_annot)[samp_annot=="LumB"]
length(lumB_samps)
# 209
lumB_samps <- lumB_samps[lumB_samps %in% colnames(exprDT)]  
length(lumB_samps)
# 167

lumA_dt <- exprDT[, colnames(exprDT) %in% lumA_samps]
stopifnot(dim(lumA_dt) > 0)

lumB_dt <- exprDT[, colnames(exprDT) %in% lumB_samps]
stopifnot(dim(lumB_dt) > 0)

####################################
### Master Regulator Analysis performed by msVIPER
####################################

#### Generating the gene expression signatures

# viper:: rowTtest that efficiently performs Student’s t-test for each row of a dataset. 
# can take an ‘ExpressionSet’ object as argument 
# can also take two matrixes as arguments, the first one containing the ‘test’ samples and the second the‘reference’ samples
# by default a 2-tail test.

# here we compare lumA as test and lumB as reference

lumAB_signature <- rowTtest(x=lumA_dt, y=lumB_dt)
# output: List of Student-t-statistic (statistic) and p-values (p.value)
# statistics : # genes x 1 matrix
# pvalue : # genes x 1 matrix


# While we could define the Gene Expression Signature (GES) by using the t-statistic, to be consistent
# with the z-score based null model for msVIPER (see section 6.2), we will estimate z-score values for the GES:
  
lumABz_signature <- (qnorm(lumAB_signature$p.value/2, lower.tail = FALSE) *
                     sign(lumAB_signature$statistic))[, 1]

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
  lumAB_nullmodel <- ttestNull(x=lumA_dt, y=lumB_dt, per = 1000,
                               repos = TRUE, verbose = FALSE)
  outFile <- file.path(outFolder, "lumAB_nullmodel.Rdata")
  save(lumAB_nullmodel, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "lumAB_nullmodel.Rdata")
  lumAB_nullmodel <- get(load(outFile))
}

# we will use the an apropriate cell context-specific regulatory network provided for brca

#### msViper analysis
# we have everything for the master regulatory analysis
if(runMSviper) {
  outFile <- file.path(outFolder, "lumAB_mrs.Rdata")
  lumAB_mrs <- msviper(lumABz_signature, brca_regul, lumAB_nullmodel, verbose = FALSE)
  save(lumAB_mrs, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "lumAB_mrs.Rdata")
  lumAB_mrs <-   lumAB_nullmodel <- get(load(outFile))
}


#### result summary and vizulization
# summary of the results
nSum <- 25
mrs_summary <- summary(lumAB_mrs)
mrs_summary <- mrs_summary[order(mrs_summary$FDR, mrs_summary$p.value, decreasing=FALSE),]
mrs_summary$Regulon_symb <- all_symbs[paste0(mrs_summary$Regulon)]
  
# results vizualization
plot(lumAB_mrs, cex = .7)


#### leading edge analysis
# msVIPER infers the relative activity of a regulatory gene based on the enrichment of its most closely-
#   regulated targets on a given GES, but does not identify which are the target genes enriched in the GES.
# a method called leading-edge analysis to identify the genes driving the
# enrichment of a gene-set on a GES based on Gene Set Enrichment Analysis (GSEA)

ledge_lumAB_mrs <- ledge(lumAB_mrs)
summary(ledge_lumAB_mrs)


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

lumAB_mrshadow <- shadow(lumAB_mrs, regulators = 25, verbose = FALSE)
summary(lumAB_mrshadow)


#### synergy analysis

# to predict synergistic interactions between regulators we first compute the enrichment of co-regulons, defined
# as the intersection between regulons. We expect that a combination of regulators will synergistically regulate
# a gene expression signature if their co-regulon show a significantly higher enrichment on the signature than
# the union of the corresponding regulons

# msviperCombinatorial = Co-regulon analysis; computes the enrichment of all
# co-regulons, generated from a selected number of MRs (indicated by the regulators parameter), on the GES
# as an example, we compute the enrichment of the co-regulons for the top 25 regulators,
lumAB_mrscomb <- msviperCombinatorial(lumAB_mrs, regulators = 25, verbose = FALSE)


# synergy analysis = comparison between the enrichment of the co-regulon versus the union of the corresponding regulons
# implemented by the function msviperSynergy
# input: msviperCombinatorial and the number of permutations used to compute the p-values (default 1000)

lumAB_synmrs <- msviperSynergy(lumAB_mrscomb, verbose = FALSE)
summary(lumAB_synmrs)
plot(lumAB_synmrs, 25, cex = .7)

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
lumAB_vpres <- viper(cbind(lumA_dt, lumB_dt), brca_regul, verbose = FALSE)

dim(lumAB_vpres)
# # regulatory genes x sample

# The differential activity of regulatory proteins between groups of samples, can be obtained
# by any hypothesis testing statistical method, like for example the Student’s t-test:
  

tmp <- rowTtest(x=lumAB_vpres[, colnames(lumA_dt)], y=lumAB_vpres[,colnames(lumB_dt)])
data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2),
          "p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:10], ]

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
lumAB_vpsig <- viperSignature(lumA_dt, lumB_dt, verbose = FALSE)
lumAB_vpres <- viper(lumAB_vpsig, brca_regul, verbose = FALSE)

#### sample clustering
# because VIPER expresses activity for all the regulatory proteins in the same scale – normalized enrich-
#   ment score –, euclidean distance is an appropriate measure of similarity between samples
# we can perform an unsupervised hierarchical cluster analysis of the samples in a similar way we would do
# it in the case of gene expression data
lumAB_vpres_dd <- dist(t(lumAB_vpres), method = "euclidean")
heatmap(as.matrix(lumAB_vpres_dd), Rowv = as.dendrogram(hclust(lumAB_vpres_dd, method = "average")), symm = T)

row_cols <- ifelse(labels(lumAB_vpres_dd) %in% lumA_samps, "darkblue", 
                   ifelse(labels(lumAB_vpres_dd) %in% lumB_samps, "forestgreen",NA))
stopifnot(!is.na(row_cols))
heatmap(as.matrix(lumAB_vpres_dd), 
        Rowv = as.dendrogram(hclust(lumAB_vpres_dd, method = "average")),
        ColSideColors = row_cols,
        symm = T)

### => this shows only lumA samples because I did the comparison with lumA as test over lumB as reference


# We have developed, and included in the viper package, a function to compute the similarity between the
# columns of a gene expression or VIPER-predicted activity matrix. It follows the same concept as the two-tail
# Gene Set Enrichment Analysis (GSEA)[7], but it is based on the aREA algorithm[12]. The viperSimilarity
# function takes an expression or activity matrix as input, and generates a matrix of similarity scores between
# sample pairs, in the form of a ‘similarityDistance’ class object.

lumAB_dd <- viperSimilarity(lumAB_vpres)

# We can use the generic function scale to ‘scale’ the similary matrix in the rage [-1; 1], and the resulting
# matrix will be analogous to a correlation matrix. In this case, identical signatures will produce a similarity
# score equal to 1, while perfectly reversed signatures will produce similarity scores equal to -1. Orthogonal
# signatures will be characterized by similarity scores close to zero. As for other matrixes of similarity, the
# ‘signatureDistance’ class object can be transformed into a ‘distance’ class object with the method as.dist,
# which in turn can be used to perform, for example, cluster analysis of the samples (Fig. 5).

heatmap(as.matrix(as.dist(lumAB_dd)), Rowv = as.dendrogram(hclust(as.dist(lumAB_dd),
   method = "average")), symm = T)

stop("---ok \n")


################################# TRASH

expr_dt <- dset
expr_samps <- colnames(expr_dt)
expr_samps_short <- substr(expr_samps, 1, 12)
sum(annot_dt$cgc_case_id %in% expr_samps_short)



annot_dt <- get(load("../ovarian_project/tcga_data/DOWNLOAD_TCGA_BRCA_LUMALUMB_RECOUNT2/tcga_sampleAnnot.Rdata"))

lumA_samps <- annot_dt$cgc_case_id[annot_dt$PAM50 == "lumA"]
stopifnot(length(lumA_samps) > 0)

lumB_samps <- annot_dt$cgc_case_id[annot_dt$PAM50 == "lumB"]
stopifnot(length(lumB_samps) > 0)


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
