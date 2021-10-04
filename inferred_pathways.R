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
regulfile <- file.path("brca_tcga_rnaseq851_signalomeregulon.rda")
regul_ <- load(regulfile)

vpres <- viper(dset, regul, verbose = FALSE)

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
