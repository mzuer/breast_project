source("TPQCI/getlv1net.r")
source("TPQCI/calcModule.R")

# Rscript tpqci_cond1_cond2_proteo_Johansson.R

plotType <- "png"
myWidth <- 400
myHeight <- 400

cond1 <- "LumB"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"




outFolder <- file.path("TPQCI_COND1_COND2_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

myHeightGG <- myWidthGG <- 7

library(ggplot2)
library(igraph)
library(foreach)
library(STRINGdb)
string_db <- STRINGdb$new( species=9606)

de_file <- file.path("STRINGDB_COND1_COND2_PROTEO_JOHANSSON/test_LumB_vs_ref_LumA/DE_topTable.Rdata")
de_dt <- get(load(de_file))
de_dt$gene <- rownames(de_dt)
nrow(de_dt)
# 9995
de_dt_tmp <- de_dt
de_dt_tmp$gene_init <- de_dt_tmp$gene
de_dt_mpd <- string_db$map( de_dt_tmp, "gene", removeUnmappedRows = TRUE )
nrow(de_dt_mpd)
# 9957
# remove ambiguity
dup_genes <- de_dt_mpd$gene[duplicated(de_dt_mpd$gene)]
dup_prots <- de_dt_mpd$STRING_id[duplicated(de_dt_mpd$STRING_id)]
de_dt_mpd_nodup <- de_dt_mpd[!de_dt_mpd$gene %in% dup_genes & !de_dt_mpd$STRING_id %in% dup_prots,]
nrow(de_dt_mpd_nodup)
# 9916
stopifnot(!duplicated(de_dt_mpd_nodup$gene))
rownames(de_dt_mpd_nodup) <-  de_dt_mpd_nodup$gene
# the string_db$map() converts lowercase to uppercase
# stopifnot(de_dt_mpd_nodup$gene %in% rownames(proteo_dt)) # not true
stopifnot(de_dt_mpd_nodup$gene_init %in% rownames(proteo_dt))
stopifnot(!duplicated(de_dt_mpd_nodup$gene_init))
rownames(de_dt_mpd_nodup) <- de_dt_mpd_nodup$gene_init

# get the graph
ppi_net <- string_db$get_subnetwork(de_dt_mpd_nodup$STRING_id)

stopifnot(!duplicated(de_dt_mpd_nodup$STRING_id))
stopifnot(!duplicated(de_dt_mpd_nodup$gene))

prot2gene <- setNames(de_dt_mpd_nodup$gene_init, de_dt_mpd_nodup$STRING_id)

stopifnot(V(ppi_net)$name %in% names(prot2gene))
length(V(ppi_net)$name)
# 9911 -> 5 disappear !
nrow(de_dt_mpd_nodup)
# 9916

V(ppi_net)$name <- prot2gene[paste0(V(ppi_net)$name)]

V(ppi_net)$indeg <- degree(ppi_net, mode = "in")
V(ppi_net)[V(ppi_net)$indeg == max(V(ppi_net)$indeg)] # TP53
V(ppi_net)[V(ppi_net)$indeg == 104] # TP53


# plot(induced.subgraph(graph=ppi_net,vids=neighbors(ppi_net, "TP53")))
plot(induced.subgraph(graph=ppi_net,vids=neighbors(ppi_net, "SEMA3F")))
# plot(induced.subgraph(graph=ppi_net,vids=sample(de_dt_mpd_nodup$gene, 1000)))

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

####################################
### select sub DT
####################################

exprDT <- proteo_dt

stopifnot(any(colnames(exprDT) %in% names(samp_annot)))
dim(exprDT)
stopifnot(names(samp_annot) %in% colnames(proteo_dt))

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

# reminder:
# cond1 <- "LumB"
# cond2 <- "LumA"

# emulate a CNV 
# for each gene of each lumB sample -> take the ratio over mean lumA for this sample
cond2_means <- rowMeans(cond2_dt)

stopifnot(names(cond2_means) == rownames(cond1_dt))

exprRatioCond2_dt <- apply(cond1_dt, 2, function(x) x/cond2_means)
stopifnot(exprRatioCond2_dt >= 0)

####################################

#######
#Input#
#######

# for CV -> for the moment, take ratio of gene.sample expression over mean of all other cond

#Cautions: All this part are virtual declarations, need to be replaced to actual data as need

#PPI is a iGraph class object
ppi <- ppi_net
#CNV is a data frame that each row is a gene, each column is a sample
#the first column is gene name with a name "Gene"
cnv <- as.data.frame(exprRatioCond2_dt)
cnv$Gene <- rownames(cnv)
cnv <- cnv[,c("Gene", colnames(exprRatioCond2_dt))]
stopifnot(prot2gene %in% cnv$Gene)
cnv <- cnv[cnv$Gene %in% prot2gene,]
#RNA is also a data frame that each row is a gene, each column is a sample
#the first column is gene name with a name "Gene"
rna <- as.data.frame(cond1_dt)
rna$Gene <- rownames(rna)
rna <- rna[,c("Gene", colnames(cond1_dt))]
rna <- rna[rna$Gene %in% prot2gene,]

stopifnot(V(ppi)$name %in% rownames(rna))
stopifnot(V(ppi)$name %in% rownames(cnv))

rna <- rna[rownames(rna) %in% V(ppi)$name,]
cnv <- cnv[rownames(cnv) %in% V(ppi)$name,]


stopifnot(dim(cnv) == dim(rna))
stopifnot(rownames(cnv) == rownames(rna))
stopifnot(colnames(cnv) == colnames(rna))
stopifnot(setequal(V(ppi)$name, rownames(rna)))

############
#Parameters#
############

#if RNA-seq is logarithmic, set this flag to T
is.RNA.log2ed <- F

#parameters while calculating TPQCI
cnv.mass.threshold<-0.3
pcc.threshold<-0.3

#parameters while module detecting

#Obsoleted parameters
#high.degree.ratio.threshold<-0.5
#lm.min.neighbors<-0

tau.g<-0.15
tau.l<-0.3
beta<-0.5
mu<-20


path.output<-"modules.csv"

#######################
#Run:Calculating TPQCI#
#######################

if(!is.RNA.log2ed)
  rna<-rnaLog2Trans(rna)

cr.cor<-getLV1cor(cnv=cnv, rna=rna, ppn=ppi)
mass<-getCNVMass(cnv,cnv.mass.threshold)

tpqci<-calTopoPot(mass,cr,cor,pcc.threshold)


######################
#Run:Module detection#
######################
tp.weight<-tpqci$topot
names(tp.weight)<-tpqci$node_c

tw.ppi<-getWeightedNetwork(ppi,tp.weight)

#Remove high degree nodes if uncommented (Obsoleted)
#tw.ppi<-removeHighDegNodes(tw.ppi,high.degree.ratio.threshold)

lm.genes<-findLocalMax(tw.ppi)
#switch to this if need limit n.neighbors of local maximal nodes(Obsoleted)
#lm.genes<-findLocalMax(tw.ppi, lm.min.neighbors)

modules<-findModules(tw.ppi,lm.genes,tau.g,tau.l,mu)
modules<-moduleMerge(modules,beta)

#if need to output the modules, uncomment next line
#writeModToCSV(modules, path.output)

