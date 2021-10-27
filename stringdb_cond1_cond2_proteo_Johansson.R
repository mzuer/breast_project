
# Rscript stringdb_cond1_cond2_proteo_Johansson.R

cond1 <- "HER2"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"

cond1 <- "LumBHer2"
cond2 <- "LumA"
annotCol <- "PAM50.subtype_merged"

plotType <- "png"
myWidth <- 400
myHeight <- 400
myHeightGG <- myWidthGG <- 7

cond1 <- "LumB"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"

nTop <- 200
signifThresh <- 0.05

outFolder <- file.path("STRINGDB_COND1_COND2_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

# gmt_file <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_file <- "c5.go.bp.v7.4.symbols.gmt"

nTop_connect <- 5


library(STRINGdb)
library(ggplot2)
library(igraph)
library(foreach)
source("breast_utils.R")

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

### DE ANALYSIS

require(edgeR)
require(limma)

cond12_dt <- cond12_dt[,c(colnames(cond1_dt), colnames(cond2_dt))]

sub_labs <- factor(c(rep(cond1, ncol(cond1_dt)),rep(cond2, ncol(cond2_dt))), levels=c(cond1, cond2))
design <- model.matrix(~ sub_labs)



v <- voom(cond12_dt, design, plot=FALSE)
fit <- lmFit(v, design)
eb_fit <- eBayes(fit)
DE_topTable <- topTable(eb_fit, coef=ncol(v$design), number=Inf, sort.by="p") ## if not 0+ in design -> coef=2

dge <- DGEList(cond12_dt)
dge <- calcNormFactors(dge)
v2 <- voom(dge, design, plot=FALSE)
fit2 <- lmFit(v2, design)
eb_fit2 <- eBayes(fit2)
DE_topTable2 <- topTable(eb_fit2, coef=ncol(v2$design), number=Inf, sort.by="p") ## if not 0+ in design -> coef=2


# geneName = "COL14A1"
# plot_dt <- data.frame(
#   cond=c(rep(cond1, ncol(cond1_dt)),rep(cond2, ncol(cond2_dt))),
#   expr = c(cond1_dt[geneName,], cond2_dt[geneName,]),
#   stringsAsFactors = FALSE
# )
# require(ggpubr)
# ggboxplot(data=plot_dt, x="cond", y="expr")

### initialization
# WARNING: You didn't specify a species. Hence we will set 9606 (Homo Sapiens) as your species.
# WARNING: Score threshold is not specified. We will be using medium stringency cut-off of 400.
# WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING  11.0 
string_db <- STRINGdb$new( species=9606)


#STRINGdb$methods()              # To list all the methods available.
#STRINGdb$help("get_graph")      # To visualize their documentation.

##### map stringdb IDs
DE_topTable <- as.data.frame(DE_topTable)
DE_topTable$gene <- rownames(DE_topTable)

outFile <- file.path(outFolder, "DE_topTable.Rdata")
save(DE_topTable, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

DE_topTable2 <- as.data.frame(DE_topTable2)
DE_topTable2$gene <- rownames(DE_topTable2)

outFile <- file.path(outFolder, "DE_topTable_withCNF.Rdata")
save(DE_topTable2, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


m_dt <- merge(DE_topTable, DE_topTable2, by="gene")

plot(m_dt$adj.P.Val.x~m_dt$adj.P.Val.y)
curve(x*1, add=T, col="red")
cor(m_dt$adj.P.Val.x,m_dt$adj.P.Val.y) # 0.992
m_dt$gene[order(m_dt$adj.P.Val.x)] == m_dt$gene[order(m_dt$adj.P.Val.y)]
plot(rank(m_dt$adj.P.Val.x)~rank(m_dt$adj.P.Val.y))
plot(m_dt$logFC.x~m_dt$logFC.y)
curve(x*1, add=T, col="red")
cor(m_dt$logFC.y,m_dt$logFC.x) # 0.99

top50_a <- DE_topTable$gene[1:50]
top50_b <- DE_topTable2$gene[1:50]
length(intersect(top50_a, top50_b))
# 46
stop("--ok\n")

###################################################
DE_topTable_mpd <- string_db$map( DE_topTable, "gene", removeUnmappedRows = TRUE )
# this adds a column with STRING_id
# 0% unmapped

sdb_genes <- string_db$mp(DE_topTable_mpd$STRING_id)

outFile <- file.path(outFolder, "DE_topTable_mpd.Rdata")
save(DE_topTable_mpd, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


##### hits
# extract the most significant 200 genes and we produce an image of the STRING network
# for those
hits <- DE_topTable_mpd$STRING_id[1:nTop]
signifs <- DE_topTable_mpd$STRING_id[DE_topTable_mpd$adj.P.Val <= signifThresh]

outFile <- file.path(outFolder, paste0("DE_top", nTop, "_prot_ppi.", plotType))
string_db$plot_network( hits )
mtext(side=3, text=paste0("prot. DE nTop=", nTop, " (", cond1, " vs. ", cond2, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))  

### add_diff_exp_color
# filter by p-value and add a color column 
# (i.e. green down-regulated gened and red for up-regulated genes)
DE_topTable_mpd_pval05 <- string_db$add_diff_exp_color( subset(DE_topTable_mpd, adj.P.Val<= signifThresh), 
                                                            logFcColStr="logFC" )    

# post payload information to the STRING server
payload_id <- string_db$post_payload( DE_topTable_mpd_pval05$STRING_id, 
                                        colors=DE_topTable_mpd_pval05$color )

##### plot_halo_network
# display a STRING network png with the "halo"

outFile <- file.path(outFolder, paste0("DE_signif", signifThresh, "_prot_ppi_with_halo.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
string_db$plot_network( signifs, payload_id=payload_id )
mtext(side=3, text=paste0("prot. DE signif. <= ", signifThresh, " (", cond1, " vs. ", cond2, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))  


### enrichment analysis on selected genes
enrichment <- string_db$get_enrichment( hits )
head(enrichment, n=20)

### to get interactions between list of proteins
# string_db$get_interactions( hits )

### annotation without the enrichment
annotations <- string_db$get_annotations( hits )
head(annotations, n=20)


###################################################
### code chunk number 15: clustering1
###################################################
# get clusters
clustersList <- string_db$get_clusters(DE_topTable_mpd_pval05$STRING_id)

cat(paste0("... ", nrow(DE_topTable_mpd_pval05), " proteins (signif. only)\n"))
cat(paste0("... clustered in ", length(clustersList), " clusters\n"))

stopifnot(!duplicated( unlist(clustersList)))
cl_l <- unlist(lapply(clustersList, length))

outFile <- file.path(outFolder, paste0("nbrGenes_byModules_clustersSignif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cl_l))
mtext(side=3, text=paste0("# modules = ", length(clustersList)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

### code chunk number 15: clustering on full network

clustersList_full <- string_db$get_clusters(DE_topTable_mpd$STRING_id)

cat(paste0("... ", nrow(DE_topTable_mpd), " proteins (all)\n"))
cat(paste0("... clustered in ", length(clustersList_full), " clusters\n"))

stopifnot(!duplicated( unlist(clustersList_full)))
cl_l_full <- unlist(lapply(clustersList_full, length))

outFile <- file.path(outFolder, paste0("nbrGenes_byModules_clustersAll.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(cl_l_full))
mtext(side=3, text=paste0("# modules = ", length(clustersList_full)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



###################################################
### code chunk number 17: clustering2
###################################################

clustersList_to_plot <- clustersList[which(cl_l > 2)]

for(i in 1:length(clustersList_to_plot)){
  outFile <- file.path(outFolder, paste0("ppi_plot_cluster", i, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  string_db$plot_network(clustersList_to_plot[[i]], payload_id=payload_id )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
}


### code chunk number 18: proteins
# string_proteins <- string_db$get_proteins()


### code chunk number 19: atmtp
# tp53 = string_db$mp( "tp53" )
# atm = string_db$mp( "atm" )


### code chunk number 20: neighbors (eval = FALSE)
## string_db$get_neighbors( c(tp53, atm) )


### code chunk number 21: interactions
# string_db$get_interactions( c(tp53, atm) )


###################################################
### code chunk number 22: paralogs (eval = FALSE)
###################################################
## # Get all homologs of TP53 in human.
## string_db$get_paralogs(tp53)


###################################################
### code chunk number 23: Closest homologs from other species (eval = FALSE)
###################################################
## # get the best hits of the following protein in all the STRING species
## string_db$get_homologs_besthits(tp53)


###################################################
### code chunk number 24: homologs_besthits in target species (eval = FALSE)
###################################################
## # get the homologs of the following two proteins in the mouse (i.e. species_id=10090)
## string_db$get_homologs_besthits(c(tp53, atm), target_species_id=10090, bitscore_threshold=60)

####################### TRASH

# res_de_1a <- topTable(eb_fit,  adjust.method="BH", coef=ncol(v$design), number=Inf, sort.by="p")
# res_de_1b <- topTable(eb_fit,  adjust.method="BH", number=Inf, sort.by="p")
# res_de_1c <- topTable(eb_fit,  adjust.method="BH", coef=1, number=Inf, sort.by="p")
# 
# sub_labs <- factor(c(rep(cond1, ncol(cond1_dt)),rep(cond2, ncol(cond2_dt))), levels=c(cond1, cond2))
# design <- model.matrix(~ 0+sub_labs)
# v <- voom(cond12_dt, design, plot=FALSE)
# fit <- lmFit(v, design)
# eb_fit <- eBayes(fit)
# res_de_2a <- topTable(eb_fit,  adjust.method="BH", coef=ncol(v$design), number=Inf, sort.by="p")
# res_de_2b <- topTable(eb_fit,  adjust.method="BH", number=Inf)
# res_de_2c <- topTable(eb_fit,  adjust.method="BH", coef=1, number=Inf)
# 
# 
# res_de_1b = res_de_1a
# res_de_1c similar but not equal res_de_2a
# res_de_1c similar = res_de_2b (with 2 cols)
# res_de_1c = res_de_2c
