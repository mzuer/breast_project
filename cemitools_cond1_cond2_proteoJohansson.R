
# Rscript cemitools_cond1_cond2_proteo_Johansson.R

cond1 <- "HER2"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"

cond1 <- "LumBHer2"
cond2 <- "LumA"
annotCol <- "PAM50.subtype_merged"

plotType <- "png"
myWidth <- 400
myHeight <- 400
myWidthGG <- myHeightGG <- 7

cond1 <- "LumB"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"

outFolder <- file.path("CEMITOOLS_COND1_COND2_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

# gmt_file <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_file <- "c5.go.bp.v7.4.symbols.gmt"

nTop_connect <- 5

myHeightGG <- myWidthGG <- 7

library("CEMiTool")
library(aracne.networks)
data("regulonbrca")

source("breast_utils.R")

# look at the size of aracne regulon

min_mod_size <- min(unlist(lapply(regulonbrca, function(x)length(x[["tfmode"]]))))
# 2 -> 47 modules
min_mod_size <- 10

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


cemi_annot_dt <- data.frame(SampleName=c(cond1_samps, cond2_samps),
                            Class = c(rep(cond1, length(cond1_samps)), rep(cond2, length(cond2_samps))), stringsAsFactors = FALSE)

cond12_dt <- cbind(cond1_dt, cond2_dt)

####################################
### run cemitools
####################################
# cemitool function receives the expression data, performs the co-expression modules analysis and returns a CEMiTool object

cem_cond12 <- cemitool(as.data.frame(cond12_dt), cemi_annot_dt)
# CEMiTool Object
# - Number of modules: 6 
# - Modules (data.frame: 546x2): 
#   genes        modules
# 1    CPB1             M4
# 2   IGSF1 Not.Correlated
# 3 CEACAM6             M4
# - Expression file: data.frame with 9995 genes and 18 samples
# - Selected data: 546 genes selected
# - Gene Set Enrichment Analysis: 
#   List containing 3 data.frames:
#   - $ es   : Enrichment Scores (6x3) 
# - $ nes  : Normalized Enrichment Scores (6x3) 
# - $ pval : p-value (x) 
# - Over Representation Analysis: null
# - Profile plot: ok
# - Enrichment plot: null
# - ORA barplot: null
# - Beta x R2 plot: null
# - Mean connectivity plot: null



####################################
###   Module inspection
####################################

# The module analysis has produced 7 modules and the allocation of genes to modules can be seen with the module_genes function:
  
# inspect modules

# nmodules(cem_cond12)
# 6
# head(module_genes(cem_cond12))
# genes        modules
# 1    CPB1             M4
# 2   IGSF1 Not.Correlated
# 3 CEACAM6             M4
# 4 SCGB1D2 Not.Correlated
# 5   SMCO3             M3
# 6    STC2             M3

cat(paste0("... # modules, default parameters:\t", nmodules(cem_cond12), "\n"))

# Genes that are allocated to Not.Correlated are genes that are not clustered into any module.

# If you wish to adjust the module definition parameters of your CEMiTool object, use find_modules(cem).
cem_cond12 <- find_modules(cem_cond12, min_ngen = min_mod_size)
cat(paste0("... # modules (min_ngen=", min_mod_size, "):\t", nmodules(cem_cond12), "\n"))


modules_dt <- module_genes(cem_cond12)
modulesFreq_dt <- as.data.frame(table(modules_dt$modules))
colnames(modulesFreq_dt) <- c("module", "nbr_genes")

outFile <- file.path(outFolder, paste0("nbrGenes_byModules.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(modulesFreq_dt$nbr_genes))
mtext(side=3, text=paste0("# modules = ", nrow(modulesFreq_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# get_hubs function to identnTop_connectify the top n genes with the highest connectivity in each module:
cond12_hubs <- get_hubs(cem_cond12,nTop_connect) 

# A summary statistic of the expression data within each module (either the module mean or eigengene) 
# can be obtained using: 
summary_stat <- mod_summary(cem_cond12)

generate_report(cem_cond12)

####################################
### Module enrichment
####################################
# 
# When sample annotation is provided, the cemitool function will automatically evaluate how the modules are up or
# down regulated between classes. 
# This is performed using the gene set enrichment analysis function from the fgsea package.
# 
# You can generate a plot of how the enrichment of the modules varies across classes with the plot_gsea function. 
# The size and intensity of the circles in the figure correspond to the Normalised Enrichment Score (NES),
# which is the enrichment score for a module in each class normalised by the number of genes in the module. 
# 
# This analysis is automatically run by the cemitool function, but it can be independently run with the function mod_gsea(cem).

# generate heatmap of gene set enrichment analysis
cem_cond12 <- mod_gsea(cem_cond12)
cem_cond12 <- plot_gsea(cem_cond12)
p <- show_plot(cem_cond12, "gsea")
outFile <- file.path(outFolder, paste0(cond1, "_", cond2,"_gsea_enrichment.svg"))
ggsave(p[[1]], filename = outFile, height=myHeightGG*1.3, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


####################################
#### Expression patterns in modules
####################################

# plot that displays the expression of each gene within a module using the plot_profile function:
  
# plot gene expression within each module
cem_cond12 <- plot_profile(cem_cond12)
prof_plots <- show_plot(cem_cond12, "profile")

outFile <- file.path(outFolder, paste0(cond1, "_", cond2,"_gsea_enrichment.svg"))
ggsave(p[[1]], filename = outFile, height=myHeightGG*1.3, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

mTopPlot <- 3
tmpgsea_dt <-cem_cond12@enrichment[["padj"]]
tmpgsea_dt <- tmpgsea_dt[order(tmpgsea_dt[,cond1], decreasing = FALSE),]
top3m_gseaPval_cond1 <- tmpgsea_dt$pathway[1:mTopPlot]
# tmpgsea_dt <- tmpgsea_dt[order(tmpgsea_dt[,cond2], decreasing = FALSE),]
# top3m_gseaPval_cond2 <- tmpgsea_dt$pathway[1:mTopPlot]



####################################
### Adding ORA analysis
####################################

# CEMiTool can determine which biological functions are associated with the modules by performing an over 
# representation analysis (ORA). To do this you must provide a pathway list in the form of GMT file. 
# CEMiTool will then analyze how these pathways are represented in the modules.
# 
# You can read in a pathway list formatted as a GMT file using the read_gmt function. 
# This example uses a GMT file that comes as part of the CEMiTool example data:

# read GMT file
gmt_in <- read_gmt(gmt_file)

# You can then perform ORA analysis on the modules in your CEMiTool object with the mod_ora function:
# perform over representation analysis
cem_cond12 <- mod_ora(cem_cond12, gmt_in)

# The numerical results of the analysis can be accessed with the ora_data function. 
# In order to visualise this, use plot_ora to add ORA plots to your CEMiTool object. 
# The plots can be accessed with the show_plot function.

# plot ora results
cem_cond12 <- plot_ora(cem_cond12)
ora_plots <- show_plot(cem_cond12, "ora")

for(i in 1:length(ora_plots)) {
  outFile <- file.path(outFolder, paste0(names(ora_plots)[i], "_GO_ora.", plotType))
  ggsave(ora_plots[[i]][["pl"]], filename = outFile, height=myHeightGG, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  ### if some signif GO -> to plot
  if(ora_plots[[i]][["numsig"]] > 0) {
    
    
    
    mg_dt <- module_genes(cem_cond12)
    moi <- names(ora_plots)[i]
    g_moi <- mg_dt$genes[mg_dt$modules == moi]
    
    p <- plot_mymodule(module_genes=g_moi, prot_dt=cond12_dt,
                  cond1_s=cond1_samps, cond2_s=cond2_samps,
                  meanExprThresh=1, absLog2fcThresh=0.5)
    
    outFile <- file.path(outFolder, paste0(moi, "_network_fc_corr.", plotType))
    ggsave(p, filename=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
  
}
# create report as html document
generate_report(cem_cond12, directory=file.path(outFolder, paste0("Report_", cond1, "_", cond2)))

# write analysis results into files
write_files(cem_cond12, directory=file.path(outFolder, paste0("Tables_", cond1, "_", cond2)))

# save all plots
save_plots(cem_cond12, "all", directory=file.path(outFolder, paste0("Plots_", cond1, "_", cond2)))




stop("---ok\n")



## $M1
## $M1$pl
####################################
### Adding interactions
####################################
# Interaction data, such as protein-protein interactions can be added to the CEMiTool object in order 
# to generate annotated module graphs. Interaction files are formatted as a data.frame or matrix 
# containing two columns for interacting pairs of genes.

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

##   gene1symbol gene2symbol
## 1         DBH      REPIN1
## 2      RBFOX2       HERC5
## 3      ZNF460      CCDC22
## 4     SNRNP40        OAZ3
## 5       SRSF6        OAZ3
## 6      SPTAN1       ARL8A

# You can add the interaction data to your CEMiTool object using the interactions_data function and 
# generate the plots with plot_interactions. Once again, the plots can be seen with the show_plot function:
  
  # plot interactions
  library(ggplot2)
interactions_data(cem_cond12) <- int_df # add interactions
cem_cond12 <- plot_interactions(cem_cond12) # generate plot
inter_plots <- show_plot(cem_cond12, "interaction") # view the plot for the first module
inter_plots[1]

####################################
### all in one
####################################

# , a CEMiTool object with all of the components mentioned above can also be constructed using just the cemitool function:
  # 
  # run cemitool
  library(ggplot2)
cem12_all <-  cemitool(as.data.frame(cond12_dt), cemi_annot_dt, gmt_in, interactions=int_df, 
                filter=TRUE, plot=TRUE, verbose=TRUE)
# create report as html document
generate_report(cem12_all, directory=file.path(outFolder, paste0("Report_", cond1, "_", cond2)))

# write analysis results into files
write_files(cem12_all, directory=file.path(outFolder, paste0("Tables_", cond1, "_", cond2)))

# save all plots
save_plots(cem12_all, "all", directory=file.path(outFolder, paste0("Plots_", cond1, "_", cond2)))


mg_dt <- module_genes(cem12_all)
moi <- "M1"
g_moi <- mg_dt$genes[mg_dt$modules == moi]
g_moi_10 <- g_moi[1:100]

gmoi_pairs <- combn(g_moi_10, 2)

moi_10_coexpr_dt <- foreach(i = 1:ncol(gmoi_pairs), .combine='rbind') %dopar% {
  g1 <- gmoi_pairs[1,i]
  g2 <- gmoi_pairs[2,i]
  g12corr <- cor(cond12_dt[g1,], cond12_dt[g2,])
  data.frame(g1 =g1, g2=g2, corr=g12corr,stringsAsFactors = FALSE )
  
}

# discussion here about signed or unsigned
# https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
#   
#   signed : correlations in the [− 1, 1] interval are scaled into the [0, 1] interval, 
#    unsigned.: while in the latter, negative correlations are made posi- tive.

###### added mz:
# draw the plot for each module, only gene from the module
# color of the edge = positive or negative correlation
# color of the node = fc between lumA and lumB
moi <- "M5"
moi_genes <- cem_cond12@module$genes[cem_cond12@module$modules== moi]

stopifnot(length(moi_genes) > 0)

moi_genes <- moi_genes[1:10]

moi_dt <- as.data.frame(t(combn(moi_genes, 2)))
colnames(moi_dt) <- c("gene1", "gene2")

stopifnot(moi_dt$gene1 %in% rownames(cond12_dt))
stopifnot(moi_dt$gene2 %in% rownames(cond12_dt))

stopifnot(moi_dt$gene1 %in% rownames(cond1_dt))
stopifnot(moi_dt$gene2 %in% rownames(cond1_dt))

stopifnot(moi_dt$gene1 %in% rownames(cond2_dt))
stopifnot(moi_dt$gene2 %in% rownames(cond2_dt))

i=1
moi_dt$corr <- foreach(i = 1:nrow(moi_dt), .combine='c') %dopar% {
  g1 <- moi_dt$gene1[i]
  g2 <- moi_dt$gene2[i]
  cor(cond12_dt[g1,], cond12_dt[g2,])
}
moi_dt$pair_label <- paste0(moi_dt$gene1, "_", moi_dt$gene2)
pair2corr <- setNames(moi_dt$corr, moi_dt$pair_label)

# average fold change between the 2 conditions
moi_genes_dt <- data.frame(gene = moi_genes, stringsAsFactors = FALSE)
moi_genes_dt$log2_fc <- foreach(i = 1:nrow(moi_genes_dt), .combine='c') %dopar% {
  g <- moi_genes_dt$gene[i]
  log2(median(cond1_dt[g,])/median(cond2_dt[g,]))
}
g2fc <- setNames(moi_genes_dt$log2_fc, moi_genes_dt$gene)

ig_obj <- graph_from_data_frame(moi_dt, directed=FALSE)
net_obj <- intergraph::asNetwork(ig_obj)
m <- network::as.matrix.network.adjacency(net_obj)
plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m,NULL))
colnames(plotcord) <- c("X1", "X2")
plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj,"vertex.names"))


edglist <- network::as.matrix.network.edgelist(net_obj)
edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[,  2], ])
colnames(edges)[colnames(edges) == "X1.1"] <- "Y1"
colnames(edges)[colnames(edges) == "X2.1"] <- "Y2"

stopifnot(edges$vertex.names.1 %in% moi_dt$gene2)
stopifnot(edges$vertex.names %in% moi_dt$gene1)

edges$pair_label <- paste0(edges$vertex.names, "_", edges$vertex.names.1)
stopifnot(edges$pair_label %in% names(pair2corr))
edges$corr <- pair2corr[paste0(edges$pair_label)]
stopifnot(!is.na(edges$corr))


plotcord$Hub <- 1
plotcord$shouldLabel <- TRUE
plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
stopifnot(plotcord$vertex.names %in% names(g2fc))
plotcord$fc <- g2fc[paste0(plotcord$vertex.names)]
plotcord$Degree <- 1
plotcord$color <- "blue"

edges$gradient_col <- seq_gradient_pal("red", "blue")(rescale(edges$corr))
plotcord$gradient_col <- seq_gradient_pal("red", "blue")(rescale(plotcord$fc))

pl <- ggplot(plotcord) + 
  geom_segment(data = edges, aes_(x = ~X1, y = ~X2, xend = ~Y1, yend = ~Y2, color=~corr), 
               size = 0.5, alpha = 0.5) + 
  geom_point(aes_(x = ~X1, y = ~X2,size = ~Degree, fill=~fc), shape=21) +
  scale_color_gradient2(  low = muted("blue"),
                           mid = "white",
                           high = muted("red")) +
  scale_fill_gradient2(  low = muted("blue"),
                          mid = "white",
                          high = muted("red"))+
  geom_label_repel(aes_(x = ~X1, 
                        y = ~X2, 
                        label = ~vertex.names, color = ~fc), 
                   box.padding = unit(1,"lines"))+
  theme_bw(base_size = 12, base_family = "") + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.key = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_blank(), 
        panel.grid = element_blank())





#######################3 TRASH
cem_cond1 <- cemitool(as.data.frame(cond1_dt))
cem_cond1
# CEMiTool Object
# - Number of modules: 7 
# - Modules (data.frame: 567x2): 
#   genes modules
# 1   STC2      M3
# 2   CHGA      M2
# 3 SH3BGR      M3
# - Expression file: data.frame with 9995 genes and 9 samples
# - Selected data: 567 genes selected
# - Gene Set Enrichment Analysis: null
# - Over Representation Analysis: null
# - Profile plot: ok
# - Enrichment plot: null
# - ORA barplot: null
# - Beta x R2 plot: null
# - Mean connectivity plot: null

# As a default, the cemitool function first performs a filtering of the gene expression data before running the remaining analyses.
# This filtering is done in accordance to gene variance.
# In this example the filtering step has reduced the gene number to 567.

nmodules(cem_cond1)
## [1] 7

head(module_genes(cem_cond1))
# genes modules
# 1    STC2      M3
# 2    CHGA      M2
# 3  SH3BGR      M3
# 4     MGP      M3
# 5    CPB1      M1
# 6 CEACAM6      M1








