
# Rscript viper_and_ppi.R

cond1="LumB"
cond2="LumA"

require(doMC)
registerDoMC(2)

outFolder <- file.path("VIPER_AND_PPI", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

# gmt_file <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_file <- "c5.go.bp.v7.4.symbols.gmt"

nTop_connect <- 5

myHeightGG <- myWidthGG <- 7

runPPIenrich <- FALSE

library(STRINGdb)
library(ggplot2)
library(igraph)
library(foreach)
source("breast_utils.R")

inFolder <- file.path("VIPER_COND1COND2_NETWORKS_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
inFile <- file.path(inFolder, "mrs_summary_all.Rdata")
mrs_summary_all <- get(load(file=inFile))
mrs_summary_all <- mrs_summary_all[order(mrs_summary_all$p.value, decreasing = FALSE),]


inFile <- file.path(inFolder, "mrs_brca_regul.Rdata")
brca_regul <- get(load( file=inFile))

inFolder2 <- file.path("STRINGDB_COND1_COND2_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
inFile <- file.path(inFolder2, "DE_topTable_mpd.Rdata")
DE_topTable_mpd <- get(load(inFile))
## unmapped: 9957/9995
# stopifnot(!duplicated(DE_topTable_mpd$gene))
# symb2ppi <- setNames(DE_topTable_mpd$STRING_id, DE_topTable_mpd$gene)
# cannot do that, can be duplicated because mulitple PPI for a given symbol

# stopifnot(!duplicated(DE_topTable_mpd$STRING_id))

### initialization
# WARNING: You didn't specify a species. Hence we will set 9606 (Homo Sapiens) as your species.
# WARNING: Score threshold is not specified. We will be using medium stringency cut-off of 400.
# WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING  11.0 
string_db <- STRINGdb$new( species=9606)

# look at top 10 regulons

nTop <- 10

plotType <- "png"
myWidth <- 1000
myHeight <- 1000


# post payload information to the STRING server
DE_topTable_mpd <- string_db$add_diff_exp_color( DE_topTable_mpd,
                                                        logFcColStr="logFC" )    

payload_id <- string_db$post_payload( DE_topTable_mpd$STRING_id, 
                                      colors=DE_topTable_mpd$color )

if(runPPIenrich) {
  out_dt <- foreach(i = 1:nrow(mrs_summary_all), .combine='rbind') %dopar% {
    
    cat(paste0("... ", i, "/", nrow(mrs_summary_all), "\n"))
    
    i_reg <- mrs_summary_all$Regulon[i]
    i_targets <- names(brca_regul[[paste0(i_reg)]][["tfmode"]])
    
    mod_genes <- c(i_reg, i_targets)
    mod_genes <- mod_genes[mod_genes %in% DE_topTable_mpd$gene]
    
    mod_dt <- DE_topTable_mpd[DE_topTable_mpd$gene %in% mod_genes,c("gene", "STRING_id")]
    
    outrow <- as.data.frame(string_db$get_ppi_enrichment(mod_dt$STRING_id))
    outrow$reg <- i_reg
    outrow$genes <- paste0(i_targets, collapse=",")
    outrow$genes_ppi <- paste0(mod_dt$gene, collapse=",")
    outrow$mod_size <- length(i_targets) + 1
    outrow$mod_size_ppi <- nrow(mod_dt)
    outrow <- outrow[,c("reg", "mod_size","mod_size_ppi", "enrichment", "edges", "lambda", "genes", "genes_ppi")]
    outrow
  }
  outFile <- file.path(outFolder, paste0("out_dt.RData"))
  save(out_dt, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))  
  
  
} else {
  outFile <- file.path(outFolder, paste0("out_dt.RData"))
  out_dt <- get(load(outFile))
  
  plot(density(out_dt$enrichment))
}
stopifnot(!duplicated(mrs_summary_all$Regulon))
reg_viperP <- setNames(mrs_summary_all$FDR, mrs_summary_all$Regulon)

stopifnot(!duplicated(out_dt$reg))
reg_ppiP <- setNames(out_dt$enrichment, out_dt$reg)

stopifnot(setequal(names(reg_viperP), names(reg_ppiP)))
myregs <- names(reg_viperP)

plot(
  x = reg_viperP[myregs],
  y = reg_ppiP[myregs],
  cex=0.7,
  pch=16,
  xlab="viper FDR",
  ylab="PPI p-val"
)
abline(h=0.05, col="red")
abline(v=0.05, col="red")
mtext(side=3, text=paste0("# reg = ", length(myregs)))


plot(
  x = -log10(reg_viperP[myregs+0.01]),
  y = -log10(reg_ppiP[myregs]+0.01),
  cex=0.7,
  pch=16,
  xlab="viper FDR [-log10]",
  ylab="PPI p-val [-log10]"
)
abline(h=-log10(0.05), col="red")
abline(v=-log10(0.05), col="red")
mtext(side=3, text=paste0("# reg = ", length(myregs)))



for(i in 1:nTop) {
  i_reg <- mrs_summary_all$Regulon[i]
  i_targets <- names(brca_regul[[paste0(i_reg)]][["tfmode"]])
  
  mod_genes <- c(i_reg, i_targets)
  mod_genes <- mod_genes[mod_genes %in% DE_topTable_mpd$gene]
  
  mod_dt <- DE_topTable_mpd[DE_topTable_mpd$gene %in% mod_genes,c("gene", "STRING_id")]
  
  mod_interactions_dt <- string_db$get_interactions( mod_dt$STRING_id )
  
  outFile <- file.path(outFolder, paste0("ppi_plot_cluster", i, "_withColFC.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  string_db$plot_network(mod_dt$STRING_id, payload_id=payload_id )
  mtext(side=3, paste0(i_reg, " regulon (viper top ", i, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
}



### ppi enrichment



