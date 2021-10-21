
# Rscript viper_and_ppi.R

cond1="LumB"
cond2="LumA"

outFolder <- file.path("VIPER_AND_PPI", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

# gmt_file <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_file <- "c5.go.bp.v7.4.symbols.gmt"

nTop_connect <- 5

myHeightGG <- myWidthGG <- 7

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

### initialization
# WARNING: You didn't specify a species. Hence we will set 9606 (Homo Sapiens) as your species.
# WARNING: Score threshold is not specified. We will be using medium stringency cut-off of 400.
# WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING  11.0 
string_db <- STRINGdb$new( species=9606)

# look at top 10 regulons

nTop <- 10

for(i in 1:nTop) {
  i_reg <- mrs_summary_all$Regulon[i]
  i_targets <- names(brca_regul[[paste0(i_reg)]][["tfmode"]])
  
  mod_genes <- c(i_reg, i_targets)
  mod_genes <- mod_genes[mod_genes %in% DE_topTable_mpd$gene]
  
  mod_dt <- DE_topTable_mpd[DE_topTable_mpd$gene %in% mod_genes,c("gene", "STRING_id")]
  
  mod_interactions_dt <- string_db$get_interactions( mod_dt$STRING_id )
  
  outFile <- file.path(outFolder, paste0("ppi_plot_cluster", i, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  string_db$plot_network(mod_dt$STRING_id)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
  
  
}






