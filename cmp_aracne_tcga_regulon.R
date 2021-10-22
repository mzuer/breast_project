# Rscript cmp_aracne_tcga_regulon.R

outFolder <- file.path("CMP_ARACNE_TCGA_REGULON")
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myWidth <- myHeight <- 400

###################### CHECKOUT COEXPR IN TCGA DATA

require(foreach)


prepTCGAdata <- FALSE
computCorrTCGA <- TRUE
computCorrProteo <- TRUE

############################################## 


if(prepTCGAdata) {
  
  library(TCGAbiolinks)
  library(recount)
  library(org.Hs.eg.db)
  
  # breast_rec2_tcga <- TCGAquery_recount2(project="tcga", tissue = "breast")
  # save(breast_rec2_tcga, file=file.path(outFolder, "breast_rec2_tcga.Rdata"))
  inFolder <- file.path("../ovarian_project/tcga_data/DOWNLOAD_TCGA_BRCA_LUMALUMB_RECOUNT2/")
  breast_rec2_tcga <- get(load(file.path(inFolder, "breast_rec2_tcga.Rdata")))
  breast_rec2_tcga_scaled <- scale_counts(breast_rec2_tcga$tcga_breast)
  tcga_dt <- assays(breast_rec2_tcga_scaled)$counts
  
  stopifnot(breast_rec2_tcga_scaled@colData[,"cgc_case_primary_site"] == "Breast")
  annot_cols <- c(grep("^cgc_", names(breast_rec2_tcga_scaled@colData)), grep("^gdc_", names(breast_rec2_tcga_scaled@colData)))
  tcgacols <- names(breast_rec2_tcga_scaled@colData)[annot_cols]
  tcga_sampleAnnot <- data.frame(breast_rec2_tcga_scaled@colData[, tcgacols])
  ####~~~ 1st FILTER HERE -> TAKE ONLY PRIMARY TUMOR !!!
  stopifnot(rownames(tcga_sampleAnnot) == colnames(tcga_dt))
  #*** update here:
  cat(paste0("... filter 1 - primary tumor: ",sum(tcga_sampleAnnot$cgc_sample_sample_type == "Primary Tumor"), "/",
             nrow(tcga_sampleAnnot), "\n"))
  tcga_sampleAnnot <- tcga_sampleAnnot[tcga_sampleAnnot$cgc_sample_sample_type == "Primary Tumor",]
  
  stopifnot(rownames(tcga_sampleAnnot) %in% colnames(tcga_dt))
  
  # any(duplicated(tcga_sampleAnnot$cgc_sample_id))
  dim(tcga_dt)
  tcga_dt <- tcga_dt[,rownames(tcga_sampleAnnot)]
  dim(tcga_dt)
  
  # remove the duplicates -> para on x chromosomes
  all_ensembl <- breast_rec2_tcga_scaled@rowRanges@elementMetadata$gene_id
  all_ensembl_short <- gsub("\\..+", "", all_ensembl)
  
  rown_short <- gsub("\\..+", "", rownames(tcga_dt))
  stopifnot(which(duplicated(all_ensembl_short)) == which(duplicated(rown_short)))
  
  torm <- all_ensembl_short[duplicated(all_ensembl_short)]
  
  tcga_nd_dt <- tcga_dt[!rown_short %in% torm,]
  stopifnot(!duplicated(gsub("\\..+", "", rownames(tcga_nd_dt))))
  rownames(tcga_nd_dt) <- gsub("\\..+", "", rownames(tcga_nd_dt))
  
  
  all_entrez <- mapIds(org.Hs.eg.db, rownames(tcga_nd_dt), 'ENTREZID', 'ENSEMBL')
  
  all(rownames(tcga_nd_dt) %in% names(all_entrez))
  
  tcga_nd_dt <- as.data.frame(tcga_nd_dt)
  tcga_nd_dt$entrez <- all_entrez[paste0(rownames(tcga_nd_dt))]
  dim(tcga_nd_dt)
  tcga_nd_dt <- tcga_nd_dt[!is.na(tcga_nd_dt$entrez),]
  dim(tcga_nd_dt)
  stopifnot(!is.na(tcga_nd_dt$entrez))
  any(duplicated(tcga_nd_dt$entrez))
  
  sum(duplicated(tcga_nd_dt$entrez))
  
  de <- unique(tcga_nd_dt$entrez[duplicated(tcga_nd_dt$entrez)])
  
  dup_dt <- tcga_nd_dt[tcga_nd_dt$entrez %in% de,]
  nodup_dt <- tcga_nd_dt[!tcga_nd_dt$entrez %in% de,]
  stopifnot(nrow(dup_dt) + nrow(nodup_dt) == nrow(tcga_nd_dt))
  
  # aggregate to take the mean value
  x=dup_dt[dup_dt$entrez==dup_dt$entrez[1],]
  cat(paste0("... agregate counts to entrez\n"))
  new_dup_dt <- do.call(rbind, by(dup_dt,dup_dt$entrez, function(x){
    tmp <- x
    tmp$entrez <- NULL
    out_dt <- as.data.frame(t(as.data.frame(colMeans(tmp))))
    stopifnot(nrow(out_dt) == 1)
    out_dt$entrez <- unique(x$entrez)
    stopifnot(ncol(out_dt) == ncol(x))
    out_dt
  }))
  stopifnot(!duplicated(new_dup_dt$entrez))
  
  tcga_uniq_dt <- rbind(new_dup_dt, nodup_dt)
  stopifnot(!duplicated(tcga_uniq_dt$entrez))
  
  tcga_entrez_dt <- tcga_uniq_dt
  rownames(tcga_entrez_dt) <- tcga_entrez_dt$entrez
  tcga_entrez_dt$entrez <- NULL
  
  
  outFile <- file.path(outFolder,"tcga_entrez_dt.Rdata" )
  save(tcga_entrez_dt, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder,"tcga_entrez_dt.Rdata" )
  tcga_entrez_dt <- get(load(outFile))
}

if(computCorrTCGA) {
  library(aracne.networks)
  data("regulonbrca")
  
  
  ### check how many go in the same direction in prot data
  all_tfms_dt <- do.call(rbind, lapply(1:length(regulonbrca), function(x) { 
    tfm <- regulonbrca[[x]][["tfmode"]]
    data.frame(g1=names(regulonbrca)[x],
               g2=names(tfm),
               tfm = as.numeric(tfm),
               stringsAsFactors = FALSE)
  }))
  all_tfms_dt$g1 <- as.character(all_tfms_dt$g1)
  all_tfms_dt$g2 <- as.character(all_tfms_dt$g2)
  
  stopifnot(any(all_tfms_dt$g1 %in% rownames(tcga_entrez_dt)))
  stopifnot(any(all_tfms_dt$g2 %in% rownames(tcga_entrez_dt)))
  nrow(all_tfms_dt)
  # 331919
  all_tfms_dt <- all_tfms_dt[all_tfms_dt$g1 %in% rownames(tcga_entrez_dt) & all_tfms_dt$g2 %in% rownames(tcga_entrez_dt),]
  nrow(all_tfms_dt)
  # 322597
  cat(paste0("... start computing corr TCGA"))
  all_tfms_dt_tcga <- all_tfms_dt
  all_tfms_dt_tcga$prot_corr <-foreach(i = 1:nrow(all_tfms_dt_tcga), .combine='c') %dopar% {
    g1 <- all_tfms_dt_tcga$g1[i]
    g2 <- all_tfms_dt_tcga$g2[i]
    if(g1 %in% rownames(tcga_entrez_dt) & g2 %in% rownames(tcga_entrez_dt)) {
      cor(as.numeric(tcga_entrez_dt[paste0(g1),]), as.numeric(tcga_entrez_dt[paste0(g2),]))
    } else {
      NA
    }
  }  
  
  stopifnot(!is.na(all_tfms_dt_tcga))
  
  outFile <- file.path(outFolder,"all_tfms_dt_tcga.Rdata" )
  save(all_tfms_dt_tcga, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder,"all_tfms_dt_tcga.Rdata" )
  all_tfms_dt_tcga <- get(load(outFile))
  
}



###################### COMPARE PROTEO DATA

if(computCorrProteo) {
  cond1="LumB"
  cond2="LumA"
  inFolder <- file.path("VIPER_COND1COND2_NETWORKS_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
  outFile <- file.path(inFolder, "mrs_brca_regul.Rdata")
  brca_regul <- get(load(file=outFile))
  brca_regul_filt3 <- brca_regul
  cat(paste0("... written: ", outFile, "\n"))
  
  proteo_dt <- read.delim("data/johansson_data_relative_ratios_to_pool.csv", sep=",")
  stopifnot(!duplicated(proteo_dt$gene_symbol))
  rownames(proteo_dt) <- proteo_dt$gene_symbol
  proteo_dt$gene_symbol <- proteo_dt$gene_symbol <- proteo_dt$ensembl_id <- NULL 
  
  
  ### check how many go in the same direction in prot data
  all_tfms_dt <- do.call(rbind, lapply(1:length(brca_regul_filt3), function(x) { 
    tfm <- brca_regul_filt3[[x]][["tfmode"]]
    data.frame(g1=names(brca_regul_filt3)[x],
               g2=names(tfm),
               tfm = as.numeric(tfm),
               stringsAsFactors = FALSE)
  }))
  stopifnot(all_tfms_dt$g1 %in% rownames(proteo_dt))
  stopifnot(all_tfms_dt$g2 %in% rownames(proteo_dt))
  
  all_tfms_dt_proteo <- all_tfms_dt
  all_tfms_dt_proteo$prot_corr <-foreach(i = 1:nrow(all_tfms_dt_proteo), .combine='c') %dopar% {
    g1 <- all_tfms_dt_proteo$g1[i]
    g2 <- all_tfms_dt_proteo$g2[i]
    if(g1 %in% rownames(proteo_dt) & g2 %in% rownames(proteo_dt)) {
      cor(as.numeric(proteo_dt[paste0(g1),]), as.numeric(proteo_dt[paste0(g2),]))
    } else {
      NA
    }
  }
  stopifnot(!is.na(all_tfms_dt_proteo))
  
  outFile <- file.path(outFolder,"all_tfms_dt_proteo.Rdata" )
  save(all_tfms_dt_proteo, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder,"all_tfms_dt_proteo.Rdata" )
  all_tfms_dt_proteo <- get(load(outFile))
}

stop("---ok\n")



outFile <- file.path(outFolder, paste0("cmp_tfmodeARACNe_proteoCorr.", plotType))
do.call(plotType, list(height=myHeight, width=myWidth))
plot(x=all_tfms_dt$tfm,
     y=all_tfms_dt$prot_corr,
     xlab="ARACNe TFm",
     ylab="proteo corr.")
mtext(side=3, text=paste0("proteo corr: ", cond1, "+", cond2, " data"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
