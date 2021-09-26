# Rscript vae_result_exploration.R

outFolder <- file.path("VAE_RESULT_EXPLORATION_JOHANSSON_MOSTVAR")
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- 400
myWidth <- 400

indata_file <- file.path("../data/johansson_data_relative_ratios_to_pool_1000mostVar.txt")
vae_latent_repr_file <- file.path("BREAST_PROTEO_JOHANSSONMOSTVAR_VAE/latent_repr_100LD_256DS_150_128_0.2_1_mmd.csv")
corr_file <- file.path("BREAST_PROTEO_JOHANSSONMOSTVAR_VAE/correlations_all_100LD_256DS_150_128_0.2_1_mmd.csv")
corrpval_file <- file.path("BREAST_PROTEO_JOHANSSONMOSTVAR_VAE/p_values_all_100LD_256DS_150_128_0.2_1_mmd.csv")
annot_file <- file.path("../data/johansson_tumor_annot.csv")

annot_dt <- read.delim(annot_file, header=T, stringsAsFactors = FALSE, sep=",")
stopifnot(!duplicated(annot_dt$Tumor.ID))
samp2pam50 <- setNames(annot_dt$PAM50.subtype, annot_dt$Tumor.ID)

in_dt <- read.delim(indata_file, header=T, stringsAsFactors = FALSE, sep=",")
dim(in_dt)
expr_dt <- in_dt
expr_dt$gene_symbol <- expr_dt$ensembl_id <- NULL
# 1000x47
stopifnot(!duplicated(in_dt$gene_symbol))
stopifnot(!duplicated(in_dt$ensembl_id))
stopifnot(setequal(names(samp2pam50), colnames(expr_dt)))
exprT_dt <- t(expr_dt)
nSamp <- ncol(expr_dt)


rownames(annot_dt) <- annot_dt$Tumor.ID
stopifnot(setequal(rownames(annot_dt), colnames(expr_dt)))
annot_dt <- annot_dt[colnames(expr_dt),]
annot_dt$PAM50.subtype <- as.factor(annot_dt$PAM50.subtype)
col_lab <- "PAM50.subtype"


vae_lr_dt <- read.delim(vae_latent_repr_file, header=F, stringsAsFactors = FALSE, sep=",")
dim(vae_lr_dt)
rownames(vae_lr_dt) <- colnames(expr_dt)
# 45x100

corr_dt <- read.delim(corrpval_file, header=F, stringsAsFactors = FALSE, sep=",")
dim(corr_dt)
rownames(corr_dt) <- in_dt$gene_symbol
# 1000  100

pval_dt <- read.delim(corrpval_file, header=F, stringsAsFactors = FALSE, sep=",")
dim(pval_dt)
rownames(pval_dt) <- in_dt$gene_symbol
# 1000  100

# PCA on the input data




data_in <- "vae_lr_dt"
for(data_in in c("vae_lr_dt","exprT_dt")) {
  data_dt <- eval(parse(text=data_in))
  pca_data <- prcomp(data_dt)
  pca_data_lowrepr <- pca_data$x
  stopifnot(nrow(pca_data_lowrepr) == nSamp)
  
  stopifnot(colnames(expr_dt) == rownames(annot_dt))
  
  # dev.off()
  
  outFile <- file.path(outFolder, paste0("pca_", data_in, "_johanssonMostVar.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(x = pca_data_lowrepr[,1],
       y = pca_data_lowrepr[,2],
       xlab="PC1", ylab="PC2",
       main=paste0(data_in),
       pch=16,
       col=as.numeric(annot_dt[,paste0(col_lab)]))
  mtext(text=paste0("# samp=",nSamp), side=3)
  legend("topleft", legend=as.character(levels(annot_dt[,paste0(col_lab)])),
         pch=16, 
         col=as.numeric(unique(annot_dt[,paste0(col_lab)])), 
         bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}



lds_var <- apply(vae_lr_dt, 2, var)
stopifnot(!is.na(lds_var))
nTopVarLDs <- 5
mostVar_lds <- order(lds_var, decreasing = TRUE)[1:nTopVarLDs]

all_ld_pairs <- combn(x=c(mostVar_lds), m=2)

i=1
for(i in c(1:ncol(all_ld_pairs))) {
  ld1 <- all_ld_pairs[1,i]
  ld2 <- all_ld_pairs[2,i]
  
  ld1_rank <- which(mostVar_lds == ld1)
  ld2_rank <- which(mostVar_lds == ld2)
  
  # dev.off()
  outFile <- file.path(outFolder, paste0("LD", ld1, "var", ld1_rank, "_vs_LD", ld2, "var", ld2_rank, ".", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(x = vae_lr_dt[,ld1],
       y = vae_lr_dt[,ld2],
       xlab=paste0("LD", ld1, " (var rank: ", ld1_rank, ")"), 
       ylab=paste0("LD", ld2, " (var rank: ", ld2_rank, ")"),
       main=paste0("Samp. in LD"),
       pch=16,
       col=as.numeric(annot_dt[,paste0(col_lab)]))
  mtext(text=paste0("# samp=",nSamp, "; # features=", ncol(vae_lr_dt)), side=3)
  legend("topleft", legend=as.character(levels(annot_dt[,paste0(col_lab)])),
         pch=16, 
         col=as.numeric(unique(annot_dt[,paste0(col_lab)])), 
         bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}





