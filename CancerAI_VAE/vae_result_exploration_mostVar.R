# Rscript vae_result_exploration.R
library(ComplexHeatmap)


outFolder <- file.path("VAE_RESULT_EXPLORATION_JOHANSSON_MOSTVAR")
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- 400
myWidth <- 400
plotTypeGG <- "svg"
myHeightGG <- 7
myWidthGG <- 7

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
lds_var_ord <- order(lds_var, decreasing = TRUE)
mostVar_lds <- lds_var_ord[1:nTopVarLDs]

all_ld_pairs <- combn(x=c(mostVar_lds), m=2)

i=1

# outFile <- file.path(outFolder, paste0("LD", ld1, "var", ld1_rank, "_vs_LD", ld2, "var", ld2_rank, ".", plotType ))
outFile <- file.path(outFolder, paste0("LDi_vsLDj_all_topVarRanks.", plotTypeGG))
# do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))
do.call(plotTypeGG, list(outFile, height=myHeightGG*3, width=myWidthGG*2.5))

par(mfrow=c(4,3))

for(i in c(1:ncol(all_ld_pairs))) {
  ld1 <- all_ld_pairs[1,i]
  ld2 <- all_ld_pairs[2,i]
  
  ld1_rank <- which(mostVar_lds == ld1)
  ld2_rank <- which(mostVar_lds == ld2)
  
  # dev.off()
  # outFile <- file.path(outFolder, paste0("LD", ld1, "var", ld1_rank, "_vs_LD", ld2, "var", ld2_rank, ".", plotType ))
  # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(x = vae_lr_dt[,ld1],
       y = vae_lr_dt[,ld2],
       xlab=paste0("LD", ld1, " (var rank: ", ld1_rank, ")"), 
       ylab=paste0("LD", ld2, " (var rank: ", ld2_rank, ")"),
       main=paste0("Samp. in LD"),
       cex=1.2,
       cex.axis=2,
       cex.lab=2,
       pch=16,
       col=as.numeric(annot_dt[,paste0(col_lab)]))
  mtext(text=paste0("# samp=",nSamp, "; # features=", ncol(vae_lr_dt)), side=3)
  legend("topleft", legend=as.character(levels(annot_dt[,paste0(col_lab)])),
         pch=16, 
         col=as.numeric(unique(annot_dt[,paste0(col_lab)])), 
         bty="n")
  # foo <- dev.off()
  # cat(paste0("... written: ", outFile, "\n"))
}
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################ least var for each subtype

# retrieve for each subtype the least var

tmp_dt <- vae_lr_dt
stopifnot(rownames(tmp_dt) %in% annot_dt$Tumor.ID)
tmp_dt$PAM50 <- samp2pam50[rownames(tmp_dt)]
stopifnot(!is.na(tmp_dt$PAM50))  

nLDs=100

x=tmp_dt[tmp_dt$PAM50==tmp_dt$PAM50[1],]
leastVar <- c(by(tmp_dt, tmp_dt$PAM50, function(x){
  subt <- unique(x$PAM50)
  x$PAM50 <- NULL
  samp_vars <- apply(x, 2, var)
  stopifnot(length(samp_vars) == nLDs)
  setNames(which.min(samp_vars), subt)
}))

stopifnot(!duplicated(leastVar))

all_ld_pairs <- combn(x=c(leastVar), m=2)

i=1

# outFile <- file.path(outFolder, paste0("LD", ld1, "var", ld1_rank, "_vs_LD", ld2, "var", ld2_rank, ".", plotType ))
outFile <- file.path(outFolder, paste0("LDi_vsLDj_all_subtypeLeastVarRanks.", plotTypeGG))
# do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))
do.call(plotTypeGG, list(outFile, height=myHeightGG*3, width=myWidthGG*2.5))

par(mfrow=c(4,3))

for(i in c(1:ncol(all_ld_pairs))) {
  ld1 <- all_ld_pairs[1,i]
  ld2 <- all_ld_pairs[2,i]
  
  ld1_rank <- names(leastVar)[leastVar == ld1]
  ld2_rank <- names(leastVar)[leastVar == ld2]
  
  # dev.off()
  # outFile <- file.path(outFolder, paste0("LD", ld1, "var", ld1_rank, "_vs_LD", ld2, "var", ld2_rank, ".", plotType ))
  # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(x = vae_lr_dt[,ld1],
       y = vae_lr_dt[,ld2],
       xlab=paste0("LD", ld1, " (least var. for: ", ld1_rank, ")"), 
       ylab=paste0("LD", ld2, " (least var. for: ", ld2_rank, ")"),
       main=paste0("Samp. in LD"),
       cex=1.2,
       cex.axis=2,
       cex.lab=2,
       pch=16,
       col=as.numeric(annot_dt[,paste0(col_lab)]))
  mtext(text=paste0("# samp=",nSamp, "; # features=", ncol(vae_lr_dt)), side=3)
  legend("topleft", legend=as.character(levels(annot_dt[,paste0(col_lab)])),
         pch=16, 
         col=as.numeric(unique(annot_dt[,paste0(col_lab)])), 
         bty="n")
  # foo <- dev.off()
  # cat(paste0("... written: ", outFile, "\n"))
}
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################ complex heatmap

annotation_colors <- list("PAM50"= c("Basal" = "red",
                                     "HER2" = "pink",
                                     "LumA" = "darkblue",
                                     "LumB" = "lightblue",
                                     "Normal" = "forestgreen"))



mat <- vae_lr_dt
colnames(mat) <- paste0("LD", 1:ncol(mat))

stopifnot(rownames(mat) %in% names(samp2pam50))
pam50 <- samp2pam50[rownames(mat)]
    
outFile <- file.path(outFolder, paste0("complex_heatmap_LDs.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*1.6))
df = data.frame(PAM50=pam50)
ha = HeatmapAnnotation(df = df, col=annotation_colors)
Heatmap(t(mat), name = "Prot. LS", top_annotation = ha, row_title = "LDs", row_labels = rep("", nrow(t(mat))))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


                  



stop("--ok\n")

mat = cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),
            rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)

df = data.frame(type = c(rep("a", 5), rep("b", 5)))
value = rnorm(10)
ha = HeatmapAnnotation(df = df)
Heatmap(mat, name = "foo", top_annotation = ha)

########################################################
library(ComplexHeatmap)
library(circlize)

expr = readRDS(system.file(package = "ComplexHeatmap", "extdata", "gene_expression.rds"))
mat = as.matrix(expr[, grep("cell", colnames(expr))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))

type = gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(type = type, annotation_name_side = "left")

ht_list = Heatmap(mat_scaled, name = "expression", row_km = 5, 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) +
  Heatmap(base_mean, name = "base mean", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), 
                                                                    height = unit(2, "cm"))),
          width = unit(15, "mm")) +
  rowAnnotation(length = anno_points(expr$length, pch = 16, size = unit(1, "mm"), 
                                     axis_param = list(at = c(0, 2e5, 4e5, 6e5), 
                                                       labels = c("0kb", "200kb", "400kb", "600kb")),
                                     width = unit(2, "cm"))) +
  Heatmap(expr$type, name = "gene type", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
          width = unit(15, "mm"))

ht_list = rowAnnotation(block = anno_block(gp = gpar(fill = 2:6, col = NA)), 
                        width = unit(2, "mm")) + ht_list

draw(ht_list)