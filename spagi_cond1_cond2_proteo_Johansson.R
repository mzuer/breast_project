

library(spagi)
require(ggplot2)
require(ggrepel)


# Rscript spagi_cond1_cond2_proteo_Johansson.R

plotType <- "png"
myWidth <- 400
myHeight <- 400

cond1 <- "LumB"
cond2 <- "LumA"
annotCol <- "PAM50.subtype"

de_file <- file.path("STRINGDB_COND1_COND2_PROTEO_JOHANSSON/test_LumB_vs_ref_LumA/DE_topTable.Rdata")

outFolder <- file.path("SPAGI_COND1_COND2_PROTEO_JOHANSSON", paste0("test_", cond1, "_vs_ref_", cond2))
dir.create(outFolder, recursive=TRUE)

myHeightGG <- myWidthGG <- 7

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


# Pre-process the query data

stopifnot(cond1_dt >= 0) # there are raw ratios

log2_cond1_dt <- log2(cond1_dt+0.01) # need the data in log2
# Here, we will use the expression cutoff as 1.8 of the query data.

plot(density(log2_cond1_dt))
# mean(ROR1.data>=1.8) = 0.25 => they keep ~25% of the data
exp_cutoff_cond1 <- quantile(log2_cond1_dt, probs = 0.75)
abline(h=exp_cutoff_cond1, col="red") # 0.25

cond1_pd_dt <- preprocess_querydata(cell.tissue.data = log2_cond1_dt, exp.cutoff.th = exp_cutoff_cond1)

# Identify active pathway paths of the processed query data

# Here we will use the background pathway.path data
# to get the active pathway paths of the processed query data.

cond1_activePath <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = cond1_pd_dt)

# Get active pathway ranking metric (i.e., activity score and number of downstream transcription factors)

# Here we will use the ROR1.active.pathway and ROR1.processed.data data sets to get the active pathway ranking metric. 
# Also we will use a high expression threshold (here we will use 7) for the processed query data.
# mean(ROR1.data > 7) # 0.25
exp_cutoff_high_cond1 <- quantile(log2_cond1_dt, probs=0.975)
cond1_activePath_rankmetric <- get_pathway_ranking_metric(active.pathway.path = cond1_activePath,
                                                          processed.query.data = cond1_pd_dt, high.exp.th = exp_cutoff_high_cond1)

# stopifnot(length(cond1_activePath_rankmetric[["activity.score"]]) == length(cond1_activePath_rankmetric[["downstream.tf.count"]]))
# not TRUE one is NULL and the other is absent if not present
cond1_activ_dt <- do.call(rbind, lapply(colnames(cond1_dt), function(x) {
  tmp <- cond1_activePath_rankmetric[["activity.score"]][[x]]
  if(is.null(tmp)) {
    data.frame(
      cond = cond1,
      samp = x,
      tf = NA,
      activity_score = NA,
      tf_count = NA,
      stringsAsFactors = FALSE
    )
  } else {
  data.frame(
    cond = cond1,
    samp = x,
    tf = names(tmp),
    activity_score = as.numeric(tmp),
    tf_count = as.numeric(cond1_activePath_rankmetric[["downstream.tf.count"]][[x]]),
    stringsAsFactors = FALSE
  )
  }}))


# Plot the ranking metric result in a 2D plane

# After getting the active pathway ranking metric result you can display them in your preferred format.
# Here we will plot them in a 2D plane where x-axis denotes the number of downstream transcription factors and 
# y-axis denotes the activity score for each pathway.

# display_pathway_ranking_metric(pathway.ranking.metric = cond1_activePath_rankmetric)
#> [1] "ROR1_LEC -- result plotting done!!"
#To separate the top ranked pathways we can do this
abline(v=45, h=0.2, lty=2, col="black")



log2_cond2_dt <- log2(cond2_dt+0.01) # need the data in log2
# Here, we will use the expression cutoff as 1.8 of the query data.

plot(density(log2_cond2_dt))
# mean(ROR1.data>=1.8) = 0.25 => they keep ~25% of the data
exp_cutoff_cond2 <- quantile(log2_cond2_dt, probs = 0.75)
abline(h=exp_cutoff_cond2, col="red") # 0.25

cond2_pd_dt <- preprocess_querydata(cell.tissue.data = log2_cond2_dt, exp.cutoff.th = exp_cutoff_cond2)

# Identify active pathway paths of the processed query data

# Here we will use the background pathway.path data
# to get the active pathway paths of the processed query data.

cond2_activePath <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = cond2_pd_dt)


# Get active pathway ranking metric (i.e., activity score and number of downstream transcription factors)

# Here we will use the ROR1.active.pathway and ROR1.processed.data data sets to get the acitve pathway ranking metric. 
# Also we will use a high expression threshold (here we will use 7) for the processed query data.
mean(ROR1.data > 7) # 0.25
exp_cutoff_high_cond2 <- quantile(log2_cond2_dt, probs=0.975)
cond2_activePath_rankmetric <- get_pathway_ranking_metric(active.pathway.path = cond2_activePath,
                                                          processed.query.data = cond2_pd_dt, high.exp.th = exp_cutoff_high_cond2)

cond2_activ_dt <- do.call(rbind, lapply(colnames(cond2_dt), function(x) {
  tmp <- cond2_activePath_rankmetric[["activity.score"]][[x]]
  if(is.null(tmp)) {
    data.frame(
      cond = cond2,
      samp = x,
      tf = NA,
      activity_score = NA,
      tf_count = NA,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      cond = cond2,
      samp = x,
      tf = names(tmp),
      activity_score = as.numeric(tmp),
      tf_count = as.numeric(cond2_activePath_rankmetric[["downstream.tf.count"]][[x]]),
      stringsAsFactors = FALSE
    )
  }}))


# Plot the ranking metric result in a 2D plane

# After getting the active pathway ranking metric result you can display them in your preferred format.
# Here we will plot them in a 2D plane where x-axis denotes the number of downstream transcription factors and 
# y-axis denotes the activity score for each pathway.

# display_pathway_ranking_metric(pathway.ranking.metric = cond2_activePath_rankmetric)
#> [1] "ROR1_LEC -- result plotting done!!"
#To separate the top ranked pathways we can do this
abline(v=45, h=0.2, lty=2, col="black")



mean_cond1_dt <- aggregate(list(cond1_activ_dt$activity_score, cond1_activ_dt$tf_count), 
                 by = list(cond1_activ_dt$cond, cond1_activ_dt$tf), mean)
colnames(mean_cond1_dt) <- c("cond", "TF", "mean_activity", "mean_TFcount")
mean_cond1_dt <- mean_cond1_dt[order(mean_cond1_dt$mean_activity, decreasing=TRUE),]
head(mean_cond1_dt)
# > head(mean_cond1_dt)
# Group.1 Group.2 mean_activity mean_TFcount
# 10    LumB    CD40     0.7500000            4
# 36    LumB   ITGAL     0.7083333            2
# 6     LumB   CD247     0.7023810            7
# 7     LumB    CD3D     0.7023810            7
# 8     LumB    CD3G     0.7023810            7
# 12    LumB    CD48     0.7023810            7

mean_cond2_dt <- aggregate(list(cond2_activ_dt$activity_score, cond2_activ_dt$tf_count), 
                           by = list(cond2_activ_dt$cond, cond2_activ_dt$tf), mean)
colnames(mean_cond2_dt) <- c("cond", "TF", "mean_activity", "mean_TFcount")
mean_cond2_dt <- mean_cond2_dt[order(mean_cond2_dt$mean_activity, decreasing=TRUE),]
head(mean_cond2_dt)
# Group.1 Group.2 mean_activity mean_TFcount
# 8     LumA    CDH2     0.6666667          2.0
# 24    LumA   GFRA1     0.5000000          1.5
# 7     LumA    CD63     0.3333333          2.0
# 15    LumA   EFNB1     0.3333333          2.0
# 19    LumA   ERBB2     0.3333333          2.0
# 25    LumA    GPC1     0.3333333          2.0


common_dt <- merge(mean_cond1_dt, mean_cond2_dt, by=c("TF"))

p <- ggplot(data=common_dt, aes(x = mean_activity.x, y = mean_activity.y, label=TF))+
     labs(x=paste0(cond1, " mean activity"),
          y=paste0(cond2, " mean activity"),
     title ="TF mean activity")+
  geom_point(color = "red")+
  geom_text_repel() + 
  theme_bw(base_size = 12, base_family = "") + 
  theme(plot.title = element_text(hjust=0.5, face = "bold"),
        plot.subtitle = element_text(hjust=0.5),
        legend.key = element_blank(), 
        axis.line=element_line(),
        axis.text=element_text(size=12),
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_blank(), 
        panel.grid = element_blank())

outFile <- file.path(outFolder, paste0(cond1, "_", cond2,"_intersectTFs_activity.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

     
p1 <- ggplot(data=mean_cond1_dt, aes(x = mean_activity, y = mean_TFcount, label=TF))+
  labs(x=paste0("mean activity"),
       y=paste0("TF count"),
       title =paste0("TF mean activity and count (", cond1, ")"))+
  geom_point(color = "red")+
  geom_text_repel() + 
  theme_bw(base_size = 12, base_family = "") + 
  theme(plot.title = element_text(hjust=0.5, face = "bold"),
        plot.subtitle = element_text(hjust=0.5),
        legend.key = element_blank(), 
        axis.line=element_line(),
        axis.text=element_text(size=12),
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_blank(), 
        panel.grid = element_blank())

outFile <- file.path(outFolder, paste0(cond1, "_TF_count_vs_activity.", plotType))
ggsave(p1, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



p2 <- ggplot(data=mean_cond2_dt, aes(x = mean_activity, y = mean_TFcount, label=TF))+
  labs(x=paste0("mean activity"),
       y=paste0("TF count"),
       title =paste0("TF mean activity and count (", cond2, ")"))+
  geom_point(color = "red")+
  geom_text_repel() + 
  theme_bw(base_size = 12, base_family = "") + 
  theme(plot.title = element_text(hjust=0.5, face = "bold"),
        plot.subtitle = element_text(hjust=0.5),
        legend.key = element_blank(), 
        axis.line=element_line(),
        axis.text=element_text(size=12),
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_blank(), 
        panel.grid = element_blank())

outFile <- file.path(outFolder, paste0(cond2, "_TF_count_vs_activity.", plotType))
ggsave(p2, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



stop("***Done\n")
stop("--ok\n")

# To save the active pathway paths in csv format
# After getting the active pathway paths, you can save the data in csv file format
# where each row will denote a path and all the paths starting from the same source comprises of a pathway of the source protein.

out_dt <- lapply(unlist(cond1_activePath$$ROR1_LEC, recursive = F, use.names = F), write, "ROR1.active.pathway.csv", append=T, ncolumns=10)









# gmt_file <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_file <- "c5.go.bp.v7.4.symbols.gmt"

nTop_connect <- 5


pvalThresh_plot <- 0.05
fcThresh_plot <- 0.1

# look at the size of aracne regulon

min_mod_size <- min(unlist(lapply(regulonbrca, function(x)length(x[["tfmode"]]))))
# 2 -> 47 modules
# min_mod_size <- 10












# Pre-process the query data

# The query data is already in CPM and log2 normalized form. 
# Here, we will use the expression cutoff as 1.8 of the query data.

ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)

# Identify active pathway paths of the processed query data

# Here we will use the background pathway.path data
# to get the active pathway paths of the processed query data.

ROR1.active.pathway <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)


ROR1.active.pathway2 <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.data)

# To save the active pathway paths in csv format

# After getting the active pathway paths, you can save the data in csv file format
# where each row will denote a path and all the paths starting from the same source comprises of a pathway of the source protein.

lapply(unlist(ROR1.active.pathway$ROR1_LEC, recursive = F, use.names = F), write, "ROR1.active.pathway.csv", append=T, ncolumns=10)

# Get active pathway ranking metric (i.e., activity score and number of downstream transcription factors)

# Here we will use the ROR1.active.pathway and ROR1.processed.data data sets to get the acitve pathway ranking metric. Also we will use a high expression threshold (here we will use 7) for the processed query data.

ROR1.active.pathway.ranking.metric <- get_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)

# Plot the ranking metric result in a 2D plane

# After getting the active pathway ranking metric result you can display them in your preferred format.
# Here we will plot them in a 2D plane where x-axis denotes the number of downstream transcription factors and 
# y-axis denotes the activity score for each pathway.

display_pathway_ranking_metric(pathway.ranking.metric = ROR1.active.pathway.ranking.metric)
#> [1] "ROR1_LEC -- result plotting done!!"
#To separate the top ranked pathways we can do this
abline(v=45, h=0.2, lty=2, col="black")




