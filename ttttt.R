
## to not have too many genes to plot, add threshold of fc

module_genes=module_genes
prot_dt=t12_dt
meanExprThresh=0.01
absLog2fcThresh=0.5
cond1_s=cond1_samps
cond2_s=cond2_samps


moi <- "M5"
moi_genes <- cem_cond12@module$genes[cem_cond12@module$modules== moi]

plot_mymodule(module_genes=moi_genes, prot_dt=cond12_dt,
              cond1_s=cond1_samps, cond2_s=cond2_samps,
              meanExprThresh=1, absLog2fcThresh=0.5 
              ) 

plot_mymodule <- function(module_genes, prot_dt, cond1_s, cond2_s, meanExprThresh, absLog2fcThresh) {
  
  stopifnot(length(module_genes) > 0)
  stopifnot(module_genes %in% rownames(prot_dt))
  stopifnot(cond1_s %in% colnames(prot_dt))  
  stopifnot(cond2_s %in% colnames(prot_dt))  
  
  # keep the ones that pass threshold
  tmp_dt <- prot_dt[paste0(module_genes),]
  t1_dt <- tmp_dt[,paste0(cond1_s)]
  t2_dt <- tmp_dt[,paste0(cond2_s)]
  stopifnot(dim(t1_dt) > 0)
  stopifnot(dim(t2_dt) > 0)
  
  t12_dt <- cbind(t1_dt, t2_dt)
  
  c1 <- rowMeans(t1_dt) >= meanExprThresh & rowMeans(t2_dt) >= meanExprThresh
  stopifnot(names(c1) == module_genes)
  
  c12 <- c1 & (abs(log2(rowMeans(t1_dt)/rowMeans(t2_dt))) >= absLog2fcThresh)
  stopifnot(names(c12) == module_genes)
  
  # get all pairs
  cat(paste0("... genes passing Expr and FC thresholds:\t", sum(c12), "/", length(c12), "\n"))
  
  module_genes <- module_genes[c12]
  stopifnot(length(module_genes) == sum(c12))
  
  module_dt <- as.data.frame(t(combn(module_genes, 2)))
  colnames(module_dt) <- c("gene1", "gene2")
  
  stopifnot(module_dt$gene1 %in% rownames(t12_dt))
  stopifnot(module_dt$gene2 %in% rownames(t12_dt))
  stopifnot(module_dt$gene1 %in% rownames(t1_dt))
  stopifnot(module_dt$gene2 %in% rownames(t1_dt))
  stopifnot(module_dt$gene1 %in% rownames(t2_dt))
  stopifnot(module_dt$gene2 %in% rownames(t2_dt))
  
  i=1
  module_dt$corr <- foreach(i = 1:nrow(module_dt), .combine='c') %dopar% {
    g1 <- module_dt$gene1[i]
    g2 <- module_dt$gene2[i]
    cor(t12_dt[g1,], t12_dt[g2,])
  }
  module_dt$pair_label <- paste0(module_dt$gene1, "_", module_dt$gene2)
  pair2corr <- setNames(module_dt$corr, module_dt$pair_label)
  
  # average fold change between the 2 conditions
  module_genes_dt <- data.frame(gene = module_genes, stringsAsFactors = FALSE)
  module_genes_dt$log2_fc <- foreach(i = 1:nrow(module_genes_dt), .combine='c') %dopar% {
    g <- module_genes_dt$gene[i]
    log2(median(t1_dt[g,])/median(t2_dt[g,]))
  }
  g2fc <- setNames(module_genes_dt$log2_fc, module_genes_dt$gene)
  
  ig_obj <- graph_from_data_frame(module_dt, directed=FALSE)
  net_obj <- intergraph::asNetwork(ig_obj)
  m <- network::as.matrix.network.adjacency(net_obj)
  plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m,NULL))
  colnames(plotcord) <- c("X1", "X2")
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj,"vertex.names"))
  
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[,  2], ])
  colnames(edges)[colnames(edges) == "X1.1"] <- "Y1"
  colnames(edges)[colnames(edges) == "X2.1"] <- "Y2"
  
  stopifnot(edges$vertex.names.1 %in% module_dt$gene2)
  stopifnot(edges$vertex.names %in% module_dt$gene1)
  
  edges$pair_label <- paste0(edges$vertex.names, "_", edges$vertex.names.1)
  stopifnot(edges$pair_label %in% names(pair2corr))
  mycorr <- pair2corr[paste0(edges$pair_label)]
  edges$corr <- mycorr
  stopifnot(!is.na(edges$corr))
  
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  stopifnot(plotcord$vertex.names %in% names(g2fc))
  plotcord$fc <- g2fc[paste0(plotcord$vertex.names)]
  plotcord$Degree <- 1
  
  
  size_min = 5
  size_max = 15
  plotcord$fc_size <- (size_max-size_min) * (abs(plotcord$fc) - min(abs(plotcord$fc)))/(max(abs(plotcord$fc))-min(abs(plotcord$fc))) + size_min

  pl <- ggplot(plotcord) + 
    geom_segment(data = edges, aes_(x = ~X1, y = ~X2, xend = ~Y1, yend = ~Y2, color=~corr), 
                 size = 0.5, alpha = 0.5) + 
    geom_point(aes_(x = ~X1, y = ~X2,
                    # size = ~Degree, 
                    size = ~fc, 
                    fill=~fc), shape=21) +
    scale_color_gradient2(  low = muted("blue"),
                            mid = "white",
                            high = muted("red")) +
    scale_fill_gradient2(  low = muted("blue"),
                           mid = "white",
                           high = muted("red"))+
    guides(size=FALSE)+
    geom_label_repel(aes_(x = ~X1, 
                          y = ~X2,
                          size = ~fc_size, 
                          # color = ~fc
                          # color="black",
                          label = ~vertex.names), 
                     box.padding = unit(1,"lines"))+
    theme_bw(base_size = 12, base_family = "") + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          legend.key = element_blank(), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_blank(), 
          panel.grid = element_blank())
  
  return(pl)
}
  
  
  
  
  
