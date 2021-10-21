require(scales)
require(ggrepel)

plot_mymodule <- function(module_genes, prot_dt, cond1_s, cond2_s, 
                          meanExprThresh, absLog2fcThresh,
                          size_min_nodes = 5,
                          size_max_nodes = 30,
                          size_min_labs = 14,
                          size_max_labs = 30,
                          size_min_edges=0.5,
                          size_max_edges=3) {
  
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
  
  if(length(module_genes) < 2) return(NA)
  
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
  
  ig_obj <- igraph::graph_from_data_frame(module_dt, directed=FALSE)
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
  
  plotcord$abs_fc <- abs(plotcord$fc)
  
  edges$corr_size_nodes <- (size_max_edges-size_min_edges) * (abs(edges$corr) - min(abs(edges$corr)))/
    (max(abs(edges$corr))-min(abs(edges$corr))) + size_min_edges
  
  if(nrow(edges) == 1) edges$corr_size_nodes[1] <-  0.5*(size_min_edges+size_max_edges)
  
  plotcord$fc_size_nodes <- (size_max_nodes-size_min_nodes) * (abs(plotcord$fc) - min(abs(plotcord$fc)))/
    (max(abs(plotcord$fc))-min(abs(plotcord$fc))) + size_min_nodes
  
  plotcord$fc_size_labs <- (size_max_labs-size_min_labs) * (abs(plotcord$fc) - min(abs(plotcord$fc)))/
    (max(abs(plotcord$fc))-min(abs(plotcord$fc))) + size_min_labs
  
  pl <- ggplot(plotcord) + 
    geom_segment(data = edges, aes_(x = ~X1, y = ~X2, xend = ~Y1, yend = ~Y2, color=~corr), #### CHANGE HERE COLOR DEFINITION
                 size = edges$corr_size_nodes, alpha = 0.5) + 
    geom_point(aes_(x = ~X1, y = ~X2,
                    # size = ~Degree, 
                    size = ~fc_size_nodes, 
                    fill=~fc), shape=21) +
    scale_color_gradient2(low="blue", mid="grey", high="red", 
                          limits = c(-1, 1))+
    scale_fill_gradient2(  low = muted("blue"),
                           mid = "white",
                           high = muted("red"))+
    guides(size=FALSE)+
    labs(color="Fold-change", fill="Corr.")+
    geom_label_repel(aes_(x = ~X1, 
                          y = ~X2,
                          size = ~fc_size_labs, 
                          fill = ~abs_fc,
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
plot_mymodule_v2 <- function(module_genes, prot_dt, de_dt, cond1_s, cond2_s, 
                             pvalThresh, absLogFCThresh,
                          size_min_nodes = 5,
                          size_max_nodes = 30,
                          size_min_labs = 14,
                          size_max_labs = 30,
                          size_min_edges=0.5,
                          size_max_edges=3) {
  
  stopifnot("adj.P.Val" %in% colnames(de_dt))
  stopifnot("logFC" %in% colnames(de_dt))
  gene2fc <- setNames(de_dt$logFC, de_dt$gene)
  gene2pval <- setNames(de_dt$adj.P.Val, de_dt$gene)
  stopifnot(module_genes %in% names(gene2fc))
  stopifnot(module_genes %in% names(gene2pval))
  c1 <- abs(gene2fc[module_genes]) >= absLogFCThresh
  c2 <- gene2pval[module_genes] <= pvalThresh
  c12 <- c1 & c2
  stopifnot(names(c12) == module_genes)
  cat(paste0("... genes passing Expr and FC thresholds:\t", sum(c12), "/", length(c12), "\n"))
  
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
  
  # get all pairs
  module_genes <- module_genes[c12]
  stopifnot(length(module_genes) == sum(c12))
  
  if(length(module_genes) < 2) return(NA)
  
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
  
  ig_obj <- igraph::graph_from_data_frame(module_dt, directed=FALSE)
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
  plotcord$abs_fc <- abs(plotcord$fc)
  
  edges$corr_size_nodes <- (size_max_edges-size_min_edges) * (abs(edges$corr) - min(abs(edges$corr)))/
    (max(abs(edges$corr))-min(abs(edges$corr))) + size_min_edges
  
  if(nrow(edges) == 1) edges$corr_size_nodes[1] <-  0.5*(size_min_edges+size_max_edges)
  
  plotcord$fc_size_nodes <- (size_max_nodes-size_min_nodes) * (abs(plotcord$fc) - min(abs(plotcord$fc)))/
    (max(abs(plotcord$fc))-min(abs(plotcord$fc))) + size_min_nodes
  
  plotcord$fc_size_labs <- (size_max_labs-size_min_labs) * (abs(plotcord$fc) - min(abs(plotcord$fc)))/
    (max(abs(plotcord$fc))-min(abs(plotcord$fc))) + size_min_labs
  
  pl <- ggplot(plotcord) + 
    geom_segment(data = edges, aes_(x = ~X1, y = ~X2, xend = ~Y1, yend = ~Y2, color=~corr), #### CHANGE HERE COLOR DEFINITION
                 size = edges$corr_size_nodes, alpha = 0.5) + 
    geom_point(aes_(x = ~X1, y = ~X2,
                    # size = ~Degree, 
                    size = ~fc_size_nodes, 
                    fill=~fc), shape=21) +
    scale_color_gradient2(low="blue", mid="grey", high="red", 
                          limits = c(-1, 1))+
    scale_fill_gradient2(  low = muted("blue"),
                           mid = "white",
                           high = muted("red"))+
    guides(size=FALSE)+
    labs(color="Fold-change", fill="Corr.")+
    geom_label_repel(aes_(x = ~X1, 
                          y = ~X2,
                          size = ~fc_size_labs, 
                          fill = ~abs_fc,
                          # color="black",
                          label = ~vertex.names), 
                     box.padding = unit(1,"lines"))+
    theme_bw(base_size = 12, base_family = "") + 
    theme(plot.title = element_text(hjust=0.5, face = "bold"),
          plot.subtitle = element_text(hjust=0.5),
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          legend.key = element_blank(), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_blank(), 
          panel.grid = element_blank())
  
  return(pl)
}




