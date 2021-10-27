require(igraph)
require(dplyr)
require(foreach)
require(doParallel)

getWeightedNetwork<-function(ppi,weight)
{
  if(length(V(ppi))!=length(weight))
    stop("Length of weights and nodes is not matching")
  
  V(ppi)$weight=weight[V(ppi)$name];
  
  return(ppi)
}

#remove nodes who have especially high degree(default remove nodes whose neighbor counts more than half of nodes)
removeHighDegNodes<-function(ppi,degratio=0.5)
{
  
  degrees<-degree(graph = ppi);
  nodes<-length(degrees)
  keep<-names(degrees[degrees<nodes*degratio])
  ppi.ret<-induced_subgraph(ppi,V(ppi)[keep]);
  
  return(ppi.ret)
}

findLocalMax<-function(weightedPPI, minNeighbors=0)
{
  nodes<-length(V(weightedPPI))
  isChecked<-rep(F,nodes);
  names(isChecked)<-names(nodes);
  seed<-round(runif(1,min=1,max=nodes))
  local.max<-numeric();
  
  while(sum(isChecked)!=nodes)
  {
    cat("Now have checked",sum(isChecked),"of",nodes, "nodes,", length(local.max),"local max node(s) founded.\r")
    while(isChecked[seed])
      seed<-round(runif(1,min=1,max=nodes))
    #isChecked[seed]<-T;
    if(length(neighbors(weightedPPI,seed))<minNeighbors) #if node have too little neighbors, skip it.
    {
      isChecked[seed]<-T
      next
    }
    
    maxneighbor<-findMaxNeighbor(weightedPPI,seed);
    if(isChecked[maxneighbor])# if node's max neighbor is checked, means the local max of the node has been found, skip it.
    {
      isChecked[seed]<-T
      next
    }
    isChecked[maxneighbor]<-T;

    while(seed != maxneighbor) #find local max node
    {
      isChecked[seed]<-T
      seed <- maxneighbor;
      maxneighbor<-findMaxNeighbor(weightedPPI,seed)
      if(isChecked[maxneighbor])
        {break}
      else
        {isChecked[maxneighbor]=T}
    }
    if(seed == maxneighbor)
    {
      local.max<-c(local.max,seed)
      isChecked[neighbors(weightedPPI,seed)]<-T
    }

    
  }
   local.max<-V(weightedPPI)[local.max]
   
   lm<-data.frame(g=local.max$name,w=local.max$weight)
   lm<-as.character(arrange(lm,desc(w))$g)
   
   return(V(weightedPPI)[lm]$name)
  
  #return(V(weightedPPI)[local.max])
}

findMaxNeighbor<-function(weightedPPI,seed)
{
  neis<-c(neighbors(weightedPPI,seed),V(weightedPPI)[seed]);
  
  return(neis[which.max(V(weightedPPI)[neis]$weight)])
}

findModules<-function(weightedPPI, localMaxNodes, lm.topot.least=0.05, topot.least=0.2, size.least=20)
{
  lm.topot.least<-max(V(weightedPPI)[localMaxNodes]$weight)*lm.topot.least
  localMaxNodes<-localMaxNodes[V(weightedPPI)[localMaxNodes]$weight>lm.topot.least]

  modules<-foreach(i=1:length(localMaxNodes)) %do%
    {
      #cat("Now processing:",localMaxNodes[i]$name,"(",i,"/",length(localMaxNodes),")\n")
      cat("Now processing:",localMaxNodes[i],"(",i,"/",length(localMaxNodes),")\n")
      modulExpand(weightedPPI,localMaxNodes[i], wl=topot.least)
    }
  
  return(modules[sapply(modules, function(x) length(x)>size.least)])
}

modulExpand<-function(weightedPPI,seed,wl)
{
  seed<-V(weightedPPI)[seed];
  visited<-rep(F,length(V(weightedPPI)));
  queue<-numeric();
  
  
  bfs<-function(weightedPPI,seed,wl)
  {
    visited[seed]<-T;
    queue<-seed;
    module<-seed;
    while(length(queue)!=0)
    {
      cat("Module expanding,",length(module),"nodes detected.\r")
      neis<-neighbors(weightedPPI,queue[1])
      this<-queue[1];
      queue<-queue[-1];
      neis<-neis[!visited[neis]]; #nodes not searched
      neis<-neis[neis$weight<this$weight]; #topology potential decrease
      visited[neis]<-T;
      neis<-neis[neis$weight>(seed$weight)*wl]; #topology potential not too small
      queue<-c(queue,neis);
      module<-c(module,neis);
    }
    cat("\n")
    return(module$name);
  }
  
  return(bfs(weightedPPI,seed,wl))
}



# modulExpand<-function(weightedPPI, seed, wl=1)
# {
#   seed<-V(weightedPPI)[seed]
#   neis<-neighbors(weightedPPI,seed)
#   seed.weight<-seed$weight;
#   keeps<-neis[V(weightedPPI)[neis]$weight<seed.weight]
#   keeps<-keeps[V(weightedPPI)[keeps]$weight>wl]
# 
#   if(length(keeps)>0)
#     keeps<-foreach(i=1:length(keeps),.combine=c,.inorder = F) %do%
#     {
#       modulExpand(weightedPPI,keeps[i])
#     }
#   return(c(seed,keeps))
#   
# }



moduleMerge<-function(mods,maxOverlap=0.5)
{


  isFinished<-F
  
  while (!isFinished) 
  {
    mods.ord<-order(sapply(mods,length),decreasing = T)
    mods<-mods[mods.ord]
    
    mods.count<-length(mods);
    
    mergeMat<-foreach(i=1:mods.count,.combine = cbind) %do%
      {
        sapply(mods,FUN = function(x,this){if(length(intersect(x,this))>maxOverlap*min(length(x),length(this)))
                                                return(T)
                                           else
                                                return(F)},this=mods[[i]])
      }
  
    mergeMat[upper.tri(mergeMat,diag = T)]<-F
    
    isMerged<-!apply(mergeMat, 1, any)
    isKept<-(apply(mergeMat,2,any) | isMerged)
    diag(mergeMat)<-isKept
    isFinished<-all(isMerged)
    mergeMat<-mergeMat[,isKept]
    
    mods_n<-foreach(i=1:sum(isKept)) %do%
      {
        Reduce(union,mods[mergeMat[,i]])
      }
    
    mods<-mods_n
    
  
  }

  return(mods)
}
  



writeModToCSV<-function(mods,file="")
{
  mod.len<-sapply(mods,function(x) length(x))
  maxl<-max(mod.len)
  mod.corrected<-foreach(i=1:length(mods),.combine = cbind) %do%
    {
      c(mods[[i]],rep(NA,maxl-length(mods[[i]])))
    }
  
  if(nchar(file)!=0)
    write.csv(x=mod.corrected,file=file,na = "")
  
  
  return(as.data.frame(mod.corrected))
}
