#' tSNE dimension reduction(2d) for cell-topic or region-topic distribution
#' @param input Data frame of cell-topic or region-topic distribution with cell ids/region names as rows, topics as columns
#' @param seed Make results reproducible
#' @param out_path_tsne File path for exported tSNE result
#' @param perplexity Param in Rtsne. Numeric, the perplexity parameter of tSNE algorithm  (should not be bigger than perplexity < (nrow(X) - 1)/3).(default:30)
#' @check_duplicates Param in Rtsne. Logical, checks whether duplicates are present. It is best to make sure there are no duplicates present and set to FALSE,
#' especially for large datasets. Recommended to set to FALSE when running Rtsne over regions. (default: TRUE)
#' pca Parameter in Rtsne. Logical; Whether an initial PCA step should be performed. Recommended to set to FALSE for our use. (default: TRUE)
run_tsne <- function(input,seed,...){
  library(Rtsne)
  set.seed(seed)
  TSNE <- Rtsne(input, pca=FALSE,...)
  rownames(TSNE$Y) <- rownames(input)
  colnames(TSNE$Y) <- paste0('tSNE', 1:ncol(TSNE$Y))
  #write.table(TSNE$Y,row.names=T,col.names=T,sep= "\t",quote=FALSE,file=out_path_tsne)
  #write_result(TSNE$Y, out_path_tsne)
  return(TSNE$Y)
}

#' Umap dimension reduction(2d) for cell-topic or region-topic distribution
#' @param input Data frame of cell-topic or region-topic distribution with cell ids/region names as rows, topics as columns
#' @param seed Make results reproducible
#' @param out_path_umap File path for exported Umap result
run_umap <- function(input,seed){
  library(umap)
  set.seed(seed)
  UMAP <- umap(input)
  rownames(UMAP$layout) <- rownames(input)
  colnames(UMAP$layout) <- paste0('UMAP', 1:ncol(UMAP$layout))
  #write.table(UMAP$layout,row.names=T,col.names=T,sep= "\t",quote=FALSE,file=out_path_umap)
  #write_result(UMAP$layout, out_path_umap)
  return(UMAP$layout)
}

#' density clustering for cell-topic or region-topic distribution ---- not suitable for extremely seperated clusters
#' @param input Data frame of cell-topic or region-topic distribution with cell ids/region names as rows, topics as columns
#' @param seed Make results reproducible
#' @param out_path_dclust File path for exported density clustering result

density_clustering <- function (input,seed){
  library(densityClust)
  set.seed(seed)
  library(Rtsne)
  DR <- Rtsne(input, check_duplicates = F, pca=F)
  DRdist <- dist(DR$Y)
  loginfo("Running density clustering:")
  dclust <- densityClust(DRdist,gaussian=T)
  #distance cut off at .......
  # rho is set to be 0.25 quantile; delta is set to be 0.997 quantile
  dclust <- findClusters(dclust, rho=quantile(dclust$rho, 0.25)[[1]], delta=quantile(dclust$delta, 0.997)[[1]])
  # Check thresholds
  #options(repr.plot.width=6, repr.plot.height=6)
  #plot(dclust$rho,dclust$delta,pch=20,cex=0.6,xlab='rho', ylab='delta')
  #points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
  #text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+1.5,labels=dclust$clusters[dclust$peaks])
  #abline(v=mean(dclust$rho))
  #abline(h=quantile(dclust$delta,0.95)[[1]])
  # Add cluster information
  densityClust <- dclust$clusters
  densityClust <- as.data.frame(densityClust)
  rownames(densityClust) <- rownames(input)
  colnames(densityClust) <- 'densityClust'
  densityClust[,1] <- as.factor(densityClust[,1])
  #write.table(densityClust,row.names=T,col.names=T,sep= "\t",quote=FALSE,file=out_path_dclust)
  #write_result(densityClust, out_path_dclust)
  return(densityClust)
}


#' louvain clustering (based on PhenoGraph algorithm) for cell-topic or region-topic distribution
#' @param input Data frame of cell-topic or region-topic distribution with cell ids/region names as rows, topics as columns
#' @param seed Make results reproducible
#' @param k Numeric, the number of nearest neighbours used in calculating the weightMatrix to represent distances between objects. (Default:15)
#' @param out_path_lclust File path for exported density clustering result

louvain_clustering <- function(input,seed,k){
  library(Rphenograph)
  loginfo("Running louvain clustering:")
  # use PhenoGraph algorithm to calculate the weightMatrix
  # and louvain clustering on the weighted graph
  graph_out <- Rphenograph(input, k = k)
  louvainClust <- as.data.frame(factor(membership(graph_out[[2]])))
  rownames(louvainClust) <- rownames(input)
  colnames(louvainClust) <- 'louvainClust'
  louvainClust[,1] <- as.factor(louvainClust[,1])
  #write.table(louvainClust,row.names=T,col.names=T,sep= "\t",quote=FALSE,file=out_path_lclust)
  #write_result(louvainClust, out_path_lclust)
  return(louvainClust)
}



#' helper function to export results
write_result <- function(dataframe, out_path){
  write.table(dataframe, row.names=T, col.names=T, sep="\t", quote=FALSE, file=out_path)
}