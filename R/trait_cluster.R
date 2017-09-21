#' @title Trait based Hierarchical clustering of species with varying proportional diversity
#'
#' @description Performs interpolation of 0 counts (treated as missing data) using
#' kriging over an ordered feature information \code{order}.
#'
#' @param counts The abundance counts data matrix with samples along the rows and the
#'               species along the columns. The sample names are provided as row names
#'               and the species names represent the column names.
#' @param traits  A matrix with species along the rows. Teh species names along the rows
#'               must match with the column names of the counts matrix input.
#' @param prop_div The proportion of the original counts matrix diversity at which
#'               to cut the dendrogram.
#'
#' @return The output is a counts matrix with the columns corresponding to the clusters and
#' the column name corresponding to the first of the species forming the cluster. This
#' matrix is taken as input to the CountClust:FitGoM function.
#'
#' @examples
#'
#' data <- get(load(system.file("extdata", "HimalayanBirdsData.rda", package = "ecostructure")))
#' species_metadata <- pData(featureData(data))
#' taxonomic_counts <- t(exprs(data))
#' bill_traits <- as.matrix(dist(scale(species_metadata[,c(1:3)])))
#' bill_trait_clust <- trait_cluster(counts = taxonomic_counts, traits = bill_traits, prop_div=0.3)
#'
#' @export




trait_cluster <- function(counts, traits, prop_div = 0.3){
  z <- round(dim(traits)[[1]]*prop_div)
  res.hc<-hclust(dist(traits))
  memb <- cutree(res.hc, k = z)
  cluster_counts <- matrix(0, dim(counts)[[1]], z)
  temp_counts <- matrix(0, dim(counts)[[1]], dim(counts)[[2]])
  rownames(cluster_counts) <- rownames(counts)
  rownames(temp_counts) <- rownames(counts)
  colnames(temp_counts) <- names(memb)[order(memb, decreasing=FALSE)]
  colnames(cluster_counts) <- colnames(counts)[1:z]
  for(m in 1:dim(counts)[[1]]){
    temp_counts[m,] <- as.numeric(counts[m,names(memb)[order(memb, decreasing=FALSE)]])
    for(i in 1:z){
      names <- names(memb)[memb==i]
      cluster_counts[m,i]<-sum(temp_counts[m,names])
      colnames(cluster_counts)[i]<-names[[1]]
    }
  }
  return(cluster_counts)
}
