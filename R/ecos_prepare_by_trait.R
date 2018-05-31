#' @title Prepare counts data for trait-collapsed species
#' 
#' @description Starting from a taxonomic species abundance matrix, produces
#' a new matrix, with features now collpased vesions of the species where the
#' collapsing is done by trait level information. Binning is done on the species
#' based on their trait level similarity- where similarity is determined by
#' a dendogram produced by hierarchical clustering. 
#'
#' @param counts The abundance counts matrix with samples/sites along the rows 
#'               and the species/features along the columns.
#' @param traits A distance metric between species - a square matrix with
#'               species along rows and columns and the diagonal element equal to
#'               0.  
#' @param prop_div The proportion of the original counts matrix diversity at which
#'               to cut the dendrogram.
#'
#' @return The output is a counts matrix with the columns corresponding to the 
#' species bins - each bin designated here by the first species entry into the bin. 
#' 
#' @examples
#'
#'data(himalayan_birds)
#'species_metadata <- pData(featureData(himalayan_birds))
#'taxonomic_counts <- t(exprs(himalayan_birds))
#'bill_traits <- as.matrix(dist(scale(species_metadata[,c(1:3)])))
#'counts_bill_traits <- ecos_prepare_by_trait(counts = taxonomic_counts, 
#'                                             traits = bill_traits, 
#'                                             prop_div=0.3)
#'                                                    
#' @export


ecos_prepare_by_trait <- function(counts, traits, prop_div = 0.3){
  
  if(any(traits - t(traits) != 0)){
    stop("traits matrix must be symmetric")
  }
  
  if(dim(counts)[2] != dim(traits)[1]){
    stop("the number of species in counts do not match with that in traits")
  }
  
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
