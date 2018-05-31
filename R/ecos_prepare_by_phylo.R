#' @title Prepare counts data for phylogenetically collapsed species
#'
#' @description Starting from a taxonomic species abundance matrix, produces
#' a new matrix, with features now collpased vesions of the species where the
#' collapsing is done by a phylogenetic tree structure. 
#'
#' @param counts The abundance counts matrix with samples/sites along the rows 
#'               and the species/features along the columns.
#' @param tree The phylogenetic tree structure on the species or the 
#'             column features of the counts matrix. Should be of the form
#'             \code{.tre}. 
#' @param collapse_at The point on the tree at which to collapse the tree
#'
#' @return Returns a counts matrix with the features being the 
#' phylogenetically collapsed species from the original data based on the 
#' tree structure.
#'
#' @importFrom ape branching.times
#' @importFrom phytools treeSlice
#' @export

ecos_prepare_by_phylo <- function(dat, tree, collapse_at){
  root_node <- length(tree$tip.label) + 1
  root_age <- ape::branching.times(tree)[names(ape::branching.times(tree)) == root_node]
  trees_at_slice <- phytools::treeSlice(tree, root_age - collapse_at)
  counts_at_slice <- as.data.frame(dat)
  num_in_groups <- c()
  for( i in 1:length(trees_at_slice)){
    slice_labels <- trees_at_slice[[i]]$tip.label
    slice_labels_filtered <- intersect(colnames(dat), slice_labels)
    if(length(slice_labels_filtered) == 0){
      next
    }else if(length(slice_labels_filtered) == 1){
      new.column <- as.data.frame(dat[, slice_labels_filtered])
      colnames(new.column) <- slice_labels_filtered
      drops <- trees_at_slice[[i]]$tip.label
      counts_at_slice <- cbind(counts_at_slice, new.column)
      num_in_groups <- c(num_in_groups, 1)
    }else{
      new.column <- as.data.frame(rowSums(dat[,slice_labels_filtered]))
      colnames(new.column) <- slice_labels_filtered[1]
      drops <- trees_at_slice[[i]]$tip.label
      counts_at_slice <- cbind(counts_at_slice, new.column)
      num_in_groups <- c(num_in_groups, length(slice_labels_filtered))
    }
    counts_at_slice <- as.matrix(counts_at_slice)
  }
  counts_at_slice_2 <- counts_at_slice[, (dim(dat)[2]+1):dim(counts_at_slice)[2]]
  return(list("outdat" = counts_at_slice_2, "num_groups" = num_in_groups))
}
