#' @title Collapse counts data using phylogenetic tree hierarchy of features
#'
#' @description To club the species in the columns of the
#'  counts data based on a phylogenetic tree hierarchy
#'
#' @param counts The abundance counts matrix with samples/sites along the rows and the
#'               species/features along the columns.
#' @param tree The phylogenetic tree structure on the species or the column features
#'             of the counts matrix.
#' @param collapse_at The point at which to collapse the tree
#'
#' @return Returns an counts matrix with the features being the collapsed species
#' from the original data based on the tree structure.
#'
#' @import ape
#' @import phytools
#' @export

collapse_counts_by_phylo <- function(counts, tree, collapse_at){
  root_node <- length(tree$tip.label) + 1
  root_age <- ape::branching.times(tree)[names(ape::branching.times(tree)) == root_node]
  trees_at_slice <- phytools::treeSlice(tree, root_age - collapse_at)
  counts_at_slice <- as.data.frame(counts)
  for( i in 1:length(trees_at_slice)){
    new.column <- as.data.frame(rowSums(counts_at_slice[,trees_at_slice[[i]]$tip.label]))
    colnames(new.column) <- trees_at_slice[[i]]$tip.label[1]
    drops <- trees_at_slice[[i]]$tip.label
    counts_at_slice <- counts_at_slice[,!(names(counts_at_slice) %in% drops)]
    counts_at_slice <- cbind(counts_at_slice, new.column)
  }
  counts_at_slice <- as.matrix(counts_at_slice)
  return(counts_at_slice)
}
