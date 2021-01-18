#' @title Null models for ecostructure_fit clustering validation
#'
#' @description Generate randomized matrix of counts 
#' given the observed data matrix for a particular null model, 
#' run \code{ecos_fit} on these null matrices and compare the fit
#' on null model data with that on the observed data. Used for validating 
#' the clustering. 
#'
#' @param counts The counts matrix (N x G): N- the number of samples (sites), 
#'               G- number of features (bird species)
#' @param K  The number of clusters to fit
#' @param tol The tolerance limit of the \code{ecos_fit}.
#' @param null.model The type of nullmodel used (similar to the
#'                   picante::randomizeMatrix() function argument in 
#'                   picante package)
#' @param iter_fill The number of swaps/fills in each randomized matrix build
#' @param iter_randomized The number of randomization matrices generated
#' @param plot If TRUE, plots density function of log Bayes factor over
#'             the randomized iterations
#'
#' @return  Returns a list with
#'        \item{BF.obs}{log Bayes Factor for the observed counts with 
#'        K=2 against the null with no clusters}
#'        \item{BF.rand}{a vector of log BF for each randomized count 
#'        matrix with K=2 against the null with no clusters}
#'        \item{pval}{the p-value of the observed log Bayes factor against the
#'        ones from randomized matrices}
#'
#' @importFrom  picante randomizeMatrix
#' @import slam
#' @importFrom CountClust FitGoM
#' @importFrom stats dmultinom density
#'
#' @examples
#'
#' data("himalayan_birds")
#' species_abundance_counts <- t(exprs(himalayan_birds));
#' out <- ecos_nullmodel(species_abundance_counts, K=2, 
#'                 iter_randomized=5, option = "BF")
#' out2 <- ecos_nullmodel(species_abundance_counts, K=2, 
#'                 iter_randomized=5, ption = "BIC")              
#'   
#' @export
              
ecos_nullmodel <- function(counts,
                           K,
                           tol=0.1,
                           null.model=c("frequency", "richness",
                                        "independentswap", "trialswap"),
                           iter_fill=100,
                           iter_randomized=30,
                           option,
                           plot=TRUE)
{
    counts = as.matrix(counts)
    if(missing(option)){
      option = "BIC"
    }
    if(option != "BF" & option != "BIC"){
      stop("The input argument option can only be BF (Bayes factor) or BIC")
    }
    counts <- counts[which(rowSums(counts) > 0), ];
    bf_gom_rand <-
        unlist(lapply(1:iter_randomized,
                      function(n)
                      {
                          rand_counts <-
                              picante::randomizeMatrix(counts,
                                                       null.model=null.model,
                                                       iterations=iter_fill);
                          suppressMessages(
                              topics_rand <- CountClust::FitGoM(rand_counts, 
                                                                K=K, 
                                                                tol=tol))
                              if(option == "BF"){
                                bf_gom  <- topics_rand$fit$BF
                              }else{
                                bf_gom <- topics_rand$BIC
                              }
                          return(bf_gom)
                      }))

    topics_obs <- suppressMessages(CountClust::FitGoM(counts, K=K, tol=tol));
    
    if(option == "BF"){
        bf_gom_obs  <- topics_obs$fit$BF
    }else{
        bf_gom_obs <- topics_obs$BIC
    }

    pval_bf_gom <- length(which(bf_gom_obs < bf_gom_rand))/iter_randomized;

    if(option == "BIC"){
      ll <- list("BIC.obs"=bf_gom_obs,
                 "BIC.rand"=bf_gom_rand,
                 "pval"=pval_bf_gom)
      if(plot){
        plot(density(bf_gom_rand), col="blue",
             main="BIC density over null matrices",
             xlab="BIC", ylab="density",
             xlim=c(min(bf_gom_obs,bf_gom_rand)-100, max(bf_gom_obs,bf_gom_rand)+100))
        abline(v=bf_gom_obs, col="red")
      }
      return(ll)
    }
    if(option == "BF"){
      ll <- list("BF.obs"=bf_gom_obs,
                 "BF.rand"=bf_gom_rand,
                 "pval"=pval_bf_gom)
      if(plot){
        plot(density(bf_gom_rand), col="blue",
             main="log-BF density over null matrices",
             xlab="log-BF", ylab="density",
             xlim=c(min(bf_gom_obs,bf_gom_rand)-100, max(bf_gom_obs,bf_gom_rand)+100))
        abline(v=bf_gom_obs, col="red")
      }
      return(ll)
    }
}


CheckCounts <- function(counts){
    if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
    if(is.null(dimnames(counts)[[1]]))
    {dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
    if(is.null(dimnames(counts)[[2]]))
    { dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
    empty <- slam::row_sums(counts) == 0
    if(sum(empty) != 0){
        counts <- counts[!empty,]
        cat(paste("Removed", sum(empty), "blank documents.\n")) }
    return(slam::as.simple_triplet_matrix(counts))
}
