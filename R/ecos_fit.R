#' @title Fitting Grade of Membership model for clustering into communities
#'
#' @description A grade of membership model for clustering a site by features
#' data - which could be presence-absence or counts abundance data of species
#' in the sites.
#' 
#' @param dat input data matrix of samples along rows and features along sites,
#'            with each entry a 0/1 signififying presence/absence or counts of
#'            abundances. 
#' @param K The number of clusters to fit
#' @param tol The tolerance level of the model.
#' @param num_trials Number of EM runs from different starting points. This is'
#'                   key for picking the best fit model across multiple runs.
#' @param fit_control The control parameters for the model.
#' 
#' @return Returns a model fit with \code{omega} as cluster membership 
#'         probabilities and \code{theta} as cluster feature matrix.
#'         
#' @importFrom maptpx topics
#' @importFrom methClust meth_topics
#' 
#' @examples 
#' 
#' data("himalayan_birds")
#' species_abundance_counts <- t(exprs(himalayan_birds));
#' fit <- ecostructure_fit(species_abundance_counts, K = 2, tol = 0.1)
#' species_pa_counts <- species_abundance_counts
#' species_pa_counts[species_pa_counts >=1] = 1
#' fi2 <- ecos_fit(species_pa_counts, K = 2, tol = 0.1)
#' 
#' @export



ecos_fit <- function(dat,
                     max_dat = NULL,
                     K,
                     tol = 0.1,
                     num_trials = 1,
                     fit_control = list()){
  
  row_names <- rownames(dat)
  if(length(row_names) != dim(dat)[1]){
    warning("row names not provided, or not proper, using fake rownames")
    rownames(dat) <- 1:dim(dat)[1]
  }
  if(all(dat - floor(dat) != 0)){
    stop("The matrix dat must be a matrix of integers - a binary or a counts matrix")
  }
  
  ids_na <- which(is.na(dat))
  if(length(ids_na) > 0){
    warning("NAs in dat matrix: replaced by 0")
    dat[is.na(dat)] = 0
  }
  
  dat_by_2 <- dat %% 2
  if(all(dat %% 2 - dat == 0)){
    if(max(dat) == 1 & min(dat) == 0){
      message("Binary matrix input: Fitting the Binomial Grade of Membership 
              model.")
      max_dat <- matrix(1, dim(dat)[1], dim(dat)[2])
    }
  }
  
  if(!is.null(max_dat)){
    if(dim(max_dat)[1] != dim(dat)[1] | dim(max_dat)[2] != dim(dat)[2]){
      stop("dimensions of max_dat should match with dat")
    }
    unmeth <- max_dat - dat
    if(all(unmeth - floor(unmeth) != 0) | min(unmeth) < 0){
      stop("max_dat is present: but data entries of max_dat must be integers 
           bigger than the corresponding entries in the dat matrix.")
    }
    meth <- dat
    
    fit_control_default <-  list(shape=NULL,
                                 initopics=NULL,
                                 ord=TRUE,
                                 verb=1,
                                 sample_init = TRUE,
                                 NUM_INDICES_START = floor(dim(dat)[1]*0.75),
                                 use_squarem=FALSE)
    
    fit_control <- modifyList(fit_control_default, fit_control)
    
    fits_list <- list()
    L_array <- c()
    
    for(m in 1:num_trials){
      counter = 0
      while(counter != 1){
        tmp <- try(do.call(methClust::meth_topics,
                           append(list(meth = meth, 
                                       unmeth = unmeth, 
                                       K = K, 
                                       tol = tol), 
                                  fit_control)), TRUE)
        if(!inherits(tmp, "try-error")){
          counter = 1
        }else{
          counter = 0
        }
      }
      fits_list[[m]] <- tmp
      cat("We are at iteration", m, "\n")
    }
    
    loglik <- unlist(lapply(fits_list, function(x) return(x$L)))
    ids <- which.max(loglik)
    topic_clus <- fits_list[[ids]]
    ll <- list("omega" = topic_clus$omega,
               "theta" = topic_clus$freq,
               "L" = topic_clus$L)
  }
  
  if(is.null(max_dat) && max(dat) > 1){
    
    fit_control_default <-  list(shape = NULL, initopics = NULL, bf = TRUE, 
                                 kill = 2, ord = TRUE, verb = 1, admix = TRUE, 
                                 tmax = 1000)
    fit_control <- modifyList(fit_control_default, fit_control)
    
    
    topic_clus <- do.call(CountClust::FitGoM,
                          list(data = dat,
                               K = K,
                               tol = tol,
                               num_trials = num_trials,
                               options = "BIC",
                               path_rda = NULL,
                               control = fit_control))

    ll <- list("omega" = topic_clus$fit$omega,
               "theta" = topic_clus$fit$theta,
               "BIC" = topic_clus$BIC)
  }
  return(ll)
}