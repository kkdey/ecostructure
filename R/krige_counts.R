#' @title Kriging smoothed counts data preparation
#'
#' @description Performs interpolation of 0 counts (treated as missing data) using
#' kriging over an ordered feature information \code{order}.
#'
#' @param counts The abundance counts data matrix with samples along the rows and the
#'               species along the columns.
#' @param order  The ordered feature information used for computing distance information
#'               in kriging.
#' @param krige.control The control parameters for kriging.
#' @param jitter A jitter value to be added in case of any singularity problem due to lack of
#'               uniqueness of ordering metadata. Default is 10^-6. The user may tune it in case
#'               of a singularity problem found.
#'
#' @return Returns a counts data with the zero count cells being interpolated by
#' a simple kriging predictor (using \code{SpatialExtremes} package).
#'
#' @importFrom SpatialExtremes kriging
#' @export


krige_counts <- function (counts, order, krige.control = list(), jitter = 1e-06){

  order <- order + abs(runif(length(order), 0, jitter))
  krige.control.default <- list(cov.mod = "whitmat",
                                sill=0.01, range=sd(order)/50, smooth=0.3)
  krige.control <- modifyList(krige.control, krige.control.default)

  z <- matrix(0, nrow(counts), ncol(counts))
  colnames(z)<- colnames(counts)
  rownames(z)<- rownames(counts)
  for(i in 1:nrow(counts)){
    index1 <- which(counts[i,] == 0)
    index2 <- setdiff(1:ncol(counts), index1)
    out <- do.call(SpatialExtremes::kriging, append(list(data =counts[i, index2 ],
                                                         data.coord = order[index2],
                                                         krig.coord = order[index1]),
                                                    krige.control))
    z[i, index1] <- round(as.vector(out$krig.est))
    z[i, index2] <- counts[i, index2]
  }
  colnames(z)<- colnames(counts)
  rownames(z)<- rownames(counts)
  return(z)
}
