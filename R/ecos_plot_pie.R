#' @title Produces ecological pie chart representations of ecostructure_fit model
#'
#' @description Make a pie chart representation of the membership probabilities
#' \code{omega} output from the \code{ecostructure_fit} model. 
#'
#' @param omega Matrix of Cluster membership probabilities of each sample
#'              obtained from \code{ecostructure_fit()}. The row sums 
#'              sum to 1 for each sample.
#' @param coords a matrix of latitude and longitude entries for each site (row)
#'               of the \code{omega} matrix.This matrix has as many rows as 
#'               the \code{omega} matrix and two columns for the latitude and 
#'               longitude.
#' @param bgmap_path  The path to the shapefile to be used for plotting the 
#'                    map data. One can choose a shapefile of their own 
#'                    preference of resolution. If missing, we use a default 
#'                    background map - which is a coarse one.
#' @param adjust A Boolean, indicating whether to adjust for the richness cluster.
#'               Defaults to TRUE.
#' @param thresh a numeric value between 0 and 1 indicating the level at which
#'               to threshold the richness cluster. Sites represented by member-
#'               ship probability above this threshold of the richness cluster
#'               are removed.
#' @param long_lim The limits for the longitude in the map plot.
#'                 The default is \code{c(-180,180)} and removes Antarctica 
#'                 from the map.
#' @param lat_lim The limits for the latitude in the map plot.
#'                The default is \code{c(-60,90)} and removes Antarctica 
#'                from the map.
#' @param coastline_lwd The line width for the coastline plotted in the map.
#'                      Defaults to 10.
#' @param intensity The intensity of the colors for the pie chart representation.
#'                  Lies between 0 and 1. Defaults to 1.
#' @param radius The radius of the pie charts. Defaults to 0.5 
#' @param color the vector of colors - of same size or more than the number 
#'              of clusters. If missing, we use a default set of colors- 
#'              chosen so that they are distinct from one another.
#' @param pie_control The list of control parameters to be passed into the 
#'                    \code{add.pie} function of the package maptools. 
#' @param image_width The width of the image output. Defaults to 1000.
#' @param image_height The height of the image output. Defaults to 800. 
#' @param path The path where the output image is saved.
#' 
#' @return Returns a pie chart representation of the memberships obtained 
#'         from \code{ecostructure_fit} on a global species 
#'         presence-absence/counts data
#'         
#' @import sf
#' @importFrom mapplots add.pie
#' @importFrom scales alpha
#' 
#' @examples 
#'  data("australia_birds")
#'  data("australia_model")
#'  ecos_plot_pie(omega = australia_model$omega,
#'                      coords = australia_birds$latlong, 
#'                      long_lim = c(110,160),
#'                      lat_lim = c(-50,-10),
#'                      color= c("orange", "red", "yellow", "deepskyblue", 
#'                               "chartreuse", "blue"))
#'
#' 
#' 
#' @export


ecos_plot_pie = function(omega = NULL,
                         coords = NULL,
                         bgmap_path = NULL,
                         adjust = FALSE,
                         thresh = 0.7,
                         long_lim = c(-180,180),
                         lat_lim = c(-60,90),
                         coastline_lwd = 10,
                         intensity = 1,
                         radius = 0.5,
                         color = c("dodgerblue2","#E31A1C", "green4", "#6A3D9A","#FF7F00", "black","gold1","skyblue2","#FB9A99",
                                   "palegreen2", "#CAB2D6", "#FDBF6F",  "gray70", "khaki2", "maroon","orchid1","deeppink1",
                                   "blue1","steelblue4", "darkturquoise","green1","yellow4","yellow3", "darkorange4","brown",
                                   "red", "cornflowerblue", "cyan", "brown4", "burlywood", "darkgoldenrod1",
                                   "azure4", "green","deepskyblue","yellow", "azure1"),
                         pie_control = list(),
                         image_width = 1000,
                         image_height = 800,
                         path = "geostructure_plot.tiff"){
  
  if(is.null(coords)){
    if(is.null(rownames(omega))){
      stop("coords not provided, omega rownames do not have latitude longitude
           information either")
    }
    latlong_chars <- rownames(omega)
    coords <- cbind.data.frame(
      as.numeric(sapply(latlong_chars, function(x) strsplit(x, "_")[[1]][1])),
      as.numeric(sapply(latlong_chars, function(x) strsplit(x, "_")[[1]][2])))
    colnames(coords) <- c("lat", "long")
  }else{
    if(dim(coords)[1] != dim(omega)[1]){
      stop("coords provided, but the number of rows in coords data does not
           match the number of rows in omega matrix")
    }
  }
                                     
  pie_control_default <- list(edges = 200, clockwise = TRUE, 
                      init.angle = 90, density = NULL, 
                      angle = 45, border = NULL,
                      lty = NULL, label.dist = 1.1)
  
  pie_control <- modifyList(pie_control_default, pie_control)
  
  if(is.null(bgmap_path)){
    message("reading background map shapefile from inst/extdata/ne_110m_coastline 
            folder")
    GlobalCoast <- sf::st_read(system.file("extdata","ne_110m_coastline",
                     "ne_110m_coastline.shp",package = "ecostructure"), quiet=T)
  }else{
    GlobalCoast <- sf::st_read(bgmap_path, quiet = T)
  }
  
  glob <- c(xmin=long_lim[1], xmax=long_lim[2], ymin=lat_lim[1], ymax=lat_lim[2])
  glob <- sf::st_bbox(glob)
  glob <- structure(glob, crs = sf::st_crs(GlobalCoast))
  GlobalCoast <-suppressWarnings(suppressMessages(sf::st_intersection(GlobalCoast,
                                                        sf::st_as_sfc(glob))))
  
  if(adjust){
    idx <- which(omega[,1] > thresh)
    omega <- omega[-idx,]
    coords <- coords[-idx,]
    omega <- omega[,-1]
    omega <- t(apply(omega, 1, function(x) return(x/sum(x))))
  }
  
  output_type <- strsplit(path, "[.]")[[1]][2]
  
  if(output_type == "tiff"){
    tiff(path, width = image_width, height = image_height)
  }else if(output_type == "png"){
    png(path, width = image_width, height = image_height)
  }else if(output_type == "pdf"){
    pdf(path, width = image_width, height = image_height)
  }else{
    stop("the output image may either be of  tiff, png or pdf extension")
  }
  
  plot(sf::st_geometry(GlobalCoast), axes = T,main="",reset=F,xaxs = "i", 
       yaxs = "i", lwd=coastline_lwd)
  par(lwd =.01)
  invisible(lapply(1:dim(omega)[1], function(r)
    do.call(mapplots::add.pie, append(list(
                      z=as.integer(100*omega[r,]),
                      x=coords[r,1], 
                      y=coords[r,2], 
                      labels=c("","",""),
                      radius = radius,
                      col=sapply(color, scales::alpha, intensity))
                      , pie_control))))
  invisible(dev.off())
}