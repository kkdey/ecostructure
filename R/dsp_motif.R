#' @title Creates heatmap for a geographic motif 
#'
#' @description Takes the motif fitted from \code{ecos_fit} and plots the theta
#' component as a heatmap across space using suitable raster information. 
#'
#' @param theta_mat The theta output from the model fit from 
#'                  \code{ecos_fit}
#' @param dsp_fld_res The resolution of the dispersion field.
#' @param color_ramp Color palette to be used for the heat map. 
#' @param raster_latlim If disp_fld_list is a list of matrices, the latitudinal range 
#'                      of the dispersion fields must be provided as a vector. 
#' @param raster_longlim If disp_fld_list is a list of matrices, the longitudinal range 
#'                      of the dispersion fields must be provided as a vector. 
#' @param outline_col The outline color of the projected map.  
#' @param proj A character string of projection arguments; the arguments must be entered 
#'            exactly as in the PROJ.4 documentation. To be evaluated by sp::CRS(). 
#'            default '+proj=longlat +ellps=WGS84'. 
#'            
#' @return Returns a heatmap plot, similar to the \code{dsp_plot_map} function
#'         for the theta matrix using the specified raster,
#'         
#' @keywords dispersion field, motif, maps
#'
#' @examples
#' 
#' data("example_theta")
#' himalayan_geo_motifs <- dsp_motif(example_theta, 
#'                                  color_ramp = c("black", "darkseagreen3",
#'                                                 "orange","red"),
#'                                 dsp_fld_res = 8)
#' himalayan_geo_motifs$motif_maps[[1]]
#' 
#' @export


dsp_motif = function(theta_mat,
                     dsp_fld_res = 8, 
                     color_ramp = c("black", "darkseagreen3",
                                    "orange","red"),
                     raster_latlim = c(5,50),
                     raster_longlim = c(50,120),
                     outline_col = "white",
                     proj = '+proj=longlat +ellps=WGS84'){
  
  global_boundaries <- sf::st_read(system.file("extdata","ne_110m_admin_0_boundary_lines_land",
                                               "ne_110m_admin_0_boundary_lines_land.shp",
                                               package = "ecostructure"),quiet=T)
  
  ext_gbl_bnd <-raster::extent(global_boundaries)
  
  global_coast <- sf::st_read(system.file("extdata","ne_110m_coastline",
                                          "ne_110m_coastline.shp",
                                          package = "ecostructure"),quiet=T)
  
  ext_gbl_cst <-raster::extent(global_coast)
  
  if(!proj=='+proj=longlat +ellps=WGS84'){
    global_boundaries <- st_transform(global_boundaries, proj)
    global_coast <- st_transform(global_coast, proj)
    
    global_boundaries <- as(global_boundaries,'Spatial')
    global_boundaries@bbox <- as.matrix(ext_gbl_bnd)
    global_coast <- as(global_coast,'Spatial')
    global_coast@bbox <- as.matrix(ext_gbl_cst)
  } else {
    global_boundaries <- as(global_boundaries,'Spatial')
    global_coast <- as(global_coast,'Spatial')
  }
  
  newtheme <- rasterVis::rasterTheme(region = grDevices::colorRampPalette(color_ramp)( 100 ))
  
  motif_thetas <- list()
  for (i in 1:dim(theta_mat)[2]){
    motif_thetas[[i]] <- matrix(theta_mat[,i],
                                nrow = (raster_latlim[2]-raster_latlim[1])/(1/dsp_fld_res),
                                ncol = (raster_longlim[2]-raster_longlim[1])/(1/dsp_fld_res))
  }
  
  theta_map <- list()
  motif_rasters <- list()
  
  for (i in 1:length(motif_thetas)){
    
    r <- raster::raster(motif_thetas[[i]],
                        xmn=raster_longlim[1],
                        xmx=raster_longlim[2],
                        ymn=raster_latlim[1],
                        ymx=raster_latlim[2],
                        crs=sp::CRS(proj))
    
    theta_map[[i]] <- rasterVis::levelplot(r, par.settings = newtheme,
                                           contour=F, margin=FALSE,
                                           at=seq(0, raster::maxValue(r), length.out=100)) +
      latticeExtra::layer(sp::sp.lines(global_boundaries,col= outline_col, lwd=0.5),
                          data=list(outline_col = outline_col,
                                    global_boundaries = global_boundaries)) +
      latticeExtra::layer(sp::sp.lines(global_coast,col = outline_col, lwd=0.5),
                          data=list(outline_col = outline_col,
                                    global_coast = global_coast))
    
    motif_rasters[[i]] <- r
    
  }
  names(theta_map) <- colnames(theta_mat)
  names(motif_rasters) <- colnames(theta_mat)
  ll <- list("motif_maps" = theta_map, "motif_rasters" = motif_rasters)
  return(ll)
}

#test2 <- dispersion_fields_motif(test$topic_fit$theta)