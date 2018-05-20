#' @title Create global dispersion field matrix from local data frame
#' and shapefiles
#'
#' @description Reads a data frame or a matrix of local site-species abundance,
#' with rows being sites and columns being species, along with shapefiles for
#' each of the species (may be more too) ocurring in the data file and then
#' produces a dispersion field of the regional presence pattern for a user-
#' specified latitudinal and longitudinal range for each site, driven by 
#' the species present in that site.
#'
#' @param local_data     abundance counts site-species matrix/data frame.
#' @param shapefiles_dir The directory that contains shapefiles for all the
#'                       species in the data
#' @param raster_latlim The latitudinal range of the dispersion field
#' @param raster_longlim the longitudinal range of the dispersion field
#' @param raster_resolution The scale.resolution of the raster.
#'                          The higher the reason, the finer the raster is
#' @param base_local If the count of a species in a site is less than base_local,
#'                    it is assumed to be absent from the site. The default is 0.
#'                    Non zero values may result from obsrvational or data
#'                    recording errors
#' @param optional_species_names If the column name corresponding to some species
#'                      does not match with the shapefiles, the user might need
#'                      to replace it by another species name. In such a
#'                      scenario, instead of using column names, the user can
#'                      enter a new vector of names, which may be slight
#'                      modification on the species from the column of the data,
#'                      taking account of issues with non-matching species names
#'                      between shape files and the local data frame.
#' @param include_shp_attributes Shp files may include a number of attributes, e.g.
#'                      breeding and non breeding layers, and building 
#'                      dispersion fields
#'                      may only require a subset of these attributes. The user 
#'                      can provide a vector of attributes in order to build 
#'                      the dispersion field based only on those species 
#'                      attributes. The vector must be in a format
#'                      that is read in by rgdal and must include an "x" in 
#'                      place of where the species name will be inserted in the
#'                      function.  Default includes breeding and resident ranges
#'                      for global birds based on reading .shp files from 
#'                      BirdLife International and NatureServe.
#'
#' @return Returns a list with the names of the elements of the list given by 
#' the row names of the input data matrix (grids/sites). Each element of the 
#' list is the dispersion field matrix generated corresponding to the species 
#' present in the specific site to which it corresponds.
#'
#' @keywords counts data, local site-species data, global dispersion field
#'
#' @import sf
#' @importFrom fasterize fasterize
#' @importFrom raster raster setValues mask calc stack
#' @export

dispersion_fields_create = function(local_data,
                                    gis_data_type = "shp",
                                    shp_dir = NULL, 
                                    gdb_dir = NULL,
                                    gdb_object = NULL,
                                    gdb_layer = "All_Species",
                                    gdb_species_feature = "SCINAME",
                                    raster_latlim = c(5,50),
                                    raster_longlim = c(50,120),
                                    raster_resolution = 5,
                                    base_local = 0,
                                    quiet = T,
                                    optional_species_names = NULL,
                                    gdb_species_names = NULL, 
                                    include_shp_attributes = "x$SEASONAL==1 | x$SEASONAL==2"){

  if(is.null(optional_species_names)){
    species_names <- colnames(local_data)
  } else {
    species_names <- optional_species_names
  }
  
  if(is.null(gdb_species_names)){
    gdb_species_names <- colnames(local_data)
  } else {
    gdb_species_names <- gdb_species_names
  }
  
  if(length(species_names)!=length(colnames(local_data))){
    stop("The length of the optional species names vector must match 
         with number of columns in local data")
  }
  
  local.species<-list()
  
  for (i in 1:length(row.names(local_data))){
    local.list <- species_names[which(local_data[i,]>base_local)]
    local.species[[i]]<-local.list
  }
  
  dispersion.field.matrix <- list()
  dispersion.field.raster <- list()
  
  
  if(gis_data_type == "shp"){
    
    if(is.null(shp_dir)){
      shp_dir <- getwd();
    }
    
    shapefile_names<-list.files(shp_dir, pattern=".shp")
    
    shapefile_namenos<- unlist(sapply(species_names,
                                      function(x) grep(x, shapefile_names)))
    
    ranges<-list()
    
    for (j in 1:length(shapefile_namenos)){
      cat("Reading shapefile for species", j, "\n")
      suppressMessages(ranges[[j]]<-sf::st_read(dsn=path.expand(paste0(shp_dir,
                   "/", shapefile_names[shapefile_namenos[j]])), quiet=quiet))
    }
    
    for (z in 1:length(local.species)){
      species_temp <- local.species[[z]]
      matched_indices <- match(species_temp, species_names)
      
      local_ranges <- lapply(1:length(matched_indices),
                      function(l) return(ranges[[matched_indices[l]]]))
      
      if(is.null(include_shp_attributes)){
        breedingranges <- local_ranges
      }else{
        breedingranges<-lapply(local_ranges, 
        function(x) eval(parse(text=paste0("x[",include_shp_attributes,",]"))))
      }
      
      xmin<-raster_longlim[1]
      xmax<-raster_longlim[2]
      ymin<-raster_latlim[1]
      ymax<-raster_latlim[2]
      
      r <- raster::raster(nrows=(ymax-ymin)*raster_resolution,
                          ncols=(xmax-xmin)*raster_resolution,
                          xmn=xmin,
                          xmx=xmax,
                          ymn=ymin,
                          ymx=ymax)
      
      r_sf <-fasterize::fasterize(breedingranges[[1]], r, fun="sum")
      r_sf[r_sf<0] <- 1
      r_sf[is.na(r_sf==0)] <- 0
      for (i in 2:length(breedingranges)){
        r_sf_2 <-fasterize::fasterize(breedingranges[[i]], r, fun="sum")
        r_sf_2[r_sf_2<0] <- 1
        r_sf_2[is.na(r_sf_2)] <- 0
        r_sf <- raster::stack(r_sf,r_sf_2)}
      sum<-sum(r_sf)
      dispersion.field.matrix[[z]] <- raster::as.matrix(sum)
      dispersion.field.raster[[z]] <- sum
      cat("Dispersion Field Map Complete for site:", z, "\n")
    }
  }
  
  if(gis_data_type == "gdb"){
    
    if(is.null(gdb_object)){
      if(is.null(gdb_dir)){
        gdb_dir <- getwd()
      } else geo_db_ranges <- sf::st_read(dsn=path.expand(paste0(gdb_dir)), 
                                          layer = gdb_layer, quiet=quiet)
    } else geo_db_ranges <- gdb_object
    
    for (z in 1:length(local.species)){
      species_temp <- gdb_species_names[which(local_data[z,]>base_local)]
      
      local_geo_db_ranges <- eval(parse(text=paste0("geo_db_ranges[geo_db_ranges$",
                                  gdb_species_feature, "%in% species_temp,]")))
      
      if(is.null(include_shp_attributes)){
        breedingranges <- local_geo_db_ranges
      }else{
        x <- local_geo_db_ranges
        breedingranges <- eval(parse(text=paste0("x[",include_shp_attributes,",]")))
      }
      
      xmin<-raster_longlim[1]
      xmax<-raster_longlim[2]
      ymin<-raster_latlim[1]
      ymax<-raster_latlim[2]
      
      r <- raster::raster(nrows=(ymax-ymin)*raster_resolution,
                          ncols=(xmax-xmin)*raster_resolution,
                          xmn=xmin,
                          xmx=xmax,
                          ymn=ymin,
                          ymx=ymax)
      
      r_sf <-fasterize::fasterize(breedingranges, r, fun="sum")
      ####NEED TO MAKE SURE EACH SPECIES IS COUNTED AS 1 
      dispersion.field.matrix[[z]] <- raster::as.matrix(r_sf)
      dispersion.field.raster[[z]] <- r_sf
      cat("Dispersion Field Map Complete for site:", z, "\n")
    }
  }
  
  names(dispersion.field.matrix) <- rownames(local_data)
  names(dispersion.field.raster) <- rownames(local_data)
  
  dispersion.field <- list("matrix" = dispersion.field.matrix, 
                           "raster" = dispersion.field.raster)
  
  return(dispersion.field)
}
