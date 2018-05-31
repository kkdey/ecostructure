#' @title Create dispersion fields from local species abundance data 
#' and associated GIS data source (shp files or geodatabase file).
#'
#' @description Reads a data frame or a matrix of local site-species abundances,
#' with rows being sites and columns being species, along with GIS data files for
#' each of the species (may be more too) ocurring in the local data matrix. The function
#' produces a dispersion field of the regional presence pattern for each site in 
#' a user-specified latitudinal and longitudinal range.
#'
#' @param local_data abundance site (row) x species (column) matrix/data frame.
#' @param gis_data_type One of c("shp","gdb"). Are the GIS data in seperate shapefiles 
#'                      or in a geodatabase file?
#' @param shp_dir The directory that contains shapefiles for all the
#'                       species in the local data matrix.
#' @param gdb_dir The directory that contains the geodatabase. 
#' @param gdb_object The geodatabase can be read in ahead of time and provided as an sf object. 
#'                    This will dramatically reduce run time, especially useful if the function 
#'                    is to be run multiple times.
#' @param gdb_layer If layer name of the gdb file that contains the species distributions. 
#' @param gdb_species_feature The name of the feature that contains the species names 
#'                             in the geodatabse.
#' @param raster_latlim The latitudinal range of the desired dispersion field.
#' @param raster_longlim The longitudinal range of the desired dispersion field.
#' @param raster_resolution The scale resolution of the raster in degrees.
#' @param raster_resolution The precision of the initial estimation provided to fasterize. 
#'                          fasterize tends to undercount small species ranges, so is required 
#'                          to keep at a low value in order detect species presence in a square degree cell.
#'                          Most users will not need to change, but increaseing this value will 
#'                          speed up the function but with potential for error. 
#' @param base_local If the count of a species in a local site is less than base_local,
#'                    it is assumed to be absent from the site. The default is 0.
#'                    Non zero values may result from obsrvational or data
#'                    recording errors
#' @param optional_species_names If the column name corresponding to some species
#'                      does not match with the GIS data, the user might need
#'                      to replace it by another species name. In such a
#'                      scenario, instead of using column names, the user can
#'                      enter a new vector of names, which may be slight
#'                      modification on the species from the column of the data,
#'                      taking account of issues with non-matching species names
#'                      between shape files and the local data frame.
#' @param include_shp_attributes Each species range may be divided according to a number of 
#'                      attributes, e.g. breeding and non breeding layers, and building dispersion fields
#'                      may require only a subset of these attributes. The user can provide a 
#'                      string of attributes in order to build the dispersion field
#'                      based only on those species attributes. The vector must be in a format
#'                      that is read in by sf and must include an "x" in place of where the 
#'                      sf object will be manipulated in the function.  Default includes 
#'                      breeding and resident ranges for global birds based on reading 
#'                      GIS data from BirdLife International.  
#'
#' @return Returns 3 lists. In each list, each element is a
#' dispersion field corresponding to the species present in the 
#' specific site to which it corresponds. The first list "matrix" 
#' contains each dispersion field as a matrix. The second list "raster" returns 
#' each dispersion field as a raster.  The third list "precise" returns 
#' each dispersion field at the high resolution of precision.  The names of the elements 
#' of each list are assigned as the row names of the input data matrix (sites). 
#'
#' @keywords counts data, local site-species data, dispersion field
#'
#' @import sf
#' @importFrom fasterize fasterize
#' @importFrom raster raster setValues stack aggregate as.matrix
#' @export

dsp_create_from_survey = function(local_data,
                                  gis_data_type = "shp",
                                  shp_dir = NULL, 
                                  gdb_dir = NULL,
                                  gdb_object = NULL,
                                  gdb_layer = "All_Species",
                                  gdb_species_feature = "SCINAME",
                                  raster_latlim = c(5,50),
                                  raster_longlim = c(50,120),
                                  raster_resolution = 1,
                                  precision = 0.05,
                                  base_local = 0,
                                  quiet = T,
                                  optional_species_names = NULL,
                                  include_shp_attributes = "x$SEASONAL==1 | x$SEASONAL==2") {
  
  
  if(is.null(optional_species_names)){
    species_names <- colnames(local_data)
  } else {
    species_names <- optional_species_names
  }
  
  if(length(species_names)!=length(colnames(local_data))){
    stop("The length of the optional species names vector must match with number of columns in local data")
  }
  
  local.species<-list()
  for (i in 1:length(row.names(local_data))){
    local.species[[i]] <- species_names[which(local_data[i,]>base_local)]
  }
  
  dispersion.field.matrix <- list()
  dispersion.field.raster <- list()
  dispersion.field.precise <- list()
  
  if(gis_data_type == "shp"){
    
    if(is.null(shp_dir)){
      shp_dir <- getwd();
    }
    
    shapefile_names <- list.files(shp_dir, pattern=".shp")
    shapefile_namenos <- unlist(sapply(species_names,function(x) grep(x, shapefile_names)))
    
    ranges<-list()
    for (j in 1:length(shapefile_namenos)){
      
      cat("Reading shapefile for species", j, "\n")
      range <- sf::st_read(dsn=path.expand(paste0(shp_dir,"/",shapefile_names[shapefile_namenos[j]])),quiet=quiet)
      
      if(!sf::st_crs(range) == sf::st_crs(4326)){
        if(is.na(sf::st_crs(range))){
          cat("No crs for shapefile of species", j, ". Assuming st_crs(4326)\n")
          sf::st_crs(range) <- sf::st_crs(4326)
        } else {
          range <- sf::st_transform(range, sf::st_crs(4326), check=T)
        }
      }
      
      ranges[[j]] <- range
    }
    
    if(!is.null(include_shp_attributes)){
      ranges <- lapply(ranges, function(x) eval(parse(text=paste0("x[",include_shp_attributes,",]"))))
    }
    
    r_prim <- raster::raster(resolution = precision,
                             xmn=raster_longlim[1],
                             xmx=raster_longlim[2],
                             ymn=raster_latlim[1],
                             ymx=raster_latlim[2])
    
    for (z in 1:length(local.species)){
      
      matched_indices <- match(local.species[[z]], species_names)
      
      breedingranges <- lapply(1:length(matched_indices),
                               function(l) return(ranges[[matched_indices[l]]]))
      
      r_sf <- fasterize::fasterize(breedingranges[[1]],
                                   r_prim, fun="any", background=0)
      
      for (i in 2:length(breedingranges)) {
        r_sf_2 <- fasterize::fasterize(breedingranges[[i]], 
                                       r_prim, fun="any", background=0)
        r_sf <- raster::stack(r_sf, r_sf_2)
      }
      
      r_sf <- sum(r_sf)
      
      r_sf_ag <- raster::aggregate(r_sf, fact=(1/precision)*raster_resolution, fun=max, expand=F) 
      
      r_sf[r_sf == 0] <- NA
      r_sf_ag[r_sf_ag == 0] <- NA
      
      dispersion.field.matrix[[z]] <- raster::as.matrix(r_sf_ag)
      dispersion.field.raster[[z]] <- r_sf_ag
      dispersion.field.precise[[z]] <- r_sf
      
      cat("Dispersion Field Map Complete for site:", z, "\n")
    }
  }
  
  if(gis_data_type == "gdb"){
    
    if(is.null(gdb_object)){
      if(is.null(gdb_dir)){
        gdb_dir <- getwd()
      } else gdb_object <- sf::st_read(dsn=path.expand(paste0(gdb_dir)), layer = gdb_layer, quiet=quiet)
    }
    
    if(!sf::st_crs(gdb_object) == sf::st_crs(4326)){
      if(is.na(sf::st_crs(gdb_object))){
        cat("No crs for gdb. Assuming st_crs(4326)\n")
        sf::st_crs(gdb_object) <- sf::st_crs(4326)
      } else {
        gdb_object <- st_transform(gdb_object, sf::st_crs(4326), check=T)
      }
    }
    
    colnames(gdb_object)[colnames(gdb_object) == gdb_species_feature] <- "ecos_temp_name"
    
    gdb_object <- gdb_object[gdb_object$ecos_temp_name %in% species_names,]
    
    if(!is.null(include_shp_attributes)){
      include_shp_attributes <- gsub("x","gdb_object",include_shp_attributes)
      gdb_object <- eval(parse(text=paste0("gdb_object[",include_shp_attributes,",]")))
    }
    
    
    
    r_prim <- raster::raster(resolution = precision,
                             xmn=raster_longlim[1],
                             xmx=raster_longlim[2],
                             ymn=raster_latlim[1],
                             ymx=raster_latlim[2])
    
    for (z in 1:length(local.species)){
      
      breedingranges <- gdb_object[gdb_object$ecos_temp_name %in% local.species[[z]],]
      
      ### fasterize has a problem with small ranges and will under count with low res - but is fast
      r_sf <-fasterize::fasterize(breedingranges, r_prim, fun="any", by="ecos_temp_name",background=0)
      
      r_sf <- sum(r_sf)
      
      r_sf_ag <- raster::aggregate(r_sf, fact=(1/precision)*raster_resolution, fun=max, expand=F) 
      
      r_sf[r_sf == 0] <- NA
      r_sf_ag[r_sf_ag == 0] <- NA
      
      dispersion.field.matrix[[z]] <- raster::as.matrix(r_sf_ag)
      dispersion.field.raster[[z]] <- r_sf_ag
      dispersion.field.precise[[z]] <- r_sf
      
      cat("Dispersion Field Map Complete for site:", z, "\n")
    }
  }
  
  names(dispersion.field.matrix) <- rownames(local_data)
  names(dispersion.field.raster) <- rownames(local_data)
  names(dispersion.field.precise) <- rownames(local_data)
  
  dispersion.field <- list("matrix" = dispersion.field.matrix, 
                           "raster" = dispersion.field.raster,
                           "precise" = dispersion.field.precise)
  
  return(dispersion.field)
}
