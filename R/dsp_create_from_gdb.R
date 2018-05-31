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


dsp_create_from_gdb = function(gdb_object,
                               raster_resolution = 5,
                               thresh = 4,
                               raster_latlim = c(-90,90),
                               raster_longlim = c(-180,180),
                               species_feature = "SCINAME",
                               species_names = NULL){
  
  colnames(gdb_object)[colnames(gdb_object) == species_feature] <- "ecos_temp_name"
  
  r <- raster::raster(resolution=raster_resolution,
                      xmn=raster_longlim[1],
                      xmx=raster_longlim[2],
                      ymn=raster_latlim[1],
                      ymx=raster_latlim[2])
  
  idx <- list()
  for (i in 1:ncell(r)){
    idx[[i]] <- as.vector(xyFromCell(r, i, spatial=FALSE))
  }
  
  idx_coords <- c()
  for (i in 1:length(idx)){
    idx_coords[i] <- cellFromXY(r, idx[[i]])
  }
  
  r_sf <- st_as_sf(as(r,'SpatialPolygons'))
  
  if(is.na(st_crs(gdb_object))){
    stop("gdb_object should have crs set")
  }
  
  if(!is.null(species_names)){
    gdb_object <- gdb_object %>% dplyr::filter(ecos_temp_name %in% species_names)
  } 
  
  if(!st_crs(gdb_object)==st_crs(4326)){
    gdb_object <- st_transform(gdb_object, st_crs(4326), check=T) # get gdb_object in same proj as r_sf
  } 
  
  cat("Subsetting gdb_object to raster boundary - may take a while if gdb_object is large \n")
  gdb_cropped <- st_relate(gdb_object, st_as_sfc(st_bbox(r_sf)), pattern = "T********", sparse=F)
  gdb_object <- gdb_object %>% dplyr::filter(gdb_cropped)
  
  cat("Overlapping ranges to raster - may take a while if gdb_object is large \n")
  overlps <- st_relate(gdb_object, r_sf, pattern = "T********", sparse=F)
  
  overlps_sp <- apply(overlps,1,which) # which cells(rows) from r_sf does each species(row) in gdb_obj_tnfd overlap
  overlps_cell <- apply(overlps,2,which) # which species(rows) from gdb_obj_tnfd overlap each cell(row) in r_sf 
  
  cells_with_overlap <- which(unlist(lapply(overlps_cell,any))) # which r_sf rows have > 0 species
  focal_LL <- idx[cells_with_overlap] # lat_long of cells with > 0 species
  
  cells_with_asblg <- which(unlist(lapply(overlps_cell,function(x) length(x) > thresh))) # which r_sf rows have > assemblage threshold
  focal_asblg_LL <- idx[cells_with_asblg] # lat_long of cells with > threshold
  
  ### want to create dispersion fields 
  ### want to name each dispersion field with lat_long
  ### want to return a pres_ab matrix (as.raster then to matrix) of r_sf cells with > 0
  ### want to return a pres_ab matrix (as.raster then to matrix) of r_sf cells with > thresh
  
  dispersion.field.raster <- list() 
  dispersion.field.matrix <- list()
  
  lat_long_names <- c()
  
  pres_ab_mat <- matrix(as.numeric(overlps), nrow=dim(overlps)[1])
  
  r2 <- setValues(r,0)
  build_stack = function(index_to_change){
    r_df <- r2
    values(r_df)[index_to_change] <- 1
    return(r_df)
  }
  
  for (i in 1:length(cells_with_asblg)){
    asblg <- overlps_cell[cells_with_asblg[i]][[1]]
    
    sp_ranges <- lapply(c(overlps_sp[asblg]),build_stack)
    
    dsp_fld <- raster::stack(sp_ranges)
    dsp_fld <- sum(dsp_fld)
    
    dsp_fld[dsp_fld == 0] <- NA
    
    dispersion.field.raster[[i]] <- dsp_fld
    dispersion.field.matrix[[i]] <- raster::as.matrix(dsp_fld)
    lat_long_names[i] <- paste(focal_asblg_LL[[i]][2],
                               focal_asblg_LL[[i]][1],sep="_")
  }
  
  
  row.names(pres_ab_mat) <- gdb_object$ecos_temp_name
  colnames(pres_ab_mat) <- sapply(idx, function(x) paste(x[2],x[1],sep="_"))
  
  
  names(dispersion.field.matrix) <- lat_long_names
  names(dispersion.field.raster) <- lat_long_names
  
  dispersion.field <- list("matrix" = dispersion.field.matrix, 
                           "raster" = dispersion.field.raster,
                           "pres_ab" = pres_ab_mat)
  
  return(dispersion.field)
}


 