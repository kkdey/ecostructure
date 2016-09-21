#' @title Create global abundance dispersion field maps from local data frame
#' and shapefiles
#'
#' @description Reads a data frame or a matrix of local site-species abundance,
#' with rows being sites and columns being species, along with shapefiles for
#' each of the bird species (may be more too) ocurring in the data file and then
#' produces global dispersion field of abundance pattern for a user-specified
#' latitudinal and longitudinal range for each site, driven by the species
#' present in that site.
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
#'                      to replace it by another species close to it. In such a
#'                      scenario, instead of using column names, the user can
#'                      enter a new vector of names, which may be slight
#'                      modification on the species from the column of the data,
#'                      taking account of issues with non-matching bird species
#'                      between shape files and the local data frame.
#' @return Returns a list with the names of the elements of the list given by the
#' row names of the input data matrix (sites). Each element of the list is the
#' dispersion field matrix generated corresponding to the bird species present
#' in the specific site to which it corresponds.
#'
#' @keywords counts data, local site-species data, global dispersion field
#'
#' @importFrom rgdal readOGR
#' @import raster
#' @importFrom utils stack
#'
#' @export

CreateGlobalAbundanceFields = function(local_data,
                                     shapefiles_dir=NULL,
                                     raster_latlim = c(5,50),
                                     raster_longlim = c(50,120),
                                     raster_resolution = 5,
                                     base_local = 0,
                                     optional_species_names = NULL){


  if(is.null(optional_species_names)){
    species_names <- colnames(local_data)
  }else{
    species_names <- optional_species_names;
  }

  if(length(species_names)!=length(colnames(local_data))){
    stop("The length of the optional species names vector must match with number of columns in local data")
  }


  dispersion.field <- list()
  local.species<-list()
  for (i in 1:length(row.names(local_data))){
    local.list <- species_names[which(local_data[i,]>base_local)]
    local.species[[i]]<-local.list
  }

  if(is.null(shapefiles_dir)){
    shapefiles_dir <- getwd();
  }
  shapefile_names<-list.files(shapefiles_dir, pattern=".shp")


  shapefile_namenos<- unlist(sapply(species_names,
                                           function(x) grep(x, shapefile_names)))
  ranges<-list()

  for (j in 1:length(shapefile_namenos)){
    cat("Reading shapefile for species", j, "\n")
    suppressMessages(ranges[[j]]<-rgdal::readOGR(dsn=paste0(shapefiles_dir,"/",
                                    shapefile_names[shapefile_namenos[j]]),
                layer=sub(".shp",  "", shapefile_names[shapefile_namenos[j]])))}




  for (z in 1:length(local.species)){
    species_temp <- local.species[[z]]
    matched_indices <- match(species_temp, species_names);

    local_ranges <- lapply(1:length(matched_indices),
                           function(l) return(ranges[[matched_indices[l]]]))

     breedingranges<-local_ranges[sapply(local_ranges, function(x) sum(x@data$SEASONAL==1 | x@data$SEASONAL==2)>0)] #to remove shapefiles with no breeding ranges in them, where species are entered twice
     breedingranges<-sapply(breedingranges, function(x) x[x@data$SEASONAL==1 | x@data$SEASONAL==2,])  ##extracts breeding ranges

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

    r<- raster::setValues(r, 1)

    a <- raster::mask(r, breedingranges[[1]], updatevalue=0)
    for (i in 2:length(breedingranges)){
      a<- utils::stack(a,mask(r, breedingranges[[i]], updatevalue=0))}
      sum<-sum(a)

      fun <- function(x) { x[x==0] <- NA; return(x) }
      dispersion.field[[z]] <- as.matrix(raster::calc(sum, fun))
      cat("Dispersion Field Map Complete for site:", z, "\n")
  }

  names(dispersion.field) <- rownames(local_data)
  return(dispersion.field)
}

