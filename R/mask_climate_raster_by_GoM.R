#' @title Mask a climate raster by the motif map raster
#'
#' @description Using the motif map raster obtained from applying Grade of
#' Membership model on the counts data, this function can mask climate raster -
#' including temperature, precipitation etc.
#'
#' @param motif_map_raster The motif map raster obtained from fitting GoM model 
#'                         on the species abundance counts data.
#' @param climate_raster the climate raster which is to be masked.
#' @param quant A scaling factor determining how fine the raster grid is
#' @param index The number of levels to be used
#'
#' @return Returns a masked climate raster.
#'
#' @importFrom raster quantile
#' @export


mask_climate_raster_by_geo_motif <- function(motif_map_raster, climate_raster, 
                                             quant = 0.10, index = 10){

  ras <- motif_map_raster
  ras[ras<(raster::quantile(ras,probs = seq(0, 1, quant), 
                            na.rm=T)[index])] <- NA
  cr <- crop(climate_raster, extent(ras))
  lr <- mask(x=cr, mask=ras)
  return(lr)
}

