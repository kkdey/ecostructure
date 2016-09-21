


################   Example  ######################

library(HimalayanBirdsAbundance)
data("HimalayanBirdsAbundance")
library(Biobase)
new_counts <- t(exprs(HimalayanBirdsAbundance));
metadata <- pData(HimalayanBirdsAbundance);

elevation_metadata=metadata$Elevation;
east_west_dir = metadata$WorE;

topic_clus <- maptpx::topics(new_counts, K=2, tol=0.000001)
omega <- topic_clus$omega
library(gridExtra)
library(grid)
BlockStructure(omega, blocker_metadata = east_west_dir,
               order_metadata = elevation_metadata,
               yaxis_label = "Elevation",
               levels_decreasing = FALSE)



########################################################################

  dispersion.field <- get(load("dispersion_field_mat.rda"))
  global_shapefile <- readShapeLines('ne_50m_admin_0_countries.shp', proj4string=proj)
  maps <- CreateMapsFromDispersionFields(dispersion.field, global_shapefile )

#######################################################################

  dispersion.field <- get(load("dispersion_field_mat.rda"))
  counts <- DispersionFieldToCounts(dispersion.field)

#####################################################################

   annotation = data.frame(x_names = c(paste0("A",1:5), paste0("B",1:5)),
   x = c(0.5,2.0, 3.2, 4.6, 6.3,  23.5, 26.4, 28.5, 29.6, 31.8),
   y1 = c(0.4, 0.3, 0.4, 0.35, 0.4, 0.8, 0.85, 0.9, 0.8, 0.75),
   y2d =c(5, 6.6, 4, 5.2, 20, 3.4, 5.6, 4.5, 8, 10))

   TopicMetaDiversity(annotation)

#####################################################################

    annotation = data.frame(x_names = c(paste0("A",1:5), paste0("B",1:5)),
    x = c(0.5,2.0, 3.2, 4.6, 6.3,  23.5, 26.4, 28.5, 29.6, 31.8),
    y1 = c(0.4, 0.3, 0.4, 0.35, 0.4, 0.8, 0.85, 0.9, 0.8, 0.75),
    y2 =c(5, 6.6, 4, 5.2, 20, 3.4, 5.6, 4.5, 8, 10))

    TopicMetaMeta(annotation)

#####################################################################

