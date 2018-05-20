library(raster)
library(fasterize)
library(sf)
library(rasterVis)


X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 10/10), ncol = 2)
poly_1 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 15/10), ncol = 2)
poly_2 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 6/10), ncol = 2)
poly_3 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 20/10), ncol = 2)
poly_4 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 13/10), ncol = 2)
poly_5 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 10/10), ncol = 2)
poly_6 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 15/10), ncol = 2)
poly_7 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 6/10), ncol = 2)
poly_8 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 20/10), ncol = 2)
poly_9 <- list(X[c(chull(X),chull(X)[1]),])
X <- matrix(stats::rnorm(2000, mean = 60,
                         sd = 13/10), ncol = 2)
poly_10 <- list(X[c(chull(X),chull(X)[1]),])

gdb_toy <- sf::st_sf(PRESENCE = rep(1,10),
                     SEASONAL = c(1,1,1,1,2,2,2,1,2,1),
                     SCINAME = c("A","B","C","D","E","F","G","H","I","J"),
                     geometry = st_sfc(lapply(list(poly_1,poly_2,poly_3,
                                                   poly_4,poly_5,poly_6,
                                                   poly_7,poly_8,poly_9,
                                                   poly_10), st_polygon)))

site_matrix <- matrix(nrow=3,ncol=10)
site_matrix[1,] <- c(0,0,0,5,6,3,1,0,0,2)
site_matrix[2,] <- c(0,2,3,1,0,0,0,1,2,0)
site_matrix[3,] <- c(1,2,0,0,6,10,2,0,1,2)
rownames(site_matrix) <- c("Site_1","Site_2","Site_3")
colnames(site_matrix) <- c("A","B","C","D","E","F","G","H","I","J")

st_write(gdb_toy[gdb_toy$SCINAME=="A",], "A.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="B",], "B.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="C",], "C.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="D",], "D.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="E",], "E.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="F",], "F.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="G",], "G.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="H",], "H.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="I",], "I.shp")
st_write(gdb_toy[gdb_toy$SCINAME=="J",], "J.shp")

#write.csv(site_matrix,file="toy_site_matrix.csv")
#site_matrix <- read.csv(file="toy_site_matrix.csv",header=T, row.names=1)

test_1 <- dispersion_fields_create(local_data = site_matrix,
                                   gis_data_type = "shp",
                                   shp_dir = getwd(),
                                   raster_latlim = c(50,70),
                                   raster_longlim = c(50,70),
                                   raster_resolution = 5,
                                   base_local = 0,
                                   quiet = T,
                                   include_shp_attributes = "x$SEASONAL==1 | x$SEASONAL==2")

file.remove(list.files(getwd(), pattern = ".shp"))
file.remove(list.files(getwd(), pattern = ".shx"))
file.remove(list.files(getwd(), pattern = ".dbf"))
