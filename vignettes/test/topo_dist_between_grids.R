######## Landscape barriers between Grid Sites ##########

library(gdistance) 
library(raster)
library(rgdal)

setwd("~/Desktop/Projects/Assemblage Motifs /alt_18")
alt18<-raster("alt_18.bil")
setwd("~/Desktop/Projects/Assemblage Motifs /alt_28")
alt28<-raster("alt_28.bil")
setwd("~/Desktop/Projects/Assemblage Motifs /alt_29")
alt29<-raster("alt_29.bil")

ext <- extent(74,94,26,35)
worldclim.alt.30asec <- merge(alt18,alt28)
worldclim.alt.30asec <- merge(worldclim.alt.30asec,alt29)
worldclim.alt.30asec<- crop(worldclim.alt.30asec, ext)
setwd("~/Desktop/Projects/Assemblage Motifs ")


metadata <- read.csv("gridsmeta.csv",header=T)
metadata <-metadata[c(9,8,1,10,6,2,3,36,7,11,5,33,4,30,14,32,31,12,13,37,15,34,19,38,18,35,17,23,16,24,28,21,26,20,27,25,22,29),]
row.names(metadata) <- metadata$GRID
east_west_dir = metadata$WorE
forest_patch_label <- row.names(metadata)
metadata.no.MA <- metadata[-8,]

ext <- extent(74,94,26,35)
him.ext <- c(74.36682,93.51216,26.12938,34.64562)
#westext <- c(72.9,79.0,30.0,34.4)
#eastext <- c(88.0,93.1,24.0,30.0)

SRTM.long.cells <- c(74:93)
SRTM.lat.cells <- c(26:34)


deg.ras <- list()
k <- 1
setwd("~/Desktop/Projects/Him_SRTM_V3_1_arc_sec/")
for (j in SRTM.lat.cells){
  for (i in SRTM.long.cells){
    deg.ras[[k]] <- raster(paste("N",j,"E0",i,".hgt",sep=""))
    k <- k+1
  }
}

deg.ras$filename <- 'test.tif'
deg.ras$overwrite <- TRUE
m <- do.call(merge, deg.ras)
Arc_sec_V3_SRTM <- m

setwd("~/Desktop/Projects/Assemblage Motifs ")

m <- worldclim.alt.30asec

#### build transition matix for each grid ####
for( r in 1:nrow(metadata.no.MA) ) {
  elev_orig      <- metadata.no.MA[r, 'Elevation']
  elev_diff      <- abs(big.ras - elev_orig)
  ## 12 minutes
  trans_diff     <- transition(elev_diff,function(x) 1/mean(x),8)
  ## 2.5 minutes
  trans_final    <- geoCorrection(trans_diff)
}

topo.set<- combn(row.names(metadata.no.MA), 2, FUN = NULL, simplify = TRUE)
topo.dist <- matrix(nrow=2, ncol=666)
colnames(topo.dist)<-topo.set[1,]
topo.dist[1,]<-topo.set[2,]
topo.dist[2,]<-NA

topo.latlong.set <- list()
for( r in 1:ncol(topo.set)){
  latlong.set <- vector()
  latlong.set[1] <- metadata.no.MA[topo.set[1,r],"East"]
  latlong.set[2] <- metadata.no.MA[topo.set[1,r],"North"]
  latlong.set[3] <- metadata.no.MA[topo.set[2,r],"East"]
  latlong.set[4] <- metadata.no.MA[topo.set[2,r],"North"]
  topo.latlong.set[[r]] <- latlong.set
}

for( r in 49:ncol(topo.dist) ) {
  print(r)
  start.time <- Sys.time()
  if (r == 1){
    elev_orig      <- metadata.no.MA[colnames(topo.dist)[r], 'Elevation']
    elev_diff      <- abs(m - elev_orig)
    ## 12 minutes
    trans_diff     <- transition(elev_diff,function(x) 1/mean(x),8)
    ## 2.5 minutes
    trans_final    <- geoCorrection(trans_diff, type ="c")
    save(trans_final, file=paste('/Users/alexanderwhite/Desktop/Projects/Assemblage Motifs /LCD_GRIDS/transMat_', colnames(topo.dist)[r],  '.Rdata', sep=''))
    #rm(trans_final)
    #trans_final <- load(paste('/Users/alexanderwhite/Desktop/Projects/Assemblage Motifs /LCD_GRIDS/transMat_', colnames(topo.dist)[r],  '.Rdata', sep=''))
    #load(paste('~/Desktop/Research Projects/beast/data/transMat_', col.names(topo.dist)[r],  '.Rdata', sep=''))
  } else if (colnames(topo.dist)[r] == colnames(topo.dist)[r-1]){
    trans_final <- trans_final 
  } else {
    elev_orig      <- metadata.no.MA[colnames(topo.dist)[r], 'Elevation']
    elev_diff      <- abs(m - elev_orig)
    ## 12 minutes
    trans_diff     <- transition(elev_diff,function(x) 1/mean(x),8)
    ## 2.5 minutes
    trans_final    <- geoCorrection(trans_diff, type ="c")
    save(trans_final, file=paste('/Users/alexanderwhite/Desktop/Projects/Assemblage Motifs /LCD_GRIDS/transMat_', colnames(topo.dist)[r],  '.Rdata', sep=''))
    #trans_final <- load(paste('/Users/alexanderwhite/Desktop/Projects/Assemblage Motifs /LCD_GRIDS/transMat_', colnames(topo.dist)[r],  '.Rdata', sep=''))
    #load(paste('~/Desktop/Research Projects/beast/data/transMat_', rownames(actv)[r],  '.Rdata', sep=''))
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  #####################################
  ## Least Cost Distance
  cost <- costDistance(trans_final, topo.latlong.set[[r]][1:2],topo.latlong.set[[r]][3:4])
  topo.dist[2,r] <- cost[1,1]
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  #save( cost, file=paste('/home/aewhite/lcd/cost/cost_', rownames(metadata.no.MA)[r], '.Rdata', sep=''))
  
  ## Commute Distance
  #coords <- rbind(cbind(longLat2[r, 'east.long'], longLat2[r, 'east.lat']), cbind(longLat2[r, 'west.long'], longLat2[r, 'west.lat']))
  #commute <- commuteDistance( trans_final, coords)
  #save(commute, file=paste('/home/aewhite/lcd/commute/commute_', rownames(actv)[r], '.Rdata', sep=''))
  
}

topo.dist.mat <- matrix(nrow=37,ncol=37)
topo.dist.mat[lower.tri(topo.dist.mat)] <- topo.dist[2,]
rownames(topo.dist.mat) <- rownames(metadata.no.MA)
colnames(topo.dist.mat) <- rownames(metadata.no.MA)
diag(topo.dist.mat) <- 0




