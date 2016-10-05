library(raster)
library(rgdal)

#westext <- c(72.9,79.0,30.0,34.4)
#eastext <- c(88.0,93.1,24.0,30.0)

biovars.west <- c(5,7,12,15)
biovars18data <- list()
biovars.east <- c(5,7,12,15)
biovars28data <- list()

biovars29data <- list()


j <- 1
setwd("~/Desktop/Projects/Assemblage Motifs /bio_18")
for (i in biovars.west){
  biovars18data[[j]]<-raster(paste("bio",i,"_18.bil",sep=""))
  j <- j+1
}


setwd("~/Desktop/Projects/Assemblage Motifs /alt_18")
alt18<-raster("alt_18.bil")


k <- 1
setwd("~/Desktop/Projects/Assemblage Motifs /bio_28")
for (i in biovars.east){
  biovars28data[[k]]<-raster(paste("bio",i,"_28.bil",sep=""))
  k <- k+1
}

setwd("~/Desktop/Projects/Assemblage Motifs /alt_28")
alt28<-raster("alt_28.bil")

k <- 1
setwd("~/Desktop/Projects/Assemblage Motifs /bio_29")
for (i in biovars.east){
  biovars29data[[k]]<-raster(paste("bio",i,"_29.bil",sep=""))
  k <- k+1
}

setwd("~/Desktop/Projects/Assemblage Motifs /alt_29")
alt29<-raster("alt_29.bil")


setwd("~/Desktop/Projects/Assemblage Motifs ")

g18 <- stack(biovars18data)
g18 <- stack(g18,alt18)
p18 <- stack(biovars28data)
p18 <- stack(p18,alt28)
p29 <- stack(biovars29data)
p29 <- stack(p29,alt29)

large.ras <- merge(g18,p18)
large.ras <- merge(large.ras,p29)

bio.him <- data.frame(rasterToPoints(large.ras))

grids.clim <- list()
elev.grids.choice <- list()
x.grids.choice <- list()
for(i in 1:nrow(metadata)){
  x.grids.choice[[i]] <- bio.him[which(abs(bio.him[,"x"]-metadata[i,4])==min(abs(bio.him[,"x"]-metadata[i,4]))),]
  elev.grids.choice[[i]] <- x.grids.choice[[i]][which(abs(x.grids.choice[[i]][,"y"]-metadata[i,3])==min(abs(x.grids.choice[[i]][,"y"]-metadata[i,3]))),]
  grids.clim[[i]] <- elev.grids.choice[[i]][which(abs(elev.grids.choice[[i]][,"layer.5"]-metadata[i,2])==min(abs(elev.grids.choice[[i]][,"layer.5"]-metadata[i,2]))),]
}
  

him.clim.pc <- prcomp(bio.him[,3:6], center=T, scale=T)
him.clim.pc$rotation
summary(him.clim.pc)

grids.clim.scores <- matrix(data = NA, nrow=38, ncol=3)
for(i in 1:length(grids.clim)){
  grids.clim.scores[i,1:3] <- him.clim.pc$x[as.numeric(rownames(grids.clim[[i]])),1:3]
}

rownames(grids.clim.scores) <- rownames(metadata)
grids.clim.dist <- as.matrix(dist(grids.clim.scores))


