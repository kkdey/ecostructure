

###################

library(ecostructure)
library(Biobase)

data <- get(load(system.file("extdata", "HimalayanBirdsData.rda",
                             package = "ecostructure")))
taxonomic_counts <- t(exprs(data))

grid_metadata <- pData(phenoData(data))
head(grid_metadata)
elevation_metadata=grid_metadata$Elevation;
east_west_dir = grid_metadata$WorE;

voom_data <- t(limma::voom(t(taxonomic_counts))$E)
pca_data <- prcomp(voom_data)$x
plot(pca_data[,1], pca_data[,2], col = east_west_dir, xlab = "PC1", ylab = "PC2", main = "PCA",
     pch = 20, lwd = 1.5)


elevation <- elevation_metadata
df <- data.frame("PC1" = pca_data[,1],
                 "PC2" = pca_data[,2])

ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(colour = elevation)) +
  scale_fill_gradientn(colours = terrain.colors(10))

library(MASS)
d <- dist(taxonomic_counts) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

df <- data.frame("MDS1" = fit$points[,1],
                 "MDS2" = fit$points[,2])

ggplot(df, aes(MDS1, MDS2)) +
  geom_point(aes(colour = elevation)) +
  scale_fill_gradientn(colours = terrain.colors(10))

out <- tsne::tsne(voom_data, k = 2)

df <- data.frame("tSNE1" = out[,1],
                 "tSNE2" = out[,2])

ggplot(df, aes(tSNE1, tSNE2)) +
  geom_point(aes(colour = elevation)) +
  scale_fill_gradientn(colours = terrain.colors(10))


