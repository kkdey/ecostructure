
library(ecostructure)


data("himalayan_birds")
species_abundance_counts <- t(exprs(himalayan_birds));
site_metadata <- pData(himalayan_birds);
elevation_metadata=site_metadata$Elevation;
east_west_dir = site_metadata$WorE;
topic_clus <- CountClust::FitGoM(species_abundance_counts, K=2, tol=0.1)
omega <- topic_clus$fit$omega

blocker <- as.character(east_west_dir)
blocker[sample(1:30, 7, replace = FALSE)] = "N"
blocker <- as.factor(blocker)


block_structure(omega,
               blocker_metadata = blocker,
               order_metadata = elevation_metadata)

out <- ecostructure_nullmodel(species_abundance_counts, K=2, 
                              iter_randomized=5, option = "BF")
out2 <- ecostructure_nullmodel(species_abundance_counts, K=2, 
                               iter_randomized=5, ption = "BIC")

species_pa_counts <- species_abundance_counts
species_pa_counts[species_pa_counts >=1] = 1

fit <- ecostructure_fit(species_pa_counts, K = 2, tol = 0.1)
block_structure(fit$omega,
                blocker_metadata = blocker,
                order_metadata = elevation_metadata)

fit2 <- ecostructure_fit(species_abundance_counts, K = 2, tol = 0.1)
block_structure(fit2$omega,
                blocker_metadata = blocker,
                order_metadata = elevation_metadata)



