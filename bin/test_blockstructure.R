
library(ecostructure)


data("himalayan_birds")
species_abundance_counts <- t(exprs(himalayan_birds));
site_metadata <- pData(himalayan_birds);
elevation_metadata=site_metadata$Elevation;
east_west_dir = site_metadata$WorE;
topic_clus <- ecostructure_fit(species_abundance_counts, K = 2, tol = 0.1)
omega <- topic_clus$omega


data("himalayan_birds")
species_abundance_counts <- t(exprs(himalayan_birds));
site_metadata <- pData(himalayan_birds);
elevation_metadata=site_metadata$Elevation;
east_west_dir = site_metadata$WorE;
data("himalayan_model")
ecostructure_blocks(himalayan_model$omega,
               blocker_metadata = east_west_dir,
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



data(himalayan_birds)
species_metadata <- pData(featureData(himalayan_birds))
taxonomic_counts <- t(exprs(himalayan_birds))
bill_traits <- as.matrix(dist(scale(species_metadata[,c(1:3)])))
counts_bill_traits <- ecostructure_prepare_by_trait(counts = taxonomic_counts, 
                                                    traits = bill_traits, 
                                                    prop_div=0.3)


data("australia_birds")
counts <- australia_birds$pr_ab_data
topic_model <- ecostructure_fit(counts, K = 3, tol = 10)

data("australia_birds")
data("australia_model")
ecostructure_plot_pie(omega = australia_model$omega,
                      coords = australia_birds$latlong, 
                      long_lim = c(110,160),
                      lat_lim = c(-50,-10),
                      color= c("orange", "red", "yellow", "deepskyblue", 
                               "chartreuse", "blue"))




                                        
