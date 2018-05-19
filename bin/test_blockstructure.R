
data(himal)
new_counts <- t(exprs(data));
metadata <- pData(data);
elevation_metadata=metadata$Elevation;
east_west_dir = metadata$WorE;
topic_clus <- maptpx::topics(new_counts, K=2, tol=0.1)
omega <- topic_clus$omega

blocker <- as.character(east_west_dir)
blocker[sample(1:30, 7, replace = FALSE)] = "N"
blocker <- as.factor(blocker)


block_structure(omega,
               blocker_metadata = blocker,
               order_metadata = elevation_metadata)


