
############ R script for preparing the data  ######################

counts_data <- read.csv("../data/Counts_5_18_17.csv", row.names = 1)
sites_metadata <- read.csv("../data/Metadata_Him_Bird_Grids.csv", row.names = 1)
species_metadata <- read.csv("../data/Species_Morphology_Him_Counts.csv", row.names = 1)

counts_data_2 <- counts_data[match(rownames(sites_metadata), rownames(counts_data)),]

counts_data_3 <- counts_data_2[,match(rownames(species_metadata), colnames(counts_data_2))]

reads <- t(counts_data_3)

HimalayanBirdsData <- new("ExpressionSet",
                          exprs = as.matrix(reads),
                          phenoData = new("AnnotatedDataFrame",
                                         data = sites_metadata),
                          featureData = new("AnnotatedDataFrame",
                                         data = species_metadata),
                   experimentData = new("MIAME",
                                        title = "Himalayan Birds Data"))

save(HimalayanBirdsData, file = "../data/HimalayanBirdsData.rda")



exprs(HimalayanBirdsData)
pData(phenoData(HimalayanBirdsData))


data <- get(load(system.file("extdata", "HimalayanBirdsData.rda",
                             package = "ecostructure")))
assayData(HimalayanBirdsData)
