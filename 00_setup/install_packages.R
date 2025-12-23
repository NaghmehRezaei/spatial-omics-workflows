# Install required packages for spatial omics workflows

install.packages(c("Seurat", "tidyverse", "patchwork"))
BiocManager::install(c("SpatialExperiment", "GeoMxTools",
                       "fgsea", "clusterProfiler"))
