library(Seurat)
library(data.table)
library(stringr)

sample_dir <- commandArgs(trailingOnly=TRUE)

coord = read.csv(paste0(sample_dir, "/processed/coordinates.csv"))
optimal_n <- dir(paste0(sample_dir, "analysis/clustering"), pattern = "obsm_STAGATE_cluster_")
optimal_n <- str_split(str_split(optimal_n, "_cluster_")[[1]][2], '.csv')[[1]][1]

so <- readRDS(paste0(sample_dir, "processed/Seurat.RDS"))
meta <- so@meta.data
meta$row <- 1:nrow(meta)

clustering_results <- fread(paste0(sample_dir, "analysis/clustering/obsm_STAGATE_cluster_", optimal_n, ".csv"))
colnames(clustering_results) <- c("new_spot_id", "STAGATE_cluster")

meta <- merge(meta, clustering_results, by = "new_spot_id")
meta <- meta[order(meta$row),]

# keep only seurat obs where they are present in the coordinates column, this is the obs in which clustering are done on
so <- so[, so@meta.data$new_spot_id %in% coord$barcode]
so[["STAGATE_cluster"]] <- meta$STAGATE_cluster
saveRDS(so, file = paste0(sample_dir, "processed/Seurat.RDS"))


