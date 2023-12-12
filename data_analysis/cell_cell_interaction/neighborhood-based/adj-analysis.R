#!/usr/bin/env Rscript

### Author: Yiming Li
###
### Description: This script performs neighborhood-based analysis
### Usage: adj-analysis.R $sample_dir

library(Seurat)
library(plyr)
library(data.table)
library(stringr)

### Read the list of DSIDs for use in our database
args <- commandArgs(trailingOnly=TRUE)

st_dir <- "/projects/b1131/SpatialT"
# st_dir <- "/share/fsmresfiles"

sample_dir <- args[1]
# sample_dir <- "/projects/b1131/SpatialT/10x/PID5/DS5A/DS5A.06_151670"
# sample_dir <- "/share/fsmresfiles/SpatialT/10x/PID5/DS5A/DS5A.12_151676"
ds_name <- str_split(sample_dir, '/')[[1]][7]
tech <- str_split(sample_dir, '/')[[1]][5]
p_name <- str_split(sample_dir, '/')[[1]][6]
ds_dir <- paste(c(st_dir, tech, p_name, ds_name), collapse = "/")
sample_name <- str_split(sample_dir, '/')[[1]][8]

### Functions
get_dist <- function(row1, col1, row2, col2) {
	return(sqrt((row1-row2)^2 + (col1-col2)^2))
}

get_adjacency <- function(row1, col1, row2, col2) {
	return(sum(abs(row1-row2) <= 1 & abs(col1-col2) <= 1))
}

get_dist_from_cell_type <- function(row, col, cell_type, one_cell_type_dfs) {
	cell_types <- names(one_cell_type_dfs)
	distances <- rep(NA, length(cell_types) * 3)
	names(distances) <- c(paste0(cell_types, "_min"), paste0(cell_types, "_median"), paste0(cell_types, "_adjacent"))
	for (compared_cell_type in cell_types) {
		if (cell_type == compared_cell_type) {
			next
		}
		compared_cell_type_df <- one_cell_type_dfs[[compared_cell_type]]
		tmp <- get_dist(row, col, compared_cell_type_df$row, compared_cell_type_df$col)
		distances[paste0(compared_cell_type, "_min")] <- min(tmp)
		distances[paste0(compared_cell_type, "_median")] <- median(tmp)
		distances[paste0(compared_cell_type, "_adjacent")] <- get_adjacency(row, col, compared_cell_type_df$row, compared_cell_type_df$col)
	}
	return(distances)
}

### Start adjacency-based analysis
### Results will be stored in, e.g., 10x/PID1/DS1D/DS1D.1/analysis/Distance
output_dir <- paste0(sample_dir, "/analysis/Distance")
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}

### Check if the analysis has already been done
if (file.exists(paste0(output_dir, "/cci_adj_results.tsv"))) {
	tmpdf <- fread(paste0(output_dir, "/cci_adj_results.tsv"))
	if ("wilcox_pval_two_sided" %in% colnames(tmpdf)) {
		### Save old results to a separate file
		file.copy(paste0(output_dir, "/cci_adj_results.tsv"), paste0(output_dir, "/cci_adj_results_old.tsv"), overwrite = TRUE)
	} else if ("avg_log2FC" %in% colnames(tmpdf)) {
		stop("\n### Adjacency analysis already done.\n")
	}
}

### Read in the QC-ed + transformed + processed Seurat object (from process_visium_standard.R)
seurat_object_tn_path <- paste0(sample_dir, "/processed/Seurat.RDS")
seurat_object_tn <- readRDS(seurat_object_tn_path)
cat("\n\n### Processed Seurat object read.\n\n### Starting adjacency analysis...\n")
### !!! Assuming that the Seurat object has underwent SCT and clustering

### Get the table with calculated distance-based statistics
if (file.exists(paste0(output_dir, "/distance_stats_dec_max.tsv"))) {
	coords <- fread(paste0(output_dir, "/distance_stats_dec_max.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	coords$spot <- as.character(coords$spot) ### In case the spot IDs are numeric
	cell_types <- sort(unique(seurat_object_tn@meta.data$cell_type_dec_max))
	cat("\n### Distance-based statistics read from file.\n")
} else {
	### Get coordinates and cell type annotations
	annotations <- seurat_object_tn@meta.data$cell_type_dec_max
	if ("slice1" %in% names(seurat_object_tn@images)) {
		coords <- seurat_object_tn@images$slice1@coordinates ### Visium
		coords <- coords[, c("row", "col")] ### Use the spot coordinates
	} else {
		coords <- seurat_object_tn@images$image@coordinates ### Others
		if ("x" %in% colnames(coords)) {
			coords <- coords[, c("x", "y")] ### Use the spot coordinates
		} else {
			### Some prepared MERFISH datasets did not follow the naming standard
			coords <- coords[, c("xcoord", "ycoord")] ### Use the spot coordinates
		}
	}
	colnames(coords) <- c("row", "col")
	coords <- tibble::rownames_to_column(coords, "spot")
	coords$annotation <- annotations
	
	### Get distance metrics, currently supporting:
	### * Minimum distance from another cell type (the same cell type -- marked as NA)
	### * Median distance from another cell type (the same cell type -- marked as NA)
	### * Whether adjacent to another cell type (0 = FALSE, 1 = TRUE; the same cell type -- marked as NA)
	cell_types <- sort(unique(annotations))
	one_cell_type_dfs <- list()
	for (cell_type in cell_types) {
		one_cell_type_dfs[[cell_type]] <- coords[coords$annotation == cell_type,]
	}
	df_colnames <- c("spot", paste0(cell_types, "_min"), paste0(cell_types, "_median"), paste0(cell_types, "_adjacent"))
	distance_stats_df <- data.frame(matrix(NA, ncol = length(df_colnames), nrow = 0))
	colnames(distance_stats_df) <- df_colnames
	for (i in 1:nrow(coords)) {
		row <- coords$row[i]
		col <- coords$col[i]
		cell_type <- coords$annotation[i]
		distances <- get_dist_from_cell_type(row, col, cell_type, one_cell_type_dfs)
		tmpdf <- data.frame(matrix(distances, ncol = length(df_colnames)-1, nrow = 1))
		rownames(tmpdf) <- coords$spot[i]
		tmpdf <- tibble::rownames_to_column(tmpdf, "spot")
		colnames(tmpdf) <- df_colnames
		distance_stats_df <- rbind(distance_stats_df, tmpdf)
	}
	coords <- join(coords, distance_stats_df, by = "spot")
	coords$row <- NULL
	coords$col <- NULL
	fwrite(coords, paste0(output_dir, "/distance_stats_dec_max.tsv"), sep = "\t")
	cat("\n### Distance-based statistics calculated.\n")
	
	### Ensure that coords is a data table
	coords <- fread(paste0(output_dir, "/distance_stats_dec_max.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	cell_types <- sort(unique(seurat_object_tn@meta.data$cell_type_dec_max))
}

### Test if in all cells of cell_type1, the gene expression levels in those adjacent / not adjacent to cell_type2 differ significantly
cci_adj_results_df <- data.frame(gene = character(0), cell_type1 = character(0), cell_type2 = character(0), adjacency = character(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val = numeric(0), p_val_adj = numeric(0))
cat("\n### Adjacency analysis started.\n")

if (length(cell_types) == 1) {
	cat(paste0("\n# Only one cell type present in the sample. Skipping this sample"))
	next
}
cat(paste0("\n# ", as.character(length(cell_types)), " cell types in total."))

for (cell_type in cell_types) {
	# cell_type <- cell_types[1]
	coords_cell_type <- coords[coords$annotation == cell_type,]
	other_cell_types <- cell_types[!cell_types == cell_type]
	seurat_cell_type <- seurat_object_tn[, coords_cell_type$spot]
	meta <- seurat_cell_type@meta.data
	meta <- tibble::rownames_to_column(meta, "spot")
	
	for (other_cell_type in other_cell_types) {
		# other_cell_type <- other_cell_types[1]
		tmp <- paste0(other_cell_type, "_adjacent")
		adjacent_cells <- coords_cell_type$spot[coords_cell_type[,..tmp] > 0]
		not_adjacent_cells <- coords_cell_type$spot[coords_cell_type[,..tmp] == 0]
		if (length(adjacent_cells) == 0 | length(not_adjacent_cells) == 0) {
			### All cell_type cells are adjacent to or not adjacent to other_cell_type cells
			### Skip this combination
			next
		}
		adj_df <- rbind(data.frame(spot = adjacent_cells, adjacency = "Adjacent"), data.frame(spot = not_adjacent_cells, adjacency = "Not adjacent"))
		adj_df <- join(meta, adj_df, by = "spot")
		adjacency <- adj_df$adjacency
		names(adjacency) <- adj_df$spot
		Idents(seurat_cell_type) <- adjacency
		
		### Thresholds decided based on: https://www.nature.com/articles/s41467-019-12266-7
		### Genes with absolute log2 fold change threshold > 0.1 and expressed in at least 10% of the cells are considered
		DGE_adjacency <- FindAllMarkers(seurat_cell_type, assay = "SCT", logfc.threshold = 0.1, min.pct = 0.1, verbose = FALSE)
		if (nrow(DGE_adjacency) == 0) {
			### No DGE genes found
			next
		}
		DGE_adjacency$cell_type1 <- cell_type
		DGE_adjacency$cell_type2 <- other_cell_type
		DGE_adjacency$adjacency <- as.character(DGE_adjacency$cluster)
		
		DGE_adjacency <- DGE_adjacency[,c("gene", "cell_type1", "cell_type2", "adjacency", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
		cci_adj_results_df <- rbind(cci_adj_results_df, DGE_adjacency)
		gc()
	}
	tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
	cat(paste0("\n# Cell type [", cell_type, "] done; memory consumed: ", as.character(tt), "M -- [", Sys.time(), "]\n\n"))
}

fwrite(cci_adj_results_df, paste0(output_dir, "/cci_adj_results.tsv"), sep = "\t")
cat("\n\n### Adjacency analysis completed. Results written to cci_adj_results.tsv.\n")

