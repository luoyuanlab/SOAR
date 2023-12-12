### Author: Yiming Li
### Example usage:
### create_input_files.R $sample_dir

library(Seurat)
library(BayesPrism)
library(data.table)
library(dplyr)
library(stringr)

### Read the list of DSIDs for use in our database
args <- commandArgs(trailingOnly=TRUE)

st_dir <- "/projects/b1131/SpatialT"

sample_dir <- args[1]
# sample_dir <- "/projects/b1131/SpatialT/10x/PID5/DS5A/DS5A.12_151676"

### Read Seurat objects and deconvolution results
seurat_object <- readRDS(paste0(sample_dir, "/processed/Seurat.RDS"))
location <- fread(paste0(sample_dir, "/processed/coordinates.csv"))
output_dir <- paste0(sample_dir, "/analysis/deconvolution")
theta <- readRDS(paste0(output_dir, "/BayesPrism_theta.RDS"))
bp.res <- readRDS(paste0(output_dir, "/BayesPrism_results.RDS"))
meta <- seurat_object@meta.data
spot_id_mapping <- meta$new_spot_id
names(spot_id_mapping) <- rownames(meta)

gc()

### Hard-assigned labels
hard_labels <- apply(theta, 1, function(x) names(x)[which(x == max(x))])
hard_labels <- data.frame(cell_type_dec_max = as.character(hard_labels), spot_id = names(hard_labels))
meta$spot_id <- rownames(meta)
meta$row_id <- 1:nrow(meta)
meta <- merge(meta, hard_labels, by = "spot_id")
meta <- meta[order(meta$row_id),]
seurat_object[["cell_type_dec_max"]] <- meta$cell_type_dec_max
# sum(colnames(seurat_object) == meta$spot_id)
saveRDS(seurat_object, paste0(sample_dir, "/processed/Seurat.RDS"))

rm(seurat_object)
gc()

### Save deconvoluted cell-type-specific expressions
all_cell_types_s <- character(0)
all_cell_types <- character(0)
ct_i <- 1
for (cell_type in colnames(theta)) {
	# cell_type <- "CAFs"
	cell_type_s <- gsub("/", ".", cell_type)
	cell_type_s <- gsub(" ", ".", cell_type_s)
	cell_type_s <- gsub("-", ".", cell_type_s)
	cell_type_s <- gsub("\\*", ".", cell_type_s)
	cell_type_s <- gsub("\\+", ".", cell_type_s)
	
	ct_exp <- get.exp(bp = bp.res, state.or.type = "type", cell.name = cell_type)
	counts_df <- data.table(ct_exp)
	counts_df$spot <- spot_id_mapping[rownames(ct_exp)]
	counts_df <- transpose(counts_df, keep.names = "gene", make.names = "spot")
	keep_spots <- intersect(colnames(counts_df), location$barcode)
	keep_spots <- c("gene", keep_spots)
	counts_df <- counts_df[,..keep_spots]
	
	fwrite(counts_df, paste0(output_dir, "/counts_", cell_type_s, "_deconv_only.csv"), sep = ",")
	all_cell_types[ct_i] <- cell_type
	all_cell_types_s[ct_i] <- cell_type_s
	ct_i <- ct_i + 1
}
fwrite(data.table(cell_type = all_cell_types, cell_type_s = all_cell_types_s), paste0(output_dir, "/all_cell_types.txt"), sep = "\t")
