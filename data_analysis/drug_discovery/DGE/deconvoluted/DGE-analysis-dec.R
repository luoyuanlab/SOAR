#!/usr/bin/env Rscript

### Author: Yiming
###
### Description: This script performs DGE analysis (using pseudo-cell data)

library(Seurat)
library(data.table)
library(dplyr)
library(stringr)

### Define paths and variables
args <- commandArgs(trailingOnly=TRUE)
st_dir <- "/projects/b1131/SpatialT"
dt_dir <- "/projects/b1131/SpatialT/drug-target/"
sample_dir <- args[1]
# sample_dir <- "/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/"
# sample_dir <- "/projects/b1131/SpatialT/DBiT-seq/PID150/DS150A/DS150A.GSM4096261/"

ds_name <- str_split(sample_dir, '/')[[1]][7]
tech <- str_split(sample_dir, '/')[[1]][5]
p_name <- str_split(sample_dir, '/')[[1]][6]
ds_dir <- paste(c(st_dir, tech, p_name, ds_name), collapse = "/")
sample_name <- str_split(sample_dir, '/')[[1]][8]

### Read all possible cell types
all_cell_types_df <- fread(paste0(sample_dir, "/analysis/deconvolution/all_cell_types.txt"))
all_cell_types <- all_cell_types_df$cell_type_s
ct_mapping <- all_cell_types_df$cell_type
names(ct_mapping) <- all_cell_types

### Read Seurat object
seurat_object_tn_path <- paste0(sample_dir, "processed/Seurat.RDS")
seurat_object_tn <- readRDS(seurat_object_tn_path)
spot_id_mapping <- seurat_object_tn@meta.data$new_spot_id
names(spot_id_mapping) <- as.character(rownames(seurat_object_tn@meta.data))

### Read cell type fractions
theta <- readRDS(paste0(sample_dir, "/analysis/deconvolution/BayesPrism_theta.RDS"))
rownames(theta) <- as.character(spot_id_mapping[rownames(theta)])

### DGE results will be stored in, e.g., 10x/PID1/DS1D/DS1D.1/analysis/DGE
output_dir <- paste0(sample_dir, "analysis/DGE")
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}
dir.create(paste0(dt_dir, ds_name))
dir.create(paste0(dt_dir, ds_name, "/", sample_name))
dir.create(paste0(dt_dir, ds_name, "/", sample_name, "/DGE_dec"))

### Create pseudo-cell-level object
if (file.exists(paste0(sample_dir, "analysis/deconvolution/binded_exp_Seurat.RDS"))) {
# if (FALSE) {
	binded_so <- readRDS(paste0(sample_dir, "analysis/deconvolution/binded_exp_Seurat.RDS"))
} else {
	ct_exp <- list()
	total_cells_each_ct <- list()
	for (cell_type in all_cell_types) {
		# cell_type <- "Malignant"
		# cell_type <- all_cell_types[6]
		exp_mat <- fread(paste0(sample_dir, "/analysis/deconvolution/counts_", cell_type, "_deconv_only.csv"))
		genes <- exp_mat$gene
		exp_mat$gene <- NULL
		exp_mat <- as.matrix(exp_mat)
		rownames(exp_mat) <- genes
		
		# Divide by fraction and leave out the spots with zero fraction of this cell type
		fractions <- theta[,as.character(ct_mapping[cell_type])]
		non_zero_fraction_spots <- names(which(fractions != 0))
		exp_mat <- t(apply(exp_mat, 1, function(x) x / fractions))
		exp_mat[is.na(exp_mat)] <- 0
		exp_mat <- exp_mat[, non_zero_fraction_spots, drop = FALSE]
		
		ct_exp[[cell_type]] <- exp_mat
		total_cells_each_ct[[cell_type]] <- ncol(exp_mat)
	}
	ct_exp <- do.call(cbind, ct_exp)
	colnames(ct_exp) <- paste0("sp", 1:ncol(ct_exp))
	saveRDS(ct_exp, paste0(sample_dir, "analysis/deconvolution/binded_exp.RDS"))
	
	# Remove "pseudo-cells" with zero deconvoluted expression in all the genes
	# Initially this was done because CellChat does not accept objects with zero total expression cells
	col_sums <- colSums(ct_exp)
	kept_pseudo_cells <- names(which(col_sums > 0))
	all_cell_types <- all_cell_types[total_cells_each_ct != 0]
	total_cells_each_ct <- total_cells_each_ct[total_cells_each_ct != 0]
	meta <- data.frame(labels = rep(all_cell_types, times = total_cells_each_ct))
	meta$cell_id <- paste0("sp", 1:ncol(ct_exp))
	meta <- meta[meta$cell_id %in% kept_pseudo_cells,]
	row.names(meta) <- meta$cell_id
	meta$cell_id <- NULL
	ct_exp_less <- ct_exp[,kept_pseudo_cells]
	
	rm(seurat_object_tn)
	binded_so <- CreateSeuratObject(ct_exp_less)
	binded_so[["cell_type"]] <- meta$labels
	binded_so <- SCTransform(binded_so, verbose = FALSE, return.only.var.genes = FALSE)
	saveRDS(binded_so, paste0(sample_dir, "analysis/deconvolution/binded_exp_Seurat.RDS"))
}

### Perform DGE analysis on different cell types
annotations <- binded_so[["cell_type"]]$cell_type
names(annotations) <- rownames(binded_so[["cell_type"]])
Idents(binded_so) <- annotations

if (length(table(binded_so[["cell_type"]])) == 1) {
	fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types_dec.tsv"), sep = "\t")
	cat("\n\n### Only one annotated cell type -- skipping DGE analysis (cell types).")
	cat("\n# Writing empty DGE_cell_types_dec.tsv to file.")
} else {
	# https://satijalab.org/seurat/archive/v3.1/future_vignette.html
	options(future.globals.maxSize = 5000 * 1024^2)
	DGE_cell_types <- FindAllMarkers(binded_so, assay = "SCT", logfc.threshold = 0.1, min.pct = 0.1, verbose = FALSE)
	if (nrow(DGE_cell_types) == 0) {
		fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types_dec.tsv"), sep = "\t")
		cat("\n\n### DGE analysis (cell types) cannot be performed due to having too few spots in one/many of the cell types.")
		cat("\n# Writing empty DGE_cell_types_dec.tsv to file.")
	} else {
		DGE_cell_types$cluster <- as.character(DGE_cell_types$cluster)
		fwrite(DGE_cell_types, paste0(output_dir, "/DGE_cell_types_dec.tsv"), sep = "\t")
		
		for (cell_type in sort(unique(DGE_cell_types$cluster))) {
			# cell_type <- "Malignant"
			DGE_cell_types_less <- DGE_cell_types[DGE_cell_types$cluster == cell_type,]
			DGE_cell_types_less <- DGE_cell_types_less[,c("gene", "avg_log2FC", "p_val", "p_val_adj")]
			colnames(DGE_cell_types_less) <- c("gene", "stat", "pval", "qval")
			fwrite(DGE_cell_types_less, paste0(dt_dir, ds_name, "/", sample_name, "/DGE_dec/", cell_type, ".txt"), sep = "\t")
		}
		cat("\n\n### DGE analysis (cell types) results written to file.")
	}
}
tt4 <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
cat(paste0("\n### Analysis completed; max memory consumed: ", as.character(tt4), "M -- [", Sys.time(), "]\n\n"))
