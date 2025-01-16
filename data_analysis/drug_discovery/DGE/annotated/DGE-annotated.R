#!/usr/bin/env Rscript

### Author: Jenny, Yiming
###
### Description: This script performs DGE analysis (single-cell level data)

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

### Read Seurat object
seurat_object_path <- paste0(sample_dir, "processed/Seurat.RDS")
seurat_object <- readRDS(seurat_object_path)

### DGE results will be stored in, e.g., 10x/PID1/DS1D/DS1D.1/analysis/DGE
output_dir <- paste0(sample_dir, "analysis/DGE")
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}
dir.create(paste0(dt_dir, ds_name))
dir.create(paste0(dt_dir, ds_name, "/", sample_name))
dir.create(paste0(dt_dir, ds_name, "/", sample_name, "/DGE_anno"))


### this is to change '/' in cell type names to '.'
seurat_object@meta.data$cell_type_annotation_class <- gsub(
  pattern = "/", 
  replacement = ".", 
  x = seurat_object@meta.data$cell_type_annotation_class
)
annotations <- seurat_object$cell_type_annotation_class
uniq_anno = unique(annotations)
Idents(seurat_object) <- "cell_type_annotation_class"
saveRDS(seurat_object, paste0(sample_dir,'/processed/Seurat_reanno.RDS'))

### Perform DGE analysis on different cell types
if (length(uniq_anno) == 1) {
	fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types_anno.tsv"), sep = "\t")
	cat("\n\n### Only one annotated cell type -- skipping DGE analysis (cell types).")
	cat("\n# Writing empty DGE_cell_types_dec.tsv to file.")
} else {
	# https://satijalab.org/seurat/archive/v3.1/future_vignette.html
	options(future.globals.maxSize = 5000 * 1024^2)
	DGE_cell_types <- FindAllMarkers(seurat_object, assay = "SCT", logfc.threshold = 0.2, min.pct = 0.1, verbose = FALSE)
	if (nrow(DGE_cell_types) == 0) {
		fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types_anno.tsv"), sep = "\t")
		cat("\n\n### DGE analysis (cell types) cannot be performed due to having too few spots in one/many of the cell types.")
		cat("\n# Writing empty DGE_cell_types_dec.tsv to file.")
	} else {
		DGE_cell_types$cluster <- as.character(DGE_cell_types$cluster)
		fwrite(DGE_cell_types, paste0(output_dir, "/DGE_cell_types_anno.tsv"), sep = "\t")
		
		for (cell_type in sort(unique(DGE_cell_types$cluster))) {
			cat(paste0(cell_type,' being anlayzed for DEG \n'))
			# cell_type <- "Malignant"
			DGE_cell_types_less <- DGE_cell_types[DGE_cell_types$cluster == cell_type,]
			DGE_cell_types_less <- DGE_cell_types_less[,c("gene", "avg_log2FC", "p_val", "p_val_adj")]
			colnames(DGE_cell_types_less) <- c("gene", "stat", "pval", "qval")
			fwrite(DGE_cell_types_less, paste0(dt_dir, ds_name, "/", sample_name, "/DGE_anno/", cell_type, ".txt"), sep = "\t")
		}
		cat("\n\n### DGE analysis (cell types) results written to file.")
	}
}
tt4 <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
cat(paste0("\n### Analysis completed; max memory consumed: ", as.character(tt4), "M -- [", Sys.time(), "]\n\n"))
