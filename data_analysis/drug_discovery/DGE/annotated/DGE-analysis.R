#!/usr/bin/env Rscript

### Author: Yiming Li
###
### Description: This script performs cluster-based and cell-type-based DGE analysis
### 	* Currently only tested on the cancer datasets
### 
### Example:
### chmod 755 DGE-analysis.R
### ./DGE-analysis.R > DGE-analysis.log

library(Seurat)
library(dplyr)
library(data.table)

### Read the list of DSIDs for use in our database
ds_list <- fread("~/stbase/DSID_list.txt", sep = "\t")

for (i in 1:nrow(ds_list)) {
	ds_name <- as.character(ds_list[i, "DSID"])
	ds_name <- paste0("DS", ds_name)
	tech <- as.character(ds_list[i, "Technology"])
	p_name <- paste0("PID", substr(ds_name, start = 3, stop = nchar(ds_name) - 1))
	
	ds_dir <- paste(c("/share/fsmresfiles/SpatialT", tech, p_name, ds_name), collapse = "/")
	if (!dir.exists(ds_dir)) {
		stop(paste0(ds_dir, " does not exist"))
	}
	if (!file.exists(paste0(ds_dir, "/metatable.tsv"))) {
		stop(paste0(ds_dir, "/metatable.tsv does not exist"))
	}
	ds_metatable <- read.table(paste0(ds_dir, "/metatable.tsv"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
	
	cat(paste(c("\n\n>>>>>>>> ", ds_dir, " <<<<<<<<"), collapse = ""))
	
	for (sample_name in ds_metatable$SampleID) {
		# sample_name <- ds_metatable$SampleID[1] # For testing only
		sample_dir <- paste(c(ds_dir, "/", sample_name), collapse = "")
		cat(paste(c("\n\n>>>>>>>> ", sample_dir, " <<<<<<<<"), collapse = ""))
		
		### DGE results will be stored in, e.g., 10x/PID1/DS1D/DS1D.1/analysis/DGE
		output_dir <- paste0(sample_dir, "/analysis/DGE")
		if (!dir.exists(output_dir)) {
			dir.create(output_dir)
		}
		
		### Read in the QC-ed + transformed + processed Seurat object (from process_visium_standard.R)
		seurat_object_tn_path <- paste0(sample_dir, "/processed/Seurat.RDS")
		seurat_object_tn <- readRDS(seurat_object_tn_path)
		cat("\n\n### Processed Seurat object read.")
		
		### !!! Assuming that the Seurat object has underwent SCT and clustering
		
		### Perform DGE analysis on the clusters
		DGE_clusters <- FindAllMarkers(seurat_object_tn, assay = "SCT", logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
		# Using Seurat's default:
		# logfc.threshold = 0.25 (Limit testing to genes which show, on average, at least 0.25-fold difference (log-scale) between the two groups of cells)
		# min.pct = 0.1 (Only test genes that are detected in a minimum fraction of 10% cells in either of the two populations)
		DGE_clusters_more <- DGE_clusters
		DGE_clusters <- DGE_clusters[DGE_clusters$avg_log2FC > 0.25 & DGE_clusters$pct.1 > 0.1,]
		
		DGE_clusters <- DGE_clusters[,c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
		# avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.
		# pct.1 : The percentage of cells where the gene is detected in the first group
		# pct.2 : The percentage of cells where the gene is detected in the second group
		# p_val_adj : Adjusted p-value, based on bonferroni correction using all genes in the dataset.
		DGE_clusters %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> DGE_clusters_top10
		fwrite(DGE_clusters_more, paste0(output_dir, "/DGE_clusters_unfiltered.tsv"), sep = "\t")
		fwrite(DGE_clusters, paste0(output_dir, "/DGE_clusters.tsv"), sep = "\t")
		fwrite(DGE_clusters_top10, paste0(output_dir, "/DGE_clusters_top10.tsv"), sep = "\t")
		cat("\n\n### DGE analysis (clusters) results written to file.")
		
		### Perform DGE analysis on the annotated cell types
		tmp <- seurat_object_tn
		annotations <- tmp[["cell_type_annotation"]]$cell_type_annotation
		names(annotations) <- rownames(tmp[["cell_type_annotation"]])
		Idents(tmp) <- annotations
		
		if (length(table(tmp[["cell_type_annotation"]])) == 1) {
			fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types_unfiltered.tsv"), sep = "\t")
			fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types.tsv"), sep = "\t")
			fwrite(data.frame(gene = character(0), cluster = integer(0), avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)), paste0(output_dir, "/DGE_cell_types_top10.tsv"), sep = "\t")
			cat("\n\n### Only one annotated cell type -- skipping DGE analysis (cell types).")
			cat("\n# Writing empty DGE_cell_types_unfiltered.tsv, DGE_cell_types.tsv and DGE_cell_types_top10.tsv to file.")
		} else {
			DGE_cell_types <- FindAllMarkers(tmp, assay = "SCT", logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
			DGE_cell_types_more <- DGE_cell_types
			DGE_cell_types <- DGE_cell_types[DGE_cell_types$avg_log2FC > 0.25 & DGE_cell_types$pct.1 > 0.1,]
			
			DGE_cell_types <- DGE_cell_types[,c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
			colnames(DGE_cell_types)[2] <- "cell_type"
			DGE_cell_types %>% group_by(cell_type) %>% top_n(n = 10, wt = avg_log2FC) -> DGE_cell_types_top10
			fwrite(DGE_cell_types_more, paste0(output_dir, "/DGE_cell_types_unfiltered.tsv"), sep = "\t")
			fwrite(DGE_cell_types, paste0(output_dir, "/DGE_cell_types.tsv"), sep = "\t")
			fwrite(DGE_cell_types_top10, paste0(output_dir, "/DGE_cell_types_top10.tsv"), sep = "\t")
			cat("\n\n### DGE analysis (cell types) results written to file.")
		}
	}
}
