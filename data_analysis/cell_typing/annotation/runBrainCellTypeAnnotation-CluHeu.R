#!/usr/bin/env Rscript

### Author: Yiming Li
###
### Description: This script automatically annotates the cell types of brain spatial transcriptomics datasets using a cluster-based approach guided by some heuristics
### 	* It reads a table listing the DSID, technology, and species (brain_DSID_list.txt) and loop over its rows
### 
### Usage: ./runBrainCellTypeAnnotation-CluHeu.R
###
### Example: To annotate all the brain datasets in our database:
# chmod 755 runBrainCellTypeAnnotation-CluHeu.R
# ./runBrainCellTypeAnnotation-CluHeu.R > runBrainCellTypeAnnotation-CluHeu.log



library(SingleCellExperiment)
library(Seurat)
library(SingleR)
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
library(AUCell)



### Read/Create variables related to cancer annotation
human_ref_dir <- "/share/fsmresfiles/SpatialT/ref/Brain/Adult/aibs_human_ctx_smart-seq/"
ref_data_sce_neuronal_human <- readRDS(paste0(human_ref_dir, "aibs_human_ctx_smart-seq_neuronal.RDS"))
ref_data_sce_non_neuronal_human <- readRDS(paste0(human_ref_dir, "aibs_human_ctx_smart-seq_non_neuronal.RDS"))

mouse_ref_dir <- "/share/fsmresfiles/SpatialT/ref/Brain/Adult/aibs_mouse_ctx-hpf_10x/"
ref_data_sce_neuronal_mouse <- readRDS(paste0(mouse_ref_dir, "aibs_mouse_ctx-hpf_10x_neuronal.RDS"))
ref_data_sce_non_neuronal_mouse <- readRDS(paste0(mouse_ref_dir, "aibs_mouse_ctx-hpf_10x_non_neuronal.RDS"))

### Define functions
clustering_seurat <- function(data) {
	data <- RunPCA(data, assay = "SCT", verbose = FALSE)
	data <- FindNeighbors(data, reduction = "pca", dims = 1:30, verbose = FALSE)
	data <- FindClusters(data, resolution = 1.2, verbose = FALSE)
	## resolution: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html advises 0.4-1.2 for around 3K cells
	data <- RunUMAP(data, reduction = "pca", dims = 1:30, verbose = FALSE)
	return(data)
}



##################



target_list <- fread("~/stbase/brain_DSID_list.txt", sep = "\t")

for (i in 1:nrow(target_list)) {
	# i <- 1 # For testing only
	target_ds_name <- as.character(target_list[i, "DSID"])
	target_species <- as.character(target_list[i, "Species"])
	target_tech <- as.character(target_list[i, "Technology"])
	### Do not differentiate between adult and non-adult datasets
	target_p_name <- paste0("PID", substr(target_ds_name, start = 3, stop = nchar(target_ds_name) - 1))
	
	target_ds_dir <- paste(c("/share/fsmresfiles/SpatialT", target_tech, target_p_name, target_ds_name), collapse = "/")
	if (!dir.exists(target_ds_dir)) {
		cat(paste0(target_ds_dir, " does not exist"))
		next
	}
	if (!file.exists(paste0(target_ds_dir, "/metatable.tsv"))) {
		cat(paste0(target_ds_dir, "/metatable.tsv does not exist"))
		next
	}
	target_ds_metatable <- read.table(paste0(target_ds_dir, "/metatable.tsv"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
	
	cat(paste(c("\n\n>>>>>>>> ", target_ds_dir, " <<<<<<<<"), collapse = ""))
	
	for (target_sample_name in target_ds_metatable$SampleID) {
		# target_sample_name <- target_ds_metatable$SampleID[1] # For testing only
		target_sample_dir <- paste(c(target_ds_dir, "/", target_sample_name), collapse = "")
		
		### Needed since some of the datasets were incorrectly prepared
		if (!dir.exists(target_ds_dir)) {
			cat(paste0(target_sample_dir, " does not exist"))
			next
		}
		
		cat(paste(c("\n\n>>>>>>>> ", target_sample_dir, " <<<<<<<<"), collapse = ""))
		
		### Annotation results will be stored in, e.g., 10x/PID1/DS1D/DS1D.1/analysis/annotation
		output_dir_1 <- paste0(target_sample_dir, "/analysis")
		output_dir <- paste0(output_dir_1, "/annotation")
		if (!dir.exists(output_dir_1)) {
			dir.create(output_dir_1)
		}
		if (!dir.exists(output_dir)) {
			dir.create(output_dir)
		}
		
		### Read in the QC-ed + transformed + processed Seurat object (from process_visium_standard.R)
		seurat_object_tn_path <- paste0(target_sample_dir, "/processed/Seurat.RDS")
		seurat_object_tn <- readRDS(seurat_object_tn_path)
		cat("\n\n### Processed Seurat object read.")
		if ("cell_type_annotation_class" %in% colnames(seurat_object_tn@meta.data) & "cell_type_annotation" %in% colnames(seurat_object_tn@meta.data)) {
			cat("\n# Seurat object already has cell type annotation results saved, skipping this sample.")
			next
		}
	
		### Check if the Seurat object is after SCT
		if (!"SCT" %in% names(seurat_object_tn)) {
			cat("\n# [Warning] Seurat object does not have an \"SCT\" assay, renaming Spatial to SCT.")
			### Some MERFISH files prepared by Yawei do not have an SCT assay
			### However, they are already normalized so we can rename the assay to SCT for ease of subsequent analysis
			file.copy(seurat_object_tn_path, paste0(target_sample_dir, "/processed/Seurat_bk.RDS"), overwrite = TRUE)
			seurat_object_tn <- RenameAssays(object = seurat_object_tn, Spatial = "SCT")
			saveRDS(seurat_object_tn, file = seurat_object_tn_path)
		}
		
		### Check if the ST data underwent dimensionality reduction and clustering
		if (!"umap" %in% names(seurat_object_tn)) {
			cat("\n# [Warning] Seurat object does not have UMAP dimensional reduction calculated, skipping this sample.")
			next
		}
	
		cat("\n\n### Start annotation")
		coords <- seurat_object_tn@meta.data
		coords <- tibble::rownames_to_column(coords, "spot")
		coords$seurat_clusters <- as.character(coords$seurat_clusters)
		
		### Read marker gene lists and cell class/subclass mappings
		if (target_species == "Human") {
			load(paste0(human_ref_dir, "supp.RData"))
		} else {
			load(paste0(mouse_ref_dir, "supp.RData"))
		}
		
		### Use AUCell to classify clusters into neuronal vs non-neuronal
		cells_rankings <- AUCell_buildRankings(GetAssayData(seurat_object_tn, assay = "SCT"), verbose = FALSE)
		### The AUC value represents the fraction of genes, within the top 20% genes in the ranking, that are included in the signature
		runAUCell <- tryCatch(
			{
				cells_AUC <- AUCell_calcAUC(gene_lists, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.2, verbose = FALSE)
			},
			error = function(e) e
		)
		if (inherits(runAUCell, "error")){
			cat(paste(c("\n# AUCell failed, ST dataset may have a small number of genes. Creating a subsetted gene list."), collapse = ""))
			data_genes <- rownames(seurat_object_tn)
			for (gene_list_name in names(gene_lists)) {
				gene_list <- gene_lists[[gene_list_name]]
				gene_lists[[gene_list_name]] <- gene_list[gene_list %in% data_genes]
			}
			gene_lists <- gene_lists[lapply(gene_lists, length) > 0]
			cells_AUC <- AUCell_calcAUC(gene_lists, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.2, verbose = FALSE)
		}
		
		clusters <- sort(unique(coords$seurat_clusters))
		cluster_class <- data.frame(seurat_clusters = clusters, subclass = rep(NA, length(clusters)))
		for (i in 1:length(clusters)) {
			cluster <- clusters[i]
			spots <- coords[coords$seurat_clusters == cluster, "spot"]
			tmpmat <- cells_AUC[,colnames(cells_AUC) %in% spots]
			tmp <- rowSums(tmpmat@assays@data@listData$AUC)
			cluster_class$subclass[i] <- names(which(tmp == max(tmp)))
		}
		cluster_class$class <- cell_type_names[cluster_class$subclass]
		coords <- join(coords, cluster_class, by = "seurat_clusters")
		coords$class[coords$class == "glutamatergic"] <- "neuronal"
		coords$class[coords$class == "GABAergic"] <- "neuronal"
		coords <- coords[, c("spot", "class")]
		
		### Use SingleR to annotate the subclasses
		cluster_results <- seurat_object_tn[["seurat_clusters"]]$seurat_clusters
		if (target_species == "Human") {
			annotation_neuronal <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce_neuronal_human, clusters = cluster_results, labels = ref_data_sce_neuronal_human$label, de.method="wilcox")
			annotation_non_neuronal <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce_non_neuronal_human, clusters = cluster_results, labels = ref_data_sce_non_neuronal_human$label, de.method="wilcox")
		} else {
			annotation_neuronal <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce_neuronal_mouse, clusters = cluster_results, labels = ref_data_sce_neuronal_mouse$label, de.method="wilcox")
			annotation_non_neuronal <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce_non_neuronal_mouse, clusters = cluster_results, labels = ref_data_sce_non_neuronal_mouse$label, de.method="wilcox")
		}
		
		### Process SingleR annotation results
		tmpdf <- data.frame(spot = colnames(seurat_object_tn), annot_neu = annotation_neuronal$labels[cluster_results])
		coords <- join(coords, tmpdf, by = "spot")
		tmpdf2 <- data.frame(spot = colnames(seurat_object_tn), annot_non = annotation_non_neuronal$labels[cluster_results])
		coords <- join(coords, tmpdf2, by = "spot")
		coords$cell_type_annotation <- ifelse(coords$class == "neuronal", coords$annot_neu, coords$annot_non)
		coords$cell_type_annotation_class <- cell_type_names[coords$cell_type_annotation]
		
		### Save annotation results to Seurat object
		seurat_object_tn[["cell_type_annotation"]] <- coords$cell_type_annotation
		seurat_object_tn[["cell_type_annotation_class"]] <- coords$cell_type_annotation_class
		saveRDS(seurat_object_tn, file = seurat_object_tn_path)
		cat("\n\n### Seurat.RDS overwritten with annotation results.")
		
		### Visualize annotated cell types
		pdf(paste0(output_dir, "/cell_type_annotation.pdf"))
		print(SpatialDimPlot(seurat_object_tn))
		print(DimPlot(seurat_object_tn, reduction = "umap"))
		print(SpatialDimPlot(seurat_object_tn, group.by = "cell_type_annotation_class"))
		print(DimPlot(seurat_object_tn, reduction = "umap", group.by = "cell_type_annotation_class"))
		print(SpatialDimPlot(seurat_object_tn, group.by = "cell_type_annotation"))
		print(DimPlot(seurat_object_tn, reduction = "umap", group.by = "cell_type_annotation"))
		dev.off()
		cat("\n\n### Visualization of annotated cell types completed.")
	}
}
