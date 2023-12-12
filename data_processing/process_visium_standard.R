### The purpose of this script is to process the 10x Visium data stored in standard structure (HDF5 + spatial/).
### The dataset IDs to be processed are read from a preexisting file DSID_list_10x_standard.txt in the shared folder. 
### This script will loop through the datasets, read, transform, and process the datasets, and save the processed data as well as a meta table for each dataset.
### 
### Usage: Rscript --no-save process_visium_standard.R > process_visium_standard.log
### 
### Author: Saya Dennis; Yiming Li; Sanaz Ghotbaldini

library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
# library(SeuratDisk)

dn <- "/projects/b1131/SpatialT/10x"
args <- commandArgs(trailingOnly=TRUE)
dsid <- args[1]

#### Get PID
dsid <- paste0("DS", dsid)
pid <- paste0("PID", substr(dsid, start=3, stop=nchar(dsid)-1)) # e.g. "PID1"
#### Create empty data frame for meta table 
meta <- data.frame(DSID = character(0), SampleID = character(0), Nspots = integer(0), Nspots_postQC = integer(0), Ngenes = integer(0), Ngenes_postQC = integer(0), Condition = character(0))
#### Set dataset directory and get a list of samples 
dsdir <- paste0(dn, "/", pid, "/", dsid) # e.g. "/projects/b1131/SpatialT/10x/PID1/DS1A"
sampleids <- dir(dsdir, pattern = dsid) # e.g. list of elements like "DS1A.1"

cat(paste(c("\n>>>>>>>> Dataset [", dsdir, "] started <<<<<<<<\n"), collapse = ""))
for (sampleid in sampleids) {
	#### Read data
	if (file.exists(paste0(dsdir, "/", sampleid, "/original/filtered_feature_bc_matrix.h5"))) {
		if (file.exists(paste0(dsdir, "/", sampleid, "/original/spatial/tissue_lowres_image.png"))) {
			data <- Load10X_Spatial(paste0(dsdir, "/", sampleid, "/original"), assay = "Spatial")
		} else {
			cat(paste(c("### Sample: ", sampleid, " [tissue_lowres_image.png] not found, using [tissue_hires_image.png]\n"), collapse = ""))
			img <- Read10X_Image(paste0(dsdir, "/", sampleid, "/original/spatial"), image.name = 'tissue_hires_image.png')
			data <- Load10X_Spatial(paste0(dsdir, "/", sampleid, "/original"), image = img, assay = "Spatial")
			data@images$slice1@scale.factors$lowres <- data@images$slice1@scale.factors$hires
		}

	} else if (file.exists(paste0(dsdir, "/", sampleid, "/original/Seurat.RDS"))) {
	  cat(paste(c("### Sample: ", sampleid, " [filtered_feature_bc_matrix.h5] not found, using Seurat object\n"), collapse = ""))
	  data <- readRDS(paste0(dsdir, "/", sampleid, "/original/Seurat.RDS"))
	  image_key <- as.character(names(data@images))
	  names(data@images)[names(data@images) == image_key] <- "slice1"
	  
	} else {
		cat(paste(c("### Sample: ", sampleid, " [filtered_feature_bc_matrix.h5] not found, using MEX files\n"), collapse = ""))
		if (file.exists(paste0(dsdir, "/", sampleid, "/original/spatial/tissue_lowres_image.png"))) {
			counts <- Read10X(paste0(dsdir, "/", sampleid, "/original"))
			data <- CreateSeuratObject(counts, assay = "Spatial")
			img <- Read10X_Image(paste0(dsdir, "/", sampleid, "/original/spatial"))
			img <- img[Cells(data)]
			DefaultAssay(img) <- DefaultAssay(data)
			data[["slice1"]] <- img
		} else {
			cat(paste(c("### Sample: ", sampleid, " [tissue_lowres_image.png] not found, using [tissue_hires_image.png]\n"), collapse = ""))
			counts <- Read10X(paste0(dsdir, "/", sampleid, "/original"))
			data <- CreateSeuratObject(counts, assay = "Spatial")
			img <- Read10X_Image(paste0(dsdir, "/", sampleid, "/original/spatial"), image.name = 'tissue_hires_image.png')
			img <- img[Cells(data)]
			DefaultAssay(img) <- DefaultAssay(data)
			data[["slice1"]] <- img
			data@images$slice1@scale.factors$lowres <- data@images$slice1@scale.factors$hires
		}
	}
	
	nspots <- ncol(data)
	ngenes <- nrow(data)
	
	#### Spot QC
	## Step 1. Remove the spots with total UMI count < 500 / the total number of genes < 500 / >= 25% mitochondrial reads
	data[["percent_mt"]] <- PercentageFeatureSet(data, "^MT-")
	data <- data[, data$nCount_Spatial >= 500 & data$nFeature_Spatial >= 500 & data$percent_mt < 25]
	## Step 2.
	## Remove the spots with total UMI count < median(total UMI count) - 3 * SD(total UMI count).
	## Remove the spots with total number of genes < median(total number of genes) - 3 * SD(total number of genes)
	data <- data[, data$nCount_Spatial >= median(data$nCount_Spatial) - 3 * sqrt(var(data$nCount_Spatial)) & data$nFeature_Spatial >= median(data$nFeature_Spatial) - 3 * sqrt(var(data$nFeature_Spatial))]
	
	#### Gene QC
	counts <- data.frame(GetAssayData(object = data, assay = "Spatial", slot = "counts"))
	counts <- counts > 0
	n_spots_per_gene <- rowSums(counts)
	data <- data[n_spots_per_gene >= 5,]
	
	nspots_qc <- ncol(data)
	ngenes_qc <- nrow(data)

	#### Exclude sample if there are < 50 spots after QC
	if (nspots_qc < 50) {
		cat(paste(c("### [", sampleid, "] has < 50 spots after QC, excluded\n"), collapse = ""))
		next
	}
	
	#### Append to the dataset metatable
	meta <- rbind(meta, data.frame(DSID = dsid, SampleID = sampleid, Nspots = nspots, Nspots_postQC = nspots_qc, Ngenes = ngenes, Ngenes_postQC = ngenes_qc, Condition = NA))
	
	#### Transform and process data
	data <- SCTransform(data, assay = "Spatial", verbose = FALSE)
	
	### Dimensionality reduction and clustering
	data <- RunPCA(data, assay = "SCT", verbose = FALSE)
	data <- FindNeighbors(data, reduction = "pca", dims = 1:30, verbose = FALSE)
	data <- FindClusters(data, resolution = 1.2, verbose = FALSE)
	## resolution: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html advises 0.4-1.2 for around 3K cells
	data <- RunUMAP(data, reduction = "pca", dims = 1:30, verbose = FALSE)
	
	### Save processed Seurat object
	## "processed" folder should already exist for all datasets but just in case
	if (!dir.exists(paste0(dsdir, "/", sampleid, "/processed/"))) {
		dir.create(paste0(dsdir, "/", sampleid, "/processed/"))
	}
	saveRDS(data, file = paste0(dsdir, "/", sampleid, "/processed/Seurat.RDS"))

	## Rename variables for the below code that is pasted from rewrite_text_files.R ##
	sample_dir <- paste0(dsdir, "/", sampleid)

	seurat_object <- readRDS(paste0(sample_dir, "/processed/Seurat.RDS"))
	
	#### Create spot IDs following R's column name requirements
	spotmeta <- seurat_object@meta.data
	spotmeta$new_spot_id <- paste0("sp", 1:nrow(spotmeta))
	seurat_object[["new_spot_id"]] <- spotmeta$new_spot_id
	saveRDS(seurat_object, paste0(sample_dir, "/processed/Seurat.RDS"))
	
	spot_id_mapping <- spotmeta$new_spot_id
	names(spot_id_mapping) <- rownames(spotmeta)
	
	#### Prepare data for deconvolution
	if ("Spatial" %in% names(seurat_object@assays)) {
		counts <- GetAssayData(object = seurat_object, assay = "Spatial", slot = "counts")
	} else {
		### MERFISH datasets
		counts <- GetAssayData(object = seurat_object, assay = "SCT", slot = "counts")
	}
	counts <- as.matrix(counts)
	gene_names <- rownames(counts)
	counts_t <- data.table(counts)
	counts_t$gene <- gene_names
	gc()
	counts_t <- transpose(counts_t, keep.names = "cell", make.names = "gene") ### Takes long
	gc()
	cell_names <- counts_t$cell
	counts_t$cell <- NULL ### Remove the gene name column
	counts_t <- as.matrix(counts_t)
	rownames(counts_t) <- cell_names
	saveRDS(counts_t, paste0(sample_dir, "/processed/bk.dat.RDS"))
	gc()
	
	#### Retrieve QC-ed counts and coordinates
	
	## Counts
	counts <- GetAssayData(object = seurat_object, assay = "SCT", slot = "counts")
	counts <- as.matrix(counts)
	colnames(counts) <- as.character(spot_id_mapping[colnames(counts)])
	counts_df <- tibble::rownames_to_column(data.frame(counts), "gene")
	
	## Coordinates
	if ("slice1" %in% names(seurat_object@images)) {
		location <- seurat_object@images$slice1@coordinates ### Visium
		location <- location[, c("col", "row")] ### Use the spot coordinates
	} else {
		location <- seurat_object@images$image@coordinates ### Others
		if ("x" %in% colnames(location)) {
			location <- location[, c("x", "y")] ### Use the spot coordinates
		} else {
			### Some prepared MERFISH datasets did not follow the naming standard
			location <- location[, c("xcoord", "ycoord")] ### Use the spot coordinates
			colnames(location) <- c("x", "y")
		}
	}
	location$barcode <- rownames(location)
	colnames(location) <- c("x", "y", "barcode")
	location$barcode <- spot_id_mapping[location$barcode]
	
	## Sometimes (e.g. /share/fsmresfiles/SpatialT/10x/PID153/DS153A/DS153A.1), the number of spots in seurat_object@images is different from ncol(seurat_object)
	## Not sure why this is the case
	counts_spots <- colnames(counts_df)
	location_spots <- location$barcode
	keep_spots <- intersect(counts_spots, location_spots)
	counts_df <- counts_df[,c("gene", keep_spots)]
	rownames(location) <- location$barcode
	location <- location[keep_spots,]
	
	#### Write counts, coordinates, and meta_spots to file
	write.table(counts_df, file = paste0(sample_dir, "/processed/counts.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
	write.table(location, file = paste0(sample_dir, "/processed/coordinates.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
	saveRDS(list("counts" = counts_df, "coordinates" = location), file = paste0(sample_dir, "/processed/data_frames.RDS"))
	
	#### Prepare data for deconvolution (relative counts)
	if ("Spatial" %in% names(seurat_object@assays)) {
		DefaultAssay(seurat_object) <- "Spatial"
		seurat_object <- NormalizeData(seurat_object, normalization.method = "RC", scale.factor = 1e6)
		counts <- GetAssayData(object = seurat_object, assay = "Spatial")
	} else {
		### MERFISH datasets
		seurat_object <- NormalizeData(seurat_object, normalization.method = "RC", scale.factor = 1e6)
		counts <- GetAssayData(object = seurat_object, assay = "SCT")
	}
	counts <- as.matrix(counts)
	gene_names <- rownames(counts)
	counts_t <- data.table(counts)
	counts_t$gene <- gene_names
	gc()
	counts_t <- transpose(counts_t, keep.names = "cell", make.names = "gene") ### Takes long
	gc()
	cell_names <- counts_t$cell
	counts_t$cell <- NULL ### Remove the gene name column
	counts_t <- as.matrix(counts_t)
	rownames(counts_t) <- cell_names
	saveRDS(counts_t, paste0(sample_dir, "/processed/bk.dat.RC.RDS"))
	
	cat(paste(c("### [", sampleid, "] Finished.\n"), collapse = ""))
}
write.table(meta, file = paste0(dsdir, "/metatable_auto.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
if (file.exists(paste0(dsdir, "/metatable_orig.tsv"))) {
	meta$Condition <- NULL
	meta_orig <- fread(paste0(dsdir, "/metatable_orig.tsv"))
	meta <- merge(meta, meta_orig, by = c("DSID", "SampleID"))
}
write.table(meta, file = paste0(dsdir, "/metatable.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

cat(paste(c(">>>>>>>> Dataset: ", dsid, " completed <<<<<<<<\n"), collapse = ""))
