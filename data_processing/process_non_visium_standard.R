#### Description
#### * The purpose of this script is to process non-Visium data stored in standard structure (counts.csv + coordinates.csv). This also applies to Visium data with no h5 + spatial data provided for public download.
#### * This script will loop through the datasets, read, transform, and process the datasets, and save the processed data as well as a metatable for each dataset.
#### * Note that if after QC, there are < 10 spots left, we will exclude the sample from our database.
#### 
#### Author: Yiming Li, Saya Dennis (edits)

library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
# library(SeuratDisk)

dn <- "/projects/b1131/SpatialT/"
args <- commandArgs(trailingOnly=TRUE)
dsid <- args[1]
tech <- args[2]
dn <- paste0(dn, tech)

#### Get PID
dsid <- paste0("DS", dsid)
pid <- paste0("PID", substr(dsid, start=3, stop=nchar(dsid)-1)) # e.g. "PID203"
#### Create empty data frame for meta table 
meta <- data.frame(DSID = character(0), SampleID = character(0), Nspots = integer(0), Nspots_postQC = integer(0), Ngenes = integer(0), Ngenes_postQC = integer(0), Condition = character(0))
#### Set dataset directory and get a list of samples 
dsdir <- paste0(dn, "/", pid, "/", dsid) # e.g. "/projects/b1131/SpatialT/Slide-seq/PID203/DS203A"
sampleids <- dir(dsdir, pattern = dsid) # e.g. list of elements like "DS203A.1"

cat(paste(c("\n>>>>>>>> Dataset [", dsdir, "] started <<<<<<<<\n"), collapse = ""))
for (sampleid in sampleids) {
	#### Read data
	## !!! Assumes that the prepared counts.csv and coordinates.csv are ready under original/
	counts <- fread(paste0(dsdir, "/", sampleid, "/original/counts.csv"), sep = ",", header = TRUE)
	coordinates <- read.table(paste0(dsdir, "/", sampleid, "/original/coordinates.csv"), sep = ",", header = TRUE)
	counts <- counts[counts$gene != "",] ## Remove empty "genes"
	## Some datasets have duplicated values in the "gene" or "barcode" column
	## Remove these "genes" / "barcodes" since it is hard to determine which observation we should keep
	tmp <- table(counts$gene)
	tmp <- names(tmp[tmp > 1])
	if (length(tmp) > 0) {
		counts <- counts[!counts$gene %in% tmp,]
		cat(paste(c("# WARNING: Excluded ", as.character(length(tmp)), " genes with multiple associated rows.\n"), collapse = ""))
	}
	coordinates <- coordinates[!duplicated(coordinates),]
	tmp <- table(coordinates$barcode)
	tmp <- names(tmp[tmp > 1])
	if (length(tmp) > 0) {
		coordinates <- coordinates[!coordinates$barcode %in% tmp,]
		cat(paste(c("# WARNING: Excluded ", as.character(length(tmp)), " spot IDs with multiple associated rows.\n"), collapse = ""))
	}
	
	#### Generate Seurat object
	gene_names <- counts$gene
	counts$gene <- NULL
	counts <- as.matrix(counts)
	rownames(counts) <- gene_names
	seurat_object <- CreateSeuratObject(counts = counts, project = 'SlideSeq', assay = "Spatial")
	rownames(coordinates) <- coordinates$barcode
	## Using the spot coordinates, instead of pixel coordinates
	## !!! This Seurat object cannot be directly overlayed on top of the tissue image (if any)
	coordinates <- coordinates[,c("row", "col")]
	# if (sum(is.na(coordinates$imagerow)) + sum(is.na(coordinates$imagecol)) == 0) {
	# 	coordinates <- coordinates[,c("imagerow", "imagecol")]
	# } else {
	# 	coordinates <- coordinates[,c("row", "col")]
	# }
	colnames(coordinates) <- c("xcoord", "ycoord")
	seurat_object[['images']]<- new(Class = "SlideSeq", assay = "Spatial", coordinates = coordinates)

	nspots <- ncol(seurat_object)
	ngenes <- nrow(seurat_object)
	
	#### Spot QC
	
	## Step 1. Remove the spots with total UMI count < 500 / the total number of genes < 500 / >= 25% mitochondrial reads
	seurat_object[["percent_mt"]] <- PercentageFeatureSet(seurat_object, "^MT-")
	qc_step1 <- seurat_object$nCount_Spatial >= 500 & seurat_object$nFeature_Spatial >= 500 & seurat_object$percent_mt < 25
	if (sum(qc_step1) < 10) {
		cat(paste(c("# NOTE: Sample has less than 10 spots after QC. This sample will be excluded from the database and metatable.tsv.\n"), collapse = ""))
		next
	}
	seurat_object <- seurat_object[, qc_step1]

	## Step 2.
	## Remove the spots with total UMI count < median(total UMI count) - 3 * SD(total UMI count).
	## Remove the spots with total number of genes < median(total number of genes) - 3 * SD(total number of genes)
	qc_step2 <- seurat_object$nCount_Spatial >= median(seurat_object$nCount_Spatial) - 3 * sqrt(var(seurat_object$nCount_Spatial)) & seurat_object$nFeature_Spatial >= median(seurat_object$nFeature_Spatial) - 3 * sqrt(var(seurat_object$nFeature_Spatial))
	if (sum(qc_step2) < 10) {
		cat(paste(c("# NOTE: Sample has less than 10 spots after QC. This sample will be excluded from the database and metatable.tsv.\n"), collapse = ""))
		next
	}
	seurat_object <- seurat_object[, qc_step2]
	
	#### Gene QC
	counts <- data.frame(GetAssayData(object = seurat_object, assay = "Spatial", slot = "counts"))
	counts <- counts > 0
	n_spots_per_gene <- rowSums(counts)
	seurat_object <- seurat_object[n_spots_per_gene >= 5,]
	
	nspots_qc <- ncol(seurat_object)
	ngenes_qc <- nrow(seurat_object)
	
	#### Append to the dataset metatable
	meta <- rbind(meta, data.frame(DSID = dsid, SampleID = sampleid, Nspots = nspots, Nspots_postQC = nspots_qc, Ngenes = ngenes, Ngenes_postQC = ngenes_qc, Condition = NA))
	
	#### Transform and process seurat_object
	seurat_object <- SCTransform(seurat_object, assay = "Spatial", verbose = FALSE)
	
	#### Dimensionality reduction and clustering
	
	#### The number of spots in a ST dataset is often small, need to set the npcs and dims parameters
	n_pcs <- min(min(dim(seurat_object)) - 1, 50)
	seurat_object <- RunPCA(seurat_object, assay = "SCT", verbose = FALSE, npcs = n_pcs)
	n_dims <- min(min(dim(seurat_object)) - 1, 30)
	seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:n_dims, verbose = FALSE)
	
	## If the number of spots is < 50, UMAP with uwot (default) will fail, and FindClusters with resolution  = 1.2 will sometimes fail
	## https://github.com/satijalab/seurat/issues/4312#issuecomment-812938288
	if (ncol(seurat_object) < 50) {
		seurat_object <- FindClusters(seurat_object, resolution = 1, verbose = FALSE)
		## resolution 1.2 sometimes fail
		seurat_object <- RunUMAP(seurat_object, umap.method = "umap-learn", reduction = "pca", dims = 1:n_dims, verbose = FALSE)
		cat(paste(c("# WARNING: Sample has less than 50 spots after QC.\n"), collapse = ""))
	} else {
		seurat_object <- FindClusters(seurat_object, resolution = 1.2, verbose = FALSE)
		## resolution: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html advises 0.4-1.2 for around 3K cells
		seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:n_dims, verbose = FALSE)
	}
	
	### Save processed Seurat object
	## "processed" folder should already exist for all datasets but just in case
	if (!dir.exists(paste0(dsdir, "/", sampleid, "/processed/"))) {
		dir.create(paste0(dsdir, "/", sampleid, "/processed/"))
	}
	saveRDS(seurat_object, file = paste0(dsdir, "/", sampleid, "/processed/Seurat.RDS"))
	
	# #### Retrieve QC-ed + transformed counts and coordinates
	# ## Counts 
	# counts <- GetAssayData(object = seurat_object, assay = "SCT", slot = "counts")
	# counts <- tibble::rownames_to_column(data.frame(counts), "gene")
	# ## Coordinates
	# coords <- seurat_object@images$image@coordinates
	# coords <- tibble::rownames_to_column(data.frame(coords), "barcode")
	# coords <- coords[,c("barcode", "x", "y")]
	
	# #### Write counts, coordinates, and meta_spots to file
	# write.table(counts, file = paste0(dsdir, "/", sampleid, "/processed/counts.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
	# write.table(coords, file = paste0(dsdir, "/", sampleid, "/processed/coordinates.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
	# saveRDS(list("counts" = counts, "coordinates" = coords), file = paste0(dsdir, "/", sampleid, "/processed/data_frames.RDS"))
	# cat(paste(c("### [", sampleid, "] Finished.\n"), collapse = ""))
	# gc()

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
