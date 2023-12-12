### conda activate R4
library(data.table)
library(splitstackshape)
library(Seurat)
library(BayesPrism)
library(SingleCellExperiment)
library(scuttle)

### This script uses Thymus, Human as an example
### This script assumes that:
### * A Seurat object of the identified scRNA-seq reference data has already been created
### * The annotated cell types are stored in seurat.object[["label"]]
organ <- "Thymus"
species <- "Human"
seurat.object <- readRDS("/share/fsmresfiles/SpatialT/ref/Tabula_Sapiens/Thymus/Seurat.RDS")
qc_plot_output_dir <- "~/stbase" # Change this to your own directory

#### Create save path
save_path <- paste0("/share/fsmresfiles/SpatialT/ref/final/", organ)
if (!dir.exists(save_path)) {
	dir.create(save_path)
}
save_path <- paste0(save_path, "/", species)
if (!dir.exists(save_path)) {
	dir.create(save_path)
}

#### Get # cells and # genes before QC
n_cells <- ncol(seurat.object)
n_genes <- nrow(seurat.object)
gc()

#### Get metadata
meta <- seurat.object@meta.data
meta$cell_id <- rownames(meta)
meta <- meta[,c("cell_id", "label")]
meta <- meta[!is.na(meta$label),]

### Perform stratified sampling (by cell type labels) on the reference if the total number of cells is larger than 30000 to avoid out-of-memory error
### You may change 30000 to a slightly larger number or skip this step if QC drops a lot of cells
if (nrow(meta) > 30000) {
	if (organ == "Brain") {
		meta <- stratified(meta, "subclass_label", size = 30000 / n_cells)
	} else {
		meta <- stratified(meta, "label", size = 30000 / n_cells)
	}
}
seurat.object[["cell_id"]] <- colnames(seurat.object)
seurat.object <- subset(seurat.object, subset = cell_id %in% meta$cell_id)
n_cells_strat <- ncol(seurat.object)

### Get cell type proportions
cell_type_string <- table(meta$label)
percs <- paste0(as.character(round(as.numeric(cell_type_string / sum(cell_type_string) * 100), 2)), "%")
cell_type_string2 <- paste0(names(cell_type_string), " (", as.character(round(as.numeric(cell_type_string), 2)), ", ", percs, ")")
cell_type_string2 <- paste(cell_type_string2, collapse = ", ")

### Sort metadata by Seurat object's cell order
meta <- meta[match(colnames(seurat.object), meta$cell_id),]
sum(meta$cell_id == colnames(seurat.object)) == nrow(meta)

### Cell QC
### !!! Please perform cell QC case-by-case instead of using uniform thresholds
seurat.object[["percent_mt"]] <- PercentageFeatureSet(seurat.object, "^MT-")
seurat.object.bk <- seurat.object
pdf(paste0(qc_plot_output_dir, "/sc_ref_", organ, "_", species, "_beforeQC.pdf"))
print(VlnPlot(seurat.object, features = "nCount_RNA"))
print(VlnPlot(seurat.object, features = "nFeature_RNA"))
dev.off()

### * Change the nCount_RNA and nFeature_RNA thresholds based on the violin plots
seurat.object <- seurat.object.bk
seurat.object <- seurat.object[, seurat.object$nCount_RNA > 500 & seurat.object$nFeature_RNA > 250 & seurat.object$percent_mt < 20]
(n_cells_qc1 <- ncol(seurat.object))
pdf(paste0(qc_plot_output_dir, "/sc_ref_", organ, "_", species, "_QC1.pdf"))
print(VlnPlot(seurat.object, features = "nCount_RNA"))
print(VlnPlot(seurat.object, features = "nFeature_RNA"))
dev.off()

### * Change the two thresholds based on the violin plots
seurat.object2 <- seurat.object[, seurat.object$nCount_RNA < 30000 & seurat.object$nFeature_RNA < 6000]
n_cells_qc2 <- ncol(seurat.object2)
pdf(paste0(qc_plot_output_dir, "/sc_ref_", organ, "_", species, "_QC2.pdf"))
print(VlnPlot(seurat.object2, features = "nCount_RNA"))
print(VlnPlot(seurat.object2, features = "nFeature_RNA"))
dev.off()

### Overwrite the Seurat object if the second two violin plots look okay
# seurat.object
# seurat.object2
seurat.object <- seurat.object2
# seurat.object

### Write the QC-ed Seurat object and the metadata to file
meta <- seurat.object@meta.data
meta$cell_id <- rownames(meta)
meta <- meta[,c("cell_id", "label")]
fwrite(meta, paste0(save_path, "/cell_types.txt"), sep = "\t")
saveRDS(seurat.object, paste0(save_path, "/Seurat.RDS"))

### Generate SingleCellExperiment object for cell type annotation
counts <- GetAssayData(seurat.object, assay = "RNA")
meta$cell_id <- NULL
ref_data_sce <- SingleCellExperiment(list(counts = counts), colData = meta)
ref_data_sce <- logNormCounts(ref_data_sce)
saveRDS(ref_data_sce, file = paste0(save_path, "/SCE.RDS"))

### Generate data for BayesPrism (cell type deconvolution)
### * Counts data
counts <- GetAssayData(seurat.object, assay = "RNA")
counts <- as.matrix(counts)
gene_names <- rownames(counts)
counts <- data.table(counts)
counts$gene <- gene_names
counts <- transpose(counts, keep.names = "cell", make.names = "gene")
cell_names <- counts$cell
counts$cell <- NULL
counts <- as.matrix(counts)
rownames(counts) <- cell_names
min_count <- min(counts)
max_count <- max(counts)
saveRDS(counts, paste0(save_path, "/mat_transposed.RDS"))
gc()

### * BayesPrism gene filtering step 1
if (species == "Human") {
	sc.dat.filtered <- cleanup.genes(input=counts, input.type="count.matrix", species="hs", gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1","chrX","chrY"), exp.cells=5)
	sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding") # only works for human
} else {
	sc.dat.filtered.pc <- cleanup.genes(input=counts, input.type="count.matrix", species="mm", gene.group = c("Rb", "Mrp", "other_Rb", "chrM","chrX","chrY"), exp.cells=5)
}

### * BayesPrism gene filtering step 2
cell.type.labels <- meta$label
cell.state.labels <- meta$label
diff.exp.stat <- get.exp.stat(sc.dat=counts[,colSums(counts>0)>3], cell.type.labels=cell.type.labels, cell.state.labels=cell.state.labels, psuedo.count=0.1, cell.count.cutoff=50, n.cores=1)

### !!! Check that all cell types has > 50 marker genes
### This threshold can be more lenient for sparser cell types
### Change pval.max and lfc.min to get more genes
sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.pc, stat=diff.exp.stat, pval.max=0.01, lfc.min=0.1)

### You may check if ngenes_filt2 is around 5000
### If the number is too small (e.g. < 2000), we may need to use ngenes_filt1 for deconvolution instead
(ngenes_filt1 <- ncol(sc.dat.filtered.pc))
(ngenes_filt2 <- ncol(sc.dat.filtered.pc.sig))

### Save to BayesPrism-filtered references to file
saveRDS(sc.dat.filtered.pc, paste0(save_path, "/sc.dat.filtered.pc.RDS"))
saveRDS(sc.dat.filtered.pc.sig, paste0(save_path, "/sc.dat.filtered.pc.sig.RDS"))

results <- data.frame(organ = organ, species = species, save_path = save_path, ncells = n_cells, ngenes = n_genes, ncells_qc1 = n_cells_qc1, ncells_qc2 = n_cells_qc2, n_cells_strat = n_cells_strat, cell_types = cell_type_string2, min_count = min_count, max_count = max_count, ngenes_filt1 = ngenes_filt1, ngenes_filt2 = ngenes_filt2)
fwrite(results, paste0(save_path, "/summary_stats.txt"), sep = "\t")
gc()
results
cat(paste0("\n>>> ", save_path, " finished\n"))
