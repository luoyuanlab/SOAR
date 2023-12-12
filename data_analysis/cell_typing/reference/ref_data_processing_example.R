### This is a pseudo-script demonstrating the possible steps of processing the downloaded scRNA-seq datasets
### All the filenames are hard-wired, and this script is for your reference only
### 
### Author: Yiming Li

### This example uses a brain scRNA-seq dataset:
### cd /share/fsmresfiles/SpatialT/ref/Brain/Non_Adult/GSE60361



###### Read the count matrix and metadata

library(data.table) ### For fread -- faster than read.table

exprMatrix <- fread("exprMatrix.tsv", sep = "\t") ### row = gene, column = cell
str(exprMatrix)
# Classes ‘data.table’ and 'data.frame':	19972 obs. of  3006 variables:
#  $ sample        : chr  "Tspan12" "Tshz1" "Fnbp1l" "Adamts15" ...
#  $ 1772071015-C02: num  0 2 2 0 1 ...
#  $ 1772071017-G12: num  0 1 1 0 1 0 0 0 0 0 ...
#  $ 1772071017-A05: num  0 0 2.81 0 1 ...
# ......

meta <- fread("meta.tsv", sep = "\t") ### row = cell
str(meta)
# Classes ‘data.table’ and 'data.frame':	3005 obs. of  12 variables:
#  $ V1             : chr  "1772071015-C02" "1772071017-G12" "1772071017-A05" "1772071014-B06" ...
#  $ tissue         : chr  "sscortex" "sscortex" "sscortex" "sscortex" ...
#  $ group          : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ total mRNA mol : int  21580 21748 31642 32916 21531 24799 31406 20389 23022 24184 ...
#  $ well           : int  11 95 33 42 48 13 50 66 29 28 ...
#  $ sex            : int  1 -1 -1 1 1 -1 1 -1 1 1 ...
#  $ age            : int  21 20 20 21 25 20 25 23 21 21 ...
#  $ diameter       : num  0 9.56 11.1 11.7 11 11.9 11.3 10.9 12.9 11.2 ...
#  $ level1class    : chr  "interneurons" "interneurons" "interneurons" "interneurons" ...
#  $ level2class    : chr  "Int10" "Int10" "Int6" "Int10" ...
# ......

### Here we see that the columns "level1class" and "level2class" are cell type labels

### We can check if the number of cells in the count matrix and the metatable are the same
### (Assuming that the rows/columns in the metatable / count matrix are unique)
if (ncol(exprMatrix) == nrow(meta) + 1) {
	### "+ 1" because the first column of exprMatrix is not a cell
	cat("\nDimensions match\n")
	if (sum(sort(colnames(exprMatrix)[2:ncol(exprMatrix)]) == sort(meta$V1)) == nrow(meta)) {
		cat("Cell IDs match\n")
	} else {
		cat("Cell IDs do not match\n")
	}
} else {
	cat("\nDimensions do not match\n")
}

###!!!!!!! Note that you will need to filter exprMatrix and meta if the dimensions and/or cell IDs do not match!

### To perform harmonization (later with another dataset), we can check the unique cell labels in the metatable

table(meta$level1class)
# astrocytes-ependymal    endothelial-mural         interneurons
#                  224                  235                  290
#            microglia     oligodendrocytes        pyramidal CA1
#                   98                  820                  939
#         pyramidal SS
#                  399

table(meta$level2class)
#  (none)    Astro1    Astro2   CA1Pyr1   CA1Pyr2 CA1PyrInt   CA2Pyr2   Choroid
#     189        68        61       380       447        49        41        10
# ClauPyr     Epend      Int1     Int10     Int11     Int12     Int13     Int14
#       5        20        12        21        10        21        15        22
#   Int15     Int16      Int2      Int3      Int4      Int5      Int6      Int7
#      18        20        24        10        15        20        22        23
#    Int8      Int9      Mgl1      Mgl2    Oligo1    Oligo2    Oligo3    Oligo4
#      26        11        17        16        45        98        87       106
#  Oligo5    Oligo6     Peric      Pvm1      Pvm2   S1PyrDL  S1PyrL23   S1PyrL4
#     125       359        21        32        33        81        74        26
# S1PyrL5  S1PyrL5a   S1PyrL6  S1PyrL6b    SubPyr     Vend1     Vend2      Vsmc
#      16        28        39        21        22        32       105        62

### We can also write the above tables to file if you would like to work in a spreadsheet later
write.table(table(meta$level1class), "GSE60361_level1class.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(table(meta$level2class), "GSE60361_level2class.txt", quote = FALSE, sep = "\t", row.names = FALSE)

###!!!!!!! Cell type label harmonization (no generic codes for this)

### You can define a function for replacing labels during harmonization
replace_labels <- function(vector, old_label, new_label) {
	return(replace(vector, vector == old_label, new_label))
}

### For example, if you would like to change the "oligodendrocytes" labels into "OLG" in this dataset
meta$level1class <- replace_labels(meta$level1class, "oligodendrocytes", "OLG")

table(meta$level1class) ### See that the labels have changed
# astrocytes-ependymal    endothelial-mural         interneurons
#                  224                  235                  290
#            microglia                  OLG        pyramidal CA1
#                   98                  820                  939
#         pyramidal SS
#                  399






###### Save the count matrix and metadata

library(data.table)
library(mltools) ### For the sparsify() function
library(plyr) ### For the join() function

exprMatrix <- fread("exprMatrix.tsv", sep = "\t") ### Read into a data.table
gene_names <- exprMatrix$sample ### Store the gene names
exprMatrix$sample <- NULL ### Remove the gene column

exprMatrix <- sparsify(exprMatrix, sparsifyNAs = TRUE) ### Convert to dgCMatrix format
rownames(exprMatrix) <- gene_names ### The column names should be there, you can check by colnames(exprMatrix)
saveRDS(exprMatrix, "mtx.rds")

### The below assumes that we will use the meta$level2class labels, and they have been harmonized with other datasets

### Change the order of cells in the metatable
meta <- fread("meta.tsv", sep = "\t")
meta <- meta[,c("V1", "level2class")]
colnames(meta) <- c("cell", "label")
meta_sorted <- data.frame(cell = colnames(exprMatrix))
meta_sorted <- join(meta_sorted, meta, by = "cell")

fwrite(meta_sorted, file = "ident.csv")






###### Create a SingleCellExperiment object based on the count matrix and the (harmonized) metatable

library(data.table)
library(SingleCellExperiment)
library(scuttle)

### Read the prepared files
exprMatrix <- readRDS("mtx.rds") ### From the previous step
meta_sorted <- fread("ident.csv", sep = ",") ### From the previous step

meta_sorted$cell <- NULL

### Create SingleCellExperiment object for SingleR and normalize it
ref_data_sce <- SingleCellExperiment(list(counts = exprMatrix), colData = meta_sorted)

ref_data_sce <- logNormCounts(ref_data_sce)
### If your downloaded counts data is already normalized, the above command will fail.
### However, SingleR expects a "logcounts" assay in the input SCE object, so you need to run the following command.
# logcounts(ref_data_sce) <- counts(ref_data_sce)

ref_data_sce <- ref_data_sce[,ref_data_sce$label != ""] ### Remove empty cell type labels
### Need to remove the cells with labels like “doublets”, “Not Assigned”, etc., e.g.:
# ref_data_sce <- ref_data_sce[,ref_data_sce$label != "not applicable"]

saveRDS(ref_data_sce, file = "GSE60361.RDS")






######### Test cell type annotation on a dataset (3A)

library(Seurat)
library(SingleR)

ref_data_sce <- readRDS("GSE60361.RDS")
dim(ref_data_sce)
table(ref_data_sce$label)
length(table(ref_data_sce$label))

### Example: human brain Visium dataset
target_p_name <- "PID3"
target_ds_name <- "DS3A"
target_ds_dir <- paste(c("/share/fsmresfiles/SpatialT/10x", target_p_name, target_ds_name), collapse = "/")
target_ds_metatable <- read.table(paste0(target_ds_dir, "/metatable.tsv"), header = TRUE, stringsAsFactors = FALSE)
target_sample_name <- target_ds_metatable$SampleID[1] # For testing only
target_sample_dir <- paste(c(target_ds_dir, "/", target_sample_name), collapse = "")
seurat_object_tn_path <- paste0(target_sample_dir, "/processed/Seurat.RDS")
seurat_object_tn <- readRDS(seurat_object_tn_path)

### Perform spot-based cell type annotation and save to Seurat
annotation <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce, labels = ref_data_sce$label, de.method="wilcox")
seurat_object_tn[["cell_type_annotation"]] <- annotation$labels

### Perform cluster-based cell type annotation and save to Seurat
cluster_results <- seurat_object_tn[["seurat_clusters"]]$seurat_clusters
annotation <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce, clusters = cluster_results, labels = ref_data_sce$label, de.method="wilcox")
seurat_object_tn[["cell_type_annotation_clusters"]] <- annotation$labels[cluster_results]

### Visualize annotated cell types
pdf("~/stbase/DS3A_test.pdf") ### Change to your own save directory/name
print(SpatialDimPlot(seurat_object_tn))
print(DimPlot(seurat_object_tn, reduction = "umap"))
print(SpatialDimPlot(seurat_object_tn, group.by = "cell_type_annotation"))
print(DimPlot(seurat_object_tn, reduction = "umap", group.by = "cell_type_annotation"))
print(SpatialDimPlot(seurat_object_tn, group.by = "cell_type_annotation_clusters"))
print(DimPlot(seurat_object_tn, reduction = "umap", group.by = "cell_type_annotation_clusters"))
dev.off()
