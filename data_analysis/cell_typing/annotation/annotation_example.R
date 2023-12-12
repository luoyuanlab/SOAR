"""
Cell Type Annotation Example Script

Author: Saya Dennis 

Usage: Rscript --no-save annotation_example.R

Requirements: 
1. You need to first process your reference scRNA-seq dataset.
2. In the code below, edit the directory/file names of your reference dataset (saved into variable ref_data_sce) 
3. You will need to create an annotation directory under your ST sample directories 
    - This should look like this: /share/fsmresfiles/SpatialT/{tech}/PID{pid}/{dsid}/{sampleid}/analysis/annotation/
    - To automate generating this directory, refer to script ST-dataset/analysis/cell_type_annotation/01_create_anno_directory.py
4. Back up your un-annotated Seurat object 
    - Back up to a file named Seurat.RDS.bk under the same directory. 
    - Refer to ST-dataset/analysis/cell_type_annotation/saya_cell_type_annotation_examples/02_backup_seurat_before_annotation.py
5. Below, edit the PID, DSID, and technology directory (e.g. /share/fsmresfiles/SpatialT/DBiT-seq)

"""

library(data.table)
library(SingleCellExperiment)
library(scuttle)
library(Seurat)
library(SingleR)

################################################
#### Cell type annotation on a dataset DS2O ####
################################################

dref <- '/share/fsmresfiles/SpatialT/ref/Heart/Adult/heart-cell-atlas/processed/' # reference directory 
ref_data_sce <- readRDS(paste0(dref, 'sce_heart.RDS'))

target_p_name <- "PID70"
target_ds_name <- "DS70B"
target_ds_dir <- paste(c("/share/fsmresfiles/SpatialT/DBiT-seq", target_p_name, target_ds_name), collapse = "/")
target_ds_metatable <- read.table(paste0(target_ds_dir, "/metatable.tsv"), header = TRUE, stringsAsFactors = FALSE)

### Loop through samples and annotate 
for (target_sample_name in target_ds_metatable$SampleID) {
    cat(paste0("Starting annotation for sample ", target_sample_name, " -- ", Sys.time(), "\n"))
    # create annotation directory 
    # if (!dir.exists(paste0(dds, "/", sampleid, "/analysis/annotation"))) {
    #     dir.create(paste0(dds, "/", sampleid, "/analysis/annotation"))
    # }

    # load processed Seurat object 
    target_sample_dir <- paste(c(target_ds_dir, "/", target_sample_name), collapse = "")
    seurat_object_tn_path <- paste0(target_sample_dir, "/processed/Seurat.RDS.bk")
    seurat_object_tn <- readRDS(seurat_object_tn_path)

    ### Perform spot-based cell type annotation and save to Seurat
    annotation <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce, labels = ref_data_sce$label, de.method="wilcox")
    seurat_object_tn[["cell_type_annotation"]] <- annotation$labels

    ### Perform cluster-based cell type annotation and save to Seurat
    cluster_results <- seurat_object_tn[["seurat_clusters"]]$seurat_clusters
    annotation <- SingleR(test = GetAssayData(seurat_object_tn, assay = "SCT"), ref = ref_data_sce, clusters = cluster_results, labels = ref_data_sce$label, de.method="wilcox")
    seurat_object_tn[["cell_type_annotation_clusters"]] <- annotation$labels[cluster_results]

    ### Overwrite the previously saved Seurat object with cell type annotated Seurat object 
    saveRDS(seurat_object_tn, file = paste0(target_sample_dir, "/processed/Seurat.RDS"))

    ### Visualize annotated cell types
    pdf(paste0(target_sample_dir, "/analysis/annotation/cell_type_annotation.pdf")) ### Change to your own save directory/name
    # pdf("analysis/annotation/cell_type_annotation.pdf")
    print(SpatialDimPlot(seurat_object_tn))
    print(DimPlot(seurat_object_tn, reduction = "umap"))
    print(SpatialDimPlot(seurat_object_tn, group.by = "cell_type_annotation"))
    print(DimPlot(seurat_object_tn, reduction = "umap", group.by = "cell_type_annotation"))
    dev.off()
    cat(paste0("Finished annotation for sample ", target_sample_name, " -- ", Sys.time(), "\n\n"))
}
