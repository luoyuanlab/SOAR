# Spatial transcriptOmics Analysis Resource

This repository contains the data curation, processing, and analysis scripts used by the article "SOAR elucidates disease mechanisms and empowers drug discovery through spatial transcriptomics" [[bioRxiv preprint]](https://www.biorxiv.org/content/10.1101/2022.04.17.488596v2) | [[Website]](https://soar.fsm.northwestern.edu/).

## Data curation

To query the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) for potential human and mouse spatial transcriptomics datasets, please run [the Python script](https://github.com/luoyuanlab/SOAR/tree/main/data_curation/geo-query.py) using different keywords.

* `python3 geo-query.py`

The retrieved GDS list with annotated meta-information will be stored in `./<%Y%m%d>/all.csv`.

## Data processing

The data processing scripts are available under [`data_processing/`](https://github.com/luoyuanlab/SOAR/tree/main/data_processing). The scripts automatically perform spot and gene quality control, data transformation, normalization, and dimensionality reduction.

**10x Visium data** in standard format can be processed using [`process_visium_standard.R`](https://github.com/luoyuanlab/SOAR/tree/main/data_processing/process_visium_standard.R). The script assumes that the directory contains one option from the below:

1. (a) `filtered_feature_bc_matrix.h5` or [MEX files](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) and (b) the image data in a subdirectory called `spatial`
2. A Visium Seurat object with `data@images` properly added

Please note that 10x Visium data with only counts and coordinates and no `spatial/` folder data should be processed using the non-Visium scripts.

**Other types of spatial transcriptomics data** transformed into a standard format (`counts.csv` and `coordinates.csv`, please see below for the guidelines) can be processed using [`process_non_visium_standard.R`](https://github.com/luoyuanlab/SOAR/tree/main/data_processing/process_non_visium_standard.R). The script can also be used on Visium data with no h5 + spatial data provided for public download.

`counts.csv`

* Comma-delimited
* Header: `gene,<spot ID 1>,...,<spot ID n>`
* Each row = one gene
* Gene symbols should be used (not Ensembl IDs, etc.)

`coordinates.csv`

* Comma-delimited
* Header: `barcode,row,col`
	* The example file has more columns but only these three columns are required
	* Use the spot coordinates (row, col) instead of pixel coordinates (imagerow, imagecol in the example file) if possible
* Each row = one spot

The barcode column of `coordinates.csv` should be exactly the same as the `counts.csv` header (after removing "gene"), i.e. the spot IDs should match.

## Data analysis

### Overview

The overall flow of data analysis is as below.

1. Perform [spatial clustering](#spatial-clustering)
2. Perform whole-tissue [spatial variability analysis](#spatial-variability-analysis)
3. Check if the spatial transcriptomics technology is at single-cell level
	* If so (e.g. MERFISH), perform [cell type annotation](#cell-type-annotation)
	* If not (e.g. 10x Visium), perform [cell type deconvolution](#cell-type-deconvolution)
4. Perform cell-type-specific [spatial variability analyis](#spatial-variability-analysis)
5. Perform [cell-cell interaction analysis](#cell-cell-interaction-analysis)

### Spatial clustering

The scripts for performing spatial clustering are stored in [`data_analysis/spatial_clustering/`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_clustering). Please refer to the [README](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_clustering/README.md) for the details.

### Cell typing

In the [`data_analysis/cell_typing/`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing) folder, scripts are available for performing [cellular deconvolution](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/deconvolution) (spatial transcriptomics technologies with multiple cells per capture location, e.g. 10x Visium) and [cell type annotation](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/annotation) (single-cell-level spatial transcriptomics technologies, e.g. MERFISH).

#### scRNA-seq reference identification and processing

To identify scRNA-seq references for cell typing, users may utilize the [GEO query helper script](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/reference/geo-download-scRNA-seq.py). [`ref_data_processing_example.R`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/reference/ref_data_processing_example.R) is an example script for processing the downloaded scRNA-seq data. Please note that cell quality control needs to be performed case-by-case, i.e. the thresholds should be chosen manually based on the QC plots.

#### Cell type annotation

[`annotation_example.R`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/annotation/annotation_example.R) is an example script for performing cell type annotation on single-cell-level spatial transcriptomics datasets (e.g. MERFISH) using scRNA-seq reference datasets and SingleR.

<details><summary>Heuristics-guided cell type annotation for brain datasets (click me)</summary>

[`runBrainCellTypeAnnotation-CluHeu.R`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/annotation/runBrainCellTypeAnnotation-CluHeu.R)

* Usage: `./runBrainCellTypeAnnotation-CluHeu.R > runBrainCellTypeAnnotation-CluHeu.log`
* Description
    * This script automatically annotates the cell types of brain Visium datasets using a cluster-based approach guided by some heuristics.
    * Note that:
        * This script requires processed mouse and human scRNA-seq references as the input, and the file paths are currently hard-wired:
            * `/share/fsmresfiles/SpatialT/ref/Brain/Adult/aibs_human_ctx_smart-seq`
                * `aibs_human_ctx_smart-seq_neuronal.RDS`
                * `aibs_human_ctx_smart-seq_non_neuronal.RDS`
                * `supp.RData`
            * `/share/fsmresfiles/SpatialT/ref/Brain/Adult/aibs_mouse_ctx-hpf_10x`
                * `aibs_mouse_ctx-hpf_10x_neuronal.RDS`
                * `aibs_mouse_ctx-hpf_10x_non_neuronal.RDS`
                * `supp.RData`
        * This script also reads a table listing the DSID, species, and technology (`brain_DSID_list.txt`) and loops over its rows. Line 52 uses a hard-wired path to this file.
        * The annotations follow the [Common Cell Type Nomenclature (CCN)](https://portal.brain-map.org/explore/classes/nomenclature). `seurat_object[["cell_type_annotation"]]` contains the annotated subclasses, and `seurat_object[["cell_type_annotation_class"]]` contains the annotated classes (i.e., glutamatergic, GABAergic, or non_neuronal).
</details>

#### Cell type deconvolution

[Analysis scripts](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_typing/deconvolution) are available for the cell type deconvolution of spatial transcriptomics datasets using scRNA-seq reference datasets and BayesPrism.

Steps for running deconvolution

1. Sample script for processing scRNA-seq reference: `process_reference_example.R`
2. Deconvolution script: `quest_deconvolution_jobarray.R`
3. Script for preparing deconvolution results for subsequent analysis (only run this after the deconvolution is complete): `create_input_files.R`

### Spatial variability analysis

The scripts for the spatial variability (SV) analysis of spatial transcriptomics data are stored in  [`data_analysis/spatial_variability/`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_variability).

To perform whole-tissue SV analysis, use the script [`quest_SpatialDE_jobarray.py`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_variability/quest_SpatialDE_jobarray.py). To run the analysis, use the command `python quest_SpatialDE_jobarray.py $sample_directory`.

To perform cell-type-specific SV analysis, use the script [`quest_SpatialDE_ct_specific.py`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_variability/quest_SpatialDE_ct_specific.py). To run the analysis, use the command `python quest_SpatialDE_ct_specific.py $sample_directory $cell_type`.

### Cell-cell interaction analysis

To perform neighborhood-based cell-cell interaction analysis, use the script [`adj-analysis.R`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_cell_interaction/neighborhood-based/adj-analysis.R). To run the analysis, use the command `./adj-analysis.R $sample_directory`.

To perform distance-based cell-cell interaction analysis, run the bash script [`cci-analysis-COMMOT-DGE.sh`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/cell_cell_interaction/distance-based/cci-analysis-COMMOT-DGE.sh) to call different analysis scripts in the pipeline.

### Drug screen

Scripts for drug discovery analysis are stored under [`data analysis/drug_discovery`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/drug_discovery).

The four types of analysis are:

1. Differential gene expression analysis. [`Scripts for deconvoluted and annotated samples`](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/drug_discovery/DGE).
   
2. Protein-protein interaction (PPI) network for spatially variable, differentially expressed (DE-SV) genes by cell type. [`Script for generating PPI network`](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/ppi_quest.py).
   
3. CMAP L1000 drug enrichment (compounds with top overall positive/negative enrichment score on SV-DE gene sets of a cell type). [`Script for CMAP drug enrichment analysis`](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/drug_screen_perturb_quest.py).
   
4. CMAP L1000 drug perturbation (top gene targets perturbed by the top postiively/negatively enriched compounds). [`Script for CMAP drug perturbation analysis`](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/drug_screen_perturb_quest.py) (contained in the same script as above).
   
