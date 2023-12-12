import os, sys, pickle, datetime
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
from collections import Counter

### Read in data
data_dir = sys.argv[1] + "/analysis/deconvolution/"
# data_dir = '/share/fsmresfiles/SpatialT/10x/PID4/DS4A/DS4A.1/analysis/deconvolution/'
# data_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/deconvolution/'
counts_dir = data_dir + 'binded_counts.csv'
coord_dir = data_dir + 'binded_coordinates.csv'
anno_dir = data_dir + 'binded_cell_types.tsv'
counts = pd.read_csv(counts_dir, index_col = 0)
coord = pd.read_csv(coord_dir, index_col = 2)
annotation = pd.read_csv(anno_dir, index_col = 0, sep = "\t")
cell_types = annotation['cell_type'].tolist()

### Get species
ref_dir = sys.argv[2]
# ref_dir = "/projects/b1131/SpatialT/ref/final/Cancer/Breast/Human"
ref_dir = ref_dir.split("/")
species = [i for i in ref_dir if i][-1].lower() # Last non-empty string in ref path

### Get spatial distance type
thr_type = sys.argv[3]
# thr_type = "short"
thr_type_multiplier = {
    "short": 500,
    "medium": 1000,
    "long": 1500,
    "xlong": 2500,
}

### Spatial distance constraint
center_to_center_dist_techs = {
    "10x": 100,
    "ST": 200,
    "DBiT-seq": 20,
    "Slide-seq": 20,
    "MERFISH": 0.334,
    "osmFISH": 0.13,
    "seqFISH": 0.26,
    "sci-Space": 222,
}
tech =  data_dir.split("/")[4]
dis_thr = thr_type_multiplier[thr_type] / center_to_center_dist_techs[tech]

### Set up anndata
adata = sc.AnnData(counts.T)
adata.var_names_make_unique()
adata = adata[coord.index,]
coor_df = coord.loc[adata.obs_names, ["x", "y"]]
adata.obsm["spatial"] = coor_df.to_numpy()
adata.raw = adata

### Data processing
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
adata_disthr = adata.copy()
sc.pp.highly_variable_genes(adata, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, svd_solver = 'arpack')
sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = 40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.4)

### Get CellChat ligand-receptors
df_cellchat = ct.pp.ligand_receptor_database(database = 'CellChat', species = species)
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_disthr, min_cell_pct = 0.05)

if (df_cellchat_filtered.shape == (0, 0)):
    raise ValueError("ct.pp.filter_lr_database() returns an empty data frame, too few overlapping genes between reference and data")

now = datetime.datetime.now()
print("Analysis started: ")
print(now)

ct.tl.spatial_communication(adata_disthr, database_name = 'cellchat', df_ligrec = df_cellchat_filtered, dis_thr = dis_thr, heteromeric = True, pathway_sum = True)

adata_disthr.write(data_dir + "adata_disthr_" + thr_type + ".h5ad")

now = datetime.datetime.now()
print("Analysis finished: ")
print(now)
