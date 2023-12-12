import os, sys, pickle, datetime, anndata
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from collections import Counter

### Read in data
data_dir = sys.argv[1] + "/analysis/deconvolution/"
# data_dir = '/share/fsmresfiles/SpatialT/10x/PID4/DS4A/DS4A.1/analysis/deconvolution/'
# data_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/deconvolution/'
counts_dir = data_dir + 'binded_counts.csv'
counts = pd.read_csv(counts_dir, index_col = 0)

out_dir = sys.argv[1] + "/analysis/Distance/COMMOT_dec"
# out_dir = '/share/fsmresfiles/SpatialT/10x/PID4/DS4A/DS4A.1/analysis/Distance/COMMOT_dec'
# out_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT_dec'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

### Get spatial distance type
thr_type = sys.argv[3]
out_dir = out_dir + "/" + thr_type
# thr_type = "short"
# out_dir = '/share/fsmresfiles/SpatialT/10x/PID4/DS4A/DS4A.1/analysis/Distance/COMMOT_dec/thr_type'
# out_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT_dec/short'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

### Get COMMOT object and pathways
adata_disthr = sc.read_h5ad(data_dir + "adata_disthr_" + thr_type + ".h5ad")
with open(out_dir + "/pathways.txt") as f:
    pathways = [line.rstrip('\n') for line in f]

### Get cell types
anno_dir = data_dir + 'binded_cell_types.tsv'
annotation = pd.read_csv(anno_dir, index_col = 0, sep = "\t")
cell_types = annotation['cell_type'].tolist()

### Cell-type-level scores
adata_disthr.obs['cell_type'] = cell_types
for pathway in pathways:
    ct.tl.cluster_communication(adata_disthr, database_name = 'cellchat', pathway_name = pathway, clustering = 'cell_type', n_permutations = 100)

### Pull LR pairs
lrpairs = [str(i).replace("commot-cellchat-", "") for i in adata_disthr.obsp]
lrpairs = [i for i in lrpairs if "-" in i]
lrpairs.sort()
lrpairs.remove("total-total")
for lrpair in lrpairs:
    ct.tl.cluster_communication(adata_disthr, database_name = 'cellchat', pathway_name = lrpair, clustering = 'cell_type', n_permutations = 100)

### Pathway-level
rows_list = []
for pathway in pathways:
    # https://github.com/zcang/COMMOT/issues/10
    # rows (first index) represent senders and the columns (second index) represent receivers
    tmp_mat = adata_disthr.uns['commot_cluster-cell_type-cellchat-' + pathway]["communication_matrix"]
    tmp_p = adata_disthr.uns['commot_cluster-cell_type-cellchat-' + pathway]["communication_pvalue"]
    for i in tmp_mat.index:
        for j in tmp_mat.columns:
            rows_list.append({"pathway": pathway, "cell_type1": i, "cell_type2": j, "score": tmp_mat.loc[i,j], "p_val": tmp_p.loc[i,j]})

df = pd.DataFrame(rows_list)
df.to_csv(out_dir + "/communication_scores_pathway.txt", index = False)

### LR-pair-level
ccdb = adata_disthr.uns['commot-cellchat-info']["df_ligrec"]
ccdb["lrpair"] = ccdb["ligand"] + "-" + ccdb["receptor"]
ccdb_dict = dict(zip(ccdb.lrpair, ccdb.pathway))
rows_list = []
lrpairs2 = lrpairs
for lrpair in lrpairs:
    if lrpair in ccdb_dict:
        tmp_mat = adata_disthr.uns['commot_cluster-cell_type-cellchat-' + lrpair]["communication_matrix"]
        tmp_p = adata_disthr.uns['commot_cluster-cell_type-cellchat-' + lrpair]["communication_pvalue"]
        for i in tmp_mat.index:
            for j in tmp_mat.columns:
                rows_list.append({"pathway": ccdb_dict[lrpair], "lrpair": lrpair, "cell_type1": i, "cell_type2": j, "score": tmp_mat.loc[i,j], "p_val": tmp_p.loc[i,j]})
    else:
        lrpairs2.remove(lrpair)

with open(out_dir + '/lrpairs.txt', 'w') as f:
    for lrpair in lrpairs2:
        _ = f.write(f"{lrpair}\n")

df = pd.DataFrame(rows_list)
df.to_csv(out_dir + "/communication_scores_lrpair.txt", index = False)

### Overwrite result file
adata_disthr.write(data_dir + "adata_disthr_" + thr_type + "_cs.h5ad")
