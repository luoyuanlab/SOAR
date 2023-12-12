import os, sys, pickle, datetime, anndata, shutil
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from collections import Counter

### Read in data
data_dir = sys.argv[1] + "/analysis/deconvolution/"
# data_dir = '/share/fsmresfiles/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/deconvolution/'
# data_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/deconvolution/'
counts_dir = data_dir + 'binded_counts.csv'
counts = pd.read_csv(counts_dir, index_col = 0)

out_dir = sys.argv[1] + "/analysis/Distance/COMMOT_dec"
# out_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT_dec'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

### Get spatial distance type
thr_type = sys.argv[3]
out_dir = out_dir + "/" + thr_type
# out_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT_dec/short'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
# else:
#     shutil.rmtree(out_dir)
#     os.makedirs(out_dir)

adata_disthr = sc.read_h5ad(data_dir + "adata_disthr_" + thr_type + ".h5ad")
adata_disthr.layers['counts'] = scipy.sparse.csr_matrix(counts.values.T)

pathways = [str(i).replace("commot-cellchat-", "") for i in adata_disthr.obsp]
pathways = [i for i in pathways if "-" not in i]
pathways.sort()
with open(out_dir + '/pathways.txt', 'w') as f:
    for pathway in pathways:
        _ = f.write(f"{pathway}\n")

for pathway in pathways:
    # pathway = "MIF"
    
    ### rpy2 does not fully work on Quest; want to achieve the following:
    ### df_deg, df_yhat = ct.tl.communication_deg_detection(adata_disthr, database_name = 'cellchat', pathway = pathway, summary = 'receiver')
    summary = 'receiver'
    database_name = "cellchat"
    pathway_dir = out_dir + "/" + pathway
    # pathway_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT/short/MIF'
    if not os.path.exists(pathway_dir):
        os.makedirs(pathway_dir)
    
    # prepare input adata for R
    adata_deg = anndata.AnnData(X = adata_disthr.layers['counts'], var = pd.DataFrame(index=list(adata_disthr.var_names)), obs = pd.DataFrame(index=list(adata_disthr.obs_names)))
    adata_deg_var = adata_deg.copy()
    sc.pp.filter_genes(adata_deg_var, min_cells=3)
    sc.pp.filter_genes(adata_deg, min_cells=3)
    sc.pp.normalize_total(adata_deg_var, target_sum=1e4)
    sc.pp.log1p(adata_deg_var)
    sc.pp.highly_variable_genes(adata_deg_var, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_deg = adata_deg[:, adata_deg_var.var.highly_variable]
    del adata_deg_var
    
    summary_name = 'commot-'+database_name+'-sum-'+summary
    if summary == 'sender':
        summary_abrv = 's'
    else:
        summary_abrv = 'r'
    
    comm_sum = adata_disthr.obsm[summary_name][summary_abrv+'-'+pathway].values.reshape(-1,1)
    cell_weight = np.ones_like(comm_sum).reshape(-1,1)
    
    ### Save data for R
    Xmat = pd.DataFrame(adata_deg.X.toarray())
    Xmat.index = adata_deg.obs.index
    Xmat.columns = adata_deg.var.index
    Xmat.to_csv(pathway_dir + "/step1_X.csv", header = True, index = True)
    pseudoTime = pd.DataFrame(comm_sum)
    pseudoTime.to_csv(pathway_dir + "/step1_pseudoTime.csv", header = False, index = False)
    cellWeight = pd.DataFrame(cell_weight)
    cellWeight.to_csv(pathway_dir + "/step1_cellWeight.csv", header = False, index = False)
    
    nknots = 6
    
    string_fitGAM = 'sce <- fitGAM(counts=X, pseudotime=pseudoTime, cellWeights=cellWeight, nknots=%d, verbose=TRUE)' % nknots
    string_fitGAM = string_fitGAM + '\nassoRes <- data.frame( associationTest(sce, global=FALSE, lineage=TRUE) )'
    string_fitGAM = string_fitGAM + '\nassoRes[is.nan(assoRes[,"waldStat_1"]),"waldStat_1"] <- 0.0'
    string_fitGAM = string_fitGAM + '\nassoRes[is.nan(assoRes[,"df_1"]),"df_1"] <- 0.0'
    string_fitGAM = string_fitGAM + '\nassoRes[is.nan(assoRes[,"pvalue_1"]),"pvalue_1"] <- 1.0\n'
    
    with open(pathway_dir + "/step1.R", "w") as text_file:
        _ = text_file.write(string_fitGAM)
