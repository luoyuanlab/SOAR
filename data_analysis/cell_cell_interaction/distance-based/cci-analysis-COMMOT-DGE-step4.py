import os, sys, pickle, datetime, anndata
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from collections import Counter

thr_type = sys.argv[2]
out_dir = sys.argv[1] + "/analysis/Distance/COMMOT_dec/" + thr_type
# thr_type = "short"
# out_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT/short'

with open(out_dir + "/pathways_step3_success.txt") as f:
    pathways = [line.rstrip('\n') for line in f]

for pathway in pathways:
    # pathway = "MIF"
    pathway_dir = out_dir + "/" + pathway
    
    yhat_scaled = pd.read_csv(pathway_dir + "/yhatScaled.txt", sep = "\t", index_col = 0)
    df_assoRes = pd.read_csv(pathway_dir + "/assoRes.txt", sep = "\t", index_col = 0)
    
    df_deg = df_assoRes.rename(columns={'waldStat_1':'waldStat', 'df_1':'df', 'pvalue_1':'pvalue'})
    df_deg = df_deg[['waldStat', 'df', 'pvalue']]
    idx = np.argsort(-df_deg['waldStat'].values)
    df_deg = df_deg.iloc[idx]
    df_yhat = yhat_scaled
    
    deg_result = {"df_deg": df_deg, "df_yhat": df_yhat}
    with open(pathway_dir + '/DEG_pt.pkl', 'wb') as handle:
        pickle.dump(deg_result, handle, protocol = pickle.HIGHEST_PROTOCOL)
    
    df_deg_clus, df_yhat_clus = ct.tl.communication_deg_clustering(df_deg, df_yhat, deg_clustering_res=0.4)
    top_de_genes = ct.pl.plot_communication_dependent_genes(df_deg_clus, df_yhat_clus, top_ngene_per_cluster=5, filename = pathway_dir + '/DEG.pdf', font_scale=1.2, return_genes = True)
    
    deg_result = {"df_deg": df_deg, "df_yhat": df_yhat, "df_deg_clus": df_deg_clus, "df_yhat_clus": df_yhat_clus, "top_de_genes": top_de_genes}
    with open(pathway_dir + '/DEG_full.pkl', 'wb') as handle:
        pickle.dump(deg_result, handle, protocol = pickle.HIGHEST_PROTOCOL)
