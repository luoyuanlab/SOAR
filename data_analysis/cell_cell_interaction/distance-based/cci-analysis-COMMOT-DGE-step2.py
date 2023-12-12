import os, sys, pickle, datetime, anndata
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from collections import Counter

thr_type = sys.argv[2]
out_dir = sys.argv[1] + "/analysis/Distance/COMMOT_dec/" + thr_type
# out_dir = '/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1/analysis/Distance/COMMOT/short'

with open(out_dir + "/pathways_step1_success.txt") as f:
    pathways = [line.rstrip('\n') for line in f]

for pathway in pathways:
    # pathway = "MIF"
    pathway_dir = out_dir + "/" + pathway
    
    df_assoRes = pd.read_csv(pathway_dir + "/assoRes.txt", sep = "\t", index_col = 0)
    n_deg_genes = df_assoRes.shape[0]
    
    n_points = 50
    deg_pvalue_cutoff = 0.05
    
    string_step3 = 'assoRes <- assoRes[which(assoRes$pvalue_1 <= %f),]' % deg_pvalue_cutoff
    string_step3 = string_step3 + '\noAsso <- order(assoRes[,"waldStat_1"], decreasing=TRUE)'
    string_cluster = 'clusPat <- clusterExpressionPatterns(sce, nPoints = %d,' % n_points\
        + 'verbose=TRUE, genes = rownames(assoRes)[oAsso][1:min(%d,length(oAsso))],' % n_deg_genes \
        + ' k0s=4:5, alphas=c(0.1))'
    
    string_step3 = string_step3 + '\n' + string_cluster
    string_step3 = string_step3 + '\nyhatScaled <- data.frame(clusPat$yhatScaled)\n'
    
    with open(pathway_dir + "/step3.R", "w") as text_file:
        _ = text_file.write(string_step3)
