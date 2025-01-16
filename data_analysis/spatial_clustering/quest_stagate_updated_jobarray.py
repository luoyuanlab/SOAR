import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import tensorflow as tf
from sklearn.mixture import GaussianMixture
from sklearn.metrics.cluster import adjusted_rand_score

import STAGATE

    
sample_dir = sys.argv[1]
tech = sample_dir.split("/")[4]
print(">>> " + sample_dir + " started first step STAGE clustering<<<")

# directory of counts and coordinates
counts_file = os.path.join(sample_dir, 'processed/counts.csv')
coor_file = os.path.join(sample_dir, 'processed/coordinates.csv')

if os.path.isfile(counts_file) and os.path.isfile(coor_file):
	# read and format data to anndata 
	counts = pd.read_csv(counts_file, index_col=0)
	coor_df = pd.read_csv(coor_file) 
	coor_df.set_index('barcode', drop=True, inplace=True)
	adata = sc.AnnData(counts.T)
	adata.var_names_make_unique()

	# keep only obs that are in coordinatesfile
	adata = adata[coor_df.index,]
	coor_df = coor_df.loc[adata.obs_names, ['x', 'y']]
	adata.obsm["spatial"] = coor_df.to_numpy()
	adata.raw = adata


	# check if need to log1p by finding non-int values across columns
	int_only = True       
	for col in counts.columns.tolist():
		col_int = counts[col].astype(str).str.isdigit().all()
		if col_int == False:
			int_only=False 
			print('NON-INT COUNTS')
			break

	if counts.to_numpy().max() >20 and int_only==True:
		sc.pp.log1p(adata)

	# normalization
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.filter_genes(adata,min_cells=5)
	sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=3000)

	tf.compat.v1.disable_eager_execution()
	rad_cur = 2
	STAGATE.Cal_Spatial_Net(adata, rad_cutoff=rad_cur)
	neighbors = adata.uns['Spatial_Net'].shape[0]/adata.n_obs
	print('TECH:', tech, ' INIT NEIGHBORS:', neighbors)

	# add radius_cutoff based on technology until reach at least 5 neighbors 
	if neighbors < 5 :
		while neighbors < 5:
			if tech == 'ST':
				rad_add = 1
			elif tech == 'DBiT-seq':
				rad_add = 1
			elif tech == '10x':
				rad_add = 2
			elif tech == 'seqFISH':
				rad_add = 5
			elif tech == 'MERFISH':
				rad_add = 30
			elif tech == 'Slide-seq':
				rad_add = 30
			elif tech == 'osmFISH':
				rad_add = 300
			else:
				rad_add = 10
				
			rad_cur = rad_cur + rad_add
			STAGATE.Cal_Spatial_Net(adata, rad_cutoff= rad_cur)
			neighbors = adata.uns['Spatial_Net'].shape[0]/adata.n_obs

	print(' FINAL RADIUS CUTOFF:', rad_cur,  'FINAL NEIGHBORS:', neighbors)
	#print(adata.uns['Spatial_Net'])

	#### Running STAGATE ####
	adata = STAGATE.train_STAGATE(adata, alpha=0)

	sc.pp.neighbors(adata, use_rep='STAGATE')
	sc.tl.umap(adata)

	# determine cluster resolution based on cell size
	if adata.shape[0] < 100:
	    res=1.2
	elif adata.shape[0] >=100 and adata.shape[0] <500:
	    res = 0.7
	elif adata.shape[0] >=500 and adata.shape[0] <5000:
	    res = 0.5    
	elif adata.shape[0] >=5000 and adata.shape[0] <20000:
	    res = 0.3    
	elif adata.shape[0] >=20000:
	    res = 0.1

	# clustering
	sc.tl.louvain(adata, resolution=res)
	sc.pl.embedding(adata, basis="spatial", color="louvain",s=6, show=False, title='STAGATE')

	# create clustering folder if does not exist
	if not os.path.exists(f'{sample_dir}analysis/clustering'):
		os.makedirs(f'{sample_dir}analysis/clustering')

	pd.DataFrame(adata.obsm['STAGATE'], index=adata.obs.index).to_csv(f'{sample_dir}analysis/clustering/STAGATE_30dim.csv')
	pd.DataFrame(adata.obs['louvain'], index=adata.obs.index).to_csv(f'{sample_dir}analysis/clustering/STAGATE_clusters.csv')
	print(">>> " + sample_dir + " finished first step STAGATE clustering<<<")
