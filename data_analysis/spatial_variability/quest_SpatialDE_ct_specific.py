import SpatialDE
import NaiveDE
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import os
import sys

st_dir = "/projects/b1131/SpatialT/" # On Quest

sample_dir = sys.argv[1]
cell_type = sys.argv[2]
print("\n\n>>> " + sample_dir + "[" + cell_type + "] started <<<")

### Read counts and coordinates
counts = pd.read_csv(sample_dir + 'analysis/deconvolution/counts_' + cell_type + '_deconv_only.csv')
counts_num = counts._get_numeric_data()
min_count = counts_num.min().min()
if (min_count < 0):
    counts_num[counts_num < 0] = 0

coordinates = pd.read_csv(sample_dir + 'processed/coordinates.csv')
counts.loc['Total',:]= counts.sum(axis=0)

### Align counts and coordinates index
error_count = 0
for i,j in zip (counts.columns.tolist()[1:], coordinates['barcode'].tolist()):
    if i != j:
        error_count = error_count + 1

if error_count > 0:
    print("[ERROR] " + sample_dir + " has not matching spot IDs.")
    sys.exit()

### Get total counts
total_counts = counts.iloc[-1][1:].tolist()

### Process data
sample_info = pd.DataFrame()
if 'x' in coordinates.columns:
    sample_info['y'] = coordinates['y']
    sample_info['x'] = coordinates['x']
else:
    print("[ERROR] " + sample_dir + " has problematic coordinates column names.")
    sys.exit(1)

sample_info['total_counts'] = total_counts
sample_info.index = coordinates['barcode']
# sample_info
reshaped_counts = counts.set_index('gene').iloc[:-1].transpose()
reshaped_counts.index = coordinates['barcode']
reshaped_counts = reshaped_counts.T[reshaped_counts.sum(0) >= 3].T
# reshaped_counts

### Run SpatialDE
try:
    norm_expr = NaiveDE.stabilize(reshaped_counts.T).T
except:
    norm_expr = np.log(reshaped_counts.T).T

resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T
# sample_resid_expr = resid_expr.sample(n=15202, axis=1, random_state=1)
X = sample_info[['x', 'y']]

try:
    results = SpatialDE.run(X.to_numpy(), resid_expr)
except:
    print("[ERROR] " + sample_dir + " SpatialDE failed, probably because the data is too sparse / contains too few spots for this cell type.")
    sys.exit(1)

if not os.path.exists(f'{sample_dir}analysis/'):
    os.makedirs(f'{sample_dir}analysis/')

if not os.path.exists(f'{sample_dir}analysis/SVG/'):
    os.makedirs(f'{sample_dir}analysis/SVG/')

### Write results to file
results.to_csv(sample_dir + 'analysis/SVG/SpatialDE_results_' + cell_type + '.tsv', sep = '\t', index = False)
# g - The name of the gene
# pval - The P-value for spatial differential expression
# qval - Significance after correcting for multiple testing
# l - A parameter indicating the distance scale a gene changes expression over
print("# SVG analysis [" + cell_type + "] finished.")

sign_results = results.query('qval < 0.05')
n_patterns = 5 # Default, hard-wired for now

if sign_results.shape[0]>0:
    ### Get average l
    l = pd.DataFrame(sign_results['l'].value_counts()).index.tolist()
    count = pd.DataFrame(sign_results['l'].value_counts())['count'].tolist()
    total_count = sum(count)
    total = 0
    for i,j in zip(l, count):
        ij = i*j
        total += ij
    
    L = round(total/total_count)
    histology_results, patterns = SpatialDE.aeh.spatial_patterns(X.to_numpy(), resid_expr, sign_results, C = n_patterns, l = L, verbosity = 1, delta_elbo_threshold = 1)
    print("# Pattern analysis [" + cell_type + "] finished.")
else:
    patterns = pd.DataFrame(columns=['0', '1'])
    histology_results = pd.DataFrame(columns=['g', 'pattern', 'membership'])
    print("# [WARNING] Cannot perform pattern analysis [" + cell_type + "], no sig genes.")

### Write results to file
histology_results.to_csv(sample_dir + 'analysis/SVG/SpatialDE_histology_results_' + cell_type + '.tsv', sep = '\t', index = False)
patterns.to_csv(sample_dir + 'analysis/SVG/SpatialDE_patterns_' + cell_type + '.tsv', sep = '\t', index = False)
print(">>> " + sample_dir + " finished <<<")
