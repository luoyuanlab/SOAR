import pandas as pd
import os
import sys

dsid = sys.argv[1]
sampleid = sys.argv[2]
tech = sys.argv[3]
pid = 'PID'+dsid.split('DS')[1] 

dge_dir = '/projects/b1131/SpatialT/drug-target/'+dsid+'/'+sampleid+'/DGE_dec/'
svg_dir = '/projects/b1131/SpatialT/'+tech+'/'+pid+'/'+dsid+'/'+sampleid+'/analysis/SVG/'
sample_dge_dir = dge_dir.split('DGE_dec/')[0]
if not os.path.exists(sample_dge_dir+'DGE_dec_SVG'):
    os.makedirs(sample_dge_dir+'DGE_dec_SVG')

for file in os.listdir(dge_dir):
    cell_type = file.split('.txt')[0]
    print(cell_type)
    deg = pd.read_csv(dge_dir+file, sep='\t')
    if os.path.isfile(svg_dir+'SpatialDE_results_'+cell_type+'.tsv'):
        svg = pd.read_csv(svg_dir+'SpatialDE_results_'+cell_type+'.tsv', sep='\t')
        deg=deg[deg['gene'].isin(svg['g'].tolist())]
        print('SVG filtered')
    else:
        print('No SVG file')
    df_up = deg[(deg['stat']>0.5) & (deg['qval']<0.05) ]
    df_down = deg[(deg['stat']<-0.5) & (deg['qval']<0.05) ]
    print('DEG up and down shape:', df_up.shape [0],df_down.shape[0])
    
    deg.to_csv(sample_dge_dir+'DGE_dec_SVG'+'/'+cell_type+'.csv',index=False)   