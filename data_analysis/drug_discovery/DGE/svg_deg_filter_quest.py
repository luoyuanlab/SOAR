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

svg = pd.read_csv(svg_dir+'SpatialDE_results.tsv', sep='\t')
svg_filter = svg[svg['qval']<=0.1]
print('svg read in')

if not os.path.exists(sample_dge_dir+'DGE_dec_SVG'):
    os.makedirs(sample_dge_dir+'DGE_dec_SVG')

for file in os.listdir(dge_dir):
    cell_type = file.split('.txt')[0]
    print(cell_type)
    deg = pd.read_csv(dge_dir+file, sep='\t')
    deg=deg[deg['gene'].isin(svg_filter['g'].tolist())]
    print('SVG filtered')

    df_up = deg[(deg['stat']>0.5) & (deg['qval']<0.05) ]
    df_down = deg[(deg['stat']<-0.5) & (deg['qval']<0.05) ]
    print('DEG up and down shape:', df_up.shape [0],df_down.shape[0])
    
    deg.to_csv(sample_dge_dir+'DGE_dec_SVG'+'/'+cell_type+'.csv',index=False)   
