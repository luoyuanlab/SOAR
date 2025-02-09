import pandas as pd
import os
import numpy as np
from math import isnan
import sys

ds = sys.argv[1]
sample = sys.argv[2]
print('PPI analysis for: ', ds, sample)
ppi_out_dir = '/projects/b1131/SpatialT/drug-target-cmap2-svg/'+ds+'/'+sample+'/PPI/'
if not os.path.exists(ppi_out_dir):
    os.makedirs(ppi_out_dir)
# for single cellular annotated samples
#deg_dir = '/projects/b1131/SpatialT/drug-target/'+ds+'/'+sample+'/DGE_anno_SVG/'
# for deconvoluted samples
deg_dir = '/projects/b1131/SpatialT/drug-target/'+ds+'/'+sample+'/DGE_dec_SVG/'

# read in gene reference and remove unformatted genes
gene = pd.read_table('/projects/b1131/SpatialT/cmap_ppi_database/HGNC.tsv')
gene = gene[['NCBI Gene ID(supplied by NCBI)','Approved symbol']]
gene.columns = ['num','name']
lst = []
for i in gene['num'].tolist():
    if pd.isna(i) != True:
        lst.append(int(i))
    else:
        lst.append(i)
gene['num'] = lst
gene_dict = gene.set_index('num').to_dict()['name']

clean_hgnc_dict = filter(lambda k: not isnan(k), gene_dict)
clean_hgnc_dict = {k: gene_dict[k] for k in gene_dict if not isnan(k)}
clean_hgnc_dict = {int(k):v for k,v in clean_hgnc_dict.items()}

# interactome for ppi network
interactome = pd.read_csv('/projects/b1131/SpatialT/cmap_ppi_database/Interactome.tsv', sep='\t')

# mapping genes to interactome
unmappable_pro = []
Protein_A_Name = []
Protein_B_Name = []

for i in interactome['#Protein A'].tolist():
    if i in clean_hgnc_dict.keys(): 
        Protein_A_Name.append(clean_hgnc_dict[i] )
    else:
        if i not in unmappable_pro: 
            unmappable_pro.append(i)
            print(i)
        Protein_A_Name.append('drop')

for i in interactome['Protein B'].tolist():
    if i in clean_hgnc_dict.keys(): 
        Protein_B_Name.append(clean_hgnc_dict[i] )
    else:
        if i not in unmappable_pro: 
            unmappable_pro.append(i)
            print(i)
        Protein_B_Name.append('drop')

# drop unmappable interactome interactions
interactome['Protein_A_Name'] = Protein_A_Name
interactome['Protein_B_Name'] = Protein_B_Name
dropa_index = interactome[interactome['Protein_A_Name'] == 'drop'].index.tolist()
dropb_index = interactome[interactome['Protein_B_Name'] == 'drop'].index.tolist()
print(len(dropa_index), len(dropb_index))
# 538 interactions removed due to unmappable gene number 
for i in dropa_index:
    if i not in dropb_index:
        dropb_index.append(i)
len(dropb_index)
clean_interactome = interactome.drop(dropb_index)
print('interactome cleaned')

# generate ppi network using SV-DGEs
for deg_file in os.listdir(deg_dir):
    # read in deg-svg file and filter (no svg fitler placed on cell types that don't have SVG analysis )
    cell_type = deg_file.split('.csv')[0]
    log2fc = pd.read_csv(deg_dir+deg_file)
    log2fc['abs_stat'] = [abs(i) for i in log2fc['stat'].tolist()]
    log2fc_top300 = log2fc.sort_values('abs_stat', ascending=False)[:300]
    log2fc_top300_sig = log2fc_top300[log2fc_top300['qval']<0.05]
    print('sv-deg file read in')

    # map degs to inteactome
    ppi_index_interactome_top300 = []
    for i in clean_interactome.iterrows():
        row = i[1]
        if row['Protein_A_Name'] in log2fc_top300_sig['gene'].tolist() and row['Protein_B_Name'] in log2fc_top300_sig['gene'].tolist():
            ppi_index_interactome_top300.append(i[0])
    ppi_top300 = clean_interactome[clean_interactome.index.isin(ppi_index_interactome_top300)][['#Protein A','Protein B', 'Protein_A_Name','Protein_B_Name']]
    print('sv-deg mapped to interactome')

    # drop ppis where both proteins are the same
    same_pro_index_top300 = []
    for i in ppi_top300.iterrows():
        row = i[1]
        if row['#Protein A'] == row['Protein B']:
            same_pro_index_top300.append(i[0])
    ppi_top300_nored = ppi_top300.drop(same_pro_index_top300)

    # format edges list for output
    ppi_top300_nored = ppi_top300_nored[['Protein_A_Name','Protein_B_Name']]
    ppi_top300_nored.columns = ['Source','Target']

    # format nodes list for output
    ppi_uniq_pro = list(set(ppi_top300_nored['Source'].tolist()) )
    for i in ppi_top300_nored['Target'].tolist():
        if i not in ppi_top300_nored['Source'].tolist() and i not in ppi_uniq_pro :
            ppi_uniq_pro.append(i)

    ppi_nodes_top300 = pd.DataFrame()
    ppi_nodes_top300['ID'] = ppi_uniq_pro
    ppi_nodes_top300['Label'] = ppi_uniq_pro

    ppi_nodes_top300 = pd.merge(ppi_nodes_top300, log2fc[['gene','stat']], how = 'left', right_on = 'gene', left_on = 'Label')
    ppi_nodes_top300 = ppi_nodes_top300.drop('gene',axis = 1)

    ppi_nodes_top300['sign'] = [np.sign(i) for i in ppi_nodes_top300['stat'].tolist()]
    ppi_nodes_top300['abs_stat'] = [abs(i) for i in ppi_nodes_top300['stat'].tolist()]

    # save
    ppi_nodes_top300.to_csv(ppi_out_dir+cell_type+'_ppi_nodes.csv', index=False)
    ppi_top300_nored.to_csv(ppi_out_dir+cell_type+'_ppi.csv', index=False)
    print(cell_type,'ppi saved')
