import numpy as np
np.random.seed(1024)
import pandas as pd
import statsmodels
from statsmodels.stats.multitest import fdrcorrection
import os
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler(feature_range = (-2,2))
import math
import sys

dsid = sys.argv[1]
sampleid = sys.argv[2]
deg_dir = '/projects/b1131/SpatialT/drug-target/'+dsid+'/'+sampleid+'/DGE_dec_SVG'
cell_direc = [i for i in os.listdir(deg_dir)]
print('STARTING:', dsid, ' ', sampleid)

# output dir
out_dir = '/projects/b1131/SpatialT/drug-target-cmap2-svg/'
enrich_dir = out_dir+dsid+'/'+sampleid+'/GSEA/'
top_enrich_dir = out_dir+dsid+'/'+sampleid+'/Enrichment/'
perturb_dir = out_dir+dsid+'/'+sampleid+'/Perturbation/'

##### gene processing #####
# hgnc genes
hgnc_v2 = pd.read_csv('/projects/b1131/SpatialT/cmap_ppi_database/HGNC091923.txt', sep='\t')
hgnc_dict = dict(zip(hgnc_v2['Approved symbol'], hgnc_v2['NCBI gene ID']))
hgnc_v2_prev = hgnc_v2[hgnc_v2['Previous symbol'].notnull()]
hgnc_dict_prev = dict(zip(hgnc_v2_prev['Previous symbol'], hgnc_v2_prev['NCBI gene ID']))

#cmap genes
gene_info = pd.read_csv('/projects/b1131/SpatialT/cmap_ppi_database/geneinfo_beta.txt', sep='\t')

# mapping gene indices to real names 
gene_cmap2 = gene_info[['gene_id','gene_symbol']]
gene_cmap2.columns = ['num','name']
lst = []
for i in gene_cmap2['num'].tolist():
    if pd.isna(i) != True:
        lst.append(int(i))
    else:
        lst.append(i)
gene_cmap2['num'] = lst
gene_dict_cmap2 = gene_cmap2.set_index('num').to_dict()['name']
gene_dict_cmap2_rev = {v:k for k,v in gene_dict_cmap2.items()}
print('Genes cleaned for mapping')


##### cmap2 enrichment setup #####
# rank file of 6hr perturbations (from z scores)- gene x compound effects
cmap2_cmap1_perturb_rank_dur = pd.read_csv('/projects/b1131/SpatialT/cmap_ppi_database/cmap2_6hr_pert_rank.csv', engine='c')
cmap2_cmap1_perturb_rank_dur = cmap2_cmap1_perturb_rank_dur.set_index('Unnamed: 0')
cmap2_cmap1_perturb_rank_dur_arry = np.array(cmap2_cmap1_perturb_rank_dur)

# gene and compound name setup for gsea 
ROWS_dur = [str(gene_dict_cmap2_rev[i]) for i in cmap2_cmap1_perturb_rank_dur.index]
COLS_dur = cmap2_cmap1_perturb_rank_dur.columns.tolist()
ROW2IDX2_val_dur = list(range(0,len(ROWS_dur)))
COL2IDX2_val_dur = list(range(0,len(COLS_dur)))
ROW2IDX_key_dur = ROWS_dur
COL2IDX_key_dur = COLS_dur
ROW2IDX_dur = {k:v for k,v in zip(ROW2IDX_key_dur, ROW2IDX2_val_dur)}
COL2IDX_dur = {k:v for k,v in zip(COL2IDX_key_dur, COL2IDX2_val_dur)}
N_GENES, N_PERTURBATIONS = cmap2_cmap1_perturb_rank_dur_arry.shape
N_REPEATS = 1000
RANDOM_RANK = np.random.rand(N_GENES, N_REPEATS).argsort(axis=0) + 1  # gene x repeat


##### cmap2 perturbation setup #####
# cmap2 perturbation zscore
z_score = np.load('/projects/b1131/SpatialT/cmap_ppi_database/cmap2_6hr_pert_z.npy')
z_score_df = pd.DataFrame(z_score)
z_score_df.columns = cmap2_cmap1_perturb_rank_dur.columns
z_score_df.index = cmap2_cmap1_perturb_rank_dur.index
print('cmap data read in')


##### format dge results ######
def get_up_down_deg (df):
    df_up = df[(df['stat']>0.5) & (df['qval']<0.05) ]
    df_down = df[(df['stat']<-0.5) & (df['qval']<0.05) ]
    return df_up, df_down

# transform from genename to number
def hgnc_gene_up_down(up_df, down_df):
    gene_num_up = []
    for i in up_df['gene'].tolist():
        if i in hgnc_dict.keys() and math.isnan(hgnc_dict[i]) == False:
            gene_num_up.append(str(int(hgnc_dict[i])))
        elif i in hgnc_dict_prev.keys() and math.isnan(hgnc_dict_prev[i]) == False:
            gene_num_up.append(str(int(hgnc_dict_prev[i])))
        else:
            print('SKIP GENE',i)

    gene_num_down = []
    for i in down_df['gene'].tolist():
        if i in hgnc_dict.keys() and math.isnan(hgnc_dict[i]) == False:
            gene_num_down.append(str(int(hgnc_dict[i])))
        elif i in hgnc_dict_prev.keys() and math.isnan(hgnc_dict_prev[i]) == False:
            gene_num_down.append(str(int(hgnc_dict_prev[i])))
        else:
            print('SKIP GENE', i) 
    return gene_num_up, gene_num_down

##### GSEA ####
def cmap2_gsea_setup(gene_num_up, gene_num_down):
    DATA = {
        "U": [ROW2IDX_dur[g] for g in gene_num_up if g in ROW2IDX_dur],
        "D": [ROW2IDX_dur[g] for g in gene_num_down if g in ROW2IDX_dur],
    }
    N_UP = len(DATA["U"])
    N_DOWN = len(DATA["D"])
    #print(N_UP, N_DOWN)
    return DATA, N_UP, N_DOWN

# up/down gene set level es
def STEP1(v):
    # deg sorted up/down reg gene perrturb array
    v = sorted(v)
    # t = length of up/down reg gene perturb array
    t = len(v)
    
    # (i+1)/t = fraction of elements up to index i in sorted up/down gene list 
    # v[i]/n_genes = fraction of the rank of the current gene at index i out of the total number of genes
    # total upward (a) and downward (b) influence (es) for each up/down reg gene set
    a = max((i + 1) / t - v[i] / N_GENES for i in range(t)) 
    b = max(v[i] / N_GENES - i / t for i in range(t))
    es = 0
    if a > b: #check what a and b is
        es = a
    elif b > a:
        es = -b
    #print('a:',a,'b:', b)
    return es

# diff of up/down gene set es 
def STEP2(rank_u, rank_d):
    #print(rank_u)
    # calculate up genes and down genes es difference for each perturbation to determine overall es  
    es_u = STEP1(rank_u) if rank_u.size else 0
    es_d = STEP1(rank_d) if rank_d.size else 0
    #print(es_u, es_d)
    if np.sign(es_u) == np.sign(es_d):
        return 0
    return es_u - es_d

def run_gsea(DATA, N_UP, N_DOWN):
    background = np.array([STEP2(RANDOM_RANK[:N_UP, c], RANDOM_RANK[N_UP:N_UP + N_DOWN, c]) for c in range(N_REPEATS)])
    result = []
    # for each perturbation, calculate the es score difference of up genes and down genes
    for col in range(N_PERTURBATIONS):
        es = STEP2(cmap2_cmap1_perturb_rank_dur_arry[DATA["U"], col], cmap2_cmap1_perturb_rank_dur_arry[DATA["D"], col])
        #print(es)
        p = 1.0
        if es > 0:
            p = np.mean(background > es)
        elif es < 0:
            p = np.mean(background < es)
        result.append((es, p))
    return result


##### format gsea #####
# all drugs enrichment score
def format_gsea_res(result, cell):
    res_formatted = pd.DataFrame(result)
    res_formatted.columns = ['es','p']
    res_formatted['es'] = [round(i,3) for i in res_formatted['es'].tolist()]
    res_formatted['p-adj']=statsmodels.stats.multitest.fdrcorrection(res_formatted['p'].tolist(), alpha=0.05, method='indep', is_sorted=False)[1]
    res_formatted = res_formatted.reset_index()
    res_formatted.index = COLS_dur
    #res_formatted = res_formatted.drop('level_0',axis = 1)
    res_formatted['drug'] = [i.split('--')[0] for i in res_formatted.index.tolist()]
    if not os.path.exists(enrich_dir ):
        os.makedirs(enrich_dir )
    res_formatted.to_csv(enrich_dir+cell+'.csv')
    return res_formatted


# top and bottom 500 enriched drugs - for website visualization
def format_enrich(ds, cell):
    res = pd.DataFrame()
    res['CMAP_instance'] = [i.split('--')[1] for i in ds.index.tolist()]
    res['Drug'] = ds['drug'].tolist()
    res['Enrichment_score'] = ds['es'].tolist()
    res['P-value'] = ds['p'].tolist()
    res['Adj-p']= ds['p-adj'].tolist()
    res = res.sort_values('Enrichment_score', ascending=False)
    res_inv = res.head(500)
    res_pos = res.tail(500)
    if not os.path.exists(top_enrich_dir ):
        os.makedirs(top_enrich_dir )
    res_inv.to_csv(top_enrich_dir+cell+'_INV_es.csv', index=False)
    res_pos.to_csv(top_enrich_dir+cell+'_POS_es.csv', index=False)
    return res_inv, res_pos


##### format perturbation #####
# format pertrbation network to create edge and node file
def get_perturb_zs(cell, direction, es_df, dge, cell_type_deg_z):
    perturb_specific_dir = perturb_dir+cell+'/'+direction+'/'
    if not os.path.exists(perturb_specific_dir):
        os.makedirs(perturb_specific_dir)  
    for j,i in es_df.iterrows():
        # adding cmap perturbation rank and z score for each gene for top enriched compounds. 
        cmap_instance = i['CMAP_instance']
        drug = i['Drug']
        cmap_drug = drug + '--' +cmap_instance
        #print(j, drug_new_idx)
        zs = cell_type_deg_z[[cmap_drug]]
        zs.columns = ['z']
        zs = zs.sort_values('z')
        perturb_df = zs.merge(dge, how = 'left', left_index=True, right_on = 'gene')
        
        # edge of perturbation network
        deg_gene_num = perturb_df.shape[0]
        perturb_edges = pd.DataFrame()
        perturb_edges['Source'] = [drug for i in range(deg_gene_num)]
        perturb_edges['CMAP_instance'] = [cmap_instance for i in range(deg_gene_num)]
        perturb_edges['Target'] = perturb_df['gene'].tolist()
        perturb_edges['log2fc'] = perturb_df['stat'].tolist()
        perturb_edges['exp_z'] = perturb_df['z'].tolist()
        perturb_edges['exp_z_norm'] = [i[0] for i in scaler.fit_transform(np.array(perturb_edges['exp_z'].tolist()).reshape(-1, 1))]
        perturb_edges['abs_log2fc'] = [abs(i) for i in perturb_df['stat'].tolist()]
        perturb_edges['sign_log2fc'] = [np.sign(i) for i in perturb_df['stat'].tolist()]
        perturb_edges = perturb_edges.sort_values('exp_z')
        head_edges = perturb_edges.head(30)
        tail_edges = perturb_edges.tail(30)
        perturb_edges_sub = pd.concat([head_edges, tail_edges])
        perturb_edges_sub.to_csv(perturb_specific_dir+cmap_drug+'_perturb.csv',index=False)
        
        # node of perturbation network 
        perturb_nodes = pd.DataFrame()
        perturb_nodes['Id'] = perturb_edges_sub['Target'] 
        perturb_nodes['Label'] = perturb_edges_sub['Target'] 
        perturb_nodes['log2fc'] = perturb_edges_sub['log2fc']
        perturb_nodes['sign_log2fc'] = perturb_edges_sub['sign_log2fc'] 
        perturb_nodes['abs_log2fc'] = perturb_edges_sub['abs_log2fc'] 
        drug_row = pd.DataFrame({'Id':[perturb_edges_sub['Source'].unique().tolist()[0]],
                             'Label':[perturb_edges_sub['Source'].unique().tolist()[0]],
                             'log2fc':[10],
                             'sign_log2fc':[0],
                             'abs_log2fc':[10]})
        perturb_nodes = pd.concat([drug_row, perturb_nodes])
        perturb_nodes.to_csv(perturb_specific_dir+cmap_drug+'_perturb_nodes.csv',index=False)

##### run enrich + perturb ######
def run_enrich_perturb(dsid, sampleid, deg_dir, cell_direc):
    for ct in cell_direc:
        # deg
        cell = ct.split('.csv')[0]
        deg = pd.read_csv(deg_dir + '/' + ct)
        print('deg file read in for:', cell)
        deg_up, deg_down = get_up_down_deg(deg)
        print('deg shape:', deg_up.shape, deg_down.shape)
        if deg_up.shape[0] > 2000 or deg_up.shape[0] <10 or deg_down.shape[0] > 2000 or deg_down.shape[0] <10:
            print('deg too few or too many')
            print('\n')
            continue
        gene_num_up, gene_num_down =  hgnc_gene_up_down(deg_up, deg_down )
        print('degs hgnc-formatted:', dsid, sampleid)    

        # enrichment
        DATA, N_UP,N_DOWN =  cmap2_gsea_setup(gene_num_up, gene_num_down)
        print('cmap matched deg:',N_UP,N_DOWN)
        gsea_res= run_gsea(DATA, N_UP, N_DOWN)
        print('finish gsea, head:', gsea_res[:3])
        gsea_res_form = format_gsea_res(gsea_res, cell)
        enrich_inv, enrich_pos = format_enrich(gsea_res_form, cell)
        print('enrich saved, shape:', enrich_inv.shape, enrich_pos.shape, enrich_inv.head(1), enrich_inv.head(1))

        # perturbation
        dge =pd.concat([deg_up, deg_down]).sort_values('stat')
        cell_type_deg_z = z_score_df[z_score_df.index.isin(dge['gene'].tolist())]
        print('deg genes matched to cmap:',cell_type_deg_z.shape[0])
        get_perturb_zs(cell, 'INV', enrich_inv, dge, cell_type_deg_z)
        get_perturb_zs(cell, 'POS', enrich_inv, dge, cell_type_deg_z)
        print('inv and pos perturb files saved, shapes:', len(os.listdir(perturb_dir+cell+'/INV/')),len(os.listdir(perturb_dir+cell+'/POS/')) )
        print('\n')

run_enrich_perturb(dsid, sampleid, deg_dir, cell_direc)