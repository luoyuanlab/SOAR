# Drug Discovery

Drug discovery analysis aims to identify repurposable and established compounds for targeting cell types of interests in pathological sample. 

## Drug Enrichment
[CMAP2 Enrichment Calculation](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap2_drug_enrichment_perturbation_sample_level.ipynb) 
- Requires `CMAP2 ranked gene x compound perturbation results` (e.g., cmap2_6hr_pert_rank.csv). Genes (12,328) are essentially ranked at the compound level using signature-level perturbation z scores published by CMAP. Only 6hr perturbations (145,491) are included to reduce confounding factor of treatment duration and ensure consistency with CMAP1 results. To adjust filter criteria (e.g., perturbation duration, cell type, compound type), source perturbation data and metadata can be downloaded from [CMAP](https://clue.io/data/CMap2020#LINCS2020) (signature-level/level5 data). More explanation can be found from this CMAP ['article'](https://clue.io/connectopedia/replicate_collapse).
- Also requires `cell-type level SVG-filtered DGE results` (e.g., DS4A/DS4A.1/DGE_dec_SVG/Endothelial.txt), which can be generated from [DGE scripts](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/DGE) and [SVG scripts](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_variability). Filter is used to select up/down DEGs with |log2fc| > 0.5 and q-value < 0.05. For cell types that don't have enough deconvoluted spots for SVG analysis, DGE results are used directly. 
- Gene set enrichment analysis (GSEA) is conducted on the resulting set of up and down DGE sets (10 <= gene set size <= 2000) for each CMAP compound, and then GSEA score difference between up and down DEG set is calculated as the enrichment score for the compound. P-value is calculated by comparing the proportion of enrichment score calculated from random gene rankings greater/less than the compound's enrichment score (depending on the sign).
- Compounds with the 500 highest (inverse enrichment, suppress DEGs) and lowest (positive enrichment, promote DEGs) enrichment score are each saved for every cell type (e.g., drug_screen/drug_enrichment_cmap2_SVG/DS4A/DS4A.1/Malignant_INV_es.csv). 


## Drug Perturbation
[CMAP2 Perturbation Network](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap2_drug_enrichment_perturbation_sample_level.ipynb)
- Requires a list of the `top 500 INV/POS enriched compounds generated from drug enrichment` above(e.g., drug_screen/drug_enrichment_cmap2_SVG/DS4A/DS4A.1/Malignant_INV_es.csv) and `gene x compound perturbation z scores` (e.g., cmap2_6hr_pert_z.npy).
- For each positively/inversely enriched compound, the CMAP perturbation z score of that compound on each up/down DEG is retrieved as the effect a compound has on the gene target (edge color).
- Top and bottom 30 DEGs with the highest absolute value of perturbation z score are retained for plotting the perturbation network.
- Resulting node and edge files for plotting drug perturbation networked are saved (e.g., drug_screen/drug_perturbation_cmap2_SVG/DS4A/DS4A.1/Malignant/INV/zolpidem--CPC020_HT29_6H:BRD-K44876623-001-02-9:10_perturb.csv and zolpidem--CPC020_HT29_6H:BRD-K44876623-001-02-9:10_perturb_nodes.csv).
- Log2fc of the top and bottom 30 SVG-filtered DEGs is also saved in the result file as node color.


## Protein-protein Interaction
[Top 300 PPI Network](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap1_drug_enrichment_perturbation_process_300ppi_create_batch.ipynb) (second secion of code, batch level) 
- Requires `human protein interactome` (350k pairs) (e.g., interactome.tsv') and `cell-type level SVG-filtered DGE results` as mentioned above. Interactome genes unmappable to HGNC genes are removed. 
- Top 300 DEGs with the highest absolute values is used to find matching PPIs in the interactome (both receiver and sender need to be from the top 300 DEG list).
- Log2fc of the top and bottom 30 DEGs is also included as node color.
- Resulting node and edge files for plotting PPI net work are saved (e.g., drug_screen/ppi_svg/DS4A/DS4A.1/Malignant_ppi.csv and Malignant_ppi_nodes.csv)

Maaping of NCBI gene IDs and symbols can be done using HGNC genes (e.g., 'HGNC.tsv'). 


