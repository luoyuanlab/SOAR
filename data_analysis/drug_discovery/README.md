# Drug Discovery

Drug discovery analysis aims to identify repurposable and established compounds for targeting cell types of interests in pathological sample. 

## Drug Enrichment
['CMAP2 Enrichment Calculation'](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap2_drug_enrichment_perturbation_create_process_ds4a1.ipynb) (first section of code, sample level)
- Requires `CMAP2 ranked gene x compund perturbation results` (e.g., cmap2_6hr_pert_rank.csv). Genes (12,328) are essentially ranked at the compound level using signature-level perturbation z scores published by CMAP. Only 6hr perturbations (145,491) are included to reduce confounding factor of treatment duration and ensure consistency with CMAP1 results. To adjust filter criteria (e.g., perturbation duration, cell type, compound type), source perturbation data and metadata can be downloaded from ['CMAP']('https://clue.io/data/CMap2020#LINCS2020') (signature-level/level5 data). More explanation can be found from this CMAP ['article']('https://clue.io/connectopedia/replicate_collapse').
- Also requires `cell-type level DGE results` (e.g., DS4A/DS4A.1/DGE_dec/Endothelial.txt), which can be generated from ['scripts']('https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/DGE'). Filter is used to select up/down DEGs with |log2fc| > 0.5 and q-value < 0.05.
- Gene set enrichment analysis (GSEA) is conducted on the resulting set of up and down DGE sets for each CMAP compound, and then GSEA score difference between up and down DEG set is calculated as the enrichment score for the compound. P-value is calculated by comparing the proportion of enrichment score calculated from random gene rankings greater/less than the compound's enrichment score (depending on the sign).
- Compounds with the 200 highest (inverse enrichment, suppress DEGs) and lowest (positive enrichment, promote DEGs) enrichment score are saved for each cell type in a sample.
- Post-calculation formatting is done to prepare for website display. 

['CMAP1 Enrichment Processing'](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap1_drug_enrichment_perturbation_process_300ppi_create_batch.ipynb) (first secion of code, batch level) 
- Post-calculation formatting to prepare for website display. 

## Drug Perturbation
['CMAP2 Perturbation Network'](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap2_drug_enrichment_perturbation_create_process_ds4a1.ipynb) (second section of code)
- Requires `CMAP2 ranked gene x compund perturbation results`, `cell-type level DGE results` as mentioned above in addition to `list of top and bottom 200 enriched compounds from above` (e.g., DS4A/DS4A.1/Endothelial_POS_es.csv') and `gene x compound perturbation z scores` (e.g., cmap2_6hr_pert_z.npy).
- For each positively/inversely enriched compound, the CMAP perturbation z score of that compound on each up/down DEG is retrieved as the effect a compound has on the gene target (edge color).
- Top and bottom 30 DEGs with the highest absolute value of perturbation z score are retained for plotting the perturbation network.
- Log2fc of the top and bottom 30 DEGs is also saved in the result file as node color.

['CMAP1 Perturbation Processing'](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap1_drug_enrichment_perturbation_process_300ppi_create_batch.ipynb) (first secion of code, batch level) 
- Post-calculation formatting to prepare for website network illustration.

## Protein-protein Interaction
['Top 300 PPI Network'](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/cmap1_drug_enrichment_perturbation_process_300ppi_create_batch.ipynb) (second secion of code, batch level) 
- Requires `human protein interactome` (350k pairs) (e.g., interactome.tsv') and `cell-type level DGE results` as mentioned above. Interactome genes unmappable to HGNC genes are removed. 
- Top 300 DEGs with the highest absolute values is used to find matching PPIs in the interactome (both receiver and sender need to be from the top 300 DEG list).
- Log2fc of the top and bottom 30 DEGs is also saved in the result file as node color.

Maaping of NCBI gene IDs and symbols can be done using HGNC genes (e.g., 'HGNC.tsv'). 


