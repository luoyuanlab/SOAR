# Drug Discovery

Drug discovery analysis aims to identify repurposable and established compounds for targeting cell types of interests in pathological sample. Analysis is conducted on spatially variable and differentially expressed (SV-DE) genes for each deconvoluated cell type. 

## Drug Enrichment and Perturbation [Script](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/drug_screen_perturb_quest.py) 

### Input
- `CMAP L1000 perturbation profile` (e.g., cmap2_6hr_pert_z.npy), which can be downloaded from [CMAP](https://clue.io/data/CMap2020#LINCS2020) signature/level 5 data. This file contains the perturbation MODZ score calculated for each compound perturbation on each gene (12,328 total). We included 6hr perturbations (145,491) to reduce redunduncy and avoid the confounding factor of treatment duration. More explanation on MODZ score can be found from this CMAP [article](https://clue.io/connectopedia/replicate_collapse). Additional metadata can also be downloade from CMAP to filter for cell type, dosage, etc. 

- `CMAP L1000 rank profile` (e.g., cmap2_6hr_pert_rank.csv), which is created by ranking the perturbation MODZ score of genes for each of the 145,491 perturbations.

- `Gene name mappings` (e.g., geneinfo_beta.txt and HGNC091923.txt) for mapping between gene IDs and symbols

- `Cell-type level SV-DE genes` (e.g., dsid/sampleid/DGE_dec_SVG/). SVG can be generated from [script](https://github.com/luoyuanlab/SOAR/tree/main/data_analysis/spatial_variability), and DGE can be generated from the [script](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/DGE). Filter is used to select spatially variable up and down DE gene sets with |log2fc| > 0.5 and q-value < 0.05. For cell types that don't have enough deconvoluted spots for SVG analysis, DGE results are used directly.

### CMAP2 Enrichment Calculation
- Gene set enrichment analysis (GSEA) is conducted on the resulting set of up and down DGE sets (10 <= gene set size <= 2000) for each CMAP compound perturbation. GSEA score difference between up and down DEG set is calculated as the enrichment score for the compound. P-value is calculated by comparing the proportion of enrichment score calculated from random gene rankings greater/less than the compound's enrichment score (depending on the sign). More information can be found at the Methods section of [manuscript] (https://www.biorxiv.org/content/10.1101/2022.04.17.488596v3).
- Compounds with the 500 highest (inverse enrichment, suppress DEGs) and lowest (positive enrichment, promote DEGs) enrichment score are saved for every cell type as output(e.g., dsid/sampleid/Enrichment/). 

### CMAP2 Perturbation Network
- For each of the 500 positively and 500 inversely enriched compound perturbation, CMAP perturbation MODZ score of the compound on the SV-DE genes DEGs represents the effect a compound has on the gene target.
- Top and bottom 30 SV-DE genes with the highest absolute value of perturbation MODZ score are saved  along with the log2fc of the SV-DE genes for plotting the perturbation network (e.g., dsid/sampleid/Perturbation/). 


## Protein-protein Interaction [Script](https://github.com/luoyuanlab/SOAR/blob/main/data_analysis/drug_discovery/PPI_Drug_Enrichment_Perturbation/ppi_quest.py)

### Input
- `Human protein interactome` (e.g., Interactome.tsv'), which contains 350k pairs of PPIs. Maping of NCBI gene IDs and symbols can be done using HGNC reference (e.g., 'HGNC.tsv'). 

- `Cell-type level SV-DE genes` similar to that required for drug screen. 

### PPI Network
- Top 300 SV-DE genes with the highest absolute values is used to find matching PPIs in the interactome (both receiver and sender need to be from the top 300 DEG list).
- Results along with log2fc of the SV-DE genes are saved for plotting PPI network. 



