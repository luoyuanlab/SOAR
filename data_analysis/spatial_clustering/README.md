# Spatial Clustering

Implementation of spatial clustering using gene expression and spatial location via the [STAGATE package](https://stagate.readthedocs.io/en/latest/index.html)
- Requires `counts.csv` and `coordinates.csv` in the processed folder and a list of sample directories for running the shell scripts. 
- Steps of running:
    - Run [quest_stagate_jobarray.sh](https://github.com/luoyuanlab/ST-dataset/blob/main/analysis/database_utilities/clustering/quest_stagate_jobarray.sh) which compiles [quest_step01_stagate_jobarray.py](https://github.com/luoyuanlab/ST-dataset/blob/main/analysis/database_utilities/clustering/quest_step01_stagate_jobarray.py) followed by [quest_step02_stagate_jobarray.R](https://github.com/luoyuanlab/ST-dataset/blob/main/analysis/database_utilities/clustering/quest_step02_stagate_jobarray.R), and finally [quest_stagate_to_seurat_jobarray.R](https://github.com/luoyuanlab/ST-dataset/blob/main/analysis/database_utilities/clustering/quest_stagate_to_seurat_jobarray.R)
        -  quest_step01_stagate_jobarray.py:
            - **Please note that your provided argument (read into `sample_dir`) should have a "/" at the end of the path.**
            - Reads counts and coordinates into AnnData object and keeps only spots that have coordinates. 
            - If the AnnData contains only integer counts and the max count is greater than 20, it is likely not log transformed during processing, so log1p is done. 
            - Conducts normalization then finds the top 3000 highly variable genes by the CellRanger approach (expects normalized and transformed counts). 
            - Initial radius cutoff (`rad_cur`) was set to **2** for STAGATE to find neighbors for each spot. If the initial number of neighbors is less than **5**, step-wise addition of radius cutoff (`rad_add`) is done to reach at least 5 neighbors per spot. 
                - For non-visium technologies, `rad_add` might need to be increased(if too slow in reaching the optimal number)/decreased(if too many neighbors) to reach optimal neighbors between **5 - 15**. 
            - STAGATE spatial net is then trained. Neighbors are found and dimensionality reduction is done
            - Result is saved in `obsm_stagate.csv` in the clustering subfolder of the analysis folder. 
        - quest_step02_stagate_jobarray.py:
            - Reads in `obsm_stagate.csv` 
            - Uses mclust to find the optimal number of clusters with the **lowest BIC**.
                - Default number of clusters is between **2 - 30**, but if a sample has lower than 60 counts, the max number of clusters is **spots / 2**.
            - The result of cluster label assignment for each spot under the optimal number of clusters is saved in `obsm_stagate_cluster_*optimal_cluster_number*`.csv in the same folder. 
        - quest_stagate_to_seurat_jobarray.R:
            - Adds cluster assignment of each spot to the metadata of `Seurat.RDS` stored in the processed folder. 
