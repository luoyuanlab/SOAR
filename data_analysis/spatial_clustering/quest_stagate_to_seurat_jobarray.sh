#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-8
#SBATCH --mem=30G
#SBATCH --mail-user=yanyi.ding@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="stag_to_seurat%a" 
#SBATCH --output=/projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/stagate/out/stag_to_seurat_jobarray%a.out

##cd ${HOME}/SpatialT/ST-dataset/
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/cosmx_colon_case_revision_samples.txt

module purge all
module load R/4.4.0

Rscript /projects/b1131/SpatialT/ST-dataset/analysis/database_utilities/clustering/quest_stagate_to_seurat_updated_jobarray.R ${input_args[$SLURM_ARRAY_TASK_ID]}