#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-3
#SBATCH --mem=50G
#SBATCH --mail-user=yanyi.ding@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="stagate%a" 
#SBATCH --output=/projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/stagate/out/stagate_jobarray%a.out

module purge all
module load python-miniconda3/4.12.0
source activate myenv

##cd ${HOME}/SpatialT/ST-dataset/
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/cosmx_colon_case_revision_samples.txt

python /projects/b1131/SpatialT/ST-dataset/analysis/database_utilities/clustering/quest_stagate_updated_jobarray.py ${input_args[$SLURM_ARRAY_TASK_ID]}
##echo ${input_args[$SLURM_ARRAY_TASK_ID]}