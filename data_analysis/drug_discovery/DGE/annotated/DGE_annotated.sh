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
#SBATCH --job-name="DGE%a"
#SBATCH --output=/projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/deg/out/DGE_%a.out

module purge all
module load R/4.4.0

cd /projects/b1131/SpatialT/

IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/cosmx_colon_case_revision_samples.txt
IFS=$' ' read -ra split_dirs <<< ${input_args[${SLURM_ARRAY_TASK_ID}]}
sample_dir=${split_dirs[0]}

echo "Sample directory: ${sample_dir}"

Rscript --vanilla /projects/b1131/SpatialT/ref_scripts/DGE-annotated.R $sample_dir