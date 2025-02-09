#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-4
#SBATCH --mem=10G
#SBATCH --mail-user=yanyi.ding@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ppi%a"
#SBATCH --output=/projects/b1131/ydn4687/SpatialT/cosmx_colon_revision/drug/ppi%a.out

module purge all
module load python-miniconda3/4.12.0
source activate myenv

IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/SpatialT/cosmx_patho_samples.tsv
IFS=$' ' read -ra split_dirs <<< ${input_args[${SLURM_ARRAY_TASK_ID}]}
dsid=${split_dirs[0]}
sampleid=${split_dirs[1]}

echo "DSID: ${dsid}"
echo "SampleiD: ${sampleid}"
python /projects/b1131/SpatialT/ST-dataset/analysis/database_utilities/drug/ppi_quest.py $dsid $sampleid