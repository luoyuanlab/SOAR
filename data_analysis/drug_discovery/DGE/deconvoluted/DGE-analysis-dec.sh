#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=2-6
#SBATCH --mem=30G
#SBATCH --mail-user=yiming.li@northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="DGE2%a"
#SBATCH --output=/projects/b1131/ylz8811/pbs-cmds/DGE2_b4_%a.out

module purge all
module load R/4.1.1
module load geos/3.8.1

cd /projects/b1131/SpatialT/

IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/SpatialT/deconvoluted_samples.txt
IFS=$' ' read -ra split_dirs <<< ${input_args[${SLURM_ARRAY_TASK_ID}]}
sample_dir=${split_dirs[0]}

echo "Sample directory: ${sample_dir}"

Rscript --vanilla /projects/b1131/SpatialT/ref_scripts/DGE-analysis-dec.R $sample_dir
