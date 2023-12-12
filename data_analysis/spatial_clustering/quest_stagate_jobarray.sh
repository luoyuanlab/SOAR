#!/bin/bash

module purge all
module load python-miniconda3/4.12.0
source activate SpatialT

IFS=$'\n' read -d '' -r -a input_args < sample_list.tsv

if ${HOME}/.conda/envs/SpatialT/bin/python quest_step01_stagate_jobarray.py ${input_args[$SLURM_ARRAY_TASK_ID]}
then
	echo "STAGATE step 1 completed"
else
	echo "STAGATE step 1 failed"
	exit 1
fi

module purge all
module load R/4.1.1
module load geos/3.8.1
if Rscript quest_step02_stagate_jobarray.R ${input_args[$SLURM_ARRAY_TASK_ID]}
then
	echo "STAGATE step 2 completed"
else
	echo "STAGATE step 2 failed"
	exit 1
fi

if Rscript quest_stagate_to_seurat_jobarray.R ${input_args[$SLURM_ARRAY_TASK_ID]}
then
	echo "STAGATE step 3 completed"
else
	echo "STAGATE step 3 failed"
	exit 1
fi
