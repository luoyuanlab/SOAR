#!/bin/bash

IFS=$'\n' read -d '' -r -a input_args < sample_list.txt
IFS=$' ' read -ra split_dirs <<< ${input_args[${SLURM_ARRAY_TASK_ID}]}
sample_dir=${split_dirs[0]}
ref_dir=${split_dirs[1]}
distance_type="short" # Supported arguments: short, medium, long, or xlong

echo "### Sample directory: ${sample_dir}"
echo "### Reference directory: ${ref_dir}"
echo "### Distance type: ${distance_type}"

module purge all
module load python-miniconda3/4.12.0
source activate SpatialT
if ${HOME}/.conda/envs/SpatialT/bin/python cci-analysis-COMMOT.py $sample_dir $ref_dir $distance_type
then
	echo "COMMOT main analysis completed"
else
	echo "COMMOT main analysis failed"
	exit 1
fi

if ${HOME}/.conda/envs/SpatialT/bin/python cci-analysis-COMMOT-DGE-step0.py $sample_dir $ref_dir $distance_type
then
	echo "DGE - Step 0 (Python) completed"
else
	echo "DGE - Step 0 (Python) failed"
	exit 1
fi

source deactivate
module purge all
module load R/4.1.1
module load geos/3.8.1
if Rscript --vanilla cci-analysis-COMMOT-DGE-step1.R $sample_dir $distance_type
then
	echo "DGE - Step 1 (R) completed"
else
	echo "DGE - Step 1 (R) failed"
	exit 1
fi

module purge all
module load python-miniconda3/4.12.0
source activate SpatialT
if ${HOME}/.conda/envs/SpatialT/bin/python cci-analysis-COMMOT-DGE-step2.py $sample_dir $distance_type
then
	echo "DGE - Step 2 (Python) completed"
else
	echo "DGE - Step 2 (Python) failed"
	exit 1
fi

source deactivate
module purge all
module load R/4.1.1
module load geos/3.8.1
if Rscript --vanilla cci-analysis-COMMOT-DGE-step3.R $sample_dir $distance_type
then
	echo "DGE - Step 3 (R) completed"
else
	echo "DGE - Step 3 (R) failed"
	exit 1
fi

module purge all
module load python-miniconda3/4.12.0
source activate SpatialT
if ${HOME}/.conda/envs/SpatialT/bin/python cci-analysis-COMMOT-DGE-step4.py $sample_dir $distance_type
then
	echo "DGE - Step 4 (Python) completed"
else
	echo "DGE - Step 4 (Python) failed"
	exit 1
fi

if ${HOME}/.conda/envs/SpatialT/bin/python cci-analysis-COMMOT-pull-scores.py $sample_dir $ref_dir $distance_type
then
	echo "Pulling scores completed"
else
	echo "Pulling scores failed"
	exit 1
fi
