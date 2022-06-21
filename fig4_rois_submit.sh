#!/bin/bash
# rsx roi plot gametes
#SBATCH --job-name=rsx_roi
#SBATCH --ntasks=1
#SBATCH --time=2:00
#SBATCH --mem=10G
#SBATCH --partition=cpu
#SBATCH --array=0-5
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=rsx_roi_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 




# define sample condition to be worked on
x=${SLURM_ARRAY_TASK_ID}
echo $x

SCRIPTS="fig4_spot_plot_embryo_rsx.R fig3_spot_plot_gametes_rsx.R fig4_spot_plot_adult_rsx.R"

arr=($SCRIPTS})
echo ${arr[x]}

# run R script

Rscript ${arr[x]}


source deactivate

echo 'JOB ENDED'
date