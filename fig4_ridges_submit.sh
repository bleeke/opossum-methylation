#!/bin/bash
# ridgeplots figure 3 
#SBATCH --job-name=ridges
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=30G
#SBATCH --partition=cpu
#SBATCH --array=0-2
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=ridges_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 




# define sample condition to be worked on
x=${SLURM_ARRAY_TASK_ID}
echo $x

SCRIPTS="fig4_adult_x_vs_auto.R fig4_adult_allele_x_vs_auto.R fig4_x_vs_auto_noEDTE.R "

arr=(${SCRIPTS})
echo ${arr[x]}

# run R script

Rscript ${arr[x]}


source deactivate

echo 'JOB ENDED'
date