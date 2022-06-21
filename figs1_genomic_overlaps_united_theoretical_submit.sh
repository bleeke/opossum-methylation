#!/bin/bash
# dataprep dataprep
#SBATCH --job-name=overlaps
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --partition=cpu
#SBATCH --array=0-1
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=genomic_overlaps_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 

# define command line arguments provided to R script

# define sample condition to be worked on
x=${SLURM_ARRAY_TASK_ID}
echo $x

STAGES="united_pooled_destranded_filt5 all_mondom_cpgs"

arr=($STAGES)
echo ${arr[x]}


DIR="/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out /camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"

dir_arr=($DIR)
echo ${dir_arr[x]}


# run R script

echo 'Rscript genomic_overlaps_united_theoretical.R'

Rscript genomic_overlaps_united_theoretical.R ${arr[x]} ${dir_arr[x]}

source deactivate

echo 'JOB ENDED'
date