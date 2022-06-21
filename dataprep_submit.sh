#!/bin/bash
# dataprep dataprep
#SBATCH --job-name=dataprep
#SBATCH --ntasks=2
#SBATCH --time=72:00:00
#SBATCH --mem=60G
#SBATCH --partition=cpu
#SBATCH --array=0-10
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=dataprep_01_%A_%a_submit.out

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

STAGES="7.5_embryo 6.5_embryo 1.5_embryo 2.5_embryo 3.5_embryo 4.5_embryo 5.5_embryo 7.5_ED 7.5_TE oocyte sperm"

arr=($STAGES)
echo ${arr[x]}

# define coverage threshold 
cov=5

# run R script

echo 'Rscript 01_dataprep.R'

Rscript 01_dataprep.R ${arr[x]} $cov

source deactivate

echo 'JOB ENDED'
date