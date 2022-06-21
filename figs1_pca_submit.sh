#!/bin/bash
# dataprep dataprep
#SBATCH --job-name=pca
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=pca_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 

# define command line arguments provided to R script


# run R script

echo 'Rscript figs1_pca.R'

Rscript figs1_pca.R

source deactivate

echo 'JOB ENDED'
date