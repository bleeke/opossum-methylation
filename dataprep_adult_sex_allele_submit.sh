#!/bin/bash
# dataprep dataprep
#SBATCH --job-name=dataprep
#SBATCH --ntasks=2
#SBATCH --time=72:00:00
#SBATCH --mem=60G
#SBATCH --partition=cpu
#SBATCH --array=0-11
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=dataprep_adult_sex_allele_%A_%a_submit.out

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

STAGES="brain_male_genome1 brain_male_genome2 brain_female_genome1 brain_female_genome2 liver_male_genome1 liver_male_genome2 liver_female_genome1 liver_female_genome2 spleen_male_genome1 spleen_male_genome2 spleen_female_genome1 spleen_female_genome2"  

arr=($STAGES)
echo ${arr[x]}

# define coverage threshold 
cov=5
echo $cov

# run R script

echo 'Rscript 01_dataprep_adult_sex_allele.R'

Rscript 01_dataprep_adult_sex_allele.R ${arr[x]} $cov

source deactivate

echo 'JOB ENDED'
date