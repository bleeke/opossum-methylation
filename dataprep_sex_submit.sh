#!/bin/bash
# dataprep dataprep
#SBATCH --job-name=dataprep
#SBATCH --ntasks=2
#SBATCH --time=72:00:00
#SBATCH --mem=60G
#SBATCH --partition=cpu
#SBATCH --array=0-19
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=dataprep_sex_%A_%a_submit.out

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

STAGES="sperm_male oocyte_female 1.5_embryo_female 2.5_embryo_female 3.5_embryo_female 4.5_embryo_female 5.5_embryo_female 6.5_embryo_female 7.5_embryo_female 7.5_ED_female 7.5_TE_female 1.5_embryo_male 2.5_embryo_male 3.5_embryo_male 4.5_embryo_male 5.5_embryo_male 6.5_embryo_male 7.5_embryo_male 7.5_ED_male 7.5_TE_male"

arr=($STAGES)
echo ${arr[x]}

# define coverage threshold 
cov=5
echo $cov

# run R script

echo 'Rscript 01_dataprep_sex.R'

Rscript 01_dataprep_sex.R ${arr[x]} $cov $sex

source deactivate

echo 'JOB ENDED'
date