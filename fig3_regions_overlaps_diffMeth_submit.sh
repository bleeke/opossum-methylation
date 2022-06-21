#!/bin/bash
# regions diff meth
#SBATCH --job-name=regionsDiffMeth
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --partition=cpu
#SBATCH --array=1-5
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=regionsDiffMeth_overlaps_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 

# define command line arguments provided to R script

# define feature type to be worked on
x=${SLURM_ARRAY_TASK_ID}
echo $x


# run R script
echo 'Rscript fig2_regions_overlaps_diffMeth.R $x'

Rscript 20210810_fig2_regions_overlaps_diffMeth.R $x

source deactivate

echo 'JOB ENDED'
date