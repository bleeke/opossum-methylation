#!/bin/bash
# repeats line plots
#SBATCH --job-name=replineplot
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=replineplot_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 




# run R script
echo 'Rscript repeat_line_plots.R $x'

Rscript repeat_line_plots.R $x

source deactivate

echo 'JOB ENDED'
date