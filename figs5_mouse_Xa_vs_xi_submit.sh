#!/bin/bash
# mouse Xa vs Xi ridgeplot
#SBATCH --job-name=mouseridgeplot
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90 
#SBATCH --output=mouseridgeplot_%A_%a_submit.out

echo 'JOB STARTED'
date

ml purge
ml Anaconda2
source activate bRy_3_6_0

# this will run in my bRy conda R env 




# run R script
echo 'figs5_mouse_Xa_vs_xi.R'

Rscript  figs5_mouse_Xa_vs_xi.R

source deactivate

echo 'JOB ENDED'
date