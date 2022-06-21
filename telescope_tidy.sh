#!/bin/bash
# data tidy telescope output
#SBATCH --job-name=tidyscope
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=10G
#SBATCH --partition=cpu
#SBATCH --array=1-707
#SBATCH --output=tidy_telescope_output_%A_%a.out

# this script takes the telescope_report.tsv output from telescope
# tidies it and calculate cpm for each repeat
# writes a file with the format unique_repeat_locus cpm 



echo "IT HAS BEGUN"
date

INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" telescope_reports.txt)
BASEDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/repeat_expression

echo "start is" $INFILE


# load a
n R session in bRy environment
ml purge
ml Anaconda2
source activate bRy_3_6_0



# run featureCounts  R script
Rscript 20210822_telescope_tidy.R $INFILE $BASEDIR


echo "IT IS FINISHED"
date