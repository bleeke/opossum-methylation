#!/bin/bash
# extract and make cytosine report
#SBATCH --job-name=extract
#SBATCH --ntasks=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu
#SBATCH --array=1-108
#SBATCH --output=bismark_extract_%A_%a.out


# BJL 
# for processing BSseq data




echo "begin"
date

## VARIABLES AND DIRECTORIES ## 

## create variable containing library number

LIBNUM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" libraries.txt)
echo "we are working on library number" "$LIBNUM"


## directories 

TEMPDIR=/camp/lab/turnerj/scratch/bryony/deduplicated

GENOMEDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/genome # opossum genome, inc pseudoY and gap-filled Xchr at RSX locus

EXTRACTDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methyl_extract


# load Bismark 

ml purge

module use -a /camp/apps/eb/dev/modules/all

ml Bismark


## METHYLATION EXTRACTION ##

cd $TEMPDIR  

bismark_methylation_extractor --comprehensive --multicore 2 --bedGraph --cytosine_report --genome_folder $GENOMEDIR --o $EXTRACTDIR ${LIBNUM}_*deduplicated.bam

echo "end"
date