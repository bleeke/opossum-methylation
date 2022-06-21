#!/bin/bash
# BSseq opossum allele-specific
#SBATCH --job-name=extract
#SBATCH --ntasks=8
#SBATCH --time=3-00:00:0
#SBATCH --mem=60G
#SBATCH --partition=cpu
#SBATCH --array=1-54
#SBATCH --output=bismark_extract_opossum_%A_%a.out



# BJL 
# 20190902
# for methylation extraction of BS-seq bams 
# from opossum allele-specific analysis
# data has been Bismark mapped single end non-dir against n-masked genome
# and run through SNPsplit programme
# then deduplicated by Bismark


echo "begin"
date

# variables, files and directories

BAMDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/bams/split_bams

EXTRACTDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/methyl_extract


INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 20191016_bams_to_extract.txt)



## METHYLATION EXTRACTION ##  

ml purge
ml Bismark
echo "module loaded are:"
ml

cd $BAMDIR

bismark_methylation_extractor --bedGraph --o $EXTRACTDIR $INFILE

echo "end"
date
