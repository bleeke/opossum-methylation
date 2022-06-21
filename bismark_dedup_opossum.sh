#!/bin/bash
# BSseq opossum allele-specific
#SBATCH --job-name=dedup
#SBATCH --ntasks=8
#SBATCH --time=3-00:00:0
#SBATCH --mem=60G
#SBATCH --partition=cpu
#SBATCH --array=1-54
#SBATCH --output=bismark_dedup_opossum_%A_%a.out




# BJL 
# 20190902
# for deduplication of BS-seq bams 
# from opossum allele-specific analysis
# data has been Bismark mapped single end non-dir against n-masked genome
# and run through SNPsplit programme


echo "begin"
date

# variables, files and directories

BAMDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/bams/split_bams

INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 20191016_bams_to_dedup.txt)


# load bismark

ml purge
ml Bismark
echo "modules loaded are:" 
ml




## DEDUPLICATE BAMS ##

cd $BAMDIR

deduplicate_bismark --bam $INFILE


echo "end"
date

