#!/bin/bash
# BSseq mouse allele-specific
#SBATCH --job-name=extract
#SBATCH --ntasks=8
#SBATCH --time=3-00:00:0
#SBATCH --mem=60G
#SBATCH --partition=cpu
#SBATCH --array=1-36
#SBATCH --output=bismark_extract_mouse_%A_%a.out



# BJL 
# 20190902
# for methylation extraction of BS-seq bams 
# from mouse allele-specific analysis
# data has been Bismark mapped single end non-dir against n-masked genome
# and run through SNPsplit programme
# then deduplicated by Bismark


echo "begin"
date

# variables, files and directories

BAMDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/bams/split_bams

EXTRACTDIR=/camp/lab/turnerj/scratch/bryony/20190906_methyl_extract
mkdir -p $EXTRACTDIR

INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 20190903_bams_to_extract.txt)



## METHYLATION EXTRACTION ##  

ml purge
ml Bismark
echo "module loaded are:"
ml

cd $BAMDIR

bismark_methylation_extractor --bedGraph --o $EXTRACTDIR $INFILE

echo "end"
date