#!/bin/bash
# extract and make cytosine report
#SBATCH --job-name=extract
#SBATCH --ntasks=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu
#SBATCH --array=1-18
#SBATCH --output=bismark_extract_%A_%a.out


# BJL 
# for processing BSseq data




echo "begin"
date

## VARIABLES AND DIRECTORIES ## 

## create variable containing library number

LIBNUM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" adult_libraries.txt)
echo "we are working on library number" "$LIBNUM"


## directories 

BAMDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/bams/split_bams

GENOMEDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/n-masked_genome # opossum genome, inc pseudoY and gap-filled Xchr at RSX locus with n-masked parental SNPs

EXTRACTDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methyl_extract


# load Bismark 

ml purge

module use -a /camp/apps/eb/dev/modules/all

ml Bismark


## METHYLATION EXTRACTION ##

cd $BAMDIR  

bismark_methylation_extractor --comprehensive --multicore 2 --bedGraph --cytosine_report --genome_folder $GENOMEDIR --o $EXTRACTDIR ${LIBNUM}_*allele_flagged.deduplicated_sorted.bam

echo "end"
date