#!/bin/bash
# SNPsplit
#SBATCH --job-name=SNPsplit
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --partition=cpu
#SBATCH --array=1-12
#SBATCH --output=SNPsplit_%A_%a.out


####------------------------------------------------------------------####
## script for running SNPsplit on BS-seq data
## the data was mapped with Bismark
####------------------------------------------------------------------####


echo "IT HAS BEGUN"
date



## create variable containing library number
LIBNUM=${SLURM_ARRAY_TASK_ID}
echo "we are working on library number" "$LIBNUM"



# directories 

BAMDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/bams # where the mapped files to be split are

SNPDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/n_mask_genome # where the list of SNPs is

OUTDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/bams/split_bams




# run SNPsplit

ml purge
ml Anaconda2
source activate SNPsplit
echo "loaded modules are:"
ml

echo "running SNPsplit"

SNPsplit --bisulfite --conflicting -o $OUTDIR --snp_file $SNPDIR/all_SNPs_SPRET_EiJ_GRCm38.txt.gz $BAMDIR/LEE73A${LIBNUM}_*.bam



echo "IT HAS FINISHED"
date