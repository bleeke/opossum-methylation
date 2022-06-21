#!/bin/bash
# featureCounts opossum adult rna-seq JvDB 
#SBATCH --job-name=opossfeatureCounts
#SBATCH --ntasks=8
#SBATCH --time=72:00:00
#SBATCH --mem=80G
#SBATCH --partition=cpu
#SBATCH --array=1
#SBATCH --output=20191007_opossum_rnaseq_featureCounts_%A_%a.out


# this script launches an R script
# takes input sam files and uses featureCounts from the RsubRead package to do read summarization
# here it is applied to rna-seq data from opossum adult JvDB ross


echo "IT HAS BEGUN"
date

# define options 
INDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/rna-seq/bams/split_bams
ANNOT=/camp/lab/turnerj/working/shared_projects/OOPs/genome/Monodelphis_domestica.monDom5.97_mod_shifted_XY.gtf
OUTMATRIX=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/rna-seq/counts_matrix.txt
OUTFILE=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/rna-seq/myFeatureCounts.Rdata


# load an R session in bRy environment
ml purge
ml Anaconda2
source activate bRy



# run featureCounts  R script
Rscript 20191007_featureCounts.R $INDIR $ANNOT $OUTMATRIX $OUTFILE


echo "IT IS FINISHED"
date