#!/bin/bash
# HISAT2 map
#SBATCH --job-name=HISATmap
#SBATCH --ntasks=4
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --partition=cpu
#SBATCH --array=1-57
#SBATCH --output=hisat2_map_%A_%a.out


####------------------------------------------------------------------####
## script for mapping RNA-seq libraries with HISAT2 
####------------------------------------------------------------------####


echo "start"
date

#directories

TRIMDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/rna-seq/data/newmap/trimmed

BAMDIR=/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/rna-seq/bams

GENOMEDIR=/camp/lab/turnerj/working/shared_projects/OOPs/genome/n-masked/mondom5_pseudoY_X-gaps-filled_220819_JZ_masked # modified and n-masked genome made by jasmin 20190927

# read in files

cd $TRIMDIR

FILE1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" trimmed.txt)

FILE2="${FILE1/R1*val_1.fq.gz/}R2*val_2.fq.gz"

OUTFILE="${FILE1/R1*val_1.fq.gz/}.sam"

echo "we will map" $FILE1 $FILE2 

# load modules and map data

ml purge
ml HISAT2
echo "modules loaded are:" 
ml

cd $BAMDIR

hisat2 --no-softclip --no-mixed --no-discordant -x $GENOMEDIR -1 $FILE1 -2 $FILE2 -S $OUTFILE


date 
echo "end"