#!/bin/bash
# Bismark_genome_prep
#SBATCH --job-name=bisgenprep
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH --partition=cpu
#SBATCH --output=bismark_genome_prep_%A_%a.out

echo "IT HAS BEGUN"
date


###------###------###------###------###------###------###------###
###     script to prepare a bismark genome index
###------###------###------###------###------###------###------###


# directories

NMASKGENOME=/camp/lab/turnerj/working/shared_projects/OOPs/genome/n-masked



# clean and load modules, list them for posterity
ml purge
ml Bismark/0.22.1-intel-2018b
echo "modules loaded are:" ml

# perform bismark genome indexing for use with bowtie2
bismark_genome_preparation --bowtie2 $NMASKGENOME


echo "IT HAS FINISHED"
date
