#!/bin/bash
# HISAT2_genome_builder
#SBATCH --job-name=hisgenbuild
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH --partition=cpu
#SBATCH --output=hisat2_genome_builder_%A_%a.out

echo "IT HAS BEGUN"
date


###------###------###------###------###------###------###------###
###     script to prepare a HISAT2 genome index
###------###------###------###------###------###------###------###


# directories

NMASKGENOME=/camp/lab/turnerj/working/shared_projects/OOPs/genome/n-masked



# clean and load modules, list them for posterity
ml purge
ml HISAT2
echo "modules loaded are:" ml


# perform HISAT2 genome indexing 
cd $NMASKGENOME
hisat2-build mondom5_pseudoY_X-gaps-filled_220819_JZ_masked.fasta mondom5_pseudoY_X-gaps-filled_220819_JZ_masked    


echo "IT HAS FINISHED"
date
