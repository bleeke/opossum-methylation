#!/bin/bash
# SNPsplit_genome_build
#SBATCH --job-name=SNP_genome
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH --partition=cpu
#SBATCH --output=SNPsplit_genome_build_%A_%a.out



###------###------###------###------###------###------###------###
###     script to make a mouse genome reference
###     where SNP positions between ref/alt_strain 
###     are N-masked
###------###------###------###------###------###------###------###

echo "IT HAS BEGUN"
date 

# variable and directories

BASEDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/n_mask_genome

REFGENOME=/camp/lab/turnerj/working/Bryony/rrbs/bismark/genomes/mouse/
echo "using reference genome:" $REFGENOME

VCFFILE=mgp.v5.merged.snps_all.dbSNP142.vcf.gz
echo "using vcf:" $VCFFILE

REFGENOMEBUILD=GRCm38

STRAIN=SPRET_EiJ
echo "using alt strain:" $STRAIN 

# load SNPsplit conda environment
ml purge
ml Anaconda2
source activate SNPsplit

# version info
echo "get version: SNPsplit_genome_preparation --version"
SNPsplit_genome_preparation --version

 
# move to appropriate directory
cd $BASEDIR

# run genome build command
SNPsplit_genome_preparation --nmasking --genome_build $REFGENOMEBUILD --vcf_file $VCFFILE --reference_genome $REFGENOME --strain $STRAIN


source deactivate 
echo "IT HAS FINISHED"
date





 	