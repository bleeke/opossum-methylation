#!/bin/bash
# Bismark map
#SBATCH --job-name=Bismarkmap
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --partition=cpu
#SBATCH --array=1-12
#SBATCH --output=bismark_map_%A_%a.out


####------------------------------------------------------------------####
## script for trimming mapping BS-seq libraries with Bismark 
####------------------------------------------------------------------####


echo "IT HAS BEGUN"
date


## create variable containing library number
LIBNUM=${SLURM_ARRAY_TASK_ID}
echo "we are working on library number" "$LIBNUM"


# directories 

INPUTDIR=/camp/lab/turnerj/inputs/babs/bryony.leeke/robert.goldstone/PM18191

BASEDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq

FASTQDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/raw

TRIMDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/trimmed

MAPDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/bams

GENOMEDIR=/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/n_mask_genome/SPRET_EiJ_N-masked





####--------------------------------------####
## 				MERGE FASTQS 
####--------------------------------------####

#echo "MERGING FASTQS"

# merge all fastqs associated with library into one file


#cd $INPUTDIR

#find -L $PWD/* -type f -name "LEE73A${LIBNUM}_*fastq.gz" > $BASEDIR/fastqs_to_merge_$LIBNUM.txt # make file containing names of files to merge by library number prefix

#cd $BASEDIR

#{ xargs cat < fastqs_to_merge_$LIBNUM.txt ; } > $FASTQDIR/LEE73A${LIBNUM}_merged_fastq.gz # merge 




####--------------------------------------####
## 				RUN TRIMGALORE
####--------------------------------------####



#echo "running TRIMGALORE"

# clean modules and load TrimGalore!, list for posterity

#ml purge
#ml Trim_Galore/0.4.4-foss-2016b
#echo "modules loaded are:" 
#ml

#cd $FASTQDIR

#trim_galore --clip_R1 6  --three_prime_clip_r1 6  --output_dir $TRIMDIR LEE73A${LIBNUM}_merged_fastq.gz 




####--------------------------------------####
## 				RUN BISMARK
####--------------------------------------####


cd $TRIMDIR

ml purge
ml Bismark/0.22.1-intel-2018b
echo "modules loaded are:"
ml

bismark --non_directional --un --ambiguous --output_dir $MAPDIR --genome $GENOMEDIR LEE73A${LIBNUM}_*.fq.gz   