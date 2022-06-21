#!/bin/bash
# BSseq opossum embryos
#SBATCH --job-name=bsseq_op
#SBATCH --ntasks=8
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu
#SBATCH --array=1-200
#SBATCH --output=BSseq_oposs_embryos_%A_%a.out


# BJL 
# for processing BSseq data
# re-running those that failed for lack of space



echo "begin"
date

## VARIABLES AND DIRECTORIES ## 

## create variable containing library number

LIBNUM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" to_map_please.txt)
echo "we are working on library number" "$LIBNUM"


## directories 


TRIMDIR=/camp/lab/turnerj/working/Bryony/BSseq/opossum_embryos/data/trimmed # files generated previously, see commented-out sections "merge" and "trim" below

BAMDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/bams

DEDUPDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/bams/deduplicated

GENOMEDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/genome # opossum genome, inc pseudoY and gap-filled Xchr at RSX locus

EXTRACTDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methyl_extract

TEMPDIR=/camp/lab/turnerj/scratch/bryony


## MERGE FASTQS ## 

# merge all fastqs associated with library number into one large file
# merge R1 and R2 bc mapping will be single-ended against all 4 bisu converted genomes 
# already done

#cd $INPUTDIR

#find -L $PWD/* -type f -name "LEE307A${LIBNUM}_*fastq.gz" > $BASEDIR/fastqs_to_merge_$LIBNUM.txt # make file containing names of files to merge by library number prefix

#cd $BASEDIR

#{ xargs cat < fastqs_to_merge_$LIBNUM.txt ; } > $FASTQDIR/LEE307A${LIBNUM}_merged_fastq.gz # merge 




## TRIM FASTQS ##

# trimming paramenters for PBAT-type BSseq libs
# already done

#cd $FASTQDIR
#ml purge

#ml Trim_Galore/0.4.4-foss-2016b

#trim_galore --clip_R1 6  --three_prime_clip_r1 6  --output_dir $TRIMDIR LEE307A${LIBNUM}_merged_fastq.gz 


## MAP FASTQS ## 

cd $TRIMDIR

ml purge

module use -a /camp/apps/eb/dev/modules/all

ml Bismark

echo "using Bismark version" 
which Bismark

bismark --non_directional --un --ambiguous --temp_dir $TEMPDIR -output_dir $BAMDIR --genome $GENOMEDIR LEE307A${LIBNUM}_*.fq.gz   


## DEDUPLICATE BAMS ##

cd $BAMDIR

deduplicate_bismark --output_dir $DEDUPDIR --bam LEE307A${LIBNUM}_*.bam



echo "end"
date