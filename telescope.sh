#!/bin/bash
# repeat quant
#SBATCH --job-name=telescope
#SBATCH --ntasks=3
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --partition=cpu
#SBATCH --array=1-707
#SBATCH --output=oposs_telescope_%A_%a.out


# using chrX shifted repeat gtf 

echo 'start'; date


# directories and variables

INPUTDIR=/camp/lab/turnerj/working/Bryony/public_data/SKM_2020/bams

OUTPUTDIR=/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/repeat_expression

GTFFILE=/camp/lab/turnerj/working/Bryony/annotations/20210123_UCSC_monDom5_RepeatMasker_sorted_shifted_locus.gtf

BAMFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SKM_bams.txt)



cd $INPUTDIR 

ml Anaconda2

source activate telescope

telescope assign $BAMFILE $GTFFILE --exp_tag $BAMFILE --outdir $OUTPUTDIR


source deactivate 


echo 'end'; date