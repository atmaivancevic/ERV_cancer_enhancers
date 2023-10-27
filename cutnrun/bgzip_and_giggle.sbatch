#!/bin/bash

## Script to search genomic intervals using giggle
##
## Example usage:
## inDir=. outDir=. gigIdx=/Shares/CL_Shared/db/giggle/hg38/cistrome/Human_Factor/indexed  db=cistrome_human_factor sbatch --array 0-3 bgzip_and_giggle.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J giggle
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# important note: make sure the giggle executable is in your PATH

# load modules
module load samtools/1.10

# define query files
queries=($(ls $inDir/*_mergedpeaks.bed | xargs -n 1 basename))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Bgzipping bed file...")

bgzip ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Starting giggle...")

giggle search -i ${gigIdx} -q ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}.gz -s -g 3209286105 \
| sed 's#sorted/##g' \
| sed 's/.gz//g' \
| grep -v "#" \
| sort -nrk8 \
| sed '1i#file\tfile_size\toverlaps\todds_ratio\tfishers_two_tail\tfishers_left_tail\tfishers_right_tail\tcombo_score' \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bed.gz}_VERSUS_${db}.tab

echo $(date +"[%b %d %H:%M:%S] Unzipping original bed file...")

gunzip ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")
