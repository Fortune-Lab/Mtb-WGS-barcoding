#! /bin/bash

#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=24G

module load python/3.10.9-fasrc01
source activate bwa_bowtie2

prefix=$1
read1=$2
read2=$3

metaphlan --bowtie2db /n/boslfs02/LABS/sfortune_lab/Lab/mchase/Databases/mpa_vJun23_CHOCOPhlAnSGB_202307/ ${read1},${read2} --bowtie2out ${prefix}.metagenome.bowtie2.bz2 --input_type fastq --nproc 12 -o ${prefix}_profiled_metagenome.txt

touch "${array[0]}.metaphlan.done"
