#! /bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-01:30
#SBATCH -p shared
#SBATCH --mem=16000

#https://stackoverflow.com/questions/51139711/hpc-cluster-select-the-number-of-cpus-and-threads-in-slurm-sbatch

module load python/3.10.9-fasrc01

export JAVA_HOME=/n/boslfs02/LABS/sfortune_lab/Lab/software/jdk1.8.0_202/bin
export PATH=$JAVA_HOME:$PATH
export PATH=/n/boslfs02/LABS/sfortune_lab/Lab/software/bin:$PATH
export PATH=/n/boslfs02/LABS/sfortune_lab/Lab/software/bwa/:$PATH

prefix=$1
read1=$2
read2=$3

ref=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/Databases/ErdmanPolished.circ.fasta

aln=bwa
lib=ErdmanNanopore

outbam=${prefix}.${aln}.${lib}.bam
bamindex=${prefix}.${aln}.${lib}.Sd.bai
outbamSortDedup=${prefix}.${aln}.${lib}.Sd.bam
outbamRemovedup=${prefix}.${aln}.${lib}.Rd.bam
nameSort=${prefix}.${aln}.${lib}.nSort.bam
Rgroup="@RG\tID:$prefix\tSM:$prefix"
outbamSort=${prefix}.S.${aln}.${lib}.bam


set -eo pipefail

source activate fastp_qc

if [ -f "${prefix}.trimming.complete" ]; then
        echo "${prefix}.trimming.complete exists."
else

        fastp --thread $SLURM_CPUS_PER_TASK \
              --length_required 50 \
              -i ${read1} -I ${read2} --html ${prefix}.fastp.html --json ${prefix}.fastp.json \
              -o ./${prefix}.R1.trim.fastq.gz -O ./${prefix}.R2.trim.fastq.gz --failed_out ./${prefix}.fail.fastq.gz
      touch "${prefix}.trimming.complete"
fi

if [[ -f "./${prefix}.downsampling.complete" ]]; then
        echo "${prefix}.downsampling.complete exists."
else
        rasusa reads --coverage 100 --genome-size 4.4mb -o ./${prefix}.ds.R1.trim.fastq.gz -o ./${prefix}.ds.R2.trim.fastq.gz \
                ./${prefix}.R1.trim.fastq.gz ./${prefix}.R2.trim.fastq.gz
        touch "./${prefix}.downsampling.complete"
fi



mamba deactivate


if [ -f "${prefix}.bwa.mapping.done" ]; then
        echo "${prefix}.bwa.mapping.done exists."
else
        bwa mem $ref ./${prefix}.R1.trim.fastq.gz ./${prefix}.R2.trim.fastq.gz -t $SLURM_CPUS_PER_TASK -R $Rgroup -M  | samtools view -h -S -b -o $outbam -
        samtools sort -T $prefix -o $outbamSort -O bam -@ $SLURM_CPUS_PER_TASK $outbam
        samtools index $outbamSort
        touch "${prefix}.bwa.mapping.done"
fi

if [ -f "${prefix}.markDuplicates.done" ]; then
        echo "${prefix}.markDuplicate.done exists."
else
        java -Xmx12g -jar /n/boslfs02/LABS/sfortune_lab/Lab/mchase/software/picard.jar MarkDuplicates \
                INPUT=$outbamSort OUTPUT=$outbamSortDedup \
                REMOVE_DUPLICATES=false \
                METRICS_FILE=$prefix.${aln}.${lib}.dedupMets.txt \
                VALIDATION_STRINGENCY=LENIENT TMP_DIR=./
        	samtools index $outbamSortDedup $bamindex
        	touch "${prefix}.markDuplicates.done"
fi

if [ -f "${prefix}.deduplicated_fastqs.done" ]; then
        echo "${prefix}.deduplicated_fastqs.done exists."
else
	samtools view -F 1024 -h -b $outbamSortDedup -o $outbamRemovedup
	samtools index $outbamRemovedup
	samtools sort -n -@ $SLURM_CPUS_PER_TASK $outbamRemovedup -o $nameSort
	samtools fastq -1 ${prefix}.R1.RD.fastq.gz -2 ${prefix}.R2.RD.fastq.gz -s ${prefix}.R0.RD.fastq.gz $nameSort
	touch "${prefix}.deduplicated_fastqs.done"
fi

if [ -f "${prefix}.coverage.done" ]; then
	echo "${prefix} completed."
else
	samtools sort -T $prefix -o ${prefix}.dedup.sort.bam -O bam -@ $SLURM_CPUS_PER_TASK $nameSort
	samtools coverage ${prefix}.dedup.sort.bam > ${prefix}.samtools.coverage
	touch "${prefix}.coverage.done"

fi

