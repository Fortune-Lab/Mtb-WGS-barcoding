#! /bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 0-01:30
#SBATCH -p shared
#SBATCH --mem=10000

module load python/3.10.9-fasrc01

set -eo pipefail

source activate tb-profiler

sample=$1
read1=$2
read2=$3

tmp_dir=/tmp/tmp.${SLURM_JOB_ID}
base_dir=`pwd`
mkdir -p ${tmp_dir}/database
echo "${tmp_dir}"
db=$(tb-profiler list_db | cut -f5); cp $db* ${tmp_dir}/database
cd ${tmp_dir} && tb-profiler profile -1 ${read1} -2 ${read2} -p ${sample} --external_db database/tbdb --no_delly --threads 4 --csv
mv ${tmp_dir}/results/${sample}.* ${base_dir}
cd ${base_dir}
touch "${sample}.tbprofiler.done"

mamba deactivate
