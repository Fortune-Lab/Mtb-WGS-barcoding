#! /bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:20
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --array=1-14   # Replace XX with the number of accessions

module load python
eval "$(mamba shell hook --shell bash)"

#mamba activate /n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates
mamba activate mtb_isolates_cluster
SRA_LIST="amplicon_sra_list.txt"
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SRA_LIST)

if [ -z "$ACCESSION" ]; then
    echo "No accession found for index ${SLURM_ARRAY_TASK_ID}, exiting."
    exit 1
fi

COMPLETE_FLAG="${ACCESSION}.download.complete"

if [ -f "$COMPLETE_FLAG" ]; then
    echo "$COMPLETE_FLAG exists."
else
    # Download SRA file to subdirectory (default for prefetch)
    prefetch ${ACCESSION} -O ./
    
    # Find the .sra file.
    SRA_PATH=$(find . -type f -name "${ACCESSION}.sra" | head -n1)
    if [ -z "$SRA_PATH" ]; then
        echo "ERROR: SRA file for $ACCESSION not found!"
        exit 1
    fi

    # Dump fastq files to the current directory
    fastq-dump --split-files "$SRA_PATH" -O ./

    # Delete SRA file and its directory if not root dir
    SRA_DIR=$(dirname "$SRA_PATH")
    rm -f "$SRA_PATH"
    # If directory is not working dir, attempt to remove it if empty
    if [ "$SRA_DIR" != "." ]; then
        rmdir "$SRA_DIR" 2>/dev/null || true
    fi

    # success flag
    touch "$COMPLETE_FLAG"
fi

mamba deactivate
