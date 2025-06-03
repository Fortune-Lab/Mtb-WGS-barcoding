#! /bin/bash

#SBATCH -c 2
#SBATCH -N 1
#SBATCH -t 0-03:00
#SBATCH -p sapphire
#SBATCH --mem=8000

#module load seqtk/1.2-fasrc01
#module load parallel/20180522-fasrc01

export PATH=/n/boslfs02/LABS/sfortune_lab/Lab/software/bin:$PATH

now=$(date +"%m_%d_%Y")

out=barcode_counts_${now}.txt

cat $1 | parallel --gnu -j 2 perl WGSBarcodeCounter18merIndices.plx {} > $out

#perl SumBarcodes.plx $out

#perl RemoveSequencingErrors.plx
