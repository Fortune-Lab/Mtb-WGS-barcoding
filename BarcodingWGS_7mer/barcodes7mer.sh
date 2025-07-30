#! /bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-03:00
#SBATCH -p sapphire
#SBATCH --mem=8000

#require seqtk and gnu parallel
#module load seqtk/1.2-fasrc01
#module load parallel/20180522-fasrc01

export PATH=/n/boslfs02/LABS/sfortune_lab/Lab/software/bin:$PATH

now=$(date +"%m_%d_%Y")

out=barcode_counts_${now}.txt

cat $1 | parallel --gnu -j 24 --colsep '\t' perl  WGSBarcodeCounter7mer.plx {1} {2} {3} > $out

perl SumBarcodes7mer.plx $out

#perl RemoveSequencingErrors.plx
