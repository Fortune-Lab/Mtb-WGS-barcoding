



# Project Title
Data associated with the manuscript "**Project title**" by **Author** et al. (**Year, manuscript DOI if available**)

## Project Description
This project describes the computational workflow used to identify barcoded Mycobacterium tuberculosis (Mtb) in tissue extracts from animals infected with a library of barcoded strains. More specifically, the workflow is applied to FASTQ files generated as follows: Animals are infected with a genetically barcoded strain of Mycobacterium tuberculosis (Mtb), strain Erdman described in Martin CJ, Cadena AM, Leung VW, Lin PL, Maiello P, Hicks N, Chase MR, Flynn JL, Fortune SM. Digitally Barcoding Mycobacterium tuberculosis Reveals In Vivo Infection Dynamics in the Macaque Model of Tuberculosis. mBio. 2017 May 9;8(3):e00312-17. doi: 10.1128/mBio.00312-17. PMID: 28487426; PMCID: PMC5424202. Bacteria is cultured out of harvested tissues, followed by bacterial DNA extraction, and Illumina WGS sequencing.

### Directories:
- ./Scripts ##Subfolder with scripts
  - Remove_duplicates.sh 
  - SumBarcodes.plx 
  - SumBarcodes18.plx 
  - WGSBarcodeCounter18mer.plx  
  - WGSBarcodeCounter.plx
  - barcodes.sh 
  - barcodes18.sh 
  - runMetaphlan.sh 
  - run_tb-profiler
  - submit_Reads_to_Metaphlan.sh
  - submit_Reads_to_tb_profiler
  - submit_reads_to_deduplicate.sh
  - transpose.plx
- ./Fig3 ##Subfolder with scripts and data associated with Figure 3
  - File name.R ##Use a descriptive file name and provide a short 1-2 sentence description here
  </br></br> **NOTE:** Each file in this directory is a small batch script generated during analysis. Not all analyses were included in the final manuscript. </br></br>
- ./Supplemental_figures ##Subfolder with scripts and data associated with Supplemental_figures
  - FigS1.R ##Use a descriptive file name and provide a short 1-2 sentence description here
  - FigS1.py ##Use a descriptive file name and provide a short 1-2 sentence description here
  - FigS4.txt ##Use a descriptive file name and provide a short 1-2 sentence description here
  - FigS5.nw ##Use a descriptive file name and provide a short 1-2 sentence description here.

### Prerequisites:

    Requires Perl, seqtk https://github.com/lh3/seqtk, fastp, gatk, tb-profile and metaphlan4
    To download the fastq files from SRA you will also need to install sra-tools: https://github.com/ncbi/sra-tools. To speed things up a little gnu parallel: https://www.gnu.org/software/parallel/
### Library analysis
   Extract, count barcodes and remove singletons.
   perl BarcodeReader0627.plx fastq_file | awk '{print $3}' | sort | uniq -c | awk '$1 > 1 {print $2 "\t" $1}' | sort -k2 -nr > LIB065162_primary_reads_sort_count_remove_singletons.tsv
   Find inflection point and output graph.
   python sliding_window.py
### Barcoding 18mer

### Barcoding 7mer

### Sample quality control
Tb-profiler
Metaphlan

### Citations:
  - Skip if this information is already in a methods section or in the script.





