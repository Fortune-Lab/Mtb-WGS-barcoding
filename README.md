



# Project Title
Data associated with the manuscript "**Project title**" by **Author** et al. (**Year, manuscript DOI if available**)

## Project Description
This project describes the computational workflow used to identify barcoded Mycobacterium tuberculosis (Mtb) in tissue extracts from animals infected with a library of barcoded strains. More specifically, the workflow is applied to FASTQ files generated as follows: Animals are infected with a genetically barcoded strain of Mycobacterium tuberculosis (Mtb), strain Erdman described in Martin CJ, Cadena AM, Leung VW, Lin PL, Maiello P, Hicks N, Chase MR, Flynn JL, Fortune SM. Digitally Barcoding Mycobacterium tuberculosis Reveals In Vivo Infection Dynamics in the Macaque Model of Tuberculosis. mBio. 2017 May 9;8(3):e00312-17. doi: 10.1128/mBio.00312-17. PMID: 28487426; PMCID: PMC5424202. Bacteria is cultured out of harvested tissues, followed by bacterial DNA extraction, and Illumina WGS sequencing.

### Directories:
- ./BarcodeLibraryAnalysisScripts ##Subfolder with scripts
  - BarcodeReaderV.plx
  - sliding_window.py
  - twoPlot.py
    
- ./BarcodingAmplicon18mer
  - BarcodeReader18.plx
  - FindThreshold.plx
    
- ./BarcodingWGS_7mer
  - Remove_duplicates.sh
  - submit_reads_to_deduplicate.sh
  - WGSBarcodeCounter.plx
  - SumBarcodes.plx
  - RemoveSequencingErrors.plx    
  - barcodes.sh

  - ./Additional_scripts
  - runMetaphlan.sh 
  - run_tb-profiler
  - submit_Reads_to_Metaphlan.sh
  - submit_Reads_to_tb_profiler
  - transpose.plx

### Prerequisites:
    Requires Perl, seqtk https://github.com/lh3/seqtk, fastp, gatk, tb-profile and metaphlan4
    To download the fastq files from SRA you will also need to install sra-tools: https://github.com/ncbi/sra-tools. To speed things up a little gnu parallel: https://www.gnu.org/software/parallel/

### Data:
   [bioproject]: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1290335/
   [PRJNA1290335][bioproject]


### Library analysis
  - Extract, count barcodes and remove singletons. <br>
   ```
   `perl BarcodeReader0627.plx fastq_file | awk '{print $3}' | sort | uniq -c | awk '$1 > 1 {print $2 "\t" $1}' | sort -k2 -nr > LIB065162_primary_reads_sort_count_remove_singletons.tsv`
   ```
   - Find inflection point and output graph. <br>
     
   `python sliding_window.py`
### Sample preprocessing for WGS Barcoding
  - Run fastp
  - Remove duplicates
  - Extract reads <br>
    `Remove_duplicates.sh` <br>
  - Outputs fastq files that are ready to run in downstream scripts.
  
### Amplicon Barcoding 18mer
  - Extract and process 18 base barcodes.
      - BarcodeCounter18.plx
      - FindThreshold18.plx <br>
      
### WGS Barcoding 7mer
  - Extract and process 7 base barcodes.
      - WGSBarcodeCounter.plx
      - SumBarcodes.plx <br>
  `sbatch barcodes.sh` <br>
  `perl RemoveSequencingErrors.plx` <br>
   - Outputs file of barcodes and L5 integration data for each sample.

### Sample quality control
  - Run Tb-profiler
  `sbatch submit_Reads_to_Metaphlan.sh` <br>
  - Run Metaphlan
`submit_Reads_to_Metaphlan.sh`
### Final curation





