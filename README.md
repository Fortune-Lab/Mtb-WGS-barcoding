



# Project Title
Data associated with the manuscript "**Project title**" by **Author** et al. (**Year, manuscript DOI if available**)

## Project Description
Immunological memory elicited either through previous or ongoing M. tuberculosis (Mtb) infection provides a critical mechanism by which hosts protect against re-infection and disease progression upon Mtb re-exposure. Conversely, the uneven competition between distinct Mtb strains suggest certain bacterial clades have enhanced ability to spread across communities and circulate globally, potentially by evading memory responses gained by prior infection with genomically different strains. To address whether memory responses induced by one strain can protect against a genetically distinct strain, we conducted a heterologous reinfection study in cynomolgus macaques involving primary infection by a Lineage 4 Erdman Mtb strain and subsequent re-challenge by a Lineage 2 strain, HT-L2. Recent epidemiologic studies have shown that the clade to which HT-L2 belongs has been spreading successfully over the last decade in Lima, Peru. Here, through microbiologic, PET-CT imaging and sequencing of Mtb genomic barcodes, we show that reinfected animals developed fewer lung lesions and controlled both pulmonary and disseminated forms of infection better than naïve animals without prior exposure to Mtb. Our data support that protection against reinfection is not limited by Mtb lineage, providing optimism that vaccines can be effective across populations and geographic locations.

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





