# README for emm-typing-tool

The emm gene typing tool assigns emm type and subtype by querying the CDC M-type specific database (ftp://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/). Genomic reads are mapped to the latest version using bowtie2 (version 2.1.0; following options used: --fr --no-unal --minins 300 --maxins 1100 -k 99999 -D 20 -R 3 -N 0 -L 20 -I S,1,0.50) 16. At this stage the validated emm genes (emm1-125) are separated from the non-validated (emm125+, STC and STG) and both sets are analysed in parallel. In each set, alleles with 100% coverage (minimum depth of 5 reads per bp) over the length of their sequences and >90% identity were selected and the allele with the highest percentage identity was reported. Following this selection, a decision algorithm (see decision_algorithm.png) is implemented to determine (a) whether the validated or not validated emm type will be reported, (b) whether emm type or subtype will be reported and (c) whether further investigation is necessary (i.e. contamination, new type or presence of two emm subtypes).

## Table of content
---------------------------

* Dependencies
* Downloading reference database
* Running emm-typing-tool
* emm-typing-tool Output
* Troubleshooting
* Contact Information
* Licence Agreement

## Dependencies
---------------------------

emm-typing-tool  is written with Python 2.7.5 and requires the following packages installed before running:
* Bowtie2 (https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
* Samtools (https://sourceforge.net/projects/samtools/files/samtools/0.1.19/) - **version specific**. To install after downloading:

```
cd <path-to>/samtools-0.1.19
make
```

* PyYaml (http://pyyaml.org/)
* numpy (http://www.scipy.org/scipylib/download.html)
* lxml  (http://lxml.de/installation.html)
* biopython (https://github.com/biopython/biopython)
* emboss (http://emboss.open-bio.org/html/adm/ch01s01.html)
* blast+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## Running emm-typing-tool
----------------------------
```
usage: emm_typing.py [-h] [-input_directory INPUT_DIRECTORY]
                     [--fastq_1 FASTQ_1] [--fastq_2 FASTQ_2]
                     [--profile_file_directory PROFILE_FILE_DIRECTORY]
                     [--output_directory OUTPUT_DIRECTORY] [--bowtie BOWTIE]
                     [--samtools SAMTOOLS] [--log_directory LOG_DIRECTORY]
                     [--verbose]

optional arguments:
  -h, --help            show this help message and exit
  -input_directory INPUT_DIRECTORY, -i INPUT_DIRECTORY
                        please provide an input directory
  --fastq_1 FASTQ_1, -1 FASTQ_1
                        Fastq file pair 1
  --fastq_2 FASTQ_2, -2 FASTQ_2
                        Fastq file pair 2
  --profile_file_directory PROFILE_FILE_DIRECTORY, -m PROFILE_FILE_DIRECTORY
                        emm variants directory
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        please provide an output directory
  --bowtie BOWTIE, -b BOWTIE
                        please provide the path for bowtie2
  --samtools SAMTOOLS, -sam SAMTOOLS
                        please provide the path for samtools
  --log_directory LOG_DIRECTORY, -log LOG_DIRECTORY
                        please provide the path for log directory
  --verbose, -v         if selected a summary.yml file is generated in the tmp
                        directory
```

