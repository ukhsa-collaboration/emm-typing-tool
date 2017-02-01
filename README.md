# README for emm-typing-tool

The emm gene typing tool assigns emm type and subtype by querying the CDC M-type specific database (ftp://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/). Genomic reads are mapped to the latest version using bowtie2 (version 2.1.0; following options used: --fr --no-unal --minins 300 --maxins 1100 -k 99999 -D 20 -R 3 -N 0 -L 20 -I S,1,0.50) 16. At this stage the validated emm genes (emm1-125) are separated from the non-validated (emm125+, STC and STG) and both sets are analysed in parallel. In each set, alleles with 100% coverage (minimum depth of 5 reads per bp) over the length of their sequences and >90% identity were selected and the allele with the highest percentage identity was reported. Following this selection, a decision algorithm (see decision_algorithm.png) is implemented to determine (a) whether the validated or not validated emm type will be reported, (b) whether emm type or subtype will be reported and (c) whether further investigation is necessary (i.e. contamination, new type or presence of two emm subtypes).

## Table of content
---------------------------

* Dependencies
* Before running emm-typing-tool
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

## Before running emm-typing-tool
----------------------------

Prior to running the tool for the first time the user needs to create a reference directory. The reference directory should contain an indexed reference.seq file and a multi-fasta file containing the trimmed EMM variant sequences (180 bps) available in the CDC database.

(a) **reference.seq**: Download reference sequence (AE004092.2; Streptococcus pyogenes M1 GAS)in fasta format and index the reference sequence by using makeblastdb command usually available in /usr/local/blast+/bin/
	
	<path-to>/makeblastdb -in reference.seq -dbtype nucl  -out reference
	
(b) EMM variant sequence:   fasta file saved as **emm_tsdna.fas**
Within the EMM_data directory run the following command to download the trimmed fasta sequences of the emm alleles from the CDC database:
	
	wget ftp://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/trimmed.tfa
	
Then copy the edit_allele_file.sh into the EMM_data directory and run:
	
	sh edit_allele_file.sh
	
This will rename the trimmed.tfa file to emm_tsdna_fas and make some edits to the content so that it can be parsed by the emmm typing tool.

**NOTES**:
* The emm_tsdna.fas file should be downloaded at regular intervals to keep up with the updates in the CDC database.
* The path to the reference directory must be provided using the -m option.

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

There are 2 pathways by which the script can accessed by:
	
1. Specify EMM_data dir and glob fastq files from input-dir

	`python emm_typing.py  -m profile_file_directory -i input_directory -o output_directory`
	
2. Specify the fastq_1, fastq_2 files, EMM_data dir and output_directory
	
 `python emm_typing.py  -m profile_file_directory -1 fastq_1 -2 fastq_2 -o output_directory`
 
### Input directory

Either an input directory with two fastq files (paired reads) or the paths to the two fastq files should be provided (options -i for input directory and -1 sample1.R1.fastq.gz -2 sample1.R2.fastq.gz for fastq files). Paired reads should be indicated as sampleID.R*.fastq*

### Output directory

If no output dir is provided, an output dir is created in the input directory, called emm_typing. However, if the -1 -2 option is used then an output directory must be provided.

## emm-typing-tool Output
-------------------------

Output directory should contain EMM_log.txt, id.results.xml a logs directory with stderr and stdout files (if not redirected using -log option) and a tmp directory.

(a) **EMM_log.txt**: this file is useful for troubleshooting unresolved emm type assignemnt. It contains:
	
- Variants (SNPs, indels and mixed bases) for allele sequences with >= 90% coverage (emm Validated and non-Validated alleles)
- coverage statistics (%coverage, %identity, meanDepth, minDepth and filtered coverage (depth >= 5)
	
(b) **Result XML file**: the XML file extracts the following values from the EMM log file:
  
- Id : NGS sample id (sample identifier)
- Version : software version number 
- Final_EMM_type : Predicted EMM subtype or type based on the decision algorithm. This can be reported.
- EMM validated : assigned validated variant (emm1-124)
- EMM_Nonvalidated: assigned non validate variant (emm125+)
	
For both EMM_validated and EMM_NonValidated alleles these metrics are provided:

- percent identity - Percentage identity across the length of the allele sequence.
- percent coverage - Percentage coverage across the length of the allele sequence (filter depth>=5 is implemented)
- mean consensus depth - the minimum average consensus depth
- minimum consensus depth- the minimum  consensus depth values
- SNPs - Indicates the number of positions with a different base called compared to the reference sequence of the allele.  
- Indels: Indicates the number of bases involved in insertion/deletion events.
- Mixed positions: Indicates the number of positions where a base call could not be made (base frequency < 80%).

An example result.xml file is provided (see Example.result.xml).

(c) **tmp directory**: the tmp directory contains files useful for troubleshooting:

- pileup file
- sorted bam file
- reference file with flanking sequences (100 bps upstream and downstream)
- summary.yml only available if -v option is used. Contains the information for all alleles in the form of a dictionary; useful for more indepth troubleshooting.

## Troubleshooting: 
------------------

### Formats for the Final_EMM_type:

- *“emm1.0.sds”*: A subtype is predicted if 100% identity across the full length of the locus is observed.
- *“emm1”*: A type is predicted if >= 92% identity across the region 30-120 of the locus is observed128.
- *“emm128.0.sds/emm256.0.sds”*: Two emm non-validated types can be predicted if 100% identity is observed across the full length of both loci. This is usually seen when no validated emm type is present.
- *“Investigate new type”*: A possible new type is predicted if < 92% identity across region 30-120 of the locus is observed. 
- *"Not determined"*: No emm type/subtype can be predicted in the following scenarios:

  - two emm non-validated types with > 100% identity are called (further investigation is required and de novo sssembly might be able to provide a type).
  - if accompanied by '\*\*' in the 'EMM validated' field then the isolate is tagged for low mapping coverage at 5' end and de novo assemly will be able to provide the correct emm type.
  
- *"Mixed sample: emm3.1.sds/emm28.0.sds"*: two emm validated types with 100% identity are called (contamination). Repeat NGS if this is the case.
- *"emm9: mixed subtypes"*: two subtypes with 100% identity are called (de novo assembly is required to resolve this)

### Errors while setting up

```
Traceback (most recent call last):
  File "/Users/georgiakapatai/repositories/emm-typing-tool/modules/utility_functions.py", line 54, in try_and_except
    return function(*parameters, **named_parameters)
  File "/Users/georgiakapatai/repositories/emm-typing-tool/modules/EMM_determiner_functions.py", line 131, in prep_SRST
    process = subprocess.Popen(['seqret',seq,'-firstonly','-auto','-out',output_directory+ '/' + bait], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  File "/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/subprocess.py", line 710, in __init__
    errread, errwrite)
  File "/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/subprocess.py", line 1335, in _execute_child
    raise child_exception
OSError: [Errno 2] No such file or directory
```
**SOLUTION**: The path for emboss needs to be defined in the PATH environmental variable:

```
PATH="$PATH:/usr/local/emboss/bin"
export PATH
```
If similar error is raised for the lines where subprocess is calling bowtie2 and samtools then either provide the paths with the -b and -sam options respectively or add the paths to the PATH variable as described above.

## Contact Information
----------------------
Georgia Kapatai 

Email: Georgia.Kapatai@phe.gov.uk 

Twitter: @gkapatai

## Licence Agreement
--------------------
This software is covered by GNU General Public License, version 3 (GPL-3.0).

