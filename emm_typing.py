#! /usr/bin/env python

def help_function():
	print """
	********************* Group A streptococci emm typing tool *********************
	================================================================================
	Last update:31/10/2016
	
	Three main directories:
	1. EMM_data 
	2. input
	3. output
	
	EMM_data directory
	--------------------
	
	:note created by user prior to running the tool for the first time
	
	Contains reference.seq and the EMM locus variant sequences in fasta format (*.fas)
	
	(a) reference.seq: Download reference sequence (AE004092.2; Streptococcus pyogenes M1 GAS)in fasta format
	and index the reference sequence by using makeblastdb command usually available in /usr/local/blast+/bin/
	
	<path-to>/makeblastdb -in reference.seq -dbtype nucl  -out reference
	
	(b) EMM variant sequence:   fasta file saved as emm_tsdna.fas
	Within the EMM_data directory run the following command to download the trimmed fasta sequences of the emm
	alleles from the CDC database:
	
	wget ftp://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/trimmed.tfa
	
	Then copy the edit_allele_file.sh into the EMM_data directory and run:
	
	sh edit_allele_file.sh
	
	This will rename the trimmed.tfa file to emm_tsdna_fas and make some edits to the content so that it can
	be parsed by the emmm typing tool.
	
	:note the emm_tsdna.fas should be downloaded at regular intervals to keep up with the updates in the CDC
	database
	:note the path to the EMM_data directory must be provided using the -m option.
	
	Input directory
	--------------------
	Either an input directory with two fastq files (paired reads) or the paths to the two fastq files should be
	provided (options -i for input directory and -1 sample1.R1.fastq.gz -2 sample1.R2.fastq.gz for fastq files)
	Paired reads should be indicated as sampleID.R*.fastq*
	
	Output directory 
	--------------------
	
	:note if no output dir is provided, an output dir is created in the input directory, called emm_typing
	:note if the -1 -2 option is used then an output directory must be provided.
	
	Output directory should contain EMM_log.txt, id.results.xml and stderr and stdout files in the logs directory (if not redirected using -log option)
	
	(a)EMM_log.txt details:
	
	- Variants (SNPs, indels and mixed bases) for allele sequences with >= 90% coverage (emm Validated and non-Validated alleles)
	- coverage statistics (%coverage, %identity, meanDepth, minDepth and filtered coverage (depth >= 5)
	
	:note this file is useful for troubleshooting unresolved emm type assignemnt.
	
	(b)XML output file 
	The XML file extracts the following values from the Final MLST log file:
	- Id : NGS sample id (sample identifier)
	- Version : software version number 
	- EMM validated : assigned validated variant (emm1-124)
	- EMM_Nonvalidated: assigned non validate variant (emm125+)
	For both EMM_validate and EMM_NonValidated alleles these metrics are provided:
		- percent identity - Percentage identity across the length of the allele sequence.
		- percent coverage - Percentage coverage across the length of the allele sequence (filter depth>=5 is implemented)
		- mean consensus depth - the minimum average consensus depth
		- minimum consensus depth- the minimum  consensus depth values
		- SNPs - Indicates the number of positions with a different base called compared to the reference sequence of the allele.  
		- Indels: Indicates the number of bases involved in insertion/deletion events.
		- Mixed positions: Indicates the number of positions where a base call could not be made (base frequency < 80%).

	
	tmp directory
	------------
	Should contain
	- pileup file
	- sorted bam file
	- reference file with flanking sequences (100 bps upstream and downstream)
	
	==========================================================================
	Notes
	--------------------------------------------------------------------------
	There are 2 pathways by which the script can accessed by:
	
	1. Specify EMM_data dir and glob fastq files from input-dir
	python emm_typing.py  -m profile_file_directory -i input_directory -o output_directory
	
	2. Specify the fastq_1, fastq_2 files, EMM_data dir and output_directory
	python emm_typing.py  -m profile_file_directory -1 fastq_1 -2 fastq_2 -o output_directory
	
	"""

# importing system libraries
import os, os.path, sys, subprocess, getopt, argparse, glob, yaml, inspect

# importing required modules
module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
	module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
	if module_folder not in sys.path:
		sys.path.insert(1, module_folder)

from utility_functions import *
import log_writer
import EMM_determiner_functions
from EMM_determiner_functions import create_xml_file


"""
This function checks if a file exists or not.

:param filepath: the name of the file
:param file_description: a description of the file that shown the error 

:returns: error message that the file does not exist.
"""

def check_file_exists(filepath, file_description):
	if not os.path.exists(filepath):
		print("The " + file_description + " (" + filepath + ") does not exist")
		sys.exit(1)

"""
We have set parser = argparse.ArgumentParser() and added all arguments by adding parser.add_argument.
"""
parser = argparse.ArgumentParser()
parser.add_argument('-input_directory', '-i', help='please provide an input directory')
parser.add_argument('--fastq_1', '-1', help='Fastq file pair 1')
parser.add_argument('--fastq_2', '-2', help='Fastq file pair 2')
parser.add_argument('--profile_file_directory', '-m', help='emm variants directory')
parser.add_argument('--output_directory', '-o', help='please provide an output directory')
parser.add_argument('--bowtie', '-b', help='please provide the path for bowtie2', default='bowtie2')
parser.add_argument('--samtools', '-sam', help='please provide the path for samtools', default='samtools')
parser.add_argument('--log_directory', '-log', help='please provide the path for log directory')
opts = parser.parse_args()

"""
This main function caters for two input options.

:note option 1: if user just wants to provide the path for a dir that has the fastq files then they can with just -i option.
:note option 2: if user chooses to provide forward and reverse fastq_files then they can specify them with -1 and -2 options.
"""

def main():
	fastq_files = []
	
	glob_pattern = "*.fastq*"
	ids = None
	version = "1.0"
		
	## option 1: if user just wants to provide the path for a dir that has the fastq files then they can with just -i option.
	if opts.input_directory:
		check_file_exists(opts.input_directory, 'input directory')
		if not opts.profile_file_directory:
			print("If you are providing the input directory for the fastqs, you must provide the path to the EMM reference data directory")
			check_file_exists(opts.profile_file_directory, 'profile_file_directory')
			print parser.print_help()
			sys.exit(1)
		fastq_files = glob.glob(opts.input_directory + "/" + glob_pattern)
		if fastq_files == []:
			print "No fastq files in the input directory, please provide a valid input directory"
			sys.exit(1)
		profile_file_directory = opts.profile_file_directory
		
		if not opts.output_directory:
			opts.output_directory = opts.input_directory + '/emm_typing'
			if not os.path.isdir(opts.output_directory): os.makedirs(opts.output_directory) #make output_directoryectory
		else:
			if not os.path.exists(opts.output_directory): os.makedirs(opts.output_directory)
		(seqDir,seqFileName) = os.path.split(fastq_files[0])	
		(ids,ext) = seqFileName.split('.',1)
		files = fastq_files
		
		if not opts.log_directory:
			opts.log_directory = opts.output_directory + '/logs'
			if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory)#make log_directory
		else:
			if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory)
			
   #option 2: if user chooses to provide forward and reverse fastq_files then they can specify them with -1 and -2 options.
	elif opts.fastq_1 or opts.fastq_2:
		check_file_exists(opts.fastq_1, 'Fastq 1')
		check_file_exists(opts.fastq_2, 'Fastq 2')
		if not opts.profile_file_directory:
			print("If you are using -1 and -2 options, you must provide the path to the EMM reference data directory")
			check_file_exists(opts.profile_file_directory, 'profile_file_directory')
			print parser.print_help()
			sys.exit(1)
	 
		if not opts.output_directory:
			print("If you are using -1 and -2 options, you must provide the path to the output directory")
			check_file_exists(opts.profile_file_directory, 'output_directory')
			print parser.print_help()
			sys.exit(1)

		if not os.path.isdir(opts.output_directory): os.makedirs(opts.output_directory)

		fastq_files.append(opts.fastq_1)
		fastq_files.append(opts.fastq_2)
		profile_file_directory = opts.profile_file_directory
		(SeqDir,seqFileName) = os.path.split(fastq_files[0])	
		(ids,ext) = seqFileName.split('.',1)
		files = fastq_files
		
		if not opts.log_directory:
			opts.log_directory = opts.output_directory + '/logs'
			if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directoryectory
		else:
			if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory)		
		
	
	stderr_log_output = opts.log_directory + "/" + 'emm_typing'+ ".stderr"
	stdout_log_output = opts.log_directory + "/" + 'emm_typing'+ ".stdout"
	
	logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)
	
	
	top_hits = EMM_determiner_functions.findST(fastq_files, opts.output_directory, profile_file_directory, opts.bowtie,opts.samtools,ids,opts.log_directory,version = version )
	
	#outfp = open(os.path.join(opts.output_directory,'top_hits.log'), 'w')
	#outfp.write(str(top_hits))
	#outfp.close()
	
	#if workflow_name != None:
	try_and_except(stderr_log_output,create_xml_file, top_hits, opts.output_directory,ids,version)	###################	SP
	try_and_except(stderr_log_output,write_component_complete,opts.output_directory)
	#elif workflow_name != None:
	#	xml_values = try_and_except(stderr_log_output,create_xml_file_for_failed_sample, opts.output_directory,ids,workflow_name,version)
	#	try_and_except(stderr_log_output,create_xml_file, xml_values, opts.output_directory,ids,workflow_name,version)
	#	try_and_except(stderr_log_output,write_component_complete,opts.output_directory)
	
	
try:
	main()	
except Exception:
	print "error with main function"
	sys.exit(1)	

	


