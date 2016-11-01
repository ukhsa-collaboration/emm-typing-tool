# import shutil # not used
from collections import Counter
import string, re
import os, os.path, sys, subprocess, inspect
import pickle
import yaml
#import getopt # not used
from optparse import OptionParser
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation # not used
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import Bio.Seq
# import Bio.SeqIO # already imported SeqIO in line 10
# import Bio.SeqRecord # already imported SeqRecord in line 12
import fileinput
from lxml import etree
from time import clock, time
import glob
from math import *
import traceback
from itertools import groupby
from operator import itemgetter
import numpy



module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
    module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
    if module_folder not in sys.path:
        sys.path.insert(1, module_folder)

import log_writer
from utility_functions import *
"""
Function
Main function -> Calls prepare, get_profiles, get_filesets, align_and_get_scores and write_scores functions
    
The option of the method
files[str]:  The path to where the fastq file located
output_directory[str]: The path to where output file located
profile_file_directory[str]:  The path to where reference.seq, and the EMM variant sequences (*.fas) files are located
bowtie[str]: The path to Bowtie2 command
samtools[str]: The path to SAMtools command
ids[str]: Unique identifier number
version[str]: version number


"""
def findST(files, output_directory, profile_file_directory, bowtie, samtools, ids, log_directory, version = ""):
    
    # start run time
    start = clock()
    
    #delete tmp file
    if not os.path.exists(output_directory + '/tmp'):
        #make tmp directory in output_directory
        os.makedirs(output_directory + '/tmp')
    
    
    workingDir  = output_directory + '/tmp'
    
    #create EMM_log file
    log = open(output_directory + "/EMM_log.txt", "w")
    

    #set stderr.log and stdout.log files.
    #stderr_log_output from try_and_except function and logger function are appended into  ids+ ".stderr.log file within output dir
    
    stderr_log_output = log_directory + "/" + 'emm_typing'+ ".stderr"
    stdout_log_output = log_directory + "/" + 'emm_typing'+ ".stdout"

    logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)
    
    
    #Extract flanking regions of 100bp upstream and downstream of each MLST locus by blast against a reference genome
    try_and_except(stderr_log_output, prep_SRST,profile_file_directory, output_directory, logger)
    
    #Concatenate flanking regions extracted by prep SRST function to correspondent locus variants sequence in fasta format.
    #Newly concatenated sequence are then indexed by Bowtie2
    lociHeader = try_and_except(stderr_log_output, prepare,output_directory + "/summary.txt", workingDir, bowtie,logger)
    
    #Calls two functions: Align and Score functions.
    #1. Align function: map each read set to reference sequence and creates SAM file, converts the SAM file to BAM file,Sort and index BAM file and Generate pileup
    #2. Score function:designate the correct allele and calculate coverage statistics for each locus 
    top_hits = try_and_except(stderr_log_output ,align_and_get_scores,workingDir, files, bowtie, samtools,log,logger,ids) ################### SP
    
    #files with the following extension:'.fasta', '.pkl', '.sam', '.tmp', '.bt2','.out','.unmapped','.unmap','fai' are removed from tmp and output files
    for root, dirs, files in os.walk(output_directory): 
        for currentFile in files:
            exts=('.fasta', '.pkl', '.sam', '.tmp', '.bt2','.out','.unmapped','.unmap','fai' )
            if any(currentFile.lower().endswith(ext) for ext in exts):
                os.remove(os.path.join(root, currentFile))
            elif currentFile in ['my_blast_tmp.xml', 'summary.txt', 'PHE221509-all.bam']:
                os.remove(os.path.join(root, currentFile))
    
    return top_hits
    
        
"""
Function
- Extract flanking regions of 100bp upstream and downstream of EMM by blast against a reference genome. BLAST uses the first variant sequence as a query.
NB Need to make sure BLAST, EMBOSS and Biopython are in the path.
- Create summary.txt file (a tab-delimited text file display the path to the variant sequences and flanking sequences) 

The option of the method
output_directory[str]: The path to where the summary.txt file will be created
profile_file_directory[str]: The path to the reference.seq and the EMM variant sequences (*.fas) files
logger[str]: The path to where the stderr and stdout logged
Return value
return  summary.txt file 
"""
def prep_SRST(profile_file_directory, output_directory, logger):
    
    reference_fasta_file = profile_file_directory + "/reference.seq"    
    refseq_record = SeqIO.read(reference_fasta_file, "fasta", generic_dna)
    locus_files = glob.glob(profile_file_directory + "/*.fas")
    locus_files = sorted(locus_files)
    summary_file_handle = open(output_directory + "/summary.txt", "w")
    for seq in locus_files:
        (seqDir,seqFileName) = os.path.split(seq)   
        (seqBaseName,ext) = os.path.splitext(seqFileName)
        bait = seqBaseName + "_bait.fasta"
        # extract first sequence to use as bait
        
        log_writer.info_header(logger, "create bait file")
        
        process = subprocess.Popen(['seqret',seq,'-firstonly','-auto','-out',output_directory+ '/' + bait], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        process.wait() 
        
        log_writer.log_process(logger, process, log_error_to = "info")
        
        cline = NcbiblastnCommandline(query=output_directory+ '/' + bait, db=profile_file_directory + "/reference",evalue=0.001, out=output_directory + "/my_blast_tmp.xml", outfmt=5)      
        stdout_log_output, stderr_log_output = cline()
        result_handle = open(output_directory + "/my_blast_tmp.xml")
        blast_record = NCBIXML.read(result_handle)
        query_length = blast_record.query_letters
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0] # only consider top hit
            if hsp.align_length/float(query_length) > 0.5:
                if hsp.sbjct_start > hsp.sbjct_end:
                    subject_start = hsp.sbjct_start + (hsp.query_start - 1)
                else:
                    subject_start = hsp.sbjct_start - (hsp.query_start - 1)
                if hsp.sbjct_start > hsp.sbjct_end:
                    subject_end = hsp.sbjct_end - (query_length - hsp.query_end)
                else:
                    subject_end = hsp.sbjct_end + (query_length - hsp.query_end)
                revcomp = 1 # hit is in forward strand
                if hsp.sbjct_start > hsp.sbjct_end:
                    revcomp = -1
                left_coords = [min(subject_start,subject_end)-100,min(subject_start,subject_end)-1]
                right_coords = [max(subject_start,subject_end)+1,max(subject_start,subject_end)+100]
                left_cmd = ["seqret ",reference_fasta_file," -sbegin ",str(left_coords[0])," -send ",str(left_coords[1])," -osformat fasta -auto -out " + output_directory + "/tmp_left_flank.fasta"]
                os.system(''.join(left_cmd)) # extract left flank using emboss
                
                right_cmd = ["seqret ",reference_fasta_file," -sbegin ",str(right_coords[0])," -send ",str(right_coords[1])," -osformat fasta -auto -out " + output_directory + "/tmp_right_flank.fasta"]
                os.system(''.join(right_cmd)) # extract right flank using emboss
                
                left_record = SeqIO.read(output_directory + "/tmp_left_flank.fasta", "fasta")
                
                if revcomp < 0:
                    left_record.id = "down"
                    left_record.seq = left_record.seq.reverse_complement() # reverse the sequence
                else:
                    left_record.id = "up"
                right_record = SeqIO.read(output_directory + "/tmp_right_flank.fasta", "fasta")
                if revcomp < 0:
                    right_record.id = "up"
                    right_record.seq = right_record.seq.reverse_complement() # reverse the sequence
                else:
                    right_record.id = "down"
                right_record.description = ""
                left_record.description = ""
                out_handle = open(output_directory + "/" + seqBaseName + "_flanks.fasta", "w")
                out_handle.write(right_record.format("fasta"))
                out_handle.write(left_record.format("fasta"))
                
                out_handle.close()
                # generate file list for srst
                summary_file_handle.write('\t'.join([seqBaseName,seq,output_directory + "/" + seqBaseName + "_flanks.fasta"]) + "\n")
    
    
    summary_file_handle.close()
                
"""
Function
1. Concatenate flanking regions extracted by prep SRST function to EMM variant sequences in fasta format. Newly concatenated sequence are then indexed by Bowtie2

2. Then extract and store as pickled object:
a. locus- variant names (loci.pkl)
b. start and end position of EMM variant sequences (without the flanking sequences)(ranges.pkl)
c. EMM variants sequence (refSeqs.pkl)

The option of the method
specFn[str]: A tab-delimited text file display the path to the flanking  and EMM sequences(summary.txt)
workDir[str] The path to where refSeqs.pkl, ranges.pkl and loci.pkl will be created
bowtie[str]: The command used to index the reference sequence
logger[str]: The path to where the stderr and stdout logged
return
return[list]  loci name
"""
def prepare(specFn, workDir, bowtie,logger):
    
    (specDir,summaryFileName) = os.path.split(specFn)
    spc = [] # (locus name, variants fasta, flanking sequences)
    for l in open(specFn):
        spc.append(l.split())
    refFn = os.path.join(workDir, "reference.fa")
    rf = open(refFn, "w") # file handle for reference sequence fasta file
    ranges = {}
    loci = [] # array of locus names
    refSeqs = {} # list of reference sequences (key = id, value = seq object)
    for (loc, variantsFn, flanksFn) in spc:
        loci.append(loc)
        fs = {} # flanking sequences at this locus (key = id, value = seq object)
        f = open(os.path.join(specDir, flanksFn))
        for r in SeqIO.parse(f, "fasta"):
            fs[r.id] = r.seq
        f = open(variantsFn)                                                            ################### SP
        for r in SeqIO.parse(f, "fasta"):
            s = Bio.Seq.MutableSeq('', generic_dna) 
            s += fs['up'] # add upstream seq, allele seq, downstream seq
            s += r.seq
            s += fs['down']
            SeqIO.write([SeqRecord(s, id=r.id)], rf, "fasta") # add to reference fasta file
            ranges[r.id] = (len(fs['up']), len(fs['up']) + len(r.seq)) # get range of allele sequence
            refSeqs[r.id] = s # store this reference sequence in list
    
    rf.close()
    rangesFn = os.path.join(workDir, "ranges.pkl") #start and end position of locus variant sequences (without the flanking sequences)
    f = open(rangesFn, 'w')
    pickle.dump(ranges, f)
    f.close()
    lociFn = os.path.join(workDir, "loci.pkl")
    f = open(lociFn, 'w')
    pickle.dump(loci, f)
    f.close()
    refSeqsFn = os.path.join(workDir, "refSeqs.pkl") #Locus variants sequence
    f = open(refSeqsFn, 'w')
    pickle.dump(refSeqs, f)
    f.close()
    bowtie2_index = bowtie + "-build"
    log_writer.info_header(logger, "bowtie_indexed")
    process = subprocess.Popen([bowtie2_index, refFn, refFn], stderr=subprocess.PIPE, stdout=subprocess.PIPE) # generate index of reference fasta for mapping
    process.wait()
    
    log_writer.log_process(logger, process, log_error_to = "info")
    os.system("rm -f summary.txt")
    return loci

"""
Function
Calls two functions: Align and Score functions.
Align function:
(a) Calls Bowtie2(with very senstivie options). Bowtie2 map each read set to reference sequence and creates SAM file
(b) Converts the SAM file to BAM file
(c) Sort and index BAM file
(d) Generate pileup file
Score function:
(a) From the pileup file read Depthofcoverage  and calculate probability score based on Depthofcoverage 
(b) Designate the correct allele based on probability score and degree of variability of the read from the locus variant (by identifying present or absence of SNPs/INDELs).
(c) Calculate coverage statistics for each locus (max_percenatge non consensus bases, Minimum total depth, Maximum total depth, Minimum consensus depth, 
Maximum consesnsus depth, mean consensus depth and stdDev of consensus depth) 
    
The option of the method
fileSets[dict]:  Keys are fastq filenames and value correspond to path to fastq file
bowtie[str]: The path to Bowtie2 command
samtools[str]: The path to SAMtools command
log[str]:  The path to where the EMM_log.txt will be created
logger[str]: The path to where the stderr and stdout logged

Return
Return   score[]probability and coverage statistics score value for each allele
"""

def align_and_get_scores(workingDir, files, bowtie, samtools, log,logger,ids):  ################### SP
    
    # insertSize = None
    # significance = -10 ## CUTOFF --> variable not used
    out = sys.stdout
    nameSep = "-"
    # verboseFiles = False --> variable not used
    scores = []
    paired = True
    
    pair = files
    if os.path.exists(pair[0]) and os.path.exists(pair[1]):
        align(workingDir, paired, pair, sys.stderr, bowtie, samtools,logger,ids)
        s = score(pair, workingDir, paired, out, log, nameSep, bowtie, samtools,logger,ids, log)
        out.flush()
        log.flush()
    else:
        log_writer.info_header(logger, "the paired reads are not labelled as as sampleid.R*.fastq*")
        
    return s


"""
Function
Prints EMM value and metrics for each variant with coverage > 90%
The option for method:
hits[dict] = {allele: [identity, coverage, meanDepth, minDepth, snps, indels, mixed, filteredCoverage], ...}
log[str]:  Location to where the EMM_log.txt file will be created

"""
def write_log(hits,log):
	# metrics[allele] = [identity, coverage, meanDepth, minDepth, snps, indels, mixed, filteredCoverage]
    if len(hits) > 0:
		print >> log
		print >> log, '=' * 70
		print >> log
		print >> log, "Results Summary"
		print >> log
		print >> log, "Allele\tidentity\tcoverage\tmeanDepth\tminDepth\tsnps\tindels\tmixed\tfilteredCoverage"
		sorted_hits = sorted(hits.items(), key=lambda x: x[1][0], reverse=True)
		i=0
		while sorted_hits[i][1][1] >= 90:
			print >> log, sorted_hits[i][0]+'\t'+str(sorted_hits[i][1][0])+'\t'+str(sorted_hits[i][1][1])+'\t'+str(sorted_hits[i][1][2])+'\t'+str(sorted_hits[i][1][3])+'\t'+str(len(sorted_hits[i][1][4]))+'\t'+str(sorted_hits[i][1][5])+'\t'+str(len(sorted_hits[i][1][6]))+'\t'+str(sorted_hits[i][1][7])
			i+=1

"""
Function
Calls bamify and pileupReads function

The option for method
workingDir[str]: The path to where  SAM, BAM, Pileup will be created
paired[bool]: true= paired end reads 
files[list]: The path to the fastq file location
logFile[str]: The path to the EMM_log.txt file location
bowtie[str]: The path to  Bowtie2 command 
samtools[str]: The path to SAMtools command
logger[str]: The path to where the stderr and stdout logged
ids[str]: unique identifier number

return
return pileup file
"""
def align(workDir, paired, files, logFile, bowtie, samtools,logger,ids): # removed insertSize since it wasn't used in the function
    refFn = os.path.join(workDir, "reference.fa")
    bam = bamify(workDir, 'all', files, refFn, True, logFile, bowtie, samtools,logger,ids)
    pileFn =  os.path.join(workDir, 'all.pileup')
    pileupReads(workDir, bam, refFn, open(pileFn, 'w'), logFile, samtools,logger)
"""
Function
(a) Map each read set to each of the possible EMM variants by calling Bowtie2 (with  very sensitive options) and create SAM and tmp file
(b) Convert the sam to tmp file by unsetting the secondary alignment bit score
(c) Convert the tmp to BAM file
(d) Sort BAM 

The option for method:
workDir[str]: The path to where the SAM, BAM and sorted BAM files will be created
pref[str]:  pref = "all" 
files[list]: The path to the fastq file location
refFn[str]:  The path to the  reference file location
expand[bool]:  True or false value 
logFile[str]: The path to the EMM_log.txt file location
bowtie[str]: The path to  Bowtie2 command 
samtools[str]: The path to SAMtools command
logger[str]: The path to where the stderr and stdout logged
ids[str]: unique identifier number

returns[string]: out: sorted BAM file
"""

def bamify(workDir, pref, files, refFn, expand, logFile, bowtie, samtools, logger, ids): # removed insertSize since it wasn't used in the function
    
    
    single = False
    out0 = os.path.join(workDir,ids + '-' + pref + '-all')# #pref = all
    out = os.path.join(workDir, ids + '-'+ pref + '-all.bam')
    tmp = os.path.join(workDir, ids + '-' + pref + '.tmp') # temporary sam output
    sam = os.path.join(workDir, ids + '-'+ pref + '.sam')
    bam = os.path.join(workDir, ids + '-'+ pref + '.bam')
    #un_reads = os.path.join(workDir, ids + '-'+ pref + '.unmapped')
    #un_conc_reads = os.path.join(workDir, ids + '-'+ pref + '.unmap')   
    
    if expand: # expand = true
        log_writer.info_header(logger, "Creating tmp file")
        #process = subprocess.Popen([bowtie, '--fr', '--minins', '300', '--maxins', '1100', '-x', refFn, '-1', files[0], '-2', files[1],'-S', tmp, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50', '--un', un_reads , '--un-conc',un_conc_reads], stderr=subprocess.PIPE, stdout=subprocess.PIPE) #refFn = refrence sequence, tmp = temorary sam output, -k = report up to 99999 good alignments per read, -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --very-sensitive option
        #process.wait()
        process = subprocess.Popen([bowtie, '--fr', '--no-unal', '--minins', '300', '--maxins', '1100', '-x', refFn, '-1', files[0], '-2', files[1], '-S', tmp, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50'], stderr=subprocess.PIPE, stdout=subprocess.PIPE) #refFn = refrence sequence, tmp = temorary sam output, -k = report up to 99999 good alignments per read, -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --very-sensitive option
        process.wait()
        #print ' '.join([bowtie, '--fr', '--minins', '300', '--maxins', '1100', '-x', refFn, '-1', files[0], '-2', files[1],'-S', tmp, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50', '--un', un_reads , '--un-conc',un_conc_reads])
        
        log_writer.log_process(logger, process, log_error_to = "info")
        
        log_writer.info_header(logger, "remove_secondary_mapping_bit")
        i = open(tmp)
        o = open(sam, 'w')
        remove_secondary_mapping_bit(tmp, sam)
        i.close()
        o.close()
        
    else:# expand = false, command is called within getNovelAllele function
        log_writer.info_header(logger, "Creating sam file")
        #process= subprocess.Popen([bowtie,  '--fr', '--minins', '300', '--maxins', '1100', '-x', refFn, '-1', files[0], '-2', files[1],'-S', sam, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50', '--un', un_reads , '--un-conc', un_conc_reads], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        #process.wait()
        process= subprocess.Popen([bowtie,  '--fr', '--no-unal', '--minins', '300', '--maxins', '1100', '-x', refFn, '-1', files[0], '-2', files[1],'-S', sam, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50'], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        process.wait()
        log_writer.log_process(logger, process, log_error_to = "info")
        
    
    log_writer.info_header(logger, "Converting sam to bam")
    process = subprocess.Popen([samtools, 'view', '-bhS', '-o', bam, sam], stderr=subprocess.PIPE, stdout=subprocess.PIPE)# convert to bam
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    
    log_writer.info_header(logger, "Sorting bam")
    process = subprocess.Popen([samtools, 'sort', bam, out0], stderr=subprocess.PIPE, stdout=subprocess.PIPE) # sort bam
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    return out

"""
Function:
Takes a SAM file and deducts 256 from the second column(FLAG) that unset the secondary alignment bit score 
NB: reads with bit(250) set are not reported when using Samtools pileup
The option for method:
sam[string]: SAM file 
sam_parsed[string]: parsed SAM file

Return
returns[string]: Parsed SAM file

"""
def remove_secondary_mapping_bit(sam,sam_parsed):
    lines = iter(fileinput.input([sam]))
    sam_parsed_file = open(sam_parsed, "w")
    headers = []
    body = []
    for line in lines:
        if line.startswith('@'):
            sam_parsed_file.write(line)
        else:
            # chomp line
            line = line.rstrip('\n')
            details = line.split("\t")
            flag = int(details[1])
            if flag > 256:
                details[1] = str(flag - 256)
            print >> sam_parsed_file, '\t'.join(details)
    sam_parsed_file.close()


"""
Function
Generate pileup file by using SAMtools mpileup command.
NB: use -B -A -f option to optimises coverage
--A flag count anomalous read 
 
The option for method:
workDir[str]: The path to where  pileup file will be created
bam[str]:  The path to the BAM file location
refFn[str]:  The path to the reference file location
pileupFile[str]: The path to pileup file location
logFile[str]: The path to the EMM_log.txt file location
samtools[str]: The path to SAMtools command
logger[str]: The path to where the stderr and stdout logged

return
returns: pileup file

"""
def pileupReads(workDir, bam, refFn, pileupFile, logFile, samtools,logger):
    log_writer.info_header(logger, "index bam file")
    process = subprocess.Popen([samtools, 'index', bam], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    
    log_writer.info_header(logger, "Generate pileup file")
    process = subprocess.Popen([samtools, 'mpileup', '-B', '-A', '-f', refFn, bam], stderr=subprocess.PIPE, stdout=subprocess.PIPE)# -A -count anomalous read pairs, -B - disable BAQ computation  and -f FILE - indexed reference sequence file
    for l in process.stdout:
        pileupFile.write(l)
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    pileupFile.close()
    
    
"""
Function

Parse through the pileup file and count the number of observed bases 

The option for method
pileupFile[str]: The path to pileup file location
refSeq[dict]: reference sequence for each allele 
ranges[dict]: start and end position of EMM variant sequences (without the flanking sequences)

Return
hash_alignment[dict]: hash_alignment[allele]: [pos, ref, orig_depth, filtered_match, filtered_mismatch, filtered_depth, total_indels, alt_bps, insertions_to_report] for each position
"""
def read_pileup(pileupFile, refSeq, ranges):
    

    # read pileup file
    with open(pileupFile) as pileup:
        
        hash_alignment = {}

        # Split all lines in the pileup by whitespace
        pileup_split = ( x.split() for x in pileup )
        # Group the split lines based on the first field (allele) 
        for allele, lines in groupby(pileup_split, itemgetter(0)):

            hash_alignment[allele] = []
            #if allele == 'emm75.0.sds':
            #    pass
            for fields in lines:
                alt_bps = {}
                locus = fields[0]
                nuc_num = int(fields[1])# Actual position in ref allele
                nuc_ref = fields[2]
                orig_depth = int(fields[3])
                if nuc_num <= ranges[locus][0] or nuc_num > ranges[locus][1]: 
                    continue
                elif nuc_num == ranges[locus][0]+1: # replaced +1 with + 30
                    ref_pos = nuc_num
                else:
                    try:
                        ref_pos += 1
                    except UnboundLocalError:
                        ref_pos = nuc_num
                if nuc_num != ref_pos:
                    for i in range(ref_pos, nuc_num-1): # if bps not covered
                        hash_alignment[allele].append([i+1, refSeq[allele][i], 0, 0, 0, 0, 'None', alt_bps, []])
                    ref_pos = nuc_num         # filter reads based on Phred scores - cutoff Q20 and calculate matches and mismatches
                if orig_depth != 0:
                    orig_match, orig_mismatch, nuc_match,nuc_mismatch, total_indels, alt_bps, report_insertions= pileup_extract_information(nuc_ref, fields[4], fields[5])
                    if nuc_num == ranges[allele][1] and report_insertions != []:
                        report_insertions = total_indels = []
                    elif report_insertions != []:
                        for ins in report_insertions:
                            ins_length = int(re.search('[0-9]+', ins).group())
                            ins_seq = re.search('[A-Z]+', ins).group()
                            if (ins_seq + str(refSeq[allele][nuc_num:ranges[allele][1]]))[:len(str(refSeq[allele][nuc_num:ranges[allele][1]]))] == str(refSeq[allele][nuc_num:ranges[allele][1]]):
                                report_insertions = report_insertions.remove(ins) if len(report_insertions)>1 else []
                            
                    nuc_depth = nuc_match + nuc_mismatch
                    if orig_match+orig_mismatch != orig_depth or nuc_depth > orig_depth:
                        print "Attention required!"
                        print "Line: {0}".format(fields)
                    elif nuc_num > ref_pos:
                        for i in range(ref_pos, nuc_num):
                            hash_alignment[allele].append([i+1, refSeq[allele][i], 0, 0, 0, 0, 'None', alt_bps, []])
                        hash_alignment[allele].append([nuc_num, nuc_ref, orig_depth, nuc_match, nuc_mismatch, nuc_depth, total_indels, alt_bps, report_insertions])
                        ref_pos = nuc_num
                    else:
                        # Hash for later processing in R
                        hash_alignment[allele].append([nuc_num, nuc_ref, orig_depth, nuc_match, nuc_mismatch, nuc_depth, total_indels, alt_bps, report_insertions])
                        #ref_pos += 1
                else: 
                    hash_alignment[allele].append([nuc_num, nuc_ref, 0, 0, 0, 0, 'None', alt_bps, []])
                    ref_pos += 1
    return hash_alignment

"""
Function

Parse through a single line of the pileup file and count the number of observed bases 

The option for method
ref_bp[str]: Reference base: 2nd column within the pileup file
align_bps[str]: reads bases: 4th column within the pileup file
qualities[str]: base qualities: 5th column within the pileup file

Return
filtered_match[int]: number of reads that match the reference with base quality > 20
filtered_mismatch[int]: number of reads that do not match the reference with base quality > 20
filtered_depth[int]: number of reads with base quality > 20
total_indels[list]: list of tuples (indel1, number of bases with indel1)
alt_bps[dict]: number of observed bases in this position
report_insertions[list]: insertion that occur in positions with depth > 4 and in more than half the reads

"""
def pileup_extract_information(ref_bp, align_bps, qualities):
    match = 0
    mismatch = 0
    probabilities = {}
    # remove all indels of format -1a and +1a since they do not have corresponding qualities and they refer to the following bp.
    # for example deletion  cat in positions 242-4 first appears in pos 241 with the actual deletions appearing as * in
    # the respective positions
    #ndh     241     C       79      ,$,,,,,,,,,,,,,,,,,,,,,-3catG,-3cat,,-3cat,,,,,,,,,,,,..,.,....,.........,.,.......,,...,,..... ;FFFFHGHIIHIJIIJJHII!)!*G)!!I!GGCD!!CCDDDDDDDDDCDBDDDDDDDDDCDDDJGJEDDJJJD@HFFFF
    #ndh     242     C       79      ,,,,,,,,,,,,,,,,,,,,-2at*A*,*,-2at,-2at,,-2at,,,,,-2at,-2at,,..,.,....,.........,.,.......,,...,,.....^~.       BCCFFFHIHBGIGGIIGHG!!!!J!!!H!GECC!!C>DDCDCBDDDDDDDDDBDDDDDCDDDJJJBCCJGJCDHFFFFC
    #ndh     243     A       79      ,,,,,,,,,,,,,,,,,,,**C*,***,*,,,,**,,..,.,....,.........,.,.......,,...,,...... @CCDDFHHGFGHGEHIFAH!!!!G!!!G!ACCC!!AADDDDCDDDDCDDDDDBDBDC?ABDDIIJD@CJIIACHHHGFC
    #ndh     244     T       79      ,$,$,$,,,,,,,,,,,,,,,,**C*,***,*,,,,**,,..,.,....,.........,.,.......,,...,,......      @B@DDFDEBDFEEFEI?CB!!!!H!!!B!@C>C!!:>DDDC@CDDDCDDDDDDDCCCC>DDDIGJF:>JJJ>CJHHHFC
    # The script below demonstrates that
    # >>> st1 = ',$,$.$,...,...,,.,,.,.,,,,.,,...,,.,,.,,.....,.,,,..,,.,.....,,......,,...,,,,.,.,,,,.,.,,.,,.,,,....,.,,.,.,,...+2AG,,,,,,,.,'
    # >>> st2 = '>>>CDB>HD<B@FBFCBFDJHDJBJJDDDJIDJI@H<DDDDDDDDHBDDDJDDDDDDDDDBDDDDDDDDDDDDDBDDDDDD>GDIDDJBDJDB<JIJJDIB@JDJBBJH8DBB<BDDCD'
    # >>> re.findall(r'[^,\.]', st1)
    # ['$', '$', '$', '+', '2', 'A', 'G']
    # >>> len(st1) - len(re.findall(r'[^,\.]', st1))
    # 119
    # >>> len(st2)
    # 119
    indels = re.findall(r'[\+\-]\d{1,2}[acgtnACGTN]*', align_bps)
    # it rectifies cases where the deletion is followed by a mismatch base - see above at position 241 (-3catG)
    Indels = list(set(indels))
    for index, i in enumerate(Indels):
        if len(re.search(r'[a-zA-Z]+', i).group()) != int(re.search(r'\d+', i).group()):
            dif = len(re.search(r'[a-zA-Z]+', i).group()) - int(re.search(r'\d+', i).group())
            indels = [i[:-dif] if x == i else x for x in indels]  
    for e in list(set(indels)):
        align_bps = (align_bps.replace(e,''))
    indels = [x.upper() for x in indels]
    indels_freq = Counter(indels).most_common() if indels != [] else 'None'
    # find all matches of format ^~. or .$ and remove extra symbols leaving only . and ,
    match1 = re.findall(r'\^[0-9a-zA-Z\!\ "#\$%&\'()\*\+,\.\-\/:;<>\?@\[\]\\\^_`\{\}\|~]{1}[,\.]{1}', align_bps)
    for e in list(set(match1)):
        align_bps = align_bps.replace(e, e[-1:])
    match2 = re.findall(r'[,\.]{1}\$', align_bps)
    for e in list(set(match2)):
        align_bps = align_bps.replace(e, e[:-1])
        
    # look for possible mismatches and remove extra symbols
    if set(','.join(align_bps)) != set([",","."]):
        mm1= re.findall(r'\^[0-9a-zA-Z\!\ "#\$%&\'()\*\+,\.\-\/:;<>\?@\[\]\\\^_`\{\}\|~]{1}[acgtACGT]{1}', align_bps)
        for e in list(set(mm1)):
            align_bps = align_bps.replace(e, e[-1:])
        mm2 = re.findall(r'[acgtACGT]{1}\$', align_bps)
        for e in list(set(mm2)):
            align_bps = align_bps.replace(e, e[:-1])
    
    # calculate matches and mismatches
    match = align_bps.count('.') + align_bps.count(',')
    mismatch = sum([align_bps.upper().count(x) for x in ('*', 'A', 'C', 'G', 'T', 'N')])
    
    # filter positions based on Q score  - cutoff 20
    alt_bps = {}
    filtered_mismatch = 0
    filtered_match = 0
    match_filtered = []
    # calculate filtered match and mismatch positions
    if mismatch > 0:
        mismatch_pos = {}
        for x in (r'\*{1}', r'A{1}', r'C{1}', 'G{1}', 'T{1}', 'N{1}'):
            pt = re.compile(x)
            alt = re.search(r'[\*\+ACGTN]{1}', x).group()[0]
            try:
                pt.search(align_bps.upper()).group()
            except AttributeError:
                continue
            else:
                mismatch_pos[alt] = [m.start() for m in pt.finditer(align_bps.upper())]
        mismatch_filtered = {}
        for bp in mismatch_pos.keys():
            mismatch_filtered[bp] = []
            for m in mismatch_pos[bp]:
                Q=ord(qualities[m])-33
                if Q > 20: mismatch_filtered[bp].append(m)
                
    
    match_pos = [m.start() for m in re.finditer(r'[,\.]{1}', align_bps)]
    for m in match_pos:
        Q= ord(qualities[m])-33
        if Q > 20: match_filtered.append(m)
        
    for bp in ('A', 'C', 'G', 'T', 'N', '*'):
        if bp != ref_bp:
            alt_bps[bp] = len(mismatch_filtered[bp]) if mismatch > 0 and bp in mismatch_filtered.keys() else 0
            filtered_mismatch += alt_bps[bp] 
        else:
            alt_bps[bp] = len(match_filtered)
            filtered_match = alt_bps[bp]
            
    insertions = [x for x in filter(lambda x:x[0]=="+", indels)]
    report_insertions = []
    for ins in list(set(insertions)):
        if filtered_match+filtered_mismatch > 4 and insertions.count(ins) > len(qualities)/2: # accepts a read if it occurs in more than half the reads
            report_insertions.append(ins)
    
            
    return match, mismatch, filtered_match, filtered_mismatch, indels_freq, alt_bps, report_insertions

"""
Function
Score function:
(a) Parse through the pileup file to capture DepthofCoverage
(b) Calculate probability score for all the EMM variants based on DepthofCoverage
(c) Calculate coverage statistics for all EMM variants (max_percentage non
consensus bases, Minimum total depth, Maximum total depth, Minimum consensus depth, Maximum consesnsus depth, mean consensus depth and stdDev of consensus depth) 
(d)For each locus display: locus names, variant number(allele number), number of snps different between the readset and locus,
minimum probability score , locus variant name, probability score value for each locus variant position (the probability score for each of the three bases other than the majority consensus base),
list of snps (SNP position, Reference base) and coverage statistics  for each locus variants 
(e) Then filter the correct allele based on probability score and degree of variability of the read from the locus variant (by identifying present or absence of SNPs/INDELs).
- If  number of snps different between the readset and locus variant is zero  and probability score greater than -10, the locus variant is assigned 

(f) Calculate percentage coverage (check if the reads are mapped to all locus variant position)


The option for method:
files[dict]:  Keys are fastq filenames and value correspond to path to fastq file
workDir[str]: output_directory -> The path to where logfile files will be created
sig[int]: -10 CUTOFF
paired[bool]:  assigned to False
insertSize[Nonetype]: insertSize assigned as None
outFile[str]: print result output 
logFile[str]:  The path to where the EMM_log.txt will be created
nameSep[str]: nameSep assigned as "-"
verboseFiles[bool]  False
bowtie[str]: The path to Bowtie2 command
samtools[str]: The path to SAMtools command
logger[str]: The path to where the stderr and stdout logged
ids[str]: unique identifier number

return
log file

"""
def score(files, workDir, paired, outFile, logFile, nameSep, bowtie, samtools,logger,ids, log):      # removed insertSize since it wasn't used in the function ###################   SP
    
    pileFn =  os.path.join(workDir, 'all.pileup')
    rangesFn = os.path.join(workDir, "ranges.pkl")# start and end position of locus variant sequences (without the flanking sequences)(ranges.pkl)
    refSeqsFn = os.path.join(workDir, "refSeqs.pkl")#Locus variants sequence ( refSeqs.pkl)
    
    ranges = pickle.load(open(rangesFn))
    refSeqs = pickle.load(open(refSeqsFn))

    # PILEUP FILE SIZE = ZERO: no reads mapped to emm references 
    if os.path.getsize(pileFn) == 0 :                                           ################### SP
        print >> logFile, "No reads mapped to any of the EMM reference sequences. Suggestion: check sequencing yield"   ################### SP
        print >> logFile                                                ################### SP
        return("Failed","No mapping to EMM references")

    validatedTypes = ['emm'+str(f) for f in range(1, 125)]

    hash_alignment = read_pileup(pileFn, refSeqs, ranges)

    top_hits = {}
    metrics = {}
    for allele in hash_alignment.keys():
        if hash_alignment[allele] == []: continue
        flanking = ranges[allele][0] + (len(refSeqs[allele])- ranges[allele][1])# added + 29
        matched_bps = len([f for f in hash_alignment[allele] if f[5]>4 and f[3]/float(f[5]) >= 0.8 and f[8] == []])
        unmatched_bps = [f for f in hash_alignment[allele] if f[5]<5 or f[3]/float(f[5]) < 0.8 or f[8] != [] or (f[6]!='None' and [m for m in f[6] if m[0].startswith('-') and m[1]/float(f[2]) > 0.5])]
        try:
            coverage = round(len(hash_alignment[allele])/float(len(refSeqs[allele])-flanking)*100, 4)
            identity = round(matched_bps/float(len(hash_alignment[allele]))*100, 4)
            filteredCoverage = round(len([f for f in hash_alignment[allele] if f[5]>4])/float(len(refSeqs[allele])-flanking)*100, 4)
        except ZeroDivisionError:
            identity = 0
        # select snps if mismatches more than 80% of the filtered depth (Q>20)
        snps = [(f[0],f[1], f[7]) for f in hash_alignment[allele] if f[5] > 4 and f[4]/float(f[5]) >= 0.8]
        indels = [(f[0], [m for m in f[6] if (m[0] in f[8]) or (m[1]/float(f[2]) > 0.5 and f[0]!=ranges[allele][1])]) for f in  unmatched_bps if f[8] != [] or f[6]!= 'None'] # an insertion at the very last bp of the allele would be a result of a mismatch in the flanking region and not a mutation in the allele
        indels = [f for f in indels if f[1] != []]
        mixed = [(f[0],f[1], f[7]) for f in hash_alignment[allele] if f[5] > 4 and f[3]/float(f[5]) < 0.8 and f[4]/float(f[5]) < 0.8]
        posDelEvents = []
        for f in indels:
            if f[1][0][0].startswith('-'):
                posDelEvents += range(f[0]+1, f[0]+1+int(re.search(r'\d+', f[1][0][0]).group()))
        # separate snps and deletions
        deletions = [f for f in snps if Counter(f[2]).most_common()[0][0] == '*' and f[0] not in posDelEvents] 
        snps = [f for f in snps if Counter(f[2]).most_common()[0][0] != '*']
        depths = [f[5] for f in hash_alignment[allele]]
        meanDepth = round(numpy.mean(depths), 2) if depths != [] else 0
        minDepth = min(depths) if depths != [] else 0
        if indels != [] and [f for f in indels if f[1][0][0].startswith('+')]: # only do this for insertions
            lenIndel = int(re.search('\d+', indels[0][1][0][0]).group())
            matched_bps -= lenIndel - 1
            identity = round(matched_bps/float(len(hash_alignment[allele])+lenIndel)*100, 4)
        metrics[allele] = [identity, coverage, meanDepth, minDepth, snps, indels, mixed, filteredCoverage]
    write_log(metrics, log)
    filtCovId = filter(lambda x: x[1][0]>=90.0 and x[1][-1]==100.0, metrics.items())
    with open(workDir+'/summary.yml', 'w') as out_fp:
        out_fp.write(yaml.dump(filtCovId, default_flow_style=True))
    if filtCovId == []:
        validatedMetrics = [(f, metrics[f]) for f in metrics.keys() if f.split('.')[0] in validatedTypes]
        nonValidatedMetrics = [(f, metrics[f]) for f in metrics.keys() if f.split('.')[0] not in validatedTypes]
        top_hit_failedValidated = sorted(validatedMetrics, key=lambda x:x[1][0], reverse=True) 
        top_hits['validated'] = ('Failed:'+top_hit_failedValidated[0][0], top_hit_failedValidated[0][1]) if validatedMetrics != [] else (None, ['n/a', 'n/a', 'n/a', 'n/a', [], [], [], 'n/a'])
        top_hit_failedNonValidated = sorted(nonValidatedMetrics, key=lambda x:x[1][0], reverse=True)
        top_hits['nonValidated'] = ('Failed:'+top_hit_failedNonValidated[0][0], top_hit_failedNonValidated[0][1]) if nonValidatedMetrics != [] else (None, ['n/a', 'n/a', 'n/a', 'n/a', [], [], [], 'n/a'])
    else:
        validated = [f for f in filtCovId if f[0].split('.')[0] in validatedTypes]
        nonValidated = [f for f in filtCovId if f[0].split('.')[0] not in validatedTypes]
        sortedValidated = sorted(validated, key=lambda x:x[1][0], reverse=True)
        sortedNonValidated = sorted(nonValidated, key=lambda x:x[1][0], reverse=True)
        if sortedValidated == []:
            top_hits['validated'] = (None, ['n/a', 'n/a', 'n/a', 'n/a', [], [], [], 'n/a'])
        elif len(set([f[1][0] for f in sortedValidated[:2]])) == 1:
            hits = [f for f in sortedValidated if f[1][0] == sortedValidated[0][1][0]]
            top_hits['validated'] = sortedValidated[0] if len(hits) == 1 else hits
        else:
            top_hits['validated'] = sortedValidated[0]
        if sortedNonValidated == []:
            top_hits['nonValidated'] = (None, ['n/a', 'n/a', 'n/a', 'n/a', [], [], [], 'n/a'])
        elif len(set([f[1][0] for f in sortedNonValidated[:2]])) == 1:
            hits = [f for f in sortedNonValidated if f[1][0] == sortedNonValidated[0][1][0]]
            top_hits['nonValidated'] = sortedNonValidated[0] if len(hits) == 1 else hits
        else:
            top_hits['nonValidated'] = sortedNonValidated[0]

    return top_hits


"""
Extracts the following values from scores data structure and  writes data to results.xml in the format below:

<ngs_sample id="PHE221920">
  <script value="emm typing tool" version="1-0"/>
  <results>
    <result type="EMM_validated" value="89.0">
      <result_data type="percentage_identity" value="100.00"/?
      <result_data type="percentage_coverage" value="100.00"/>
      <result_data type="mean_consensus_depth" value="109.09"/>
      <result_data type="minimum_consensus_depth" value="26"/>
    </result>
        <result type="EMM_Nonvalidated" value="232.0">
      <result_data type="percentage_identity" value="100.00"/?
      <result_data type="percentage_coverage" value="100.00"/>
      <result_data type="mean_consensus_depth" value="109.09"/>
      <result_data type="minimum_consensus_depth" value="26"/>
    </result>
  </results>
</ngs_sample>



- Id : NGS sample id (sample identifier)
- Version : software version number 
- EMM validated : assigned validated variant (emm1-124)
- ENN_Nonvalidated: assigned non validate variant (emm125+)
- mean consensus depth - the minimum average consensus depth
- number_of_reads_mapped: the number of reads mapped across allele length
- percentage_coverage: percentage coverage across allele length
- percentage_identity: percentage identity across allele length
- minimum consensus depth- the minimum  consensus depth values

The option for method:
output_directory[str]: The path to the result.xml file
xml_values: EMM value, emm score, QC mean consensus depth, QC max percentage non consensus base value, QC percentage coverage and QC minimum consensus depth
ids : NGS sample id (sample identifier)
Workflow_name: streptococcus-pyogenes-typing
Version : version number 

return
print results to "id".results.xml file
"""

def create_xml_file(top_hits,output_directory,ids,version):                            ################### SP
    #print scores
    xml_log_file = open(output_directory + "/" + ids + ".results.xml", "w") # open a file and write  xml result
    root = etree.Element("ngs_sample", id = ids)
    script = etree.SubElement(root, "script", value="emm typing tool", version = version) 
    results = etree.SubElement(root, 'results')
    if type(top_hits) == tuple:
        result = etree.SubElement(results, "result", type='Final_EMM_type', value = ':'.join(top_hits)) # Failed:No mapping to EMM references

    else:
        if type(top_hits['validated']) == tuple and top_hits['validated'][1][0] == 100:
            finalMtype = top_hits['validated'][0]
        elif top_hits['validated'][0] == None and ((type(top_hits['nonValidated']) == tuple and top_hits['nonValidated'][1][0] == 100) or (type(top_hits['nonValidated']) != tuple and top_hits['nonValidated'][0][1][0] == 100)):
            finalMtype = top_hits['nonValidated'][0] if type(top_hits['nonValidated']) == tuple else '/'.join([f[0] for f in top_hits['nonValidated']])
        elif (type(top_hits['validated']) != tuple and top_hits['validated'][0][1][0] < 100) or (type(top_hits['validated']) == tuple and top_hits['validated'][1][0]<100):
            if (type(top_hits['validated']) != tuple and len(set([f[0].split('.')[0] for f in top_hits['validated']])) == 1) or (type(top_hits['validated']) == tuple):
                top_hit = top_hits['validated'][0] if type(top_hits['validated']) != tuple else top_hits['validated']

                mutPos = 0 # investigate region 30-120 to determine whether the type can be reported or the sample should be ivestigated further before new type is reported
                for i,f in enumerate(top_hit[1][4:-1]):
                    for m in f:
                        if m[0] in range(131, 222):
                            if i == 1: # if ind
                                for indel in m[1]:
                                    mutPos += int(re.search('\d+', indel[0]).group())
                            else:
                                mutPos+=1
                newID = (90-mutPos)/float(90) * 100
                if newID >= 92:
                    finalMtype = top_hit[0].split('.')[0]
                else:
                    if type(top_hits['nonValidated']) == tuple and (type(top_hits['nonValidated']) == tuple and top_hits['nonValidated'][1][0] == 100) or (type(top_hits['nonValidated']) != tuple and top_hits['nonValidated'][0][1][0] == 100):
                        finalMtype = top_hits['nonValidated'][0]
                    else:
                        finalMtype = 'Investigate new type'
            elif (type(top_hits['validated']) != tuple and len(set([f[0].split('.')[0] for f in top_hits['validated']])) != 1):
                finalMtype = 'Not determined'
            elif [f for f in top_hits.values() if f[0].startswith('Failed')]:
                finalMtype = 'Failed'            
            else:
                finalMtype =[]
                for top_hit in top_hits['validated']:         
                    mutPos = 0 # investigate region 30-120 to determine whether the type can be reported or the sample should be investigated further before new type is reported
                    for i,f in enumerate(top_hit[1][4:-1]):
                        for m in f:
                            if m[0] in range(131, 222):
                                if i == 1: # if indels
                                    for indel in m[1]:
                                        mutPos += int(re.search('\d+', indel[0]).group())
                                else:
                                    mutPos+=1
                    newID = (90-mutPos)/float(90) * 100
                    if newID >= 92:
                        if top_hit[0].split('.')[0] not in finalMtype: finalMtype.append(top_hit[0].split('.')[0])
                finalMtype='/'.join(finalMtype)
                if finalMtype==[]:
                    if type(top_hits['nonValidated']) == tuple and (type(top_hits['nonValidated']) == tuple and top_hits['nonValidated'][1][0] == 100) or (type(top_hits['nonValidated']) != tuple and top_hits['nonValidated'][0][1][0] == 100):
                        finalMtype = top_hits['nonValidated'][0]
                    else:
                        finalMtype = 'Investigate new type'
                    
        elif top_hits['validated'][0] == None and (type(top_hits['nonValidated']) == tuple and top_hits['nonValidated'][1][0] < 100) or (type(top_hits['nonValidated']) != tuple and top_hits['nonValidated'][0][1][0] < 100):
            if type(top_hits['nonValidated']) != tuple and len(set([f[0].split('.')[0] for f in top_hits['nonValidated']])) != 1:
                finalMtype = 'Not determined' # if two different types are reported with the same coverage and covereage < 100 --> de novo assembly to resolve
            else:
                top_hit = top_hits['nonValidated'] if type(top_hits['nonValidated']) == tuple else top_hits['nonValidated'][0]
                mutPos = 0 # investigate region 30-120 to determine whether the type can be reported or the sample should be ivestigated further before new type is reported
                for i,f in enumerate(top_hit[1][4:-1]):
                    for m in f:
                        if m[0] in range(131, 232):
                            if i == 1: # if indels
                                for indel in m[1]:
                                    mutPos += int(re.search('\d+', indel[0]).group())
                            else:
                                mutPos+=1
                newID = (90-mutPos)/float(90) * 100
                if newID >= 92:
                    finalMtype = top_hit[0].split('.')[0]
                else:
                    finalMtype = 'Investigate new type'
        else:
            finalMtype = 'Not determined'
        comment = etree.Comment('(START) EMM Typing Results (START)')
        results.append(comment)
        result = etree.SubElement(results, "result", type='Final_EMM_type', value = finalMtype) 
        for key in top_hits.keys():
            emm_type = 'EMM_'+key
            EMM = str(top_hits[key][0]) if type(top_hits[key]) == tuple else '/'.join([f[0] for f in top_hits[key]])
            result = etree.SubElement(results, "result", type=emm_type, value = EMM) 
            pctID = top_hits[key][1][0] if type(top_hits[key]) == tuple else top_hits[key][0][1][0]
            pct_coverage = top_hits[key][1][7] if type(top_hits[key]) == tuple else top_hits[key][0][1][7]
            mean_depth = str(top_hits[key][1][2]) if type(top_hits[key]) == tuple else '/'.join([str(f[1][2]) for f in top_hits[key]])
            min_depth = str(top_hits[key][1][3]) if type(top_hits[key]) == tuple else '/'.join([str(f[1][3]) for f in top_hits[key]])
            snps = str(len([f for f in top_hits[key][1][4] if Counter(f[2]).most_common()[0][0] != '*'])) if type(top_hits[key]) == tuple else '/'.join([str(len([f for f in hit[1][4] if Counter(f[2]).most_common()[0][0] != '*'])) for hit in top_hits[key]])
            if type(top_hits[key]) == tuple:
                indels = re.search('\d+', top_hits[key][1][5][0][1][0][0]).group() if top_hits[key][1][5] != [] else '0'
            else:
                indels = []
                for hit in top_hits[key]:
                    if hit[1][5] != []:
                        indels.append(re.search('\d+', hit[1][5][0][1][0][0]).group())
                    else:
                        indels.append('0')
                indels = '/'.join(indels)
            mixed = str(len(top_hits[key][1][6])) if type(top_hits[key]) == tuple else '/'.join([str(len(hit[1][6])) for hit in top_hits[key]]) 
            etree.SubElement(result, "result_data", type="percentage_identity", value = str(pctID))
            etree.SubElement(result, "result_data", type="percentage_coverage", value=str(pct_coverage))
            etree.SubElement(result, "result_data", type="mean_consensus_depth", value=mean_depth)
            etree.SubElement(result, "result_data", type="minimum_consensus_depth", value=min_depth)
            etree.SubElement(result, "result_data", type="snps", value=snps)                
            etree.SubElement(result, "result_data", type="indels", value = indels)
            etree.SubElement(result, "result_data", type="mixed", value = mixed)
            print ids, emm_type, EMM, "indels:", indels, "SNPS:", snps, "mixed:", mixed
    print >> xml_log_file, etree.tostring(root, pretty_print=True)
