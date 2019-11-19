#!/usr/bin/env python

'''
Title: change_sam_queries.py
Date: 20190522
Author: Adam Nunn
Description:
	This program takes an input BAM file, and either masks base positions that are
	potentially influenced by bisulfite-seq OR all unassociated genomic variation
	observed from alignments between query and reference sequences. The script
	requires either an input FASTA sequence of the reference genome (slow) or a
	BAM file with MD tags (fast) as generated by 'samtools calmd'. Only works with
	alignments that have been generated from a directional Bisulfite-Seq protocol.

	eg. 'samtools calmd -b input.bam genome.fa 1> calmd.bam 2> calmd.log'

	++ R1 forward -> CT
	+- R1 reverse -> GA
	-+ R2 forward -> GA
	-- R2 reverse -> CT

	mask bisulfite -> use output to try to identify SNPs
	mask genomic -> convert output back to fastq, then kWIP clustering

	callMD before + after (should be no case where number of mismatches increases)
	use XB tag for non-directional libraries?


List of functions:
	main()
	build_genome()
	get_aligned_pairs_with_sequence()
	validate_MD()
	mask_genomic()
	mask_bisulfite()
	worker()
	merger()
	recursion()


Procedure:
	1. Try to build the genome dictionary if 'FASTA' is given
	2. Open pysam.AlignmentFile objects for 'BAM' and 'OUT'
	3. Iterate through each alignment in 'BAM' pysam.AlignmentFile object
	4. Determine whether CT or GA mismatches are relevant to the current alignment
	5. Mask either genomic information or bisulfite-seq information
	6. Write the new, modified alignment to the 'OUT' pysam.AlignmentFile object


Usage:
		./change_sam_queries.py [-h, --help] \
		[-f, --fasta] FASTA \
		[-G, --genomic] \
		[-Q, --quality] \
		[-T, --threads] <int> \
		[-t, --temp] TEMP \
		BAM OUT

eg. ./change_sam_queries.py -f reference.fa -G -Q -T 1 -t /var/tmp input.bam output.bam
eg. ./change_sam_queries.py input.bam output.bam
'''

###################
## INIT ENVIRONMENT

import multiprocessing as mp
import argparse, tempfile, resource, os
import pysam
from array import array


##################
## DEFINE __MAIN__
def main(BAM,OUT,GENOMIC,QUALITY,THREADS=1,TEMP=None,FASTA=None):

    ######
    # 1) DECLARE ENVIRONMENT

    # Build genome (if given)
    try: genome = build_genome(FASTA)
    except (TypeError):
        genome = None
        pass

    # Build threads pool
    pool = mp.Pool(int(THREADS))
    jobs = []

    # Determine open file limit
    rlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
    rlimit = rlimit[0]

    # Define temp directory
    TDIR = tempfile.mkdtemp(dir=TEMP)


    ######
    # 2) MODIFY ALIGNMENTS BASED ON USER CRITERIA

    # Open initial BAM file instance to parse reference contigs
    with pysam.AlignmentFile(BAM, "rb") as original:

        for ref in original.get_index_statistics():

            # Fire off workers to process reads for each scaffold 'ref' from 'BAM'
            if ref.mapped > 0:
                if genome: reference = genome[ref.contig]
                else: reference = None
                job = pool.apply_async(worker, (BAM, TDIR, ref.contig, reference, GENOMIC, QUALITY))
                jobs.append(job)


    ######
    # 3) MERGE OUTPUT FILES BASED ON OPEN FILE LIMIT

    # Merge results from the workers with a recursive merging function
    bam = recursion(jobs,rlimit,pool,TDIR,bool(genome))
    bam = bam[0].get()
    os.replace(bam, OUT)

    # close the pool
    pool.close()
    pool.join()

    ###################
    # Goodbye message
    print("\n----------------")
    print("Success!\n")


## END OF __MAIN__
##################


###################
## DEFINE FUNCTIONS

####### Function to build 'genome' dictionary from 'FASTA' file object (alternative to MD tag)
def build_genome(FASTA):

	######################################################################################################
	## FASTA = string of path to reference genome eg. "/path/to/reference.fa"                           ##
	######################################################################################################

	# declare an empty dictionary and stage the 'FASTA' file
	genome = {}
	with open(FASTA, 'r') as fasta:
		first = True

		# iterate through the lines of the 'FASTA' file
		for line in fasta:
			line = line.rstrip()

			# identify the first sequence
			if line.startswith(">") and (first == True):
				line = line.split(" ")
				chr = line[0][1:]
				sequence = ''
				first = False

			# identify all following sequences
			elif line.startswith(">") and (first == False):
				genome[chr] = sequence.upper()
				line = line.split(" ")
				chr = line[0][1:]
				sequence = ''

			# append the line to the concatenated sequence
			else: sequence += line

		# append the final sequence to the dictionary according to chr
		genome[chr] = sequence.upper()

	# return the constructed genome dictionary
	return genome



####### Function to get_aligned_pairs() with sequence from either MD tag or 'FASTA'
def get_aligned_pairs_with_sequence(read,genome):

	######################################################################################################
	## read = current iteration of pysam.AlignmentFile                                                  ##
	## genome = boolean evaluation on whether the genome has been built                                 ##
	######################################################################################################

	# test if genome has been built
	if genome:
		pairs = read.get_aligned_pairs()
		return pairs, False

    # try to use MD tags
	else:

		# error handling to report error message when MD tags are missing
		try: MD = read.get_tag("MD").upper()

		# in this case there is neither a genome or an MD tag in the current read
		except:
			print('\n{}\n'.format("ERROR: Valid FASTA file must be given if MD tags are missing in BAM file"))
			raise SystemExit(1)

		# MD tag is present, so can be used with get_aligned_pairs
		else:
			read.set_tag("MD",MD,"Z")
			pairs = read.get_aligned_pairs(with_seq=True)
			return pairs, MD



####### Function to count consecutive alignments with seemingly non-standard MD tags
def validate_MD(qseq,MD,CT,count):

	######################################################################################################
	## qseq = sequence string eg. "ACACGACTAGTCGT..."                                                   ##
	## MD = MD tag string eg. "11C17C8C16C4C12C8C6C21..."                                               ##
	## CT = boolean decision whether to read CT or GA mismatches eg. True or False                      ##
	## count = integer eg. 21                                                                           ##
	######################################################################################################

	# compare qseq with MD to determine count
	if (("T" in qseq) and (not "C" in MD) and CT) or (("A" in qseq) and (not "G" in MD) and not CT):
		count += 1

	return count



###### Function to mask matches, indels and non-CT/GA mismatches
def mask_genomic(reference,qseq,qual,pairs,CT):

	######################################################################################################
	## reference = dictionary[key] result of build_genome(), or None                                    ##
	## qseq = sequence string eg. "ACACGACTAGTCGT..."                                                   ##
	## qual = array of unsigned chars eg. array('B',[20,20,21,21,30,39,...])                            ##
	## pairs = result of read.get_aligned_pairs() eg. [(0, 78688), (1, 78689), (2, 78690), ...]         ##
	## CT = boolean decision whether to read CT or GA mismatches eg. True or False                      ##
	######################################################################################################

	new_seq, new_qual, new_cigar = "", array('B',[]), list()

	# modifying genomic data
	for bp in pairs:

		# remove insertions by ignoring them
		if bp[1] is None: continue
		else:

			# get reference base from genome or MD tag
			if reference: refbase = reference[bp[1]].upper()
			else: refbase = bp[2].upper()

		# restore deletions by incorporating the reference positions to the reads
		if bp[0] is None:

			# mask possible bisulfite positions with Ns since methylation status is unknown in restored deletions
			if (refbase == "C" and CT) or (refbase == "G" and not CT): new_seq += "N"
			else: new_seq += refbase
			new_qual.append(0)

		# all other cases are matches / mismatches
		else:

			# leave CT or GA mismatches unchanged (depending on read and strand)
			if qseq[bp[0]] == "T" and refbase == "C" and CT:
				new_seq += "T"
				new_qual.append(qual[bp[0]])
			elif qseq[bp[0]] == "A" and refbase == "G" and not CT:
				new_seq += "A"
				new_qual.append(qual[bp[0]])

			# leave matches unchanged
			elif qseq[bp[0]] == refbase:
				new_seq += refbase
				new_qual.append(qual[bp[0]])

			# mask all non-methylation mismatches with original reference base and reduce base quality
			else:
				new_seq += refbase
				new_qual.append(0)

	# generate new cigarTuple
	new_cigar.append((0, len(new_seq)))

    # return modified sequence, quality string and cigar
	return new_seq, new_qual, new_cigar




###### Function to mask CT/GA mismatches by modifying either nucleotides or base qualities
def mask_bisulfite(reference,qseq,qual,pairs,CT):

	######################################################################################################
	## reference = dictionary[key] result of build_genome(), or None                                    ##
	## qseq = a string of 'ACATCAGTCG...'                                                               ##
	## qual = an array of unsigned chars eg. array('B',[20,20,21,21,30,39,...])                         ##
	## pairs = result of read.get_aligned_pairs() eg. [(0, 78688), (1, 78689), (2, 78690), ...]         ##
	## CT = boolean decision whether to read CT or GA mismatches eg. True or False                      ##
	######################################################################################################

	new_seq, new_qual = "", array('B',[])

	# modifying bisulfite data
	for bp in pairs:

		# skip deletions
		if bp[0] is None: continue

		# positions unchanged on insertions
		elif bp[1] is None:
			new_qual.append(qual[bp[0]])
			new_seq += qseq[bp[0]]

		# all other cases are matches / mismatches
		else:

			# get reference base from genome or MD tag
			if reference: refbase = reference[bp[1]].upper()
			else: refbase = bp[2].upper()

			# check for CT or GA mismatches depending on read and strand
			if qseq[bp[0]] == "T" and refbase == "C" and CT:
				new_seq += "C"
				new_qual.append(0)
			elif qseq[bp[0]] == "A" and refbase == "G" and not CT:
				new_seq += "G"
				new_qual.append(0)

			# positions unchanged for matches and non-bisulfite mismatches
			else:
				new_seq += qseq[bp[0]]
				new_qual.append(qual[bp[0]])

    # return modified sequence and quality string
	return new_seq, new_qual



####### Function for READING the input reads, modifying, and sending them to 'q'
def worker(BAM,TDIR,ref,reference,GENOMIC,QUALITY):

	######################################################################################################
    ## BAM = path to input bam file eg. "/path/to/input.bam"                                            ##
    ## TDIR = path to temp dir for processed bam files eg. "/var/tmp/adfas7d"                           ##
    ## ref = current scaffold or chromosome reference eg. "Chr1"                                        ##
	## reference = dictionary[key] result of build_genome(), or None                                    ##
	## GENOMIC = boolean decision whether to mask genomic or bisulfite eg. True or False                ##
	## QUALITY = boolean decision whether to mask base nucleotides or base qualities eg. True or False  ##
	######################################################################################################

    # generate a name for the modified bam file
    name = TDIR + "/" + ref + ".bam"

    # Open pysam.AlignmentFile objects for reading and writing
    with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(name, "wb", header=original.header) as modified:

        # define global count variables
        rcount, count = 0, 0

        # iterate over each alignment from 'original'
        for alignment in original.fetch(ref):
            rcount += 1

            # filter out unmapped, secondary alignments, or qcfailed alignments
            if (alignment.is_unmapped) or (alignment.is_secondary) or (alignment.is_qcfail): continue

            # get alignment variables
            qseq = alignment.query_sequence
            qual = alignment.query_qualities

            # skip erroneous reads
            if (qseq == None) or (qual == None): continue

            # determine whether to read CT mismatches or GA mismatches
            CT = bool((alignment.is_read1 and not alignment.is_reverse) or (alignment.is_read2 and alignment.is_reverse))

            # get aligned pairs, with sequence from either genome (if given) or MD-tag
            pairs, MD = get_aligned_pairs_with_sequence(alignment,reference)

            # evaluate MD tag
            if MD:
                try: count = validate_MD(qseq,MD,CT,count)
                except (UnboundLocalError): pass
            
            # determine whether to mask genomic info info
            if GENOMIC:
                new = mask_genomic(reference,qseq,qual,pairs,CT)
                alignment.query_sequence = new[0]
                alignment.query_qualities = new[1]
                alignment.cigartuples = new[2]

            else:
                new = mask_bisulfite(reference,qseq,qual,pairs,CT)
                if QUALITY: alignment.query_qualities = new[1]
                else: alignment.query_sequence = new[0]

            # write read to 'modified' file
            modified.write(alignment)

    # return the filename path
    return name, rcount, count




####### Function for merging a list of bam files
def merger(TDIR,batch):

	######################################################################################################
	## TDIR = path to temp dir for temporary merge files eg. "/var/tmp/adfas7d"                         ##
	## batch = list of files to merge                                                                   ##
	######################################################################################################

	# set the tempfile name
	bam = tempfile.mktemp(dir=TDIR)

	# merge bams
	arguments = ["-f",bam]
	arguments = arguments + batch
	pysam.merge(*arguments)

	# return merged bam
	return bam



####### Function for splitting merge jobs into batches according to open file limit
def recursion(jobs,rlimit,pool,TDIR,genome):

	######################################################################################################
    ## jobs = a list of jobs submitted to the multiprocessing pool                                      ##
    ## rlimit = integer representing the system open file limit                                         ##
    ## pool = multiprocessing pool object                                                               ##
    ## TDIR = path to temp dir for temporary merge files eg. "/var/tmp/adfas7d"                         ##
    ## genome = boolean evaluation of whether genome exists                                             ##
	######################################################################################################

    # jobs contains a list of files to be merged
    if len(jobs) > 1:
        merges = []

        # set count vars for metadata
        rcount, count = 0, 0

        # build batches of files and according to open file limit
        for i in range(0, len(jobs), int(rlimit/2)):
            batch = []

            # add individual jobs to batch
            for job in jobs[i:i+int(rlimit/2)]:
                j = job.get()

                # on the first level, jobs contain metadata
                if isinstance(j, str): batch.append(j)
                else:
                    batch.append(j[0])
                    rcount += j[1]
                    count += j[2]
            
            # submit jobs for merging batchs to the pool
            merge = pool.apply_async(merger, (TDIR, batch))
            merges.append(merge)
        
        # print metadata
        if (rcount > 0):
            print("Total alignments read from BAM:\t\t{}".format(rcount))
            if not genome:
                perc = (count/rcount)*100
                perc = format(perc, '.2f')
                print("Potentially non-standard MD tags counted: {} (~{}% of total alignments)".format(count, perc))


        # now we have a new list of files so reapply the function for the next level of merging
        jobs = recursion(merges,rlimit,pool,TDIR,genome)

    # return the final list of jobs when the list contains 1
    return jobs



## END OF FUNCTIONS
###################

#############
## RUN SCRIPT

# define argparse
usage = ''' This program takes an input BAM file, and either masks base positions that are
	potentially influenced by bisulfite-seq OR all unassociated genomic variation
	observed from alignments between query and reference sequences. The script
	requires either an input FASTA sequence for the reference genome (slow) or a
	BAM file with MD tags (fast) as generated by 'samtools calmd'. Only works with
	alignments that have been generated from a directional Bisulfite-Seq protocol. '''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('inbam', metavar='<BAM>', help='[REQUIRED] Path to input BAM file')
parser.add_argument('outbam', metavar='<OUT>', help='[REQUIRED] Path to output BAM file')
parser.add_argument('-f','--fasta', metavar='<FASTA>', help='[OPTIONAL] Path to input FASTA file for reference genome. Not required if standard MD tags are present in BAM file.')
parser.add_argument('-t','--temp', metavar='<TEMP>', help='[OPTIONAL] Path to temp directory for use with --perform option. [default: /var/tmp]')
parser.add_argument('-G','--genomic', dest='GENOMIC', help='Mask all genomic differences unrelated to bisulfite sequencing [default: OFF]', default=False, action='store_true')
parser.add_argument('-Q','--quality', dest='QUALITY', help='Apply filter to base qualities instead of base nucleotides [default: OFF]', default=False, action='store_true')
parser.add_argument('-T','--threads', dest='THREADS', help='Number of threads for processing Input (NB: an additional thread is used for writing Output) [default: 1]', type=int, default=1)

args = parser.parse_args()
parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.inbam,args.outbam,args.GENOMIC,args.QUALITY,args.THREADS,args.temp,args.fasta)

## END OF SCRIPT
################