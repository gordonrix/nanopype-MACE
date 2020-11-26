"""
Demultiplexing of BAM files.

Input: BAM file, fasta file of terminal barcodes and/or internal barcodes

Output: Multiple BAM files containing demultiplexed reads, with the file name indicating the distinguishing barcodes.
    If terminal barcodes are not used, output name will be all_X.input.bam, where X is the internal barcode
    If internal barcode is not used, output name will be Y_all.input.bam, where Y is the terminal barcode pair used
"""

import datetime
import multiprocessing as mp
import os
import re
import shutil
import time
import gzip

import numpy as np
import pysam
from Bio import Align, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio import SeqIO

### Asign variables from config file
config = snakemake.config
tag = snakemake.wildcards.tag
# if 'demux_terminal_fasta' in config:
#     term_bc_fasta = config['demux_terminal_fasta']
# if 'demux_internal_fasta' in config:
#     internal_bc_fasta = config['demux_internal_fasta']
# runname = snakemake.wildcards.runname
# term_bc_pairs = config['runs'][runname]['barcodePairs']
BAMin = str(snakemake.input.alignment)
samfile = pysam.AlignmentFile(BAMin, 'rb')

refSeqfasta = str(snakemake.input.reference)
refSeqsList = list(SeqIO.parse(refSeqfasta, 'fasta'))
refDict = SeqIO.to_dict(SeqIO.parse(refSeqfasta, 'fasta'))

###

# function that creates a dictionary of barcodes from a fasta
#   barcodes must be in the form: >barcodeName
#                                 NNNNNNNN
#   function will then create a dictionary of these barcodes in the form NNNNNNNN: barcodeName
def createBarcodesDict(bcFASTA):
    bcDict = {}
    for entry in SeqIO.parse(bcFASTA, 'fasta'):
        bcName = entry.id
        bc = str(entry.seq)
        bcDict[bc]=bcName
    return bcDict

def aligned_reference(CIGAR_tuples, reference):
    """given a CIGAR tuple from pysam.AlignmentFile entry and a reference sequence as a string,
    this function will build the reference alignment string with indels accounted for"""
    index = 0
    string = ''
    for cTuple in CIGAR_tuples:
        if cTuple[0] == 0: #match
            string += reference[index:index+cTuple[1]]
            index += cTuple[1]
        elif cTuple[0] == 1: #insertion
            string += '-'*cTuple[1]
        elif cTuple[0]: #deletion
            index += cTuple[1]
    return string

def find_barcode_context(sequence, context):
    """given a sequence with barcode locations marked as Ns, and a barcode sequence context (e.g. ATCGNNNNCCGA),
    this function will identify the beginning and end of the Ns within the appropriate context if it exists only once.
    If the context does not exist or if the context appears more than once, then will return None as the first output
    value, and the reason why the barcode identification failed as the second value"""

    location = sequence.find(context)
    if location == -1:
        return None, 'barcode context not present in reference sequence'
    N_start = location + context.find('N')
    N_end = location + len(context) - context[::-1].find('N')

    if sequence[N_end:].find(context) == -1:
        return N_start, N_end
    else:
        return None, 'barcode context appears in reference more than once'

def sortByBarcodes(pairedRecords, bcDict):
    """given a list of pairs of fastq records (pairedRecords) and a barcode dictionary generated
    from the createBarcodesDict function (bcIn), this script will identify the pair
    of barcodes present at the beginning and end of each sequence in the fastq file
    and sort them into appropriately named files. If a name from the snakemake config
    file is assigned to the particular barcode pair that is identified, then that
    name will be used to name the file. If not, then the barcode names will be used
    to name the file. Sequences for which barcodes present in the bcFASTA file are
    not present at the sequence termini will be sorted into an 'unsorted' file."""

    # get length of barcode from barcodeDict
    for key in bcDict:
        bcLen = len(key)
        break

    # dictionary that will be populated with barcode pairs (or 'unsorted') as keys,
    #   and a pair of lists of fwd and rvs fastq records as values
    sortedbybarcodesDict = {'unsorted':([],[])}

    for fwd, rvs in pairedRecords:
        # assign sequence to Bio seq variable, forward and reverse barcodes to separate variables as strings
        fwdbc = str(fwd.seq[:bcLen])
        rvsbc = str(rvs.seq[:bcLen])

        # find name of barcode, move on to next sequence if either barcode is not in bcDict
        try:
            fwdbcName = bcDict[fwdbc]
            rvsbcName = bcDict[rvsbc]
        except:
            sortedbybarcodesDict['unsorted'][0].append(fwd)
            sortedbybarcodesDict['unsorted'][1].append(rvs)
            continue
        if (fwdbcName, rvsbcName) in sortedbybarcodesDict:
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][0].append(fwd)
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][1].append(rvs)
        else:
            sortedbybarcodesDict[(fwdbcName, rvsbcName)] = [[],[]]
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][0] = [fwd]
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][1] = [rvs]

    return sortedbybarcodesDict

def sortByBarcodes(BAMs, bcDict):
    """given a list of BAM file entries (pairedRecords) a barcode dictionary generated
    from the createBarcodesDict function (bcDict), and the 'barcodeInfo' subdictionary
    for the relevant tag from the config file, this script will identify the barcode(s)
    in each sequence and sort them into a nested dictionary, where global keys are the
    terminal barcode(s), values of keys are the name of the internal barcode, and
    values of values of keys are the BAM file entries belonging to that barcode
    combination in the form:

    {fwdbcName, rvsbcName):{terminalbcName:BAM file entry}}

    if a particular barcode is not specified in the config file, it is set to None.
    if a particular barcode is specified but cannot be identified according to the
    provided barcodes, it is set to 'notFound'
    """

    # dictionary that will be populated with barcode pairs (or 'unsorted') as keys,
    #   and a pair of lists of fwd and rvs fastq records as values
    sortedbybarcodesDict = {'unsorted':([],[])}

    for fwd, rvs in pairedRecords:
        # assign sequence to Bio seq variable, forward and reverse barcodes to separate variables as strings
        fwdbc = str(fwd.seq[:bcLen])
        rvsbc = str(rvs.seq[:bcLen])

        # find name of barcode, move on to next sequence if either barcode is not in bcDict
        try:
            fwdbcName = bcDict[fwdbc]
            rvsbcName = bcDict[rvsbc]
        except:
            sortedbybarcodesDict['unsorted'][0].append(fwd)
            sortedbybarcodesDict['unsorted'][1].append(rvs)
            continue
        if (fwdbcName, rvsbcName) in sortedbybarcodesDict:
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][0].append(fwd)
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][1].append(rvs)
        else:
            sortedbybarcodesDict[(fwdbcName, rvsbcName)] = [[],[]]
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][0] = [fwd]
            sortedbybarcodesDict[(fwdbcName, rvsbcName)][1] = [rvs]

    return sortedbybarcodesDict


def barcode_name_dict(barcodePairs):
    """builds dictionary used to find the name assigned to a barcode pair, if it exists"""
    outDict = {'unsorted':'unsorted'}
    for name in barcodePairs:
        bcPair = (barcodePairs[name]['fwdbc'], barcodePairs[name]['rvsbc'])
        outDict[bcPair] = name
    return outDict

bcDict_fwd = createBarcodesDict(snakemake.input.fwdBCfasta)
bcDict_rvs = createBarcodesDict(snakemake.input.rvsBCfasta)
context_fwd = config['runs'][tag]['barcodeInfo']['fwdBC']['context']
context_rvs = config['runs'][tag]['barcodeInfo']['rvsBC']['context']
# if config['runs'][tag]['barcodeInfo']['rvsBC']['reverseComplement'] == True:
    # context_rvs
print(context_fwd, context_rvs)

for reference in refDict:
    for entry in samfile.fetch('first_TrpB'):
        # print(entry)
        refAln = aligned_reference(entry.cigartuples, refDict[reference].seq)
        print(refAln)
        fwdStart, fwdStop = barcode_position(refAln, context_fwd)
        rvsStart, rvsStop = barcode_position(refAln, context_rvs)
        print(fwdStart, fwdStop, rvsStart, rvsStop)
        print(entry.query_alignment_sequence[fwdStart:fwdStop], entry.query_alignment_sequence[rvsStart:rvsStop])

def write_sorted_fqs(sortedDict):
    """given output from the sortByBarcodes function, this will loop through all fastq lists
    in the dictionary and write each to a separate new .fastq file. For barcode combinations
    that are specified in the config file, the provided name for the barcode pair
    will be used in the file name. Otherwise, the names for the individual barcodes,
    separated by '—' will be used in the file name"""
    bcNameLookup = barcode_name_dict(barcodePairs)
    for bcPair in sortedDict:
        if bcPair in bcNameLookup:
            fwdFileName = f'{sample}_{bcNameLookup[bcPair]}_fwd.fastq'
            rvsFileName = f'{sample}_{bcNameLookup[bcPair]}_rvs.fastq'
        else:
            fwdFileName = f'{sample}_{bcPair[0]}-{bcPair[1]}_fwd.fastq'
            rvsFileName = f'{sample}_{bcPair[0]}-{bcPair[1]}_rvs.fastq'
        fwdFile = os.path.join(outputDirectory, fwdFileName)
        rvsFile = os.path.join(outputDirectory, rvsFileName)
        with open(fwdFile, 'w') as fwdOut, open(rvsFile, 'w') as rvsOut:
            SeqIO.write(sortedDict[bcPair][0], fwdOut, 'fastq')
            SeqIO.write(sortedDict[bcPair][1], rvsOut, 'fastq')

def main():
    bcDict = createBarcodesDict(bcFASTA)
    pairedRecords = []
    for fwd, rvs in zip(SeqIO.parse(fwdFastq, 'fastq'), SeqIO.parse(rvsFastq, 'fastq')):
        pairedRecords.append((fwd, rvs))
    sortedDict = sortByBarcodes(pairedRecords, bcDict)
    write_sorted_fqs(sortedDict)
    return

# if __name__ == '__main__':
#     main()
