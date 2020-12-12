from Bio.Seq import Seq
from Bio import SeqIO
# from Bio.Alphabet import generic_dna
# from Bio.Alphabet import IUPAC
# from Bio import pairwise2
# from Bio.pairwise2 import format_alignment
# from Bio import Align
# import itertools as IT
# import numpy as np
# import multiprocessing as mp
# import datetime
# import pandas as pd
# import os
# import time
# import shutil
import re
import pysam
import sys

### Asign variables from config file
config = snakemake.config
tag = snakemake.wildcards.tag
barcodes = snakemake.wildcards.barcodes
BAMin = str(snakemake.input.bam)

minQualThreshold = config['mutation_analysis_quality_score_minimum']
###

outputDir = 'mutation_data'

def main():
    x = MutationAnalysis(config, tag, BAMin)
    x.test()

class MutationAnalysis:

    @staticmethod
    def longest_ORF(sequence):
        return max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', sequence), key = len) # taken from https://stackoverflow.com/questions/31757876/python-find-longest-orf-in-dna-sequence#31758161

    def __init__(self, config, tag, BAMin):
        """
        arguments:

        config      - snakemake config dictionary for all runs
        tag             - tag for which all BAM files will be demultiplexed, defined in config file
        """
        refSeqfasta = config['runs'][tag]['reference']
        self.ref = list(SeqIO.parse(refSeqfasta, 'fasta'))[0]
        self.refTrimmed = list(SeqIO.parse(refSeqfasta, 'fasta'))[1]
        self.refStr = str(self.ref.seq)
        self.refTrimmedStr = str(self.refTrimmed.seq)

        assert (str(self.ref.seq).find(str(self.refTrimmed.seq)) != -1), f'Trimmed reference sequence ({self.refTrimmed.id}) must be substring of reference sequence ({self.ref.id})'

        self.refTrimmedStart = self.refStr.find(self.refTrimmedStr)
        self.refTrimmedEnd = len(self.refTrimmedStr) + self.refTrimmedStart
        # self.alignedTrimmedRef = ' '*(self.refStr.find(self.refTrimmedStr)) + self.refTrimmedStr + ' '*(self.refStr.find(self.refTrimmedStr)+len(self.refTrimmedStr)) # pad trimmed reference with spaces for alignment of reference cigar operations

        self.QSminimum = config['mutation_analysis_quality_score_minimum']
        self.doAAanalysis = config['do_AA_analysis']

        if self.doAAanalysis: # define protein sequence
            if config['auto_detect_longest_ORF']:
                try:
                    self.refProtein = self.longest_ORF(self.refTrimmedStr)
                except ValueError:
                    raise Exception(f'No ORF identifiable within {self.refTrimmed.id}')
            else:
                assert (len(list(SeqIO.parse(refSeqfasta, 'fasta'))) >= 3), f'Reference sequence for ORF not provided. Please provide reference ORF as third sequence in reference fasta file `{refSeqfasta}` or consider setting `auto_detect_longest_ORF` to True'
                refProtein = list(SeqIO.parse(refSeqfasta, 'fasta'))[2]
                assert self.refStr.find(refProtein), f'Provided reference for protein sequence `{refProtein.id}` not found within reference sequence ({self.ref.id}). Ensure protein sequence reference is a substring of reference'
                assert len(refProtein)%3 == 0, f'length of protein sequence reference `{refProtein.id}` not a multiple of 3, and therefore cannot be used as ORF'
                self.refProtein = refProtein
            self.refProteinStart = self.refTrimmedStr.find(self.refProtein)
            self.refProteinEnd = len(self.refProtein) + self.refProtein
            if len(self.refProtein)/len(self.refTrimmedStr) < 0.7:
                print(f'[WARNING] Length of protein sequence is under 70% of the length of trimmed reference sequence `{self.refTrimmed.id}`', file=sys.stderr)

        self.AAs = "ACDEFGHIKLMNPQRSTVWY*"
        self.NTs = "ATGC"


    def clean_alignment(self, BAMentry):
        """given a pysam.AlignmentFile BAM entry,
        trims the ends off the query and reference sequences according to the trimmed reference,
        aligns these two strings as well as the quality scores, and creates an alignment string between the
        two trimmed sequences  where '|'=match, and '.'=mismatch, and returns these three strings, the list of quality
        scores, a list of all insertions, and a list of all deletions
        """

        if BAMentry.reference_name != self.ref.id:
            self.alignmentFailureReason = 'alignment uses wrong reference sequence'
            return None
        
        if BAMentry.reference_length != len(self.ref.seq):
            self.alignmentFailureReason = 'length of alignment different from reference'
            return None

        refIndex = 0
        queryIndex = 0
        refAln = ''
        queryAln = ''
        queryQualities = []
        alignStr = ''
        insertions = [] # list of tuples where first element is index within trimmed reference, second element is sequence inserted
        deletions = []  # list of tuples where first element is index within trimmed reference, second element is number of bases deleted

        for cTuple in BAMentry.cigartuples: # build pretty alignment through combining consecutive segments, which are determined by the cigar tuples, in which the first value determines the operation and the second value determines the length

            if cTuple[0] == 0: #no indel
                refSegment = self.refStr[refIndex:refIndex+cTuple[1]]
                querySegment = BAMentry.query_alignment_sequence[queryIndex:queryIndex+cTuple[1]]
                alignSegment = ''.join(['|' if r==q else '.' for r,q in zip(refSegment,querySegment)])
                refAln += refSegment
                queryAln += querySegment
                queryQualities += BAMentry.query_alignment_qualities[queryIndex:queryIndex+cTuple[1]]
                alignStr += alignSegment
                refIndex += cTuple[1]
                queryIndex += cTuple[1]

            elif cTuple[0] == 1: #insertion, not added to sequence to maintain alignment
                if self.config['AA_analysis'] and cTuple[1]%3 != 0: # frameshift, discard sequence if protein sequence analysis is being done
                    self.alignmentFailureReason = 'frameshift insertion'
                    return None
                insertions.append((refIndex-self.refTrimmedStart, BAMentry.query_alignment_sequence[queryIndex:queryIndex+cTuple[1]]))
                queryIndex += cTuple[1]

            elif cTuple[0] == 2: #deletion, 'x' added to sequence to maintain alignment
                if self.config['AA_analysis'] and cTuple[1]%3 != 0: # frameshift, discard sequence if protein sequence analysis is being done
                    self.alignmentFailureReason = 'frameshift deletion'
                    return None
                refAln += self.refStr[refIndex:refIndex+cTuple[1]]
                queryAln += 'x'*cTuple[1]
                alignStr += ' '*cTuple[1]
                deletions.append((refIndex-self.refTrimmedStart, cTuple[1]))
                refIndex += cTuple[1]

        return  [refAln[self.refTrimmedStart:self.refTrimmedEnd],
                alignStr[self.refTrimmedStart:self.refTrimmedEnd],
                queryAln[self.refTrimmedStart:self.refTrimmedEnd],
                queryQualities[self.refTrimmedStart:self.refTrimmedEnd],
                insertions,
                deletions]


    def mut_array(self, cleanAlignment):
        """ Identify mutations in an aligned sequence

        Input: output of clean_alignment(). Mutations only counted if quality score
        is above globally set threshold.

        Outputs:
            ntMutArray  - x by y array where x is the length of self.refTrimmed, representing each
                            nucleotide position within this sequence, and y is four, representing the
                            possible nucleotides that a base could be mutated to. values in the array
                            will be set to 1 if a mutation at position x is mutated to nucleotide that
                            corresponds to position y. If no mutations are found, an array of 0s is returned 
            aaMutArray  - x by y array where x is the length of self.refORF, representing each
                            amino acid position within this sequence, and y is the length of self.AAs,
                            representing the number of
                            possible nucleotides that a base could be mutated to. values in the array
                            will be set to 1 if a mutation at position x is mutated to amino acid that
                            corresponds to position y. If no mutations are found, an array of 0s is returned
            genotype    - list of strings where each string is a list of different types of mutations separated by ', '
        """
        assert ref==self.refTrimmed, 'not supposed to happen' #remove ref if this doesn't throw an error ever                      <--- LOOK
        ref, alignStr, seq, qScores, insertions, deletions = cleanAlignment

        mismatches = [i for i,a in enumerate(alignStr) if a=='.']

        if self.doAAanalysis:
            indelCodons = [] # list of amino acid positions that are affected by indel (for insertion, insertion is within a codon; for deletion, at least one base of codon deleted)
            for index, _ in insertions:
                if self.refProteinStart <= index < self.refProteinEnd:
                    protIndex = index-self.refProteinStart
                    if protIndex%3 == 0: continue # ignore if insertion occurs between codons
                    indelCodons.append( int(protIndex/3) )
            for index, length in deletions:
                if self.refProteinStart <= index < self.refProteinEnd and self.refProteinStart <= index+length < self.refProteinEnd:
                    protIndexStart = index-self.refProteinStart
                    protIndexEnd = (index+length)-self.refProteinStart
                    firstCodon = int(protIndexStart/3)
                    lastCodon = int(protIndexEnd/3)
                    indelCodons.extend([i for i in range(firstCodon,lastCodon+1)])
            AAmutArray = np.zeros((int(len(self.refProtein)/3), len(self.AAs)), dtype=int)
        else:
            AAmutArray = None

        NTmutArray = np.zeros((int(len(self.)), len(self.NTs)), dtype=int)
        codonsChecked = []
        NTsubstitutions = []
        AAsubstitutionsNonsynonymous = []
        AAsubstitutionsSynonymous = []

        for i in mismatches:

            if qScores[i] < self.QSminimum: continue

            wtNT = ref[i]
            mutNT = seq[i]
            NTmutArray[i,self.NTs.find(mutNT)] += 1
            NTsubstitutions.append(wtNT+str(i)+mutNT)

            if self.doAAanalysis and self.refProteinStart <= i < self.refProteinEnd:

                protIndex = i-self.refProteinStart

                codon = int(protIndex/3)
                if codon in codonsChecked: continue
                codonsChecked.append(codon)
                codonPosi = protIndex%3
                codonIndices = list(range(i-codonPosi, i+(3-codonPosi)))

                # check that codon doesn't contain any bases influenced by an indel
                if any(i in codonIndices for i in indelCodons): continue

                #check that all three quality scores in codon are above threshold
                QStooLow = False
                codonQS = qScores[codonIndices[0]:codonIndices[2]]
                for qs in codonQS:
                    if qs < minQualThreshold:
                        QStooLow = True
                if QStooLow: continue

                wtAA = str(Seq(ref[codonIndices[0]:codonIndices[2]]).translate())
                mutAA = str(Seq(seq[codonIndices[0]:codonIndices[2]]).translate())
                if wtAA!=mutAA:
                    AAmutArray[codon, self.AAs.find(mutAA)] += 1
                    AAsubstitutionsNonsynonymous.append(wtAA+str(codon)+mutAA)
                else:
                    AAsubstitutionsSynonymous.append(wtAA+str(codon))
            
        return NtmutArray, AAmutArray
    
    def test(self):
        count = 0
        for bam in pysam.AlignmentFile(BAMin, 'rb'):
            if count<5:
                L = self.clean_alignment(bam)
                for a,b in zip([-50],[-1]):
                    print(L[0][a:b])
                    print(L[1][a:b])
                    print(L[2][a:b])
                    print(L[3][a:b])
                print(L[4])
                print(L[5])
            count+=1




def print_aligned_record(recordAln):
    """given the output from align_with_QS, will print the alignment
    in an easily readable format"""
    for i in range(0,len(recordAln[0]),50):
        start = i-50
        if start>=0:
            for n in range(0,4):
                s = recordAln[n]
                if n == 1:
                    spacer = ' '*(3-len(str(start+1)))
                    print(spacer+str(start+1), s[start:i], i)
                else:
                    print((' '*3),s[start:i])
            print('')
    return

def process_chunk(fqList, refSeq, firstResi):
    """Accepts a list of any number of fastq records (fqList), a reference sequence (refSeq),
        and the first residue in the provided protein sequence, and outputs the NT/AA distributions,
        NT/AA arrays, and a list of records that failed to be properly processed"""
    protLength = int(len(refSeq)/3)
    NTmutArray = np.zeros((int(len(refSeq)), len(NTs)), dtype=int)
    AAmutArray = np.zeros((protLength,len(AAs)), dtype=int)
    NTmutDist = np.zeros(int(len(refSeq)), dtype=int)
    AAmutDist = np.zeros(protLength, dtype=int)
    failures = []

    for record in fqList:
        recordAligned = align_with_QS(record, refSeq)
        if recordAligned:
            # SeqNtMuts, SeqAaMuts = mut_array(recordAligned)
            try:
                SeqNtMuts, SeqAaMuts = mut_array(recordAligned)
                NTmutArray = NTmutArray + SeqNtMuts
                AAmutArray = AAmutArray + SeqAaMuts
            except:
                failures.append(record)
                continue
            NTtotalMuts = sum(sum(SeqNtMuts))
            AAtotalMuts = sum(sum(SeqAaMuts))
            NTmutDist[NTtotalMuts] += 1
            AAmutDist[AAtotalMuts] += 1
            if printMut:
                print(record.id,'\n')
                print_aligned_record(recordAligned)
                print('')
        else:
            failures.append(record)

    return NTmutArray, AAmutArray, NTmutDist, AAmutDist, failures

# count lines in a file quickly
def linecount(filename):
    f = open(filename, 'rb')
    bufgen = IT.takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in IT.repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def get_mut_data(fq, refSeq, firstResi):
    """ Generates .csv files for amino acid- and nucleotide-level mutations and distributions
    of mutations in sequences. Parallelizes alignment by grouping up to 10,000 fastq
    records and processing separately.

    input: (List of fastq file names, Seq object reference sequence as a fasta file name,
    int residue number of first codon in provided reference sequence)
    Reference sequence must be in desired reading frame, without extraneous
    nucleotides on either end
    """

    if printMut:
        printRecs = []
        recordIter = SeqIO.parse(fq, 'fastq')
        for i,record in enumerate(recordIter):
            if i<30:
                printRecs.append(record)
            else: break
        NTmutArray, AAmutArray, NTmutDist, AAmutDist, failures = process_chunk(printRecs, refSeq, firstResi)
    else:
        # determine num of seqs to process per job such that they are relatively evenly distributed
        seqcount = linecount(fq)/4
        numWorkers = mp.cpu_count()
        print('available CPU count: ', numWorkers)
        seqsPerWorker = int(seqcount/numWorkers) + 1 #add 1 to ensure no leftovers on first pass if possible
        while seqsPerWorker > 10000:
            seqsPerWorker = int(seqsPerWorker/2)
        first = True

        # parallelize work, breaking up fastq into lists of SeqIO fastq records
        recordIter = SeqIO.parse(open(fq), 'fastq')
        pool = mp.Pool(numWorkers)
        for chunk in iter(lambda: list(IT.islice(recordIter, int(seqsPerWorker*numWorkers))), []): # removes chunks from iterator until only [] remains
            chunk = iter(chunk)
            pieces = list(iter(lambda: (list(IT.islice(chunk, seqsPerWorker)),refSeq,firstResi), ([],refSeq,firstResi))) # removes pieces from iterator and combines them with args for process_chunk until iterator is empty
            result = pool.starmap(process_chunk, pieces)

            # combining results
            resultsIter = iter(result)
            if first:
                allResults = list(result[0])
                next(resultsIter)
                first = False
            for r in resultsIter:
                for i in range(0,4):
                    allResults[i] += r[i]
                allResults[4].extend(r[4])

        pool.close()
        pool.join()

        NTmutArray, AAmutArray, NTmutDist, AAmutDist, failures = allResults

    ntIDs = list(str(refSeq))
    ntPositions = [f'{str(i)}' for i in range(0, len(refSeq))]
    WTnts = [ID+ntPosi for ID,ntPosi in zip(ntIDs,ntPositions)]
    NTmutDF = pd.DataFrame(NTmutArray, columns=list(NTs))
    NTmutDF['wt_nucleotides'] = pd.Series(WTnts)
    NTmutDF.set_index('wt_nucleotides', inplace=True)
    NTmutDF = NTmutDF.transpose()

    resiIDs = list(str(refSeq.translate()))
    protLength = int(len(refSeq)/3)
    resiPositions = [f'{str(i)}' for i in range(firstResi, firstResi+protLength)]
    WTresis = [ID+posi for ID,posi in zip(resiIDs,resiPositions)]
    AAmutDF = pd.DataFrame(AAmutArray, columns=list(AAs))
    AAmutDF['wt_residues'] = pd.Series(WTresis)
    AAmutDF.set_index('wt_residues', inplace=True)
    AAmutDF = AAmutDF.transpose()

    NTdistDF = pd.DataFrame(NTmutDist, columns=['seqs_with_n_NTsubstitutions'])
    NTdistDF.index.name = 'n'
    AAdistDF = pd.DataFrame(AAmutDist, columns=['seqs_with_n_AAsubstitutions'])
    AAdistDF.index.name = 'n'
    return NTmutDF, AAmutDF, NTdistDF, AAdistDF, failures

if __name__ == '__main__':

    main()



    # def Flip_Strand_If_Needed(fq, ref):
#     """ given a SeqIO record and a reference sequence,
#     this function will return either the same sequence and quality score or the reverse complement of that sequence
#     and the reverse of the quality score such that the output is the same strandedness as the reference sequence

#     Orientation is checked by attempting to identify any 15 bp chunks present in the reference
#     sequence that are also present in the sequence. If any are identified, then the sequence is not flipped,
#     if none are identified, then the sequence will be flipped (reverse complement is returned)
#     """
#     refS = str(ref)
#     seqS = str(fq.seq)
#     seqChunkList = []
#     i1 = -1
#     for i2 in range(0,len(ref),15):
#         if i1 != -1:
#             seqChunkList.append(seqS[i1:i2])
#         i1 = i2
#     # search 15 bp chunks of seq for presence in ref. If any are present, don't flip
#     doFlip = True
#     for chunk in seqChunkList:
#         if chunk in refS:
#             doFlip = False
#             break
#     if doFlip:
#         fq = fq.reverse_complement(id=True)
#     return fq

# def align_with_QS(fastqRecord, ref):
#     """ Produce alignment of a fastqrecord as a string, with quality scores
#     for bases as list of ints, all trimmed appropriately. Outputs None if
#     sequence does not pass alignment quality checks

#     inputs: individual fastq record produced from seqIO.parse and a sequence
#     to align the record to, as seq object
#     """
#     if flippyDo:
#         fastqRecord = Flip_Strand_If_Needed(fastqRecord,ref)
#     sequence = fastqRecord.seq
#     qScore = fastqRecord.letter_annotations['phred_quality']
#     A = list(map(float,alnScoringAlg))
#     alignment = pairwise2.align.localms(ref, sequence, A[0],A[1],A[2],A[3])
#     refAln, seqAln, alnScore, begin, end = (alignment[0])
#     alnScoreThreshold = len(ref)/2
#     if alnScore < alnScoreThreshold:
#         return None
#     alignF = format_alignment(*alignment[0])
#     newLineIndices = [i for i, a in enumerate(alignF) if a == '\n']
#     alignString = alignF[newLineIndices[0]+1:newLineIndices[1]]

#     # determines number of gaps at beginning of alignment string
#     startGap = 0
#     for i in alignString:
#         if i in ['|', '.']:
#             break
#         elif i == ' ':
#             startGap += 1

#     # put three alignment elements into same reference frame, trimming ends
#     refAln = refAln[begin:end]
#     alignString = alignString[startGap:end]
#     seqAln = seqAln[begin:end]
#     qScore = qScore[begin:end]

#     if len(ref)!=len(refAln)!=len(alignString)!=len(seqAln):
#         return None

#     refAAstart = ref[:3].translate()
#     refAAend = ref[-3:].translate()

#     # quality checks
#     if discardInDelSeqs and ' ' in alignString:
#         return None
#     if Seq(refAln[:3]).translate()!=refAAstart or Seq(refAln[-3:]).translate()!=refAAend:
#         return None

#     # clean gaps unless sequences with gaps are discarded
#     if not discardInDelSeqs:
#         try:
#             newRefAln = ''
#             newAlignString = ''
#             newSeqAln = ''
#             newQS = []
#             for i in range(0, len(alignString)):
#                 if alignString[i] == ' ':
#                     if seqAln[i] == '-':
#                         newRefAln += refAln[i]
#                         newAlignString += '|'
#                         newSeqAln += refAln[i]
#                         newQS.append(-1)
#                     else: pass
#                 else:
#                     newRefAln += refAln[i]
#                     newAlignString += alignString[i]
#                     newSeqAln += seqAln[i]
#                     newQS.append(qScore[i])
#             refAln = newRefAln
#             alignString = newAlignString
#             seqAln = newSeqAln
#             qScore = newQS
#         except:
#             return None

#     return [refAln, alignString, seqAln, qScore]