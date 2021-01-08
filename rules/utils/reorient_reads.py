""" script for nanoMACE pipeline

Uses biopython pairwise2 alignment to identify the orientation of each read in a fastq file
and reorient the read if necessary. Accepts a fastq file as input, and uses the snakemake
config file to determine the length of read to align, and outputs a fastq file of the same
length as the input but with all reads reoriented properly
"""

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
import statistics
import collections

### Asign variables from config file and inputs
config = snakemake.config
inputFastq = str(snakemake.input)
outputFastq = str(snakemake.output)
tag = snakemake.wildcards.tag
###

def main():

    refSeqFasta = config['runs'][tag]['reference']
    refList = list(SeqIO.parse(refSeqFasta, 'fasta'))
    refFull = refList[0]
    
    fStrandAlignRegion = config['runs'][tag]['orientation_alignment_sequence_f'] if 'orientation_alignment_sequence_f' in config['runs'][tag] else None
    rStrandAlignRegion = config['runs'][tag]['orientation_alignment_sequence_r'] if 'orientation_alignment_sequence_r' in config['runs'][tag] else None
    assert fStrandAlignRegion != None or rStrandAlignRegion != None, 'If `consistent_orientation` == False, must assign alignment for either beginning or end of sequences'

    # beginnings and ends of alignment regions in reference sequence
    if fStrandAlignRegion:
        fStrandAlignStart = str(refFull.seq).find(fStrandAlignRegion)
        assert fStrandAlignStart != -1, 'provided `orientation_alignment_sequence_f` must be present in reference sequence'
        fStrandAlignEnd = fStrandAlignStart+len(fStrandAlignRegion)
    if rStrandAlignRegion: # use reverse indexing
        rStrandAlignStart = str(refFull.seq)[::-1].find(rStrandAlignRegion[::-1])
        assert rStrandAlignStart != -1, 'provided `orientation_alignment_sequence_r` must be present in reference sequence'
        rStrandAlignEnd = rStrandAlignStart+len(rStrandAlignRegion)

    match, mismatch, gapOpen, gapExtend = 2, -1, -1, -1

    with open(outputFastq, 'w') as outFile:
        for record in SeqIO.parse(inputFastq,'fastq'):
            seq = record.seq

            fOrientationScores = []
            rOrientationScores = []

            if fStrandAlignRegion:
                fOrientation_fPrimer_aln = pairwise2.align.localms( str(seq)[fStrandAlignStart:fStrandAlignEnd], fStrandAlignRegion, match, mismatch, gapOpen, gapExtend, one_alignment_only=1 )[0]
                fOrientationScores.append(fOrientation_fPrimer_aln.score)
                rOrientation_fPrimer_aln = pairwise2.align.localms( str(seq[-fStrandAlignEnd:-fStrandAlignStart].reverse_complement()), fStrandAlignRegion, match, mismatch, gapOpen, gapExtend, one_alignment_only=1 )[0]
                rOrientationScores.append(rOrientation_fPrimer_aln.score)
            if rStrandAlignRegion:
                fOrientation_rPrimer_aln = pairwise2.align.localms( str(seq)[-rStrandAlignEnd:-rStrandAlignStart], rStrandAlignRegion, match, mismatch, gapOpen, gapExtend, one_alignment_only=1 )[0]
                fOrientationScores.append(fOrientation_rPrimer_aln.score)
                rOrientation_rPrimer_aln = pairwise2.align.localms( str(seq[rStrandAlignStart:rStrandAlignEnd].reverse_complement()), rStrandAlignRegion, match, mismatch, gapOpen, gapExtend, one_alignment_only=1 )[0]
                rOrientationScores.append(rOrientation_rPrimer_aln.score)

            if np.mean(fOrientationScores) >= np.mean(rOrientationScores):
                SeqIO.write(record, outFile, 'fastq')
            else:
                SeqIO.write(record.reverse_complement(id='rc_'+record.id, description='reverse complement'), outFile, 'fastq')

if __name__=='__main__':
    main()