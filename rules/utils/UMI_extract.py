"""part of nanopype-MACE pipeline, written by Gordon Rix
UMI_extract.py
    identifies UMI barcodes in each sequence of a fastq based on regex pattern,
    appends this sequence to the sequence name, and writes these to a new .fastq file"""

import gzip
import re
import pandas as pd

### Asign variables from config file and inputs
config = snakemake.config
tag = snakemake.wildcards.tag
fastqIn = snakemake.input[0]
fastqOut = snakemake.output.extracted
logFile = snakemake.output.log
fwd_regex_pattern = re.compile(snakemake.params.fwd_regex_pattern)
rvs_regex_pattern = re.compile(snakemake.params.rvs_regex_pattern)

index = 0
fwdMatches = 0
rvsMatches = 0
noMatches = 0

with gzip.open(fastqOut, 'wb') as fqOut:
    for line in gzip.open(fastqIn, 'rt'):
        if index%4 == 0:
            seqID = line
        elif index%4 == 1:
            seq = line
        elif index%4 == 2:
            spacer = line
        elif index%4 == 3:
            qualities = line

            fwdMatch = fwd_regex_pattern.search(seq)
            rvsMatch = rvs_regex_pattern.search(seq)

            groupDict = None

            if fwdMatch and rvsMatch:
                assert False, 'BOTH MATCH'

            if fwdMatch:
                fwdMatches += 1
                groupDict = fwdMatch.groupdict()

            elif rvsMatch:
                rvsMatches += 1
                groupDict = rvsMatch.groupdict()
            
            if groupDict: # if match is found, append combined UMIs to sequence ID and write to compressed file

                seqIDsep = seqID.split(' ')

                seqID_UMIs = [seqIDsep[0] + '_' + groupDict['umi_1'] + groupDict['umi_2']]

                tags = seqIDsep[1:]

                seqID = ' '.join(seqID_UMIs + tags)

                for l in [seqID, seq, spacer, qualities]:
                    fqOut.write(l.encode('ascii'))
            
            else:
                noMatches += 1

        index += 1

pd.DataFrame({'tag':[tag], 'forward_umi_matches':[fwdMatches], 'reverse_umi_matches':[rvsMatches], 'total_umi_matches':[fwdMatches+rvsMatches], 'failures':[noMatches]}).to_csv(logFile, index=False)