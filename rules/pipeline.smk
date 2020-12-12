# imports
import os, sys, glob
from rules.utils.get_file import get_batch_ids_raw, get_signal_batch, get_sequence_batch
from rules.utils.storage import get_flowcell, get_kit, get_ID
# local rules
localrules: basecaller_merge_tag


# get batches
def get_batches_basecaller(wildcards):
    return expand("sequences/batches/{tag}/{runname}/{batch}.fastq.gz",
                        tag=wildcards.tag,
                        runname=config['runs'][wildcards.tag]['runname'],
                        batch=get_batch_ids_raw(config['runs'][wildcards.tag]['runname'], config=config, tag=wildcards.tag, checkpoints=checkpoints))

def get_batches_basecaller2(wildcards):
    batches = []
    for tag in config['runs']:
        batches.extend(
            expand("sequences/batches/{tag}/{runname}/{batch}.fastq.gz",
                                tag=tag,
                                runname=config['runs'][tag]['runname'],
                                batch=get_batch_ids_raw(config['runs'][wildcards.tag]['runname'], config=config, tag=wildcards.tag, checkpoints=checkpoints))
        )
    return batches

rule basecaller_merge_tag:
    input:
        lambda wildcards: get_batches_basecaller2(wildcards)
    output:
        "sequences/{tag, [^\/]*}.fastq.gz"
    run:
        with open(output[0], 'wb') as fp_out:
            for f in input:
                with open(f, 'rb') as fp_in:
                    fp_out.write(fp_in.read())

rule basecaller_stats:
    input:
        lambda wildcards: get_batches_basecaller(wildcards)
    output:
        "sequences/batches/{tag}/{runname}.hdf5"
    run:
        import gzip
        import pandas as pd
        def fastq_iter(iterable):
            while True:
                try:
                    title = next(iterable)
                    assert title[0] == '@'
                    seq = next(iterable)
                    _ = next(iterable)
                    qual = next(iterable)
                except StopIteration:
                    return
                mean_q = sum([ord(x) - 33 for x in qual]) / len(qual) if qual else 0.0
                yield len(seq), mean_q
        line_iter = (line for f in input for line in gzip.open(f, 'rb').read().decode('utf-8').split('\n') if line)
        df = pd.DataFrame(fastq_iter(line_iter), columns=['length', 'quality'])
        df.to_hdf(output[0], 'stats')


# get alignment batches
def get_batches_aligner(wildcards, config):
    r = expand("alignments/batches/{tag}/{runname}/{batch}.bam",
                        tag=wildcards.tag,
                        runname=config['runs'][wildcards.tag]['runname'],
                        batch=get_batch_ids_raw(config['runs'][wildcards.tag]['runname'], config=config, tag=wildcards.tag, checkpoints=checkpoints))
    return r

# minimap alignment
rule minimap2:
    input:
        sequence = 'sequences/{tag}.fastq.gz',
        reference = lambda wildcards: config['runs'][wildcards.tag]['reference']
    output:
        pipe("alignments/{tag, [^\/]*}.sam")
    threads: config['threads_alignment']
    group: "minimap2"
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    singularity:
        config['singularity_images']['alignment']
    shell:
        """
        {config[bin_singularity][minimap2]} -t {threads} {config[alignment_minimap2_flags]} {input.reference} {input.sequence} 1>> {output} 2> >(tee {output}.log >&2)
        if [ $(grep 'ERROR' {output}.log | wc -l) -gt 0 ]; then exit 1; else rm {output}.log; fi
        """

# sam to bam conversion and RG tag
rule aligner_sam2bam:
    input:
        sam = "alignments/{tag}.sam"
    output:
        bam = "alignments/{tag, [^\/]*}.bam",
        bai = "alignments/{tag, [^\/]*}.bam.bai"
    shadow: "minimal"
    threads: 1
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: int((1.0 + (0.2 * (attempt - 1))) * 5000)
    singularity:
        config['singularity_images']['alignment']
    shell:
        """
        {config[bin_singularity][samtools]} view -b {input.sam} | {config[bin_singularity][samtools]} sort -m 4G > {output.bam}
        {config[bin_singularity][samtools]} index {output.bam}
        """


# I think samtools install needs to be modified to do this. See https://www.biostars.org/p/389132/
# rule aligner_bam2cram:
#     input:
#         bam = "alignments/{tag}.bam",
#         bai = "alignments/{tag}.bam.bai",
#         reference = lambda wildcards: config['runs'][wildcards.tag]['reference']
#     output:
#         cram = "alignments/{tag}.cram",
#         crai = "alignments/{tag}.cram.crai"
#     shadow: "minimal"
#     threads: 1
#     resources:
#         threads = lambda wildcards, threads: threads,
#         mem_mb = lambda wildcards, attempt: int((1.0 + (0.2 * (attempt - 1))) * 5000)
#     singularity:
#         config['singularity_images']['alignment']
#     shell:
#         """
#         {config[bin_singularity][samtools]} view -T {input.reference} -C {input.bam} | {config[bin_singularity][samtools]} sort -m 4G > {output.cram}
#         {config[bin_singularity][samtools]} index {output.cram}
#         """

# mapping stats
rule aligner_stats:
    input:
        "alignments/batches/{tag}/{runname}_aligned-file-names.txt.txt"
    output:
        "alignments/batches/{tag, [^\/]*}/{runname, [^.\/]*}.hdf5"
    threads: config.get('threads_samtools') or 1
    resources:
        threads = lambda wildcards, threads: threads,
    singularity:
        config['singularity_images']['alignment']
    shell:
        """
        while IFS= read -r bam_file; do {config[bin_singularity][samtools]} view ${{bam_file}}; done < {input} | {config[bin_singularity][python]} {config[sbin_singularity][alignment_stats.py]} {output}
        """

rule demultiplex:
    input:
        'alignments/{tag}.bam'
    output:
        'demux/{tag}_{barcodes}.bam'
    script:
        'utils/demux.py'

rule index_demuxed:
    input:
        'demux/{tag}_{barcodes}.bam'
    output:
        'demux/{tag}_{barcodes}.bam.bai'
    shell:
        """
        samtools index {input}
        """

rule mutation_analysis:
    input:
        bam = 'demux/{tag}_{barcodes}.bam',
        bai = 'demux/{tag}_{barcodes}.bam.bai'
    output:
        expand('mutation_data/{{tag}}_{{barcodes}}_{datatype}', datatype = ['genotypes.csv', 'NT-muts-aggregated.csv', 'AA-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-distribution.cvs', 'failures.bam'])
    script:
        'utils/mutation_analysis.py'


# ---------------------------------------------------------------------------------
# Copyright (c) 2018-2020, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann, modified by Gordon Rix
# ---------------------------------------------------------------------------------