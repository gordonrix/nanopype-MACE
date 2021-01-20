# imports
import os, sys, glob
from rules.utils.get_file import get_sequence_batch
from rules.utils.storage import get_flowcell, get_kit, get_ID
# local rules
localrules: basecaller_merge_tag

# check if tag maps to barcode PROBABLY NOT USEFUL
def get_tag_barcode(tag, runname, config):
    bc_batches = [None]
    if '__default__' in config['barcodes']:
        bc_batches.extend([config['barcodes']['__default__'][bc] for bc in config['barcodes']['__default__'].keys() if bc in tag])
    elif runname in config['barcodes']:
        bc_batches.extend([config['barcodes'][runname][bc] for bc in config['barcodes'][runname].keys() if bc in tag])
    else:
        pass
    return bc_batches[-1]

# get batch of reads as IDs or fast5
def get_signal_batch(wildcards, config):
    raw_dir = config['storage_data_raw']
    if hasattr(wildcards, 'tag'):
        tag_barcode = get_tag_barcode(wildcards.tag, wildcards.runname, config)
        if tag_barcode:
            return os.path.join('demux', config['demux_default'], 'barcodes', wildcards.runname, tag_barcode, wildcards.batch + '.txt')
    batch_file = os.path.join(raw_dir, wildcards.runname, 'reads', wildcards.batch)
    if os.path.isfile(batch_file + '.tar'):
        return batch_file + '.tar'
    elif os.path.isfile(batch_file + '.fast5'):
        return batch_file + '.fast5'
    else:
        return []

# prefix of raw read batches
def get_batch_ids_raw(runname, config, tag=None, checkpoints=None):
    tag_barcode = get_tag_barcode(tag, runname, config) if tag else None
    if tag_barcode and checkpoints:
        if hasattr(checkpoints, config['demux_default'] + '_barcode'):
            barcode_batch_dir = getattr(checkpoints, config['demux_default'] + '_barcode').get(runname=runname).output.barcodes
        else:
            raise NotImplementedError("Demultiplexing with {} is not implemented.".format(config['demux_default']))
        barcode_batch = os.path.join(barcode_batch_dir, tag_barcode, '{id}.txt')
        batches_txt, = glob_wildcards(barcode_batch)
        return batches_txt
    else:
        batches_tar, = glob_wildcards("{datadir}/{runname}/reads/{{id}}.tar".format(datadir=config["storage_data_raw"], runname=runname))
        batches_fast5, = glob_wildcards("{datadir}/{runname}/reads/{{id}}.fast5".format(datadir=config["storage_data_raw"], runname=runname))
        return batches_tar + batches_fast5

# get batches
def get_batches_basecaller(wildcards):
    return expand("sequences/batches/{tag}/{runname}/{batch}.fastq.gz",
                        tag = wildcards.tag,
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

if config['do_basecalling']:

    rule megalodon:
        input:
            runDir = lambda wildcards: os.path.join(config['storage_data_raw'], config['runs'][wildcards.tag]['runname']),
            reference = 'ref/trpB_nanopore_reference.fasta'
        output:
            outDir = directory('megalodon/{tag, [^\/_]*}'),
            done = 'megalodon/{tag, [^\/_]*}/done.txt'
        threads: 4
        params:
            guppy_config = lambda wildcards: config.get('basecalling_guppy_config')
        resources:
            GPU=1,
            threads=4
        shell:
            """
            megalodon {input.runDir} --reference {input.reference} --output-directory {output.outDir} --guppy-config {params.guppy_config} --overwrite --guppy-server-path ~/miniconda3/envs/megalodontest/bin/guppy_basecall_server --processes {threads} --guppy-timeout 1000.0 --outputs basecalls mappings
            touch {output}
            """

    # guppy basecalling
    rule guppy:
        input:
            batch = lambda wildcards : get_signal_batch(wildcards, config),
            run = lambda wildcards : [os.path.join(config['storage_data_raw'], wildcards.runname)] + ([os.path.join(config['storage_data_raw'], wildcards.runname, 'reads.fofn')] if get_signal_batch(wildcards, config).endswith('.txt') else [])
        output:
            ["sequences/batches/{tag}/{runname}/{batch, [^.]*}.fastq.gz"] +
            ["sequences/batches/{tag}/{runname}/{batch, [^.]*}.sequencing_summary.txt"] +
            (["sequences/batches/{tag}/{runname}/{batch, [^.]*}.hdf5"] if config.get('basecalling_guppy_config') and 'modbases' in config['basecalling_guppy_config'] else [])
        shadow: "shallow"
        threads: config['threads_basecalling']
        resources:
            threads = lambda wildcards, threads: threads,
            mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.1 * (attempt - 1))) * (config['memory']['guppy_basecaller'][0] + config['memory']['guppy_basecaller'][1] * threads)),
            time_min = lambda wildcards, threads, attempt: int((1440 / threads) * attempt * config['runtime']['guppy_basecaller']), # 90 min / 16 threads
            GPU = 1
        params:
            guppy_config = lambda wildcards : '-c {cfg}{flags}'.format(
                                cfg = config.get('basecalling_guppy_config') or 'dna_r9.4.1_450bps_fast.cfg',
                                flags = ' --fast5_out' if config.get('basecalling_guppy_config') and 'modbases' in config['basecalling_guppy_config'] else ''),
            guppy_server = lambda wildcards, input : '' if (config.get('basecalling_guppy_flags') and '--port' in config['basecalling_guppy_flags']) else '--port ' + config['basecalling_guppy_server'][hash(input.batch) % len(config['basecalling_guppy_server'])] if config.get('basecalling_guppy_server') else '',
            guppy_flags = lambda wildcards : config.get('basecalling_guppy_flags') or '',
            filtering = lambda wildcards : '--qscore_filtering --min_qscore {score}'.format(score = config['basecalling_guppy_qscore_filter']) if config['basecalling_guppy_qscore_filter'] > 0 else '',
            index = lambda wildcards : '--index ' + os.path.join(config['storage_data_raw'], wildcards.runname, 'reads.fofn') if get_signal_batch(wildcards, config).endswith('.txt') else '',
            mod_table = lambda wildcards, input, output : output[2] if len(output) == 3 else ''
        singularity:
            config['singularity_images']['basecalling']
        shell:
            """
            mkdir -p raw
            {config[bin_singularity][python]} {config[sbin_singularity][storage_fast5Index.py]} extract {input.batch} raw/ {params.index} --output_format bulk
            {config[bin_singularity][guppy_basecaller]} -i raw/ --recursive --num_callers 1 --cpu_threads_per_caller {threads} -s workspace/ {params.guppy_config}  {params.filtering} {params.guppy_flags} {params.guppy_server}
            FASTQ_DIR='workspace/pass'
            if [ \'{params.filtering}\' = '' ]; then
                FASTQ_DIR='workspace'
            fi
            find ${{FASTQ_DIR}} -regextype posix-extended -regex '^.*f(ast)?q' -exec cat {{}} \; | gzip > {output[0]}
            find ${{FASTQ_DIR}} -name 'sequencing_summary.txt' -exec mv {{}} {output[1]} \;
            if [ \'{params.mod_table}\' != '' ]; then
                {config[bin_singularity][python]} {config[sbin_singularity][basecalling_guppy_mod.py]} extract `find workspace/ -name '*.fast5'` {params.mod_table}
            fi
            """
    
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

rule basecaller_merge_tag:
    input:
        lambda wildcards: get_batches_basecaller(wildcards)
    output:
        "sequences/{tag, [^\/_]*}.fastq.gz"
    run:
        with open(output[0], 'wb') as fp_out:
            for f in input:
                with open(f, 'rb') as fp_in:
                    fp_out.write(fp_in.read())

rule UMI_extract:
    input:
        'sequences/{tag}.fastq.gz'
    output:
        extracted = 'sequences/{tag, [^\/_]*}_UMI-extract.fastq.gz',
        log = 'sequences/{tag, [^\/_]*}_UMI-extract.csv'
    params:
        fwd_regex_pattern = lambda wildcards: config['fwd_regex_pattern'],
        rvs_regex_pattern = lambda wildcards: config['rvs_regex_pattern']
    script:
        'utils/UMI_extract.py'

# rule unpack_sequences:
#     input:
#         'sequences/{tag}.fastq.gz'
#     output:
#         'sequences/{tag, [^\/_]*}.fastq'
#     shell:
#         """
#         gunzip -k {input}
#         """

# rule reorient_reads:
#     input:
#         'sequences/{tag}.fastq'
#     output:
#         'sequences/{tag, [^\/_]*}_reoriented.fastq'
#     script:
#         'utils/reorient_reads.py'

# rule reverse_complement:
#     input:
#         'sequences/{tag}_reoriented.fastq'
#     output:
#         'sequences/{tag, [^\/_]*}_reoriented_RC.fastq'
#     run:
#         from Bio import SeqIO
#         inputSeqIterator = SeqIO.parse(open(input[0], 'r'), format='fastq')
#         RCseqIterator = (record.reverse_complement(id=record.id+'-revcomp') for record in inputSeqIterator)
#         with open(output[0], 'w') as out:
#             SeqIO.write(sequences=RCseqIterator, handle=out, format='fastq')

# rule UMI_cluster:
#     input:
#         F = 'sequences/{tag}_reoriented.fastq',
#         R = 'sequences/{tag}_reoriented_RC.fastq'
#     output:
#         'sequences/{tag, [^\/_]*}.cluster'
#     params:
#         outDir = lambda wildcards, output: str(output).split('/')[:-1],
#         outPrefix = lambda wildcards: wildcards.tag+'.',
#         inputFileNameF = lambda wildcards, input: str(input.F).split('/')[-1],
#         inputFileNameR = lambda wildcards, input: str(input.R).split('/')[-1]
#     shell:
#         """
#         cd {params.outDir}
#         calib -f {params.inputFileNameF} -r {params.inputFileNameR} -o {params.outPrefix} -l1  -e 2 -k 8 -t 2 -m 7
#         cd ..
#         """

# rule UMI_consensus:
#     input:
#         cluster = 'sequences/{tag}.cluster',
#         F = 'sequences/{tag}_reoriented.fastq',
#         R = 'sequences/{tag}_reoriented_RC.fastq'
#     output:
#         Fcons = 'sequences/{tag, [^\/_]*}_Fconsensus.fastq',
#         FconsGZ = 'sequences/{tag, [^\/_]*}_Fconsensus.fastq.gz',
#         Fmsa = 'sequences/{tag, [^\/_]*}_Fconsensus.msa',
#         Rcons = 'sequences/{tag, [^\/_]*}_Rconsensus.fastq',
#         Rmsa = 'sequences/{tag, [^\/_]*}_Rconsensus.msa'
#     params:
#         Fprefix = 'sequences/{tag, [^\/_]*}_Fconsensus',
#         Rprefix = 'sequences/{tag, [^\/_]*}_Rconsensus'
#     shell:
#         """
#         calib_cons -c {input.cluster} -q {input.F} {input.R} -o {params.Fprefix} {params.Rprefix}
#         gzip -k {output.Fcons}
#         """

# # get alignment batches
# def get_batches_aligner(wildcards, config):
#     r = expand("alignments/batches/{tag}/{runname}/{batch}.bam",
#                         tag=wildcards.tag,
#                         runname=config['runs'][wildcards.tag]['runname'],
#                         batch=get_batch_ids_raw(config['runs'][wildcards.tag]['runname'], config=config, tag=wildcards.tag, checkpoints=checkpoints))
#     return r

# minimap alignment

def alignment_sequence_input(wildcards):
    if config['UMI_consensus']:
        return 'sequences/{tag}_Fconsensus.fastq'
    else:
        return 'sequences/{tag}.fastq.gz'

rule minimap2:
    input:
        sequence = alignment_sequence_input,
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

# sam to bam conversion
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

if config['demux']:

    checkpoint demultiplex:
        input:
            'alignments/{tag}.bam'
        output:
            'demux/{tag}_demultiplex_complete.txt'
        script:
            'utils/demux.py'

    rule index_demuxed:
        input:
            'demux/{tag}_{barcodes}.bam'
        output:
            'demux/{tag, [^\/]*}_{barcodes}.bam.bai'
        shell:
            """
            samtools index {input}
            """

    rule mutation_analysis:
        input:
            bam = 'demux/{tag}_{barcodes}.bam',
            bai = 'demux/{tag}_{barcodes}.bam.bai'
        output:
            expand('mutation_data/{{tag, [^\/]*}}_{{barcodes}}_{datatype}', datatype = ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-aggregated.csv', 'AA-muts-distribution.csv'] if config['do_AA_analysis'] else ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv'])
        script:
            'utils/mutation_analysis.py'

else:
    rule mutation_analysis:
        input:
            bam = 'alignments/{tag}.bam',
            bai = 'alignments/{tag}.bam.bai'
        output:
            expand('mutation_data/{{tag, [^\/]*}}_{{barcodes}}_{datatype}', datatype = ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-aggregated.csv', 'AA-muts-distribution.csv'] if config['do_AA_analysis'] else ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv'])
        script:
            'utils/mutation_analysis.py'

def mut_stats_input(wildcards):
    datatypes = ['genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-aggregated.csv', 'AA-muts-distribution.csv'] if config['do_AA_analysis'] else ['genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv']
    if config['demux']:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex_complete')[0]
        checkpoint_demux_files = checkpoint_demux_prefix + '{BCs}.bam'
        return expand('mutation_data/{tag}_{barcodes}_{datatype}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, datatype=datatypes)
    else:
        return expand('mutation_data/{tag}_all_{datatype}', tag=wildcards.tag, datatype=datatypes)

rule mut_stats:
	input:
		mut_stats_input
	output:
		'{tag, [^\/]*}_mutation-stats.csv'
	script:
		'utils/mutation_statistics.py'

rule plot_mutation_spectrum:
    input:
        '{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/]*}_mutation-spectra.html'
    script:
        'utils/plot_mutation_spectrum.py'

def plot_mutations_aggregated_input(wildcards):
    if config['demux']:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex_complete')[0]
        checkpoint_demux_files = checkpoint_demux_prefix + '{BCs}.bam'
        return expand('mutation_data/{tag}_{barcodes}_{AAorNT}-muts-aggregated.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, AAorNT = wildcards.AAorNT)
    else:
        return 'mutation_data/{tag}_all_{AAorNT}-muts-aggregated.csv'

rule plot_mutations_aggregated:
    input:
        aggregated = plot_mutations_aggregated_input,
        mutStats = '{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/]*}_{AAorNT}-mutations-aggregated.html'
    script:
        'utils/plot_mutations_aggregated.py'

def plot_mutations_distribution_input(wildcards):
    if config['demux']:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex_complete')[0]
        checkpoint_demux_files = checkpoint_demux_prefix + '{BCs}.bam'
        return expand('mutation_data/{tag}_{barcodes}_{AAorNT}-muts-distribution.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, AAorNT = wildcards.AAorNT)
    else:
        return 'mutation_data/{tag}_all_{AAorNT}-muts-distribution.csv'

rule plot_mutations_distribution:
    input:
        plot_mutations_distribution_input
    output:
        'plots/{tag, [^\/]*}_{AAorNT}-mutation-distributions.html'
    script:
        'utils/plot_mutation_distribution.py'


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