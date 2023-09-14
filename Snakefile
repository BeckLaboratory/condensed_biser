"""
Condense genome assembly by removing all masked bases, run BISER, and de-condense.
"""

import gzip
import numpy as np
import os
import pandas as pd
import shutil
import sys
import tempfile

import Bio.bgzf

PIPELINE_DIR = os.path.dirname(os.path.realpath(workflow.snakefile))

sys.path.append(PIPELINE_DIR)
sys.path.append(os.path.join(PIPELINE_DIR, 'submodules/pav'))
sys.path.append(os.path.join(PIPELINE_DIR, 'submodules/pav/dep/svpop'))
sys.path.append(os.path.join(PIPELINE_DIR, 'submodules/pav/dep/svpop/dep'))
sys.path.append(os.path.join(PIPELINE_DIR, 'submodules/pav/dep/svpop/dep/ply'))
sys.path.append(os.path.join(PIPELINE_DIR, 'submodules/dupmapper'))

import svpoplib
import pavlib
import libdupmap
import cbis


#
# Configuration
#

# Read config
if 'config_file' in config:
    configfile: config['config_file']
elif os.path.isfile('config.yaml'):
    configfile: 'config.yaml'
elif os.path.isfile('config.json'):
    configfile: 'config.json'

# Read parameters
run_params = cbis.config.get_params(config)

# Shell prefix
shell.prefix(run_params['shell_prefix'])


#
# Resources
#

# Sample table
df_sample_global = None

def get_sample_table():
    """
    Get the sample table. If it has already been read in this pipeline, return the global cache table.
    """

    global df_sample_global
    global run_params

    if df_sample_global is not None:
        return df_sample_global

    df_sample_global = cbis.config.read_sample_table(run_params)

    return df_sample_global

def get_assembly_path(wildcards):
    """
    Get the path to an assembly.

    :param wildcards: Rule wildcards with "sample" entry.
    """

    df_sample = get_sample_table().set_index('SAMPLE')['ASSEMBLY']

    if wildcards.sample not in df_sample.index:
        raise RuntimeError(f'No sample table entry for sample: {wildcards.sample}')

    return df_sample.loc[wildcards.sample]


#
# Rules
#

wildcard_constraints:
    sample='[^/]+'

localrules: cbis_all

rule cbis_all:
    input:
        bed=expand('{sample}/biser.bed.gz', sample=get_sample_table()['SAMPLE'])

ruleorder: cbis_mask_bed > cbis_decondense

# Decondense BISER run.
rule cbis_decondense:
    input:
        bed_mask='{sample}/mask.bed.gz',
        bed='temp/{sample}/biser_condensed.bed'
    output:
        bed='{sample}/biser.bed.gz',
        bed_dropped='{sample}/biser_dropped.bed.gz',
        elem='{sample}/elements.txt.gz'
    threads: 1
    run:

        # Read mask, translate to a PAV alignment BED with DELs for removed loci (PAV's liftover can read this file)
        df_mask = pd.read_csv(input.bed_mask, sep='\t')

        pos = 0
        chrom_match_len = 0

        chrom_set = set()

        chrom = None
        cigar = ''

        df_map_list = list()

        chrom_index = 0

        condensed_chrom_len = dict()

        for index, row in df_mask.iterrows():

            if chrom != row['#CHROM']:
                if chrom is not None:
                    if chrom in chrom_set:
                        raise RuntimeError(f'Error at mask BED {index}: Records for {row["#CHROM"]} are not contiguous.')

                    if chrom_match_len > 0:  # No lift if whole chromosome was masked (crashes align lift)
                        df_map_list.append(pd.Series(
                            [
                                chrom, 0, pos, chrom_index,
                                chrom, 0, chrom_match_len, 0, chrom_match_len,
                                np.nan, np.nan,
                                60, False, 0x00,
                                cigar
                            ],
                            index=[
                                '#CHROM', 'POS', 'END', 'INDEX',
                                'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_TIG_POS', 'QUERY_TIG_END',
                                'RG', 'AO',
                                'MAPQ', 'REV', 'FLAGS',
                                'CIGAR'
                            ]
                        ))

                    chrom_index += 1
                    chrom_set.add(chrom)

                    cigar = ''

                    condensed_chrom_len[chrom] = chrom_match_len

                chrom = row['#CHROM']
                cigar = ''
                pos = 0
                chrom_match_len = 0

            if row['POS'] < pos:
                raise RuntimeError(f'Error at mask BED {index}: Record POS is less than the current position in the genome: {row["POS"]} < {pos}')

            mask_len = row['END'] - row['POS']

            if mask_len < 1:
                raise RuntimeError(f'Error at mask BED {index}: Masked length must be 1 or greater: {mask_len}')

            match_len = row['POS'] - pos

            if match_len > 0:
                cigar += f'{match_len}='
                chrom_match_len += match_len

            cigar += f'{mask_len}D'

            pos = row['END']

        if chrom is not None:
            df_map_list.append(pd.Series(
                [
                    chrom, 0, pos, chrom_index,
                    chrom, 0, chrom_match_len, 0, chrom_match_len,
                    np.nan, np.nan,
                    60, False, 0x00,
                    cigar
                ],
                index=[
                    '#CHROM', 'POS', 'END', 'INDEX',
                    'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_TIG_POS', 'QUERY_TIG_END',
                    'RG', 'AO',
                    'MAPQ', 'REV', 'FLAGS',
                    'CIGAR'
                ]
            ))

            condensed_chrom_len[chrom] = chrom_match_len

        df_align = pd.concat(df_map_list, axis=1).T

        df_len = pd.Series(condensed_chrom_len)
        df_len.name = 'LEN'
        df_len.index.name = 'CHROM'

        # Remove chromosomes with no sequence
        chrom_rm_set = set(df_len[df_len == 0].index)

        df_align = df_align.loc[~ df_align['#CHROM'].isin(chrom_rm_set)]
        df_len = df_len.loc[~ df_len.index.isin(chrom_rm_set)]

        ### Lift ###

        # Setup lift object
        align_lift = pavlib.align.AlignLift(df_align, df_len)

        # Process records
        df_iter = pd.read_csv(input.bed, sep='\t', iterator=True, chunksize=1000, header=None)

        record_list = list()
        dropped_record_list = list()
        dropped_record_head = None

        for df in df_iter:

            if dropped_record_head is None:
                dropped_record_head = list(df.columns) + [
                    'POS_LEFT', 'END_LEFT', 'POS_RIGHT', 'END_RIGHT'
                    'MISSING_COORD'
                ]

            for index, row in df.iterrows():

                lift_coord = [
                    align_lift.lift_to_sub(row[0], row[1]),
                    align_lift.lift_to_sub(row[0], row[2] - 1),

                    align_lift.lift_to_sub(row[3], row[4]),
                    align_lift.lift_to_sub(row[3], row[5] - 1),
                ]

                if any([coord is None for coord in lift_coord]):
                    missing_coord = ', '.join([
                        coord_name for coord, coord_name in zip(
                            lift_coord,
                            ['Left POS', 'Left END', 'Right POS', 'Right END']
                        )
                            if coord is None
                    ])

                    row = pd.concat(
                        [
                        row,
                            pd.Series(
                                [
                                    lift_coord[0], lift_coord[1], lift_coord[2], lift_coord[3], missing_coord
                                ], index=[
                                    'POS_LEFT', 'END_LEFT', 'POS_RIGHT', 'END_RIGHT',
                                    'MISSING_COORD'
                                ]
                            )
                        ]
                    )

                    dropped_record_list.append(row)

                    # raise RuntimeError(f'Error lifting record {row[0]}:{row[1]}-{row[2]} <-> {row[3]}:{row[4]}-{row[5]}: No lift for position(s): {missing_coord}')
                else:
                    row[1] = lift_coord[0][1]
                    row[2] = lift_coord[1][1] + 1
                    row[4] = lift_coord[2][1]
                    row[5] = lift_coord[3][1] + 1

                    record_list.append(row)

        df = pd.concat(record_list, axis=1).T

        if len(dropped_record_list) > 0:
            df_dropped = pd.concat(dropped_record_list, axis=1).T
        else:
            df_dropped = pd.DataFrame([], columns=dropped_record_head)

        df.sort_values([1, 2, 3, 4, 5, 6], inplace=True)

        # Write BED
        df.to_csv(output.bed, sep='\t', index=False, header=False, compression='gzip')
        df_dropped.to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')

        # Compress elements file
        elem_in_filename = input.bed + '.elem.txt'
        elem_out_filename = f'{wildcards.sample}/elements.txt.gz'

        if os.path.isfile(elem_in_filename):
            with open(elem_in_filename, 'rb') as in_file:
                with gzip.open(elem_out_filename, 'wb') as out_file:
                    shutil.copyfileobj(in_file, out_file)


# Run BISER on the condensed FASTA
rule cbis_biser:
    input:
        fa='{sample}/assembly_condensed.fa.gz'
    output:
        bed=temp('temp/{sample}/biser_condensed.bed')
    threads: 1
    run:

        temp_dir = None

        try:
            os.makedirs(run_params['temp_dir'], exist_ok=True)
            temp_dir = tempfile.mkdtemp(prefix=f'biser_{wildcards.sample}.', dir=run_params['temp_dir'])

            print(f'Using temp dir: {temp_dir}', flush=True)

            shell(
                """zcat {input.fa} > {temp_dir}/assembly.fa && """
                """biser """
                    """-o {output.bed} """
                    """-t {threads} """
                    """--temp {temp_dir} """
                    """--gc-heap 10G """
                    """{temp_dir}/assembly.fa"""
            )

        finally:
            if temp_dir is not None and os.path.exists(temp_dir) and not run_params['keep_temp']:
                shutil.rmtree(
                    temp_dir,
                    onerror=lambda function, path, excinfo: print(f'Failed cleaning temp: {path}', file=sys.stderr)
                )


# Make condensed FASTA
rule cbis_condense_fa:
    input:
        asm=get_assembly_path,
        bed_mask='{sample}/mask.bed.gz'
    output:
        fa='{sample}/assembly_condensed.fa.gz'
    threads: 1
    run:

        with Bio.bgzf.BgzfWriter(output.fa, 'wt') as out_file:
            Bio.SeqIO.write(
                cbis.fa.fa_condensed_iter(
                    input.asm,
                    pd.read_csv(input.bed_mask, sep='\t')
                ),
                out_file,
                'fasta'
            )

# Make filter BED
rule cbis_mask_bed:
    input:
        asm=get_assembly_path
    output:
        bed_mask='{sample}/mask.bed.gz'
    threads: 1
    run:

        # Get reference masked locations
        df_mask = libdupmap.ref.masked_fasta_to_bed(
            input.asm,
            soft=True,
            dist=config.get('mask_condense', run_params['mask_condense'])
        )

        # Write masked BED
        df_mask.to_csv(output.bed_mask, sep='\t', index=False, compression='gzip')
