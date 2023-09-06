"""
Condense genome assembly by removing all masked bases, run BISER, and de-condense.
"""

import numpy as np
import os
import pandas as pd
import sys

import Bio.SeqIO
import Bio.bgzf

sys.path.append('dep/pav')
sys.path.append('dep/pav/dep/svpop')
sys.path.append('dep/pav/dep/svpop/dep')
sys.path.append('dep/pav/dep/svpop/dep/ply')
sys.path.append('dep/pav/dep/dupmapper')

import svpoplib
import pavlib
import libdupmap


#
# Default values
#

DEFAULT_MASK_CONDENSE = 50  # Masked regions within this distance are merged


#
# Configuration
#

if 'config_file' in config:
    configfile: config['config_file']
elif os.path.isfile('config.yaml'):
    configfile: 'config.yaml'
elif os.path.isfile('config.json'):
    configfile: 'config.json'

has_config = False
assembly_table = None

def init_parameters():
    """
    Read and check the assembly table. If the assembly table has already been read, do nothing.
    """

    global has_config

    # Config already completed.
    if has_config:
        return

    shell_prefix = 'set -euo pipefail; '

    # Check for unknown configuration parameters
    unknown_config = set(config.keys()) - {'sample', 'assembly', 'sample_table', 'mask_condense', 'config_file'}

    if unknown_config:
        raise RuntimeError('Found {} unknown configuration parameters: {}{}'.format(
            len(unknown_config),
            ', '.join([f'"{val}"' for val in sorted(unknown_config)[:3]]),
            '...' if len(unknown_config) > 3 else ''
        ))
    
    # Read sample table
    has_sample = 'sample' in config
    has_asm = 'assembly' in config
    has_tab = 'sample_table' in config
    
    if has_sample and has_asm and has_tab:
        raise RuntimeError('Cannot use "sample_table" with "sample" and "assembly" parameters (found all three in config).')

    if has_sample and has_asm:
        
        pass

    elif has_sample or has_asm:

        # sample and assembly are not used together
        found_param = 'sample' if has_sample else 'assembly'
        nfound_param = 'sample' if not has_sample else 'assembly'
        
        tab_warn = '. Parameter "sample_table" may not be used with "sample" and "assembly".' if has_tab else ''

        raise RuntimeError('Found "{found_param}" configuration parameter without "{nfound_param}" (must be used together){tab_warn}')

    elif has_tab:
        pass

    else:
        raise RuntimeError('Pipeline requires a sample table or "sample" and "assembly" configuration paramters. See README.md for setup help.')
    
    # Use sample/assembly from config
    if 'sample' in config or 'assembly' in config:

        if 'sample' not in config:
            raise RuntimeError('If config uses ')

        # Read from config first
        asm_filename = config['assembly'].strip()

        if not asm_filename:
            raise RuntimeError('Parameter "asm_filename" is defined with an empty value')

        if not os.path.isfile(asm_filename):
            raise RuntimeError(f'Parameter "asm_filename" does not point to a regular file: {asm_filename}')

        return asm_filename

    elif 'sample' in config:
        raise RuntimeError('Defined parameter "sample" without "assembly". Use parameters together to point to a single sample, or setup "samples.tsv" with SAMPLE and ASSEMBLY columns to point to multiple')

    if 'sample_table' in config:
        sample_table_filename = config['sample_table'].strip()

        if not sample_table_filename:
            raise RuntimeError('Parameter "sample_table" is defined with an empty value')

        if not os.path.isfile(sample_table_filename):
            raise RuntimeError(f'Parameter "sample_table" does not point to a regular file: {sample_table_filename}')
        
    else:
        sample_table_filename = 'samples.tsv'

        if not os.path.isfile(sample_table_filename):
            raise RuntimeError(f'Default sample table filename {sample_table_filename} is not a file. Use "sample_table" configuration parameter to point to one or "sample" and "assembly" to name a sample and point directly to the input FASTA.')

    # Read table
    df_sample = pd.read_csv(sample_table_filename, sep='\t')

    df_sample.columns = [col.upper() for col in df_sample.columns]

    missing_cols = [col for col in ('SAMPLE', 'ASSEMBLY') if col not in df_sample.columns]

    if missing_cols:
        raise RuntimeError(f'Missing columns in "{sample_table_filename}": {", ".join(missing_cols)}')

    df_sample.set_index['SAMPLE']['ASSEMBLY']

    # Get filename
    if sample not in df_sample.index:
        raise RuntimeError(f'Missing sample "{sample}" in "{sample_table_filename}"')

    asm_filename = df_sample[sample]

    if not os.path.isfile(asm_filename):
        raise RuntimeError(f'Sample "{sample}" in table "{sample_table_filename}" does not point to a regular file: {asm_filename}')

    # Mark configuration complete
    has_config = True


#
# Definitions
#


def get_assembly(sample):
    """
    Search pipeline configuration for the input FASTA file for a sample.

    :param sample: Sample name.

    :return: Assembly FASTA path (string).

    :raises RuntimeException: If the FASTA file cannot be found in the config or
        assembly table, the path is an empty string, or the FASTA file path does
        not point to a regular file.
    """





def fa_condensed_iter(seq_file_name, df_bed, format='fasta'):
    """
    Get an iterator over condensed sequence records.

    :param seq_file_name: Input sequence file name.
    :param df_bed: A BED file in sequence space indicating the masked loci to be condensed out (removed).
    """

    with svpoplib.seq.PlainOrGzReader(seq_file_name, 'rt') as in_file:
        seq_id_set = set()

        for seq in Bio.SeqIO.parse(in_file, format):

            if seq.id in seq_id_set:
                raise RuntimeError(f'Multiple sequence entries for sequence with ID "{seq.id}"')

            seq_id_set.add(seq.id)

            subdf = df_bed.loc[df_bed['#CHROM'] == seq.id].reset_index()

            if subdf.shape[0] > 0:
                seq_list = list()

                for index in range(subdf.shape[0] - 1):
                    seq_list.append(str(seq[subdf.loc[index, 'END']:subdf.loc[index + 1, 'POS']].seq))

            else:
                seq_list = [str(seq)]
        
        seq_condensed = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(''.join(seq_list)), id=seq.id, name='', description=''
        )

        yield seq_condensed


#
# Rules
#

rule cbis_decondense:
    input:
        bed_mask='{sample}/mask.bed.gz',
        bed=temp('temp/{sample}/biser_condensed.bed'),
        elem=temp('temp/{sample}/biser_condensed.bed.elem.txt')
    output:
        bed=temp('{sample}.bed.gz'),
        elem=temp('elem/{sample}.bed.elem.txt.gz')
    threads: 1
    run:

        # Read mask, translate to a PAV alignment BED with DELs for removed loci (PAV's liftover can read this file)
        df_bed_mask = pd.read_csv(input.bed_mask, sep='\t')

        pos = 0
        chrom_match_len = 0

        chrom_set = set()

        chrom = None
        cigar = ''

        df_map_list = list()

        chrom_index = 0

        condensed_chrom_len = dict()

        for index, row in df_bed_mask.iterrows():

            if chrom != row['#CHROM']:
                if chrom is not None:
                    if chrom in chrom_set:
                        raise RuntimeError(f'Error at mask BED {index}: Records for {row["#CHROM"]} are not contiguous.')

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

        # Write
        os.makedirs('temp/mm10_chr19', exist_ok=True)
        df_align.to_csv('temp/mm10_chr19/condensed_align.bed.gz', sep='\t', index=False, compression='gzip')

        df_len.to_csv('temp/mm10_chr19/condensed_chrom_len.tsv.gz', sep='\t', index=True, header=True, compression='gzip')



        ### Lift ###

        # Setup lift object
        df_align = pd.read_csv('temp/mm10_chr19/condensed_align.bed.gz', sep='\t')
        df_len = pd.read_csv('temp/mm10_chr19/condensed_chrom_len.tsv.gz', sep='\t', usecols=('CHROM', 'LEN'), index_col='CHROM').squeeze(axis=1)

        align_lift = pavlib.align.AlignLift(df_align, df_len)


        # Process records
        df_iter = pd.read_csv(BED_BISER, sep='\t', iterator=True, chunksize=1000, header=None)

        record_list = list()

        for df in df_iter:
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

                    raise RuntimeError('Error lifting record {row[0]}:{row[1]}-{row[2]} <-> {row[3]}:{row[4]}-{row[5]}: No lift for position(s): {missing_coord}')


                row[1] = lift_coord[0][1]
                row[2] = lift_coord[1][1] + 1
                row[4] = lift_coord[2][1]
                row[5] = lift_coord[3][1] + 1

                record_list.append(row)

        df = pd.concat(record_list, axis=1).T

        # Write
        df.to_csv('mm10_chr19/biser_uncondensed.bed.gz', sep='\t', index=False, header=False, compression='gzip')



# Run BISER on the condensed FASTA
rule cbis_biser:
    input:
        fa='{sample}/assembly_condensed.fa.gz'
    output:
        bed=temp('temp/{sample}/biser_condensed.bed'),
        elem=temp('temp/{sample}/biser_condensed.bed.elem.txt')
    params:
        temp_dir=lambda wildcards: config.get('temp_dir', 'temp/{sample}/biser_temp').format(sample=wildcards.sample)
    threads: 16
    shell:
        """biser """
            """-o {output.bed} """
            """-t {threads} """
            """--temp {params.temp_dir} """
            """{input.fa}"""

# Make condensed FASTA
rule cbis_condense_fa:
    input:
        bed_mask='{sample}/mask.bed.gz'
    output:
        fa='{sample}/assembly_condensed.fa.gz'
    threads: 1
    run:

        with Bio.bgzf.BgzfWriter(output.fa, 'wt') as out_file:
            Bio.SeqIO.write(
                fa_condensed_iter(
                    get_assembly(wildcards.sample),
                    pd.read_csv(input.bed_mask, sep='\t')
                ),
                out_file,
                'fasta'
            )

# Make filter BED
rule cbis_mask_bed:
    output:
        bed_mask='{sample}/mask.bed.gz'
    threads: 1
    run:

        # Get reference masked locations
        df_mask = libdupmap.ref.masked_fasta_to_bed(
            get_assembly(wildcards.sample),
            soft=True,
            dist=config.get('mask_condense', DEFAULT_MASK_CONDENSE)
        )

        # Write masked BED
        df_mask.to_csv(output.bed_mask, sep='\t', index=False, compression='gzip')
