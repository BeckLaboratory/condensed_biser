# Pipeline configuration

import collections
import os
import pandas as pd

import svpoplib

DEFAULT_MASK_CONDENSE = 500  # Masked regions within this distance are merged

DEFAULT_TEMP_DIR = 'temp/tmp'


def get_params(config):
    """
    Get a parameters dictionary from config.

    :param config: Pipeline configuration dictionary.

    Keys:
    * sample_source [str]: "config" for single sample, "table" for a table of samples.
    * sample [str]: Sample name. Defined if sample_source is "config".
    * assembly [str]: Path to assembly. Defined if "sample_source" is "config"
    * sample_table [str]: Path to the sample table. Defined if "sample_source" is "table".
    * sample_table_format [str]: Format of the sample table ("tsv" or "xlsx"). Defined if "sample_source" is "table".
    * mask_condense [int]: Mask condense parameter. Set to DEFAULT_MASK_CONDENSE if not defined in the configuration.
    * shell_prefix [str]: String to prepend to shell commands.
    * temp_dir [str]: Path to temporary runtime files.
    * keep_temp [bool]: If True, keep BISER temporary directory.
    """

    if config is None:
        raise RuntimeError('config is None')

    params = dict()

    # Check for unknown configuration parameters
    unknown_config = set(config.keys()) - {
        'assembly', 'sample', 'sample_table',
        'mask_condense', 'config_file', 'shell_prefix',
        'temp_dir', 'keep_temp'
    }

    if unknown_config:
        raise RuntimeError('Found {} unknown configuration parameters: {}{}'.format(
            len(unknown_config),
            ', '.join([f'"{val}"' for val in sorted(unknown_config)[:3]]),
            '...' if len(unknown_config) > 3 else ''
        ))

    # Locate sample source
    has_sample = 'sample' in config
    has_asm = 'assembly' in config
    has_tab = 'sample_table' in config

    if has_sample and has_asm:

        if has_tab:
            raise RuntimeError('Cannot used configuration parameter "sample_table" with "sample" and "assembly"')

        params['sample_source'] = 'config'
        params['sample'] = config['sample'].strip()
        params['assembly'] = config['assembly'].strip()

        if params['sample'] == '':
            raise RuntimeError('Configuration parameter "sample" is empty')

        if params['assembly'] == '':
            raise RuntimeError('Configuration parameter "assembly" is empty')

    elif has_sample or has_asm:
        # Sample and assembly not used together

        tab_warn = '. Parameter "sample_table" may not be used with "sample" and "assembly".' if has_tab else ''

        if has_sample:
            raise RuntimeError(f'Found "sample" configuration parameter without "assembly" (must be used together){tab_warn}')
        else:
            raise RuntimeError(f'Found "assembly" configuration parameter without "sample" (must be used together){tab_warn}')

    else:
        params['sample_source'] = 'table'

        # Locate sample table name
        if 'sample_table' in config:
            sample_table = config['sample_table'].strip()

            if not sample_table:
                raise RuntimeError('Configured sample table file name (parameter "sample_table") is empty')

            if not os.path.isfile(sample_table):
                raise RuntimeError(f'Configured sample table file name (parameter "sample_table") does not exist or is not a regular file: {sample_table}')

        elif os.path.isfile('samples.xlsx'):
            sample_table = 'samples.xlsx'

        elif os.path.isfile('samples.tsv'):
            sample_table = 'samples.tsv'

        if sample_table.lower().endswith('.xlsx'):
            params['sample_table_format'] = 'xlsx'
        elif sample_table.lower().endswith('.tsv') or sample_table.lower().endswith('.tsv.gz'):
            params['sample_table_format'] = 'tsv'
        else:
            raise RuntimeError(f'Configured sample table file name (parameter "sample_table") does not have a ".xlsx", ".tsv", or ".tsv.gz" file extension: {sample_table}')

        params['sample_table'] = sample_table

    # Parameter shel_prefix
    params['shell_prefix'] = 'set -euo pipefail; ' + (config['shell_prefix'] if 'shell_prefix' in config else '')

    # Parameter mask_condense
    if 'mask_condense' in config:
        try:
            params['mask_condense'] = int(params['mask_condense'])
        except ValueError:
            raise RuntimeError(f'Parameter "mask_condense" is not an integer: {params["mask_condense"]}')

    else:
        params['mask_condense'] = DEFAULT_MASK_CONDENSE

    if params['mask_condense'] < 0:
        params['mask_condense'] = 0

    # Parameter temp_dir
    if 'temp_dir' in config:
        temp_dir = config['temp_dir'].strip()
    else:
        temp_dir = ''

    if temp_dir == '':
        temp_dir = DEFAULT_TEMP_DIR

    params['temp_dir'] = temp_dir

    # Parameter keep_temp
    if 'keep_temp' in config:
        if config['keep_temp'] is None:
            raise RuntimeError(f'Parameter "keep_temp" is present with a missing value')

        try:
            params['keep_temp'] = svpoplib.util.as_bool(config['keep_temp'])
        except ValueError:
            raise RuntimeError(f'Parameter "keep_temp" value is not a recognized boolean string: "{config["keep_temp"].strip()}" (try "True" or "False")')
    else:
        params['keep_temp'] = False

    # Return configured parameters
    return params


def read_sample_table(params):
    """
    Read the sample table and return as a Pandas DataFrame object with "SAMPLE" and "ASSEMBLY" fields.

    :param params: Parameters returned by get_params().
    """

    if 'sample_source' not in params:
        raise RuntimeError('Missing "sample_source" in parameters')

    if params['sample_source'] == 'table':

        # Check for required parameters
        if 'sample_table' not in params:
            raise RuntimeError('Missing "sample_table" in parameters when "sample_source" is "table"')

        if 'sample_table_format' not in params:
            raise RuntimeError('Missing "sample_table_format" in parameters when "sample_source" is "table"')

        # Read
        if params['sample_table_format'] == 'xlsx':
            df_sample = pd.read_excel(params['sample_table'])
        elif params['sample_table_format'] == 'tsv':
            df_sample = pd.read_csv(params['sample_table'], sep='\t')
        else:
            raise RuntimeError(f'Unknown "sample_table_format" value in parameters: {params["sample_table_format"]}')

        # CHeck columns
        missing_cols = [f'"{col}"' for col in ['ASSEMBLY', 'SAMPLE'] if col not in df_sample.columns]

        if missing_cols:
            raise RuntimeError(
                'Missing header columns in sample table: {}: table {}'.format(
                    ', '.join(missing_cols), params['sample_table']
                )
            )

        # Normalize missing values
        df_sample = df_sample[['SAMPLE', 'ASSEMBLY']].fillna('')
        df_sample['SAMPLE'] = df_sample['SAMPLE'].apply(lambda val: val.strip())
        df_sample['ASSEMBLY'] = df_sample['ASSEMBLY'].apply(lambda val: val.strip())

        # Check for empty table
        if df_sample.shape[0] == 0:
            raise RuntimeError(f'No records in sample table: {params["sample_table"]}')

        # Check for missing values
        df_sample = df_sample.loc[
            (df_sample['SAMPLE'] != '') | (df_sample['ASSEMBLY'] != '')
        ]

        df_sample_missing = df_sample.loc[(df_sample['SAMPLE'] == '') | (df_sample['ASSEMBLY'] == '')]

        if df_sample_missing.shape[0] > 0:
            n = df_sample_missing.shape[0]
            missing_index = ', '.join(sorted([str(val) for val in df_sample_missing.index])[:3])
            raise RuntimeError(f'Missing SAMPLE or ASSEMBLY values for {n} record(s): index={missing_index}{"..." if n > 3 else ""}')

        # Check for duplicate SAMPLE keys
        dup_sample = [sample for sample, count in collections.Counter(df_sample['SAMPLE']).items() if count > 1]

        if dup_sample:
            sample_list = ', '.join(sorted(dup_sample)[:3])
            raise RuntimeError(f'Sample table contains multiple records for the {len(dup_sample)} sample(s): {sample_list}{"..." if len(dup_sample) > 3 else ""}')

        # Warn about duplicate ASSEMBLY entries
        real_path = df_sample['ASSEMBLY'].apply(os.path.realpath)

        dup_asm = {
            asm for asm, count in collections.Counter(real_path).items() if count > 1
        }

        if dup_asm:
            df_sample_dup = df_sample.loc[list(real_path[real_path.isin(dup_asm)].index)]
            sample_list = list(df_sample_dup['SAMPLE'])

            raise RuntimeWarning(f'Found records pointing to the same {len(dup_asm)} input file(s) for sample(s): {", ".join(sample_list[:3])}{"..." if len(sample_list) > 3 else ""}')

    elif params['sample_source'] == 'config':

        # Check for required parameters
        if 'sample' not in params:
            raise RuntimeError('Missing "sample" in parameters when "sample_source" is "config"')

        if 'assembly' not in params:
            raise RuntimeError('Missing "assembly" in parameters when "sample_source" is "config"')

        # Build table
        df_sample = pd.DataFrame(
            pd.Series([params['sample'], params['assembly']], index=['SAMPLE', 'ASSEMBLY'])
        ).T

    else:
        raise RuntimeError(f'Unknown "sample_source" value in parameters: {params["sample_source"]}')

    # Return table
    return df_sample
