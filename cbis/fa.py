"""
FASTA utilities
"""

import svpoplib

import Bio.SeqIO


def fa_condensed_iter(seq_file_name, df_bed, format='fasta'):
    """
    Get an iterator over condensed sequence records.

    :param seq_file_name: Input sequence file name.
    :param df_bed: A BED file in sequence space indicating the masked loci to be condensed out (removed).
    :param format: File format.
    """

    with svpoplib.seq.PlainOrGzReader(seq_file_name, 'rt') as in_file:
        seq_id_set = set()

        for seq in Bio.SeqIO.parse(in_file, format):

            if seq.id in seq_id_set:
                raise RuntimeError(f'Multiple sequence entries for sequence with ID "{seq.id}"')

            seq_id_set.add(seq.id)

            subdf = df_bed.loc[df_bed['#CHROM'] == seq.id].reset_index()

            if subdf.shape[0] > 0:
                seq_list = [
                    str(seq[subdf.loc[index, 'END']:subdf.loc[index + 1, 'POS']].seq)
                        for index in range(subdf.shape[0] - 1)
                ]

                if subdf.iloc[-1]['END'] < len(seq):
                    seq_list.append(str(seq[subdf.iloc[-1]['END']:len(seq)].seq))

            else:
                seq_list = [str(seq)]

            seq = ''.join(seq_list)

            if seq:
                yield Bio.SeqRecord.SeqRecord(
                    Bio.Seq.Seq(seq), id=seq.id, name='', description=''
                )
