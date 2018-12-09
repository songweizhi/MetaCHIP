import sys


# https://github.com/Ecogenomics/GTDBTk/blob/cbe4146c5f019fe46b341c3863e7f077cc6a42f8/gtdbtk/markers.py
def _apply_mask(self, gtdb_msa, user_msa, msa_mask, min_perc_aa):
    """Apply canonical mask to MSA file."""

    aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)

    mask = open(msa_mask).readline().strip()

    if len(mask) != len(aligned_genomes.values()[0]):
        self.logger.error('Mask and alignment length do not match.')
        sys.exit()

    output_seqs = {}
    pruned_seqs = {}
    for seq_id, seq in aligned_genomes.iteritems():
        masked_seq = ''.join(
            [seq[i] for i in xrange(0, len(mask)) if mask[i] == '1'])

        valid_bases = len(masked_seq) - \
                      masked_seq.count('.') - masked_seq.count('-')
        if seq_id in user_msa and valid_bases < len(masked_seq) * min_perc_aa:
            pruned_seqs[seq_id] = masked_seq
            continue

        output_seqs[seq_id] = masked_seq

    return output_seqs, pruned_seqs