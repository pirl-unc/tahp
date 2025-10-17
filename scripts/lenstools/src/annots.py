"""
Functions related to annotation.
"""

import re

def extract_variant_tx_metadata(txs, gtf):
    """
    Extract relevant metadata (cds, strand, chromosome, start, stop, etc.)
    about a set of transcripts of interest.

    Args:
      txs (list): List of transcripts of interest
      gtf (string): Path to GTF annotation file.

    Returns:
      variant_txs_metadata (dict): Dictionary where keys are transcript
                                   identifiers and values are strand and
                                   coordinate information. 
    """
    # Dictionary for holding exon coding sequences coordinates for combining to
    # create transcript sequences.
    # ['cds'] contains a list with a set of exonic coordinates (chr:start-stop).
    # ['strand'] contains the strand.
    variant_txs_metadata = {}

    for expressed_tx in txs:
        variant_txs_metadata[expressed_tx] = {}
        variant_txs_metadata[expressed_tx]['cds'] = []
        variant_txs_metadata[expressed_tx]['strand'] = ''

    with open(gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            meta_entries = line.split('\t')[8].split('; ')
            for meta_entry_idx, meta_entry in enumerate(meta_entries):
                if re.search('transcript_id', meta_entry):
                    tx_id_idx = meta_entry_idx
                    break
            variant_tx = line.split('\t')[8] \
                             .split('; ')[tx_id_idx] \
                             .replace('"', '') \
                             .replace('transcript_id ', '')
            splitl = line.split('\t')
            chrom, start, stop, strand = (splitl[0], splitl[3], splitl[4], splitl[6])
            variant_txs_metadata[variant_tx]['strand'] = strand
            coords = "{}:{}-{}".format(chrom, start, stop)
            if coords not in variant_txs_metadata[variant_tx]['cds']:
                variant_txs_metadata[variant_tx]['cds'].append("{}:{}-{}".format(chrom, start, \
                                                                                 stop))
    print("Loaded variant transcripts metadata.")
    return variant_txs_metadata
