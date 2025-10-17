"""
Utilities for working with RNA-seq data within the context of LENS.
"""

import csv
import re
import numpy as np
import pysam

from Bio import SeqIO
from Bio.SeqUtils import Seq

def load_tx_abundances(quants, metric, exclude_zeros):
    """
    Given a quants file and metric of interest, return counts. Zero counts can
    be excluded if needed.

    Args:

    Returns:

    """
    count_column = ''
    counts = np.array([])
    with open(quants) as quants_fo:
        header = []
        cfo = csv.reader(quants_fo, delimiter='\t')
        for line_idx, line in enumerate(cfo):
            if line_idx == 0:
                header = line
                for i,j in enumerate(header):
                    if j == metric:
                        count_column = i
            else:
                count = np.log2(float(line[count_column]) + 1)
                if exclude_zeros:
                    if count > 0:
                        counts = np.append(counts, count)
                else:
                    counts = np.append(counts, count)
    return counts

def get_tx_threshold(percentile, counts):
    """
    Given a percentile value and a list of counts, return the corresponding
    count value.

    Args:
        percentile (int): Desired percentile
        counts (list): List of count data

    Return:
        threshold (float): Corresponding value threshold from counts
    """
    threshold = np.percentile(counts, int(percentile))
    return threshold

def get_expressed_txs(quants, metric, threshold):
    """
    Given an annotated VCF, quants file, and expression percentile (or
    threshold), emit SNV/InDel-harboring transcripts that exceed expression
    threshold. This code block is quite similar to load_tx_abundances, but
    keeps track of transcript. These two functions should be refactored.

    Args:

    Returns:
    """
    count_column = ''
    txid_column = ''
    expressed_txids = []
    with open(quants) as quants_fo:
        header = []
        cfo = csv.reader(quants_fo, delimiter='\t')
        for line_idx, line in enumerate(cfo):
            if line_idx == 0:
                header = line
                for i,j in enumerate(header):
                    if j ==  metric:
                        count_column = i
                    elif j == 'Name':
                        txid_column = i
            else:
                count = np.log2(float(line[count_column]) + 1)
                if float(count) >= float(threshold):
                    expressed_txids.append(line[txid_column])
    return expressed_txids


def write_expressed_txs(expressed_txids, outfile):
    """
    Given a list of transcript identifiers and an outfile path, write the
    transcript identifiers to that output file.
    """
    with open(outfile, 'w') as out_fo:
        for expressed_txid in expressed_txids:
            out_fo.write("{}\n".format(expressed_txid))


def filter_expressed_txs(quants, metric, outfile, exclude_zeros, percentile):
    """
    Given a quant file, metric of interest, output file, and desired
    percentile, filter quant file for transcripts with expression that meets or
    exceeds the percentile threshold.

    Args:

    Returns:
    """
    counts = load_tx_abundances(quants, metric, exclude_zeros)
    threshold = get_tx_threshold(percentile, counts)
    expressed_txids = get_expressed_txs(quants, metric, threshold)
    write_expressed_txs(expressed_txids, outfile)


def get_snv_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    md5m = {}
    tx_meta = {}

    for seq_record in SeqIO.parse(args.nt_fasta, "fasta"):
        seq_id = seq_record.id
        print("ID: {}".format(seq_id))
        print("Description: {}".format(seq_record.description))
#        if re.search('MD5', seq_id):
        seq_id = seq_record.description
        split_id = seq_id.split(' ')
        print(split_id)
        md5m[split_id[0].replace(':', '_')[:15]] = {}
        md5m[split_id[0].replace(':', '_')[:15]]['cdna_pos'] = split_id[4].replace('CDNA_POS:', '')
        md5m[split_id[0].replace(':', '_')[:15]]['tx'] = split_id[5].replace('TRANSCRIPT:', '')
        md5m[split_id[0].replace(':', '_')[:15]]['seq'] = seq_record.seq


    printable_lines = []

    txs = [md5m[x]['tx'] for x in md5m.keys()]

    # Loading tx-specific metadata from GTF.
    for transcript in txs:
        tx_meta[transcript] = {}
        tx_meta[transcript]['max_coord'] = 0
    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            print(line)
            transcript = line.split('\t')[8].split(';')[1].split(' ')[-1].replace('"', '')
            strand = line.split('\t')[6]
            max_coord = int(line.split('\t')[4])
            if transcript in tx_meta.keys():
                tx_meta[transcript]['strand'] = strand
                if max_coord > tx_meta[transcript]['max_coord']:
                    tx_meta[transcript]['max_coord'] = max_coord



#             elif re.search(transcript.split('.')[0], line):
#                    if re.search('\tCDS\t', line):
#                        strand = line.split('\t')[6]
#                        max_coord = int(line.split('\t')[4])
#                        tx_meta[transcript]['strand'] = strand
#                        if max_coord > tx_meta[transcript]['max_coord']:
#                            print("New highest coord...updating.")
#                            tx_meta[transcript]['max_coord'] = max_coord
#
    col_map = {}

    header = []
    print("Parsing pMHCs...")
    with open(args.netmhcpan) as netmhcpan_fo:
        for line_idx, line in enumerate(netmhcpan_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                header = line
                for col_idx, col in enumerate(header):
                    col_map[col] = col_idx
                header.extend(['rna_reads_covering_genomic_origin',
                               'rna_reads_covering_genomic_origin_with_peptide_cds',
                               'proportion_rna_reads_covering_genomic_origin_with_peptide_cds',
                               'primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds'])
            elif line[col_map['antigen_source']] == 'SNV':
                print(line)
                if line[col_map['antigen_source']] == 'SNV' and 'pos' in col_map.keys():
                    lo_coord = int(line[col_map['pos']]) * 3
                    hi_coord = lo_coord + (int(len(line[2]) * 3))
                    print("{}\t{}".format(lo_coord, hi_coord))
                    matched_id = md5m[line[col_map['identity']]]['tx']
                    peptide_nt_seq = str(md5m[line[col_map['identity']]]['seq']\
                                     [lo_coord-3:hi_coord-3])
                    print(peptide_nt_seq)
                    print(lo_coord)
                    print(hi_coord)
                    print(tx_meta[md5m[line[col_map['identity']]]['tx']])

                    lo_coord_gnm = max(int(md5m[line[col_map['identity']]]['cdna_pos']) - 30 + \
                                       lo_coord - 4, 1)
                    hi_coord_gnm = min(int(md5m[line[col_map['identity']]]['cdna_pos']) - 30 + \
                                       hi_coord - 4, \
                                   int(tx_meta[md5m[line[col_map['identity']]]['tx']]['max_coord']))

                    print("lo_coord_gnm:{}\thi_coord_gnm:{}".format(lo_coord_gnm, hi_coord_gnm))

                    contig_id = matched_id
                    contig_id_no_suffix = contig_id.split('.')[0]

                    bam_contig_id = contig_id
                    if (contig_id not in list(rna_bam.references) and
                        contig_id_no_suffix not in list(rna_bam.references)):
                        print("Contig not discovered in references.")
                        continue

                    if contig_id_no_suffix in list(rna_bam.references):
                        bam_contig_id = contig_id_no_suffix
                elif line[col_map['antigen_source']] == 'SNV' and 'pos' not in col_map.keys():
                    nt_seq = str(md5m[line[col_map['identity']]]['seq'])
                    peptide = line[col_map['peptide']]
                    pep_df = {}
                    for kmer_len in (24, 27, 30, 33):
                        for i in range(0,len(nt_seq)-kmer_len+3, 3):
                            coding_seq = nt_seq[i:i+kmer_len]
                            pep_seq = Seq(nt_seq[i:i+kmer_len]).translate()
                            if pep_seq in pep_df.keys():
                                print("Peptide is already in df, cannot add.")
                            pep_df[pep_seq] = {}
                            pep_df[pep_seq]['coding_seq'] = coding_seq
                            pep_df[pep_seq]['lo_coord'] = i
                            pep_df[pep_seq]['hi_coord'] = i + kmer_len

                    print(pep_df[peptide])
                    #Adding 3 to compensate subtracting 3 below
                    lo_coord = pep_df[peptide]['lo_coord'] + 3
                    hi_coord = pep_df[peptide]['hi_coord'] + 3
                    peptide_nt_seq = pep_df[peptide]['coding_seq']
                    matched_id = md5m[line[col_map['identity']]]['tx']
                    print(peptide_nt_seq)
                    print(lo_coord)
                    print(hi_coord)
                    print(tx_meta[md5m[line[col_map['identity']]]['tx']])

                    lo_coord_gnm = max(int(md5m[line[col_map['identity']]]['cdna_pos']) - 30 + \
                                       lo_coord - 4, 1)
                    hi_coord_gnm = min(int(md5m[line[col_map['identity']]]['cdna_pos']) - 30 + \
                                       hi_coord - 4, \
                                   int(tx_meta[md5m[line[col_map['identity']]]['tx']]['max_coord']))

                    print("lo_coord_gnm:{}\thi_coord_gnm:{}".format(lo_coord_gnm, hi_coord_gnm))

                covered_reads = []
                contig_id = matched_id
                contig_id_no_suffix = contig_id.split('.')[0]

                bam_contig_id = contig_id
                if (contig_id not in list(rna_bam.references) and
                    contig_id_no_suffix not in list(rna_bam.references)):
                    print("Contig not discovered in references.")
                    continue

                if contig_id_no_suffix in list(rna_bam.references):
                    bam_contig_id = contig_id_no_suffix


                print(bam_contig_id)
                all_reads = rna_bam.fetch(region=f"{bam_contig_id}:{lo_coord_gnm}-{hi_coord_gnm}")
                for read in all_reads:
                    if read.get_overlap(lo_coord_gnm, hi_coord_gnm) >= len(peptide_nt_seq):
                        covered_reads.append(read)
                print("Number of reads covering: {}".format(len(covered_reads)))

                prim_positive_hits = 0
                positive_hits = 0
                peptide_nt_seq_rc = str(Seq(peptide_nt_seq).reverse_complement())
                print(len(covered_reads))
                if len(covered_reads) > 0:
                    for read in covered_reads:
                        if (bool(re.search(peptide_nt_seq, read.query_sequence, re.IGNORECASE)) or
                            bool(re.search(peptide_nt_seq_rc, read.query_sequence, re.IGNORECASE))):
                            positive_hits += 1
                            if not(read.is_secondary):
                                prim_positive_hits += 1
                    print("Positive hits: {}".format(positive_hits))
                    print("Primary alignment positive hits: {}".format(prim_positive_hits))
                    prop_var = "N/A"
                    prop_var = positive_hits/(len(covered_reads) + 0.0)
                    if prop_var > 0.0:
                        line.extend([str(len(covered_reads)), str(positive_hits), str(prop_var), str(prim_positive_hits)])
                        printable_lines.append(line)

    if printable_lines:
        with open(args.output, 'w') as ofo:
            print(header)
            ofo.write('{}\n'.format('\t'.join(header)))
            for line in printable_lines:
                ofo.write('{}\n'.format('\t'.join(line)))


def get_indel_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    md5m = {}
    tx_meta = {}

    # Loading sequences from nucleotide FASTA
    for seq_record in SeqIO.parse(args.nt_fasta, "fasta"):
        seq_id = seq_record.id
#        if re.search('MD5', seq_id):
        seq_id = seq_record.description
        split_id = seq_id.split(' ')
        print(split_id)
        md5m[split_id[0].replace(':', '_')[:15]] = {}
        md5m[split_id[0].replace(':', '_')[:15]]['cdna_pos'] = split_id[3]\
                                                               .replace('CDNA_POS:', '')
        md5m[split_id[0].replace(':', '_')[:15]]['chr'] = split_id[2]\
                                                          .replace('VARIANT_POS:', '')\
                                                          .split(':')[0]
        md5m[split_id[0].replace(':', '_')[:15]]['var_pos'] = split_id[2]\
                                                              .replace('VARIANT_POS:', '')\
                                                              .split(':')[1]
        md5m[split_id[0].replace(':', '_')[:15]]['tx'] = split_id[4]\
                                                         .replace('TRANSCRIPT:', '')
        md5m[split_id[0].replace(':', '_')[:15]]['seq'] = seq_record.seq
        md5m[split_id[0].replace(':', '_')[:15]]['type'] = split_id[7].replace('INDEL_TYPE:', '')
        md5m[split_id[0].replace(':', '_')[:15]]['len'] = len(split_id[6].replace('ALT:', '')\
                                                          .replace('[','').replace(']',''))

    printable_lines = []

    txs = [md5m[x]['tx'] for x in md5m.keys()]

    # Loading tx-specific metadata from GTF.
    for transcript in txs:
        tx_meta[transcript] = {}
        tx_meta[transcript]['max_coord'] = 0
    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            print(line)
            transcript = line.split('\t')[8].split(';')[1].split(' ')[-1].replace('"', '')
            strand = line.split('\t')[6]
            max_coord = int(line.split('\t')[4])
            if transcript in tx_meta.keys():
                tx_meta[transcript]['strand'] = strand
                if max_coord > tx_meta[transcript]['max_coord']:
                    tx_meta[transcript]['max_coord'] = max_coord

#    # Loading tx-specific metadata from GTF.
#    for transcript in txs:
#        tx_meta[transcript] = {}
#        tx_meta[transcript]['max_coord'] = 0
#    with open(args.gtf) as gtf_fo:
#        for line in gtf_fo.readlines():
#            for transcript in txs:
#                if re.search(transcript, line):
#                    if re.search('\tCDS\t', line):
#                        strand = line.split('\t')[6]
#                        max_coord = int(line.split('\t')[4])
#                        tx_meta[transcript]['strand'] = strand
#                        if max_coord > tx_meta[transcript]['max_coord']:
#                            print("New highest coord...updating.")
#                            tx_meta[transcript]['max_coord'] = max_coord
#                elif re.search(transcript.split('.')[0], line):
#                    if re.search('\tCDS\t', line):
#                        strand = line.split('\t')[6]
#                        max_coord = int(line.split('\t')[4])
#                        tx_meta[transcript]['strand'] = strand
#                        if max_coord > tx_meta[transcript]['max_coord']:
#                            print("New highest coord...updating.")
#                            tx_meta[transcript]['max_coord'] = max_coord
#
#    print(tx_meta)

    col_map = {}

    header = []
    with open(args.netmhcpan) as netmhcpan_fo:
        for line_idx, line in enumerate(netmhcpan_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    col_map[col] = col_idx
                header = line
            elif line[col_map['antigen_source']] == 'INDEL':
                matched_id = md5m[line[col_map['identity']]]['chr']
                print("cdna pos: {}".format(md5m[line[col_map['identity']]]['cdna_pos']))
                print("indel len: {}".format(md5m[line[col_map['identity']]]['len']))
                print("contig: {}".format(matched_id))
                print("Seq: {}".format(md5m[line[col_map['identity']]]['seq']))
                print("peptide_nt_seq upperbound: {}"\
                       .format(34 + int(md5m[line[col_map['identity']]]['len'])))
                peptide_nt_seq = str(md5m[line[col_map['identity']]]['seq']\
                                 [28:(34 + int(md5m[line[col_map['identity']]]['len']))])

                print("Peptide NT Seq: {}".format(peptide_nt_seq))
                print("type: {}".format(md5m[line[col_map['identity']]]['type']))
                # Originally looking for a larger region than needed. Going
                # to look specifically around the variant position.
                lo_coord_gnm = max(int(md5m[line[col_map['identity']]]['var_pos']) - 3 - 4, 1)
                hi_coord_gnm = min(int(md5m[line[col_map['identity']]]['var_pos']) + 3 - 4, \
                                   int(tx_meta[md5m[line[col_map['identity']]]['tx']]['max_coord']))

                print("genome coords: {}-{}".format(lo_coord_gnm, hi_coord_gnm))

                covered_reads = []
                prim_covered_reads = []

                contig_id = matched_id
                contig_id_no_suffix = contig_id.split('.')[0]
                bam_contig_id = contig_id
                if (contig_id not in list(rna_bam.references) and
                    contig_id_no_suffix not in list(rna_bam.references)):
                    print("Contig not discovered in references.")
                    continue
                if contig_id_no_suffix in list(rna_bam.references):
                    bam_contig_id = contig_id_no_suffix
                for read in rna_bam.fetch(bam_contig_id, lo_coord_gnm, hi_coord_gnm):
                    covered_reads.append(read)
                for read in rna_bam.fetch(bam_contig_id, lo_coord_gnm, hi_coord_gnm):
                    if not(read.is_secondary):
                        prim_covered_reads.append(read)

                print("Number of full overlap reads: {}".format(len(covered_reads)))
                print("Number of primary full overlap reads: {}".format(len(prim_covered_reads)))

                prim_positive_hits = 0
                positive_hits = 0
                peptide_nt_seq_rc = str(Seq(peptide_nt_seq).reverse_complement())

                if len(covered_reads) > 0:
                    indel_reads = []
                    for read in covered_reads:
                        cigar = read.cigartuples
                        indel_cigar = [(x,y) for x,y in cigar if x in [1, 2]]
                        if indel_cigar:
                            indel_reads.append(read)
                    for indel_read in indel_reads:
                        if (bool(re.search(peptide_nt_seq, indel_read.query_sequence, \
                                           re.IGNORECASE))
                         or bool(re.search(peptide_nt_seq_rc, indel_read.query_sequence, \
                                           re.IGNORECASE))):
                            print("Read: {}".format(indel_read))
                            print("Read Seq: {}".format(indel_read.query_sequence))
                            print("peptide Seq: {}".format(peptide_nt_seq))
                            print("peptide Seq RC: {}".format(peptide_nt_seq_rc))
                            print("Read cigar: {}".format(indel_read.cigarstring))
                            positive_hits += 1
                            if not(indel_read.is_secondary):
                                prim_positive_hits += 1
                    print("InDel reads: {}".format(len(indel_reads)))
                    print("Positive hits: {}".format(positive_hits))
                    print("Primary alignment positive hits: {}".format(prim_positive_hits))
                    prop_var = positive_hits/(len(covered_reads) + 0.0)
                    if prop_var > 0:
                        line.extend([str(len(covered_reads)), str(positive_hits), str(prop_var), str(prim_positive_hits)])
                        printable_lines.append(line)
#                if len(prim_covered_reads) > 0:
#                    prim_indel_reads = []
#                    for read in prim_covered_reads:
#                        cigar = read.cigartuples
#                        indel_cigar = [(x,y) for x,y in cigar if x in [1, 2]]
#                        if indel_cigar:
#                            prim_indel_reads.append(read)
#                    for indel_read in prim_indel_reads:
#                        if (bool(re.search(peptide_nt_seq, indel_read.query_sequence, \
#                                           re.IGNORECASE))
#                         or bool(re.search(peptide_nt_seq_rc, indel_read.query_sequence, \
#                                           re.IGNORECASE))):
#                            print("Read: {}".format(indel_read))
#                            print("Read Seq: {}".format(indel_read.query_sequence))
#                            print("peptide Seq: {}".format(peptide_nt_seq))
#                            print("peptide Seq RC: {}".format(peptide_nt_seq_rc))
#                            print("Read cigar: {}".format(indel_read.cigarstring))
#                            prim_positive_hits += 1
#                    print("Primary InDel reads: {}".format(len(prim_indel_reads)))
#                    print("Primary positive hits: {}".format(prim_positive_hits))
#                    prop_var = prim_positive_hits/(len(covered_reads) + 0.0)
#                    if prop_var > 0:
#                        line.extend([str(prim_positive_hits)])
#                        printable_lines.append(line)
#    print(printable_lines)

    header.extend(['rna_reads_covering_genomic_origin',
                   'rna_reads_covering_genomic_origin_with_peptide_cds',
                   'proportion_rna_reads_covering_genomic_origin_with_peptide_cds',
                   'primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds'])

    if printable_lines:
        with open(args.output, 'w') as output_fo:
            print(header)
            output_fo.write('{}\n'.format('\t'.join(header)))
            for line in printable_lines:
                output_fo.write('{}\n'.format('\t'.join(line)))



def get_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    contig_seqs = {}
    md5_to_contig = {}
    tx_to_utr_buffer = {}

    for seq_record in SeqIO.parse(args.cds_fasta, "fasta"):
        seq_id = str(seq_record.id)
        print(seq_id)
        # Handling checksum identifiers
        seq_id = seq_record.description
        split_id = seq_id.split(' ')
        print(split_id)
        if len(split_id) > 1 and re.search('NAME', split_id[1]):
            md5_to_contig[split_id[0].replace(':', '_')[:15]] = split_id[1].replace('NAME:', '')
        elif len(split_id) > 1 and re.search('NAME', split_id[2]):
            md5_to_contig[split_id[0].replace(':', '_')[:15]] = split_id[2].replace('NAME:', '')
        else:
            # This is to compensate for viruses.
            viral_contig_id = '_'.join(split_id[0].split('_')[3:5]).replace('_1', '.1')\
                                       .replace('P.', 'P_')
            md5_to_contig[viral_contig_id] = split_id[0]
#            md5_to_contig[viral_contig_id] = '_'.join(split_id[0].split(':')[0].split('_')[:5])
        if re.search('UTR_BUFFER', seq_id):
            tx_to_utr_buffer[split_id[2].replace('NAME:','')] = split_id[3]\
                                                                .replace('UTR_BUFFER:', '')
            seq_id = split_id[0].replace(':', '_')
            contig_seqs[split_id[2].replace('NAME:', '')] = seq_record.seq
        # Handling human and mouse Ensembl transcript identifiers
        elif re.search('ENST|ENSMUST', seq_id):
            seq_id = seq_record.description
            split_id = seq_id.split()
            tx_to_utr_buffer[split_id[0]] = split_id[1].replace('UTR_BUFFER:', '')
            contig_seqs[split_id[0]] = seq_record.seq
        else:
            seq_id = seq_id.split(' ')[0]
#            seq_id = seq_id.split(' ')[0].split(':')[0]
#            seq_id = '_'.join(split_id[0].split(':')[0].split('_')[:5])
            contig_seqs[seq_id] = seq_record.seq

    printable_lines = []

    header = []
    header_map = {}
    with open(args.pmhcs) as pmhcs_fo:
        for line_idx, line in enumerate(pmhcs_fo.readlines()):
            line = line.rstrip().split()
            covered_reads = []
            prim_covered_reads = []
            if line_idx == 0:
                header = line
                for col_idx, col in enumerate(header):
                    header_map[col] = col_idx
            else:
                contig_id = ''
                if line[header_map['antigen_source']] == args.antigen_source and \
                        'pos' in header_map.keys():
                    print("Matched antigen source. 'pos' in header.")
                    pos = int(line[header_map['pos']])
                    peptide = line[header_map['peptide']]
                    lo_coord = pos * 3
                    hi_coord = lo_coord + (int(len(line[header_map['peptide']]) * 3))
                    matched_id = line[header_map['identity']]
                    contig_id = ''
                    if line[header_map['identity']] in md5_to_contig.keys():
                        contig_id = md5_to_contig[line[header_map['identity']]]
                    elif args.antigen_source == 'VIRUS':
                        print(line[header_map['virus_id']])
                        print(md5_to_contig.keys())
                        try:
                            contig_id = md5_to_contig[line[header_map['virus_id']]\
                                                                      .partition(':')[0]]
                        except:
                            contig_id = md5_to_contig['_'.join(line[header_map['virus_id']]\
                                                                                .split('_')[:2])]
                    print(f"Contig ID: {contig_id}")
                    print("{}\t{}".format(lo_coord, hi_coord))
                    print(peptide)
                    print(f"Matched ID: {matched_id}")
                    peptide_nt_seq = ''
                    if contig_id in contig_seqs.keys():
                        peptide_nt_seq = str(contig_seqs[contig_id][lo_coord-3:hi_coord-3])
                    elif matched_id in contig_seqs.keys():
                        peptide_nt_seq = str(contig_seqs[matched_id][lo_coord-3:hi_coord-3])
                    print(contig_id)
                    print("{}-{}".format(lo_coord-3, hi_coord-3))
                    print(f"Peptide NT sequence: {peptide_nt_seq}")

                    print(f"Contig ID: {contig_id}")

                    # Adjusting coordinates based on upstream UTRs. UTRs are
                    # included in the coordinates of transcript alignments, but
                    # not in the coordinates of peptide generation (as the
                    # first base of the first start coding is the first
                    # position from netMHCpan results.
                    if re.search('Mmus|Hsap', contig_id) and contig_id in tx_to_utr_buffer.keys():
                        # This buffer shouldn't be needed for gEVE references
                        # since peptides are generated downstream of first
                        # Methionine.
                        print("Adding UTR buffer...")
                        lo_coord = lo_coord + int(tx_to_utr_buffer[matched_id])
                        hi_coord = hi_coord + int(tx_to_utr_buffer[matched_id])
                        print("new coords: {}\t{}".format(lo_coord, hi_coord))
                    if re.match('ENST|ENSMUST', contig_id):
                        lo_coord = lo_coord + int(tx_to_utr_buffer[contig_id])
                        hi_coord = hi_coord + int(tx_to_utr_buffer[contig_id])

                    print("orig:")
                    print(peptide)
                    print(lo_coord)
                    print(hi_coord)
                    print(peptide_nt_seq)

                elif line[header_map['antigen_source']] == args.antigen_source and \
                        'pos' not in header_map.keys():
                    print("Matched antigen source. No 'pos' in header.")
                    peptide = line[header_map['peptide']]
                    matched_id = line[header_map['identity']]
                    contig_id = ''
                    if line[header_map['identity']] in md5_to_contig.keys():
                        contig_id = md5_to_contig[line[header_map['identity']]]
                    elif args.antigen_source == 'VIRUS':
                        print(line[header_map['virus_id']])
                        print(md5_to_contig.keys())
                        try:
                            contig_id = md5_to_contig[line[header_map['virus_id']]\
                                                                      .partition(':')[0]]
                        except:
                            contig_id = md5_to_contig['_'.join(line[header_map['virus_id']]\
                                                                               .split('_')[:2])]
                    print(f"Contig ID: {contig_id}")
                    print(peptide)
                    print(f"Matched ID: {matched_id}")
                    if contig_id in contig_seqs.keys():
                        nt_seq = str(contig_seqs[contig_id])
                    elif matched_id in contig_seqs.keys():
                        nt_seq = str(contig_seqs[matched_id])
                    print(nt_seq)

                    pep_df = {}
                    for kmer_len in (24, 27, 30, 33):
                        for i in range(0,len(nt_seq)-kmer_len+3, 3):
                            coding_seq = nt_seq[i:i+kmer_len]
                            pep_seq = Seq(nt_seq[i:i+kmer_len]).translate()
                            if pep_seq in pep_df.keys():
                                print("Peptide is already in df, cannot add.")
                            pep_df[pep_seq] = {}
                            pep_df[pep_seq]['coding_seq'] = coding_seq
                            pep_df[pep_seq]['lo_coord'] = i
                            pep_df[pep_seq]['hi_coord'] = i + kmer_len

                    print(pep_df[peptide])
                    #Adding 3 to compensate subtracting 3 below
                    lo_coord = pep_df[peptide]['lo_coord'] + 3
                    hi_coord = pep_df[peptide]['hi_coord'] + 3
                    peptide_nt_seq = pep_df[peptide]['coding_seq']

                    if re.search('Mmus|Hsap', contig_id) and contig_id in tx_to_utr_buffer.keys():
                        # This buffer shouldn't be needed for gEVE references
                        # since peptides are generated downstream of first
                        # Methionine.
                        print("Adding UTR buffer...")
                        lo_coord = lo_coord + int(tx_to_utr_buffer[matched_id])
                        hi_coord = hi_coord + int(tx_to_utr_buffer[matched_id])
                        print("new coords: {}\t{}".format(lo_coord, hi_coord))
                    if re.match('ENST|ENSMUST', contig_id):
                        lo_coord = lo_coord + int(tx_to_utr_buffer[contig_id])
                        hi_coord = hi_coord + int(tx_to_utr_buffer[contig_id])

                    print("no_pos:")
                    print(peptide)
                    print(lo_coord)
                    print(hi_coord)
                    print(peptide_nt_seq)


                # Getting reads that map to coordinates of peptide's genomic origin.
                if contig_id:
                    for read in rna_bam.fetch(contig_id.split(':')[0], lo_coord-3, hi_coord-3):
                        if read.get_overlap(lo_coord-3, hi_coord-3) == len(peptide_nt_seq):
                            covered_reads.append(read.query_sequence)
                        if read.get_overlap(lo_coord-3, hi_coord-3) == len(peptide_nt_seq) and not(read.is_secondary):
                            prim_covered_reads.append(read.query_sequence) 
                    print("Number of full overlap reads: {}".format(len(covered_reads)))
                    print("Number of primary full overlap reads: {}".format(len(prim_covered_reads)))

                    positive_hits = 0
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            if bool(re.search(peptide_nt_seq, read, re.IGNORECASE)):
                                positive_hits += 1
                        print("Positive hits: {}".format(positive_hits))
                        prop_var = positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0:
                            line.extend([str(len(covered_reads)), str(positive_hits), \
                                         str(prop_var)])
                    prim_positive_hits = 0
                    if len(prim_covered_reads) > 0:
                        for read in prim_covered_reads:
                            if bool(re.search(peptide_nt_seq, read, re.IGNORECASE)):
                                prim_positive_hits += 1
                        print("Primary positive hits: {}".format(prim_positive_hits))
                        prop_var = prim_positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0:
                            line.extend([str(prim_positive_hits)])
                            printable_lines.append(line)
    print(printable_lines)

    header.extend(['rna_reads_covering_genomic_origin',
                   'rna_reads_covering_genomic_origin_with_peptide_cds',
                   'proportion_rna_reads_covering_genomic_origin_with_peptide_cds',
                   'primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds'])

    if printable_lines:
        with open(args.output, 'w') as output_fo:
            print(header)
            output_fo.write('{}\n'.format('\t'.join(header)))
            for line in printable_lines:
                output_fo.write('{}\n'.format('\t'.join(line)))


def get_fusion_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """

    contig_seqs = {}

    for seq_record in SeqIO.parse(args.nt_fasta, "fasta"):
        # The _NUCLEOTIDE addition is required to match to the NETMHCpan outputs.
        seq_id = '{}'.format(seq_record.description[:15].replace(' ', '_').replace('.', '_'))
        contig_seqs[seq_id] = seq_record.seq

    reads_and_seqs = {}
    for seq_record in SeqIO.parse(args.fusion_reads, "fastq"):
        new_id = seq_record.id

        # FASTQ read label sanitization
        if seq_record.id.endswith('/1'):
            new_id = seq_record.id.strip('/1')
        elif seq_record.id.endswith('/2'):
            new_id = seq_record.id.strip('/2')

        if new_id not in reads_and_seqs.keys():
            reads_and_seqs[new_id] = [str(seq_record.seq)]
        else:
            reads_and_seqs[new_id].append(str(seq_record.seq))


    reads_and_seqs_keys = sorted(list(reads_and_seqs.keys()))

    print("reads_and_seqs head: {}".format(reads_and_seqs_keys[:10]))

    fusion_col_to_idx = {}

    fusion_reads = {}
    with open(args.fusions) as fusions_fo:
        for line_idx, line in enumerate(fusions_fo.readlines()):
            if line_idx == 0:
                line = line.split()
                for idx, i in enumerate(line):
                    fusion_col_to_idx[i] = idx
            elif line_idx > 0:
                print(line)
                line = line.split()
                fusion_id = line[0].replace(' ', '_').replace('.', '_')
                reads = []
                reads_col_idxs = []
                for reads_col_names in ["JunctionReads", "SpanningFrags"]:
                    reads_col_idxs.append(fusion_col_to_idx["JunctionReads"])
                for reads_col_idx in reads_col_idxs:
                    reads.extend(['{}'.format(x) for x in line[reads_col_idx].split(',') if \
                                 line[reads_col_idx] != '.'])
                fixed_read_names = []
                for read in reads:
                    fixed_read = read
                    if '@' in read:
                        fixed_read = '{}'.format(read.partition('@')[2])
                    fixed_read_names.append(fixed_read)
                fusion_reads[fusion_id] = fixed_read_names

    printable_lines = []

    header = []
    header_map = {}
    with open(args.netmhcpan) as netmhcpan_fo:
        for line_idx, line in enumerate(netmhcpan_fo.readlines()):
            line = line.rstrip().split()
            if line_idx == 0:
                header = line
                for col_idx, col in enumerate(header):
                    header_map[col] = col_idx
            else:
                matched_id = ''
                if line[header_map['antigen_source']] == 'FUSION' and 'pos' in header_map.keys():
                    lo_coord = int(line[header_map['pos']]) * 3 - 1
                    hi_coord = lo_coord + (int(len(line[header_map['peptide']])) * 3)
                    matched_id = line[header_map['fusion_id']].replace('.', '_')
                    contig_id = line[header_map['identity']].replace('.', '_')

                    peptide_nt_seq = str(contig_seqs[contig_id][lo_coord-3:hi_coord-3]).upper()
                elif line[header_map['antigen_source']] == 'FUSION' and \
                        'pos' not in header_map.keys():
                    peptide = line[header_map['peptide']]
                    matched_id = line[header_map['fusion_id']].replace('.', '_')
                    contig_id = line[header_map['identity']].replace('.', '_')
                    print(f"Contig ID: {contig_id}")
                    print(peptide)
                    print(f"Matched ID: {matched_id}")
                    nt_seq = line[header_map['nt_context']]
                    print(nt_seq)

                    pep_df = {}
                    for kmer_len in (24, 27, 30, 33):
                        for start in (0, 1, 2):
                            for i in range(start,len(nt_seq)-kmer_len+3, 3):
                                coding_seq = nt_seq[i:i+kmer_len]
                                pep_seq = Seq(nt_seq[i:i+kmer_len]).translate()
                                if pep_seq in pep_df.keys():
                                    print("Peptide is already in df, cannot add.")
                                pep_df[pep_seq] = {}
                                pep_df[pep_seq]['coding_seq'] = coding_seq
                                pep_df[pep_seq]['lo_coord'] = i
                                pep_df[pep_seq]['hi_coord'] = i + kmer_len
 
#                    print(pep_df.keys())
                    print(nt_seq)

                    print(pep_df[peptide])
                    #Adding 3 to compensate subtracting 3 below
                    lo_coord = pep_df[peptide]['lo_coord'] + 3
                    hi_coord = pep_df[peptide]['hi_coord'] + 3
                    peptide_nt_seq = pep_df[peptide]['coding_seq'].upper()

                if matched_id:
                    reads_with_peptide_nt = 0
                    potential_reads = fusion_reads[matched_id]
                    print(len(potential_reads))

                    for read in potential_reads:
                        reads = list(set(reads_and_seqs[read]))
                        for subread in reads:
                            print(peptide_nt_seq)
                            print(subread)
                            if (re.search(peptide_nt_seq, subread) or
                                re.search(peptide_nt_seq, str(Seq(subread).reverse_complement()))):
                                reads_with_peptide_nt += 1
                    if reads_with_peptide_nt >0:
                        line.extend(['NA', str(reads_with_peptide_nt), 'NA', str(reads_with_peptide_nt)])
                        printable_lines.append(line)

    header.extend(['rna_reads_covering_genomic_origin',
                   'rna_reads_covering_genomic_origin_with_peptide_cds',
                   'proportion_rna_reads_covering_genomic_origin_with_peptide_cds',
                   'primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds'])
    if printable_lines:
        print("printable lines")
        with open(args.output, 'w') as output_fo:
            output_fo.write('{}\n'.format('\t'.join(header)))
            for line in printable_lines:
                output_fo.write('{}\n'.format('\t'.join(line)))
    else:
        print("No printable lines.")


def get_splice_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    header = ''
    col_map = {}

    printable_lines = []

    with open(args.neosplice_summary) as neosplice_summ_fo:
        for line_idx, line in enumerate(neosplice_summ_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                header = line
                for col_idx, col in enumerate(header):
                    col_map[col] = col_idx
                header.extend(['rna_reads_covering_genomic_origin',
                               'rna_reads_covering_genomic_origin_with_peptide_cds',
                               'proportion_rna_reads_covering_genomic_origin_with_peptide_cds',
                               'primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds'])
            elif line[col_map['antigen_source']] == 'SPLICE':
                dna_seq = str(line[col_map['coding_sequence']])
                if not dna_seq:
                    continue
                chrom = line[col_map['splice_coords']].split(':')[0]
                if chrom in ['chrMT']:
                    continue
                strand = line[col_map['splice_coords']].split('(')[1].replace(')', '')
                coord_lo = int(line[col_map['splice_coords']].split(':')[1].split('-')[0])
                coord_hi = int(line[col_map['splice_coords']].split(':')[1].split('(')[0].split('-')[1])
                if coord_lo > coord_hi:
#                if strand == '-':
                    coord_lo, coord_hi = coord_hi, coord_lo
                prim_covered_reads = []
                covered_reads = []
                prim_positive_hits = 0
                positive_hits = 0
                print("{} {} {} {} {}".format(dna_seq, chrom, strand, coord_lo, coord_hi))

                for read in rna_bam.fetch(chrom, coord_lo-150, coord_lo):
                    # greater than or equal to instead of equal to due to
                    # looking at a larger range than actual DNA sequence
                    # Note: This condition breaks splice counting since the
                    # splice region contains the entire splice region (e.g.
                    # intron) and is _more_ than the peptide's coding sequence.
                    # Can try reducing this range to just the peptide in the
                    # future and see if that improves quantification.
#                    if read.get_overlap(coord_lo, coord_hi) >= len(dna_seq):
                    covered_reads.append(read.query_sequence)
                for read in rna_bam.fetch(chrom, coord_hi, coord_hi+150):
                    # greater than or equal to instead of equal to due to
                    # looking at a larger range than actual DNA sequence
                    # Note: This condition breaks splice counting since the
                    # splice region contains the entire splice region (e.g.
                    # intron) and is _more_ than the peptide's coding sequence.
                    # Can try reducing this range to just the peptide in the
                    # future and see if that improves quantification.
#                    if read.get_overlap(coord_lo, coord_hi) >= len(dna_seq):
                    covered_reads.append(read.query_sequence)

                for read in rna_bam.fetch(chrom, coord_lo-150, coord_lo):
                    # greater than or equal to instead of equal to due to
                    # looking at a larger range than actual DNA sequence
                    # Note: This condition breaks splice counting since the
                    # splice region contains the entire splice region (e.g.
                    # intron) and is _more_ than the peptide's coding sequence.
                    # Can try reducing this range to just the peptide in the
                    # future and see if that improves quantification.
#                    if read.get_overlap(coord_lo, coord_hi) >= len(dna_seq):
                    if not(read.is_secondary):
                        prim_covered_reads.append(read.query_sequence)
                for read in rna_bam.fetch(chrom, coord_hi, coord_hi+150):
                    # greater than or equal to instead of equal to due to
                    # looking at a larger range than actual DNA sequence
                    # Note: This condition breaks splice counting since the
                    # splice region contains the entire splice region (e.g.
                    # intron) and is _more_ than the peptide's coding sequence.
                    # Can try reducing this range to just the peptide in the
                    # future and see if that improves quantification.
#                    if read.get_overlap(coord_lo, coord_hi) >= len(dna_seq):
                    if not(read.is_secondary):
                        prim_covered_reads.append(read.query_sequence)

                print("Covered reads: {}".format(len(covered_reads)))
                if len(covered_reads) > 0:
                    for read in covered_reads:
                        if (bool(re.search(dna_seq, read)) or
                            bool(re.search(dna_seq, str(Seq(read).reverse_complement())))):
                            print("Peptide NT: {}".format(dna_seq))
                            print("Read: {}".format(read))
                            print(re.search(dna_seq, read))
                            positive_hits += 1
                    print("Positive hits: {}".format(positive_hits))
                    prop_var = positive_hits/(len(covered_reads) + 0.0)
                    if prop_var > 0.01:
                        print(positive_hits)
                        print(prop_var)
                if positive_hits > 0:
                    line.extend(['NA', str(positive_hits), 'NA'])
                print("Primary covered reads: {}".format(len(prim_covered_reads)))
                if len(prim_covered_reads) > 0:
                    for read in covered_reads:
                        if (bool(re.search(dna_seq, read)) or
                            bool(re.search(dna_seq, str(Seq(read).reverse_complement())))):
                            print("Peptide NT: {}".format(dna_seq))
                            print("Read: {}".format(read))
                            print(re.search(dna_seq, read))
                            prim_positive_hits += 1
                    print("Primary positive hits: {}".format(prim_positive_hits))
                    prop_var = prim_positive_hits/(len(covered_reads) + 0.0)
                    if prop_var > 0.01:
                        print(prim_positive_hits)
                if prim_positive_hits > 0:
                    line.extend([str(positive_hits)])
                    printable_lines.append(line)

    if printable_lines:
        with open(args.output, 'w') as output_fo:
            print(header)
            output_fo.write('{}\n'.format('\t'.join(header)))
            for line in printable_lines:
                output_fo.write('{}\n'.format('\t'.join(line)))


def filter_expressed_ervs(args):
    """
    Filters quant.sf file for expressed ERVs. Useful for cases where
    tissue-matched sample for differential ERV expression is not available.

    Currently, ERVs are filtered by having a 'Hsap' prefix (for human ERVs) or
    'Mmus' (for mouse ERVs). This naming convention is derived from the gEVE
    database ERV names.

    Args:

    Returns:
    """
    tpms = {}
    quant_map = {}
    all_quant_meta = {}

    with open(args.quants) as quants_fo:
        for line_idx, line in enumerate(quants_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    quant_map[col] = col_idx
            else:
                name = line[quant_map['Name']]
                if re.search('[Mmus|Hsap].*[\+\-]$', name):
                    tpms[name] = line[quant_map['TPM']]
                    all_quant_meta[name] = ','.join(line[1:])

    tx_abundances = [float(x) for x in tpms.values()]

    if not args.abundance_threshold:
        print("No abundance threshold!")
        tx_threshold = np.mean(tx_abundances) + (5*np.std(tx_abundances))
        print("Transcript threshold (mean + 5*std): {}".format(tx_threshold))
    else:
        print("Abundance threshold: {}".format(args.abundance_threshold))
        tx_threshold = float(args.abundance_threshold)
    print("tx_thresshold: {}".format(tx_threshold))

    output_header = 'Name,Tumor_CPM,Norm_CPM,log2(Tumor_CPM+1)-log2(Norm_CPM+1),Length,'\
                    'EffectiveLength,TPM,NumReads\n'
    with open(args.output, 'w') as ofo:
        ofo.write(output_header)
        for herv, tpm in tpms.items():
            if float(tpm) > float(tx_threshold):
                ofo.write("{},NA,NA,NA,{}\n".format(herv, all_quant_meta[herv]))


def expressed_self_genes(args):
    """
    Filters expressed self-antigen/CTA transcripts. Given a gene list and
    either a expression percentile or specific threshold, map genes to their
    corresponding transcripts, and filter those transcripts by expression.

    Args:

    Returns:
    """
    expressed_selfs = []
    gene_list = []
    with open(args.gene_list) as glo:
        for line in glo.readlines():
            gene_list.append(line.rstrip().split('\t')[0])
    tx_abundances = load_tx_abundances(args.quants, args.metric, args.exclude_zeros)
    tx_threshold = 0
    if not args.abundance_threshold:
        tx_threshold = get_tx_threshold(args.percentile, tx_abundances)
    else:
        tx_threshold = float(args.abundance_threshold)
    expressed_txids = get_expressed_txs(args.quants, args.metric, tx_threshold)
    print(expressed_txids[:10])
    print("tx_thresshold: {}".format(tx_threshold))
    print("# of expressed transcripts: {}".format(len(expressed_txids)))
    print("Some expressed transcripts: {}".format(expressed_txids[:10]))
    print("Some genes: {}".format(gene_list[:10]))


    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if (len(line) > 3 and
                line[2] == 'transcript' and
                'gene_name' in str(line) and
                (re.search('transcript_type "protein_coding"', '\t'.join(line)) or
                 re.search('transcript_biotype "protein_coding"', '\t'.join(line)))):
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0]
                tx_to_gene[tx_id] = gene_name

    for tx, gene in tx_to_gene.items():
        if gene in gene_list and tx in expressed_txids:
            print(gene, tx)
            expressed_selfs.append("{}:{}".format(gene, tx))

    if expressed_selfs:
        with open(args.output, 'w') as ofo:
            for expressed_self in expressed_selfs:
                ofo.write("{}\n".format(expressed_self))


def filter_virdetect_by_counts(args):
    """
    Filters VirDetect output by counts. Fitlering occurs by either returning
    viruses that exceed a specific threshold (if min_threshold provided) or by
    returning viruses that exceed a mean + (20 * standard deviation) threshold.

    Args:

    Returns:

    TODO: Need a check here to ensure the raw_counts and raw_ref lists contain the
    same number of elements.
    """
    raw_counts = []
    raw_refs = []
    with open(args.viral_quants) as vco:
        raw_counts = vco.readlines()[0].rstrip().split()[1:]
    with open(args.viral_ref) as vro:
        for line in vro.readlines():
            if line.startswith('>'):
                line = line.rstrip('\n').lstrip('>')
                raw_refs.append(line.split(' ')[0])
    all_viral_counts = {raw_refs[i]:int(raw_counts[i]) for i in range(len(raw_refs))}

    threshold = ''
    if not args.min_threshold:
        counts = [float(x) for x in all_viral_counts.values()]
        threshold = np.mean(counts) + (20*np.std(counts))
        print("Count threshold (mean + 20*std): {}".format(threshold))
    else:
        threshold= args.min_threshold
    expressed_viruses = {k:v for k, v in all_viral_counts.items() if int(v) > int(threshold)}
    if expressed_viruses:
        with open(args.output, 'w') as ofo:
            ofo.write("virus_ref\tread_count\n")
            for k, v in expressed_viruses.items():
                ofo.write("{}\t{}\n".format(k, v))


def check_depth_and_coverage(targets_dict, cov_file, coverage_thresh, mean_depth_thresh, cov_type):
    """
    """
    rna_covered = []
    col_map = {}
    with open(cov_file) as cov_fo:
        for line_idx, line in enumerate(cov_fo.readlines()):
            if line_idx == 0:
                line = line.strip().split('\t')
                for col_idx, col in enumerate(line):
                    col_map[col] = col_idx
            else:
                line = line.strip().split('\t')
                target = line[col_map['#rname']]
                if target in targets_dict.keys():
                    print(line)
                    mean_depth = line[col_map['meandepth']]
                    coverage = line[col_map['coverage']]
                    print("depth {}	{}".format(mean_depth, mean_depth_thresh)) # delete
                    print("coverage {}	{}".format(coverage, coverage_thresh)) # delete
                    if (float(coverage) > float(coverage_thresh) and
                        float(mean_depth) > float(mean_depth_thresh)):
                        if cov_type == 'viral':
                            rna_covered.append("{}\t{}".format(target, line[col_map['numreads']]))
                        elif cov_type == 'erv':
                            rna_covered.append("{},{}".format(target, targets_dict[target]))
    return rna_covered


def check_erv_rna_coverage(args):
    """
    """
    expd_ervs = {}
    out_header = ''

    # Loading set of expressed ERVs
    with open(args.expressed_ervs) as expd_ervs_fo:
        for line_idx, line in enumerate(expd_ervs_fo.readlines()):
            if line_idx == 0:
                out_header = line
            else:
                line = line.rstrip().split(',')
                expd_ervs[line[0]] = ','.join(line[1:])

    rna_covered_ervs = check_depth_and_coverage(expd_ervs, args.coverage_file,
                                                args.coverage, args.mean_depth, 'erv')

    if rna_covered_ervs:
        with open(args.output, 'w') as output_fo:
            output_fo.write("{}".format(out_header))
            for i in rna_covered_ervs:
                output_fo.write("{}\n".format(i))


def check_virus_rna_coverage(args):
    """
    """
    expressed_viruses = {}
    with open(args.expressed_viruses) as exp_vir_fo:
        for line in exp_vir_fo.readlines():
            expressed_viruses[line.rstrip().split('\t')[0]] = ''

    rna_covered_viruses = check_depth_and_coverage(expressed_viruses, args.coverage_file,
                                                   args.coverage, args.mean_depth, 'viral')

    if rna_covered_viruses:
        with open(args.output, 'w') as output_fo:
            for i in rna_covered_viruses:
                output_fo.write("{}\n".format(i))
