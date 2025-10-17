import sys

import re
from glob import glob
from math import log2
import os
import hashlib
from pprint import pprint

from Bio import SeqIO
from Bio import Align
from Bio.SeqUtils import Seq
import pandas as pd

from src import var_utils
from src import parsers
from src import seq_utils
from src import annots


def make_snv_peptides_context(args):
    """
    Creates tumor and normal sequences from missense SNVs

    - Extract somatic missense SNVS from annotated VCF
    - Load transcripts harboring missense SNVs and gather metadata from GTF

    Args:

    Returns:

    """
    transcript_protein_seqs = parsers.parse_fasta(args.pep_ref)

    mutant_peptides = {}
    mutant_seqs = {}
    wildtype_peptides = {}

    print("Loading missense SNVs...")
    missense_snvs = var_utils.extract_missense_snvs(args.somatic_vcf)
    print("Loaded missense SNVs.")

    somatic_txs = []

    print("Loading transcripts harboring somatic variants...")
    with open(args.somatic_txs) as somatic_txs_fo:
        for line in somatic_txs_fo.readlines():
            somatic_txs.append(line.strip().replace('"', ''))
    print("Loaded transcripts harboring somatic variants.")

    print("Loading transcript metadata...")
    variant_txs_metadata = {}

    for variant_tx in somatic_txs:
        variant_txs_metadata[variant_tx] = {}
        variant_txs_metadata[variant_tx]['cds'] = []
        variant_txs_metadata[variant_tx]['strand'] = ''

    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            tx_id_idx = ''
            meta_entries = line.split('\t')[8].split('; ')
            for meta_entry_idx, meta_entry in enumerate(meta_entries):
                if re.search('transcript_id ', meta_entry):
                    tx_id_idx = meta_entry_idx
                    break
            variant_tx = line.split('\t')[8].split('; ')[tx_id_idx]\
                             .replace('"', '')\
                             .replace('transcript_id ', '')
            chrom = line.split('\t')[0]
            start = line.split('\t')[3]
            stop = line.split('\t')[4]
            strand = line.split('\t')[6]
            variant_txs_metadata[variant_tx]['strand'] = strand
            coords = "{}:{}-{}".format(chrom, start, stop)
            if coords not in variant_txs_metadata[variant_tx]['cds']:
                variant_txs_metadata[variant_tx]['cds']\
                    .append("{}:{}-{}".format(chrom, start, stop))
    print("Loaded variant transcript metadata.")


    for entry in missense_snvs.keys():
        for record in missense_snvs[entry]:
            variant_tx = record['transcript']
            if variant_tx not in somatic_txs or variant_tx not in transcript_protein_seqs.keys():
                continue
            norm_nts = []
            norm_aas = []
            tumor_nts = []
            tumor_aas = []
            print("Creating peptides for record: {}".format(record))
            record_coords = "{}_{}".format(record['meta'].CHROM, record['meta'].POS)
            print("Record coordinates: {}".format(record_coords))
            print("Record transcript: {}".format(record['transcript']))

            exon_seqs = {}
            exon_seqs['norm'] = {}
            exon_seqs['tumor'] = {}
            globbed_normal_fa = glob(os.path.join(args.var_tx_seqs,
                                     '*{}_{}.normal.fa'.format(variant_tx, record_coords)))
            globbed_tumor_fa = glob(os.path.join(args.var_tx_seqs,
                                    '*{}_{}.tumor.fa'.format(variant_tx, record_coords)))
            if globbed_normal_fa and globbed_tumor_fa:
                print("Loading normal exonic sequences...")
                norm_exon_count = 0
                for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs,
                                              '*{}_{}.normal.fa'.format(variant_tx, record_coords)
                                              ))[0],
                                              "fasta"):
                    exon_seqs['norm'][seq_record.description] = seq_record.seq
                    norm_exon_count += 1
                print("Loaded {} normal exonic sequences for transcript {}."\
                      .format(norm_exon_count, variant_tx))

                print("Loading tumor exonic sequences...")
                tum_exon_count = 0
                for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs,
                                              '*{}_{}.tumor.fa'.format(variant_tx, record_coords)
                                              ))[0],
                                              "fasta"):
                    exon_seqs['tumor'][seq_record.description] = seq_record.seq
                    tum_exon_count += 1
                print("Loaded {} tumor exonic sequences for transcript {}."\
                      .format(tum_exon_count, variant_tx))

            else:
                print(("Either normal or tumor exonic FASTAs for "
                       "transcript {} variant {} could not be found."
                       .format(variant_tx, record_coords)))
                continue

            variant_tx_seqs = {}

            samp_types = ['norm', 'tumor']
            for i in samp_types:
                variant_tx_seqs[i] = ''
            if variant_txs_metadata[variant_tx]['strand'] == '+':
                print("Transcript {} is on positive strand. Creating transcript sequences now..."\
                      .format(variant_tx))
                for i in samp_types:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
                        variant_tx_seqs[i] += exon_seqs[i][cds]
                    print("Created {} nucleotide sequence for transcript {}".format(i, variant_tx))
            elif variant_txs_metadata[variant_tx]['strand'] == '-':
                print("Transcript {} is on negative strand. Creating transcript sequences now..."\
                      .format(variant_tx))
                for i in samp_types:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
                        variant_tx_seqs[i] += exon_seqs[i][cds].reverse_complement()
                    print("Created {} nucleotide sequence for transcript {}".format(i, variant_tx))


            # Both tumor and normal transcript sequences created at this point.

            snv_pos = min(record['codon_pos'], ((int(args.length) - 1) * 3))
            u_offset = 0
            d_offset = 2
            # Mutations can occur in any position of an amino acid's codon.
            # This code is attempting to determine the correct reading frame
            # prior to translation.
            if (float(record['nt_pos']) - 1) % 3 == 0:
                print("Direct translate")
            if (float(record['nt_pos']) - 1) % 3 == 1:
                print("Minus one translate")
                snv_pos += 1
            elif (float(record['nt_pos']) - 1) % 3 == 2:
                print("Minus two translate")
                snv_pos += 2

            for i in samp_types:
                target_seq = ''
                # Walking through to find position
                seq = variant_tx_seqs[i]
                print("Sample: {}".format(i))

                steps = int(record['codon_pos'] - 1)
                # Add bases upstream or downstream to compensate for neighboring deletions.
                suff_aa_len = False
                to_add_l = 0
                to_add_r = 0
                while not suff_aa_len:
                    #print(seq[steps:steps+3])
                    print("No offset: {}"
                          .format(seq[min(0, steps-((int(args.length)-1)*3))
                                      :max(len(seq), steps+((int(args.length)-1)*3))]))
                    lower_bound = max(0,
                                      steps-((int(args.length)-1)*3)+u_offset-to_add_l)
                    upper_bound = min(steps+((int(args.length)-1)*3)+d_offset+to_add_r+1,
                                      len(seq))
                    # This is to get the context around the peptide of interest.
                    context_lower_bound = max(0,
                             steps-((int(args.length)-1)*3)+u_offset-to_add_l - 75)
                    context_upper_bound = min(steps+((int(args.length)-1)*3)+d_offset+to_add_r+1+75,
                             len(seq))
                    print("SNV position: {}".format(snv_pos))
                    target_seq = seq[lower_bound:upper_bound]
                    context_seq = seq[context_lower_bound:context_upper_bound]
                    print("Peptide lower bound index: {}".format(lower_bound))
                    print("Peptide upper bound index: {}".format(upper_bound))
                    print("Target seq (given bounds): {}".format(str(target_seq)))
                    reduced_target_seq = str(target_seq).replace('X', '')

                    if (len(reduced_target_seq) >= 6*(int(args.length) - 1)+3 or
                        int(record['nt_pos']) < (int(args.length)-1)*3 or
                        int(record['nt_len']) - int(record['nt_pos']) < (int(args.length)-1)*3):
                        print("Sufficient AA length!")
                        suff_aa_len = True
                    else:
                        print(len(reduced_target_seq))
                        print("Insufficient AA length!")
                        to_add_l = target_seq[:snv_pos].count('X')
                        to_add_r = target_seq[snv_pos:].count('X')
                        print("To add left: {}".format(to_add_l))
                        print("To add right: {}".format(to_add_r))
                        if to_add_l == 0 and to_add_r == 0:
                            print("Insufficient length, but nothing to add. Investigate.")
                            suff_aa_len = True

                converted_context_seq = ''
                if re.search('tumor', i):
                    # Ensure that the SNV position within the target sequence
                    # indicates the variant base. This may be either a direct
                    # match _or_ be a case where the target sequence position
                    # is the iupac code for a wt/mt het call (more likely).
                    print("Tumor snv: {}".format(target_seq[snv_pos]))
                    print("Record alt: {}".format(record['alt_nt']))
                    if (target_seq[snv_pos] == record['alt_nt'] or
                        (record['ref_nt'] in seq_utils.iupac_conversion(target_seq[snv_pos]) and
                         record['alt_nt'] in seq_utils.iupac_conversion(target_seq[snv_pos]))):
                        print("Tumor transcript sequence contains expected variant.")
                        print("Before conversion: {}".format(''.join(target_seq)))
                        list_target_seq = list(target_seq)
                        list_target_seq[snv_pos] = record['alt_nt']
                        converted_seq = ''.join(list_target_seq)
                        # Conversion for the context
                        list_context_seq = list(context_seq)
                        print(len(list_context_seq))
                        print(snv_pos)
                        print(min(75, (len(list_context_seq) - snv_pos - 1)))
                        list_context_seq[snv_pos + min(75, (len(list_context_seq) - snv_pos - 1))] = record['alt_nt']
                        converted_context_seq = ''.join(list_context_seq)
                        print("After  conversion: {}".format(converted_seq))
                        tumor_nts.append(''.join(converted_seq).replace('X', ''))
                    else:
                        print(("Cannot add this transcript."
                               "The variant position differs between "
                               "annotations and transcript sequence. Check your references.\n"
                               "Skipping to next entry."))
                else:
                    print("Normal transcript sequence contains the expected wildtype sequence.")
                    if target_seq[snv_pos] == record['ref_nt']:
                        norm_nts.append(''.join(target_seq))
                    else:
                        print(target_seq[snv_pos])
                        print(target_seq[snv_pos-10:snv_pos+10])
                        print(record['ref_nt'])

            for seq in norm_nts:
                norm_aas.append(Seq(''.join(seq).replace('X', '')).translate())
            for seq in tumor_nts:
                tumor_aas.append(Seq(''.join(seq).replace('X', '')).translate())
            if not converted_context_seq:
                continue
            tumor_context_aa = Seq(''.join(converted_context_seq).replace('X', '')).translate()


            if len(list(set(tumor_aas))) != 1 or len(list(norm_aas)) == 0 or \
                    record['transcript'] not in transcript_protein_seqs.keys():
                print("Something is wrong. Check your tumor sequences.")
                print("Something is wrong. Check your transcript versions.")
            else:
                protein_acceptable = 1
                # Making checksum identifier
                print("Aligning normal and mutant amino acids to find variant sequences.")

                aligner = Align.PairwiseAligner()
                aligner.open_gap_score = -5
                print("Norm: {} ".format(norm_aas[0]))
                print("Tum: {} ".format(tumor_aas[0]))
                mt_wt_align = aligner.align(norm_aas[0], tumor_aas[0])[0]
                wt_pep_seq = str(transcript_protein_seqs[record['transcript']])
                print("tumor aas: {}".format(tumor_aas[0]))
                if variant_tx_seqs['tumor'][:3] != 'ATG':
                    # If the transcripts coding sequence isn't properly
                    # annotated then we can't handle it for now.
                    continue
                mt_known_align = aligner.align(wt_pep_seq, tumor_aas[0])[0]
                if args.debug_output:
                    with open(args.debug_output, 'a') as deofo:
                        deofo.write("-------------------\n")
                        for i,j in record.items():
                            deofo.write("{}: {}\n".format(i, j))
                        deofo.write("Derived wt. vs. mt. alignment:\n")
                        deofo.write("{}".format(str(mt_wt_align)))
                        deofo.write("Derived known ref. vs. mt. alignment:\n")
                        deofo.write("{}".format(str(mt_known_align)))


                print("Test line1: {}".format(str(mt_wt_align)))
                print("Test line2: {}".format(str(mt_wt_align).split()[5]))
                if '.' not in str(mt_wt_align).split()[5]:
                    protein_acceptable = 0
                else:
                    mut_aa_pos = "{}".format(str(mt_wt_align).split()[5].index('.'))
                    hash_components = ':'.join([str(record['meta'].CHROM),
                                                str(record['meta'].POS),
                                                str(record['transcript']),
                                                str(record['meta'].REF),
                                                str(record['meta'].ALT)])
                    print(hash_components)
                    var_md5 = hashlib.md5("{}".format(hash_components).encode('utf-8'))\
                                                                      .hexdigest()[:15]
                # Peptide entry
                if protein_acceptable:
                    meta_line_comp1 = "{}\n;antigen_source:SNV variant_coords:{}:{} "\
                                      .format(var_md5,
                                              record['meta'].CHROM,
                                              record['meta'].POS)
                    meta_line_comp2 = "variant_position_in_cds:{} transcript_id:{} "\
                                      .format(record['nt_pos'],
                                              record['transcript'])
                    meta_line_comp3 = "snv_ref_allele:{} snv_alt_allele:{} snv_type:{} "\
                                      .format(record['meta'].REF,
                                              record['meta'].ALT,
                                              'missense')
                    meta_line_comp4 = "mut_aa_pos:{} pep_context:{} nt_context:{}"\
                                    .format(mut_aa_pos,
                                            tumor_context_aa,
                                            converted_context_seq)
                    meta_line = meta_line_comp1 + meta_line_comp2 + meta_line_comp3 + \
                                meta_line_comp4

                    print("Adding tumor sequence to mutant peptides...")
                    alts = record['meta'].ALT
                    if len(alts) == 1:
                        new_meta_line = "{}:{}:{}:{}".format(record['meta'].CHROM, record['meta'].POS, record['meta'].REF, alts[0])
                        mutant_peptides[new_meta_line] = tumor_aas[0]

                    # Nucleotide entry

                    nt_meta_line_comp1 = "{} SRC:SNV VARIANT_POS:{}:{} TX_POS:{} "\
                                          .format(var_md5,
                                                  record['meta'].CHROM,
                                                  record['meta'].POS,
                                                  record['nt_pos'])
                    nt_meta_line_comp2 = "CDNA_POS:{} TRANSCRIPT:{} REF:{} ALT:{} "\
                                          .format(record['cdna_pos'],
                                                  record['transcript'],
                                                  record['meta'].REF,
                                                  record['meta'].ALT)
                    nt_meta_line_comp3 = "SNV_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA"\
                                          .format('missense')

                    nt_meta_line = nt_meta_line_comp1 + nt_meta_line_comp2 + nt_meta_line_comp3
                    mutant_seqs[nt_meta_line] = tumor_nts[0]

                    wt_meta_line_comp1 = "{}\n;SRC:SNV VARIANT_POS:{}:{} TRANSCRIPT:{} "\
                                         .format(var_md5,
                                                 record['meta'].CHROM,
                                                 record['meta'].POS,
                                                 record['transcript'])
                    wt_meta_line_comp2 = "REF:{} ALT:{} SNV_TYPE:{} "\
                                         .format(record['meta'].REF,
                                                 record['meta'].ALT,
                                                 'missense')
                    wt_meta_line = wt_meta_line_comp1 + wt_meta_line_comp2 + \
                                   "PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA"

                    if len(list(set(norm_aas))) > 1:
                        for norm_aa in norm_aas:
                            #wt entry (multiple alleles)
                            wildtype_peptides[wt_meta_line] = norm_aa
                    else:
                        #wt entry (single allele)
                        wildtype_peptides[wt_meta_line] = norm_aas[0]


    with open(args.mt_output, 'w') as ofo:
        for mutant_peptide_id, mutant_peptide_seq in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(mutant_peptide_id, mutant_peptide_seq))
    with open(args.nt_output, 'w') as ofo:
        for mutant_seq_id, mutant_seq in mutant_seqs.items():
            ofo.write('>{}\n{}\n'.format(mutant_seq_id, mutant_seq))
    with open(args.wt_output, 'w') as ofo:
        for wt_peptide_id, wt_peptide_seq in wildtype_peptides.items():
            ofo.write('>{}\n{}\n'.format(wt_peptide_id, wt_peptide_seq))


def make_indel_targeted_peptides(indel_start, indel_type, variant_tx_seqs, pep_len, indel_len=0):
    """
    """
    print("variant_tx_seqs: {}".format(variant_tx_seqs))
    seq = variant_tx_seqs['tumor']
    upstream_offset = 0
    downstream_offset = 0
    print("Indel start: {}".format(indel_start))
    print("Indel start % 3: {}".format(float(indel_start) % 3))

    # Want to shift "focus" to the first nucleotide of the variant
    # codon for the purposes of translation and extracting appropriate
    # sequence.
    indel_start = float(indel_start)
    if (float(indel_start)) % 3 == 0:
        print("Direct translate")
        upstream_offset += 0
        downstream_offset += 0
    elif (float(indel_start)) % 3 == 1:
        print("Minus one")
        upstream_offset += -1
        downstream_offset += -1
    elif (float(indel_start)) % 3 == 2:
        print("Minus two")
        upstream_offset += -2
        downstream_offset += -2

    indel_start = int(indel_start)
    # This is for handling cases where variants occur at the end of a transcript.
    # Internal note: caused by sample ni 06JUN2025 notes.
#    indel_start = min(int(indel_start), len(seq) - 1)
    target_seq = ''
    # Walking through to find variant position. The bcftools consensus
    # approach includes all neighboring variants and these may shift
    # the location of the somatic variant of interest for the purpose
    # of generating the peptide. The idea is to count only positions
    # that were in in the consensus position and skip counting any
    # inserted nucleotides. This is performed by having bcftools
    # consensus emit inserted variants as an alternate case (e.g. lower
    # case where reference sequence is upper case, etc.).

    # snv_pos is the position that you want the variant to be in your resulting
    # pre-translation nucleotide sequence.



    snv_pos = (int(pep_len)-1) * 3
    print("InDel position (where it starts in targeted seq) is: {}".format(snv_pos))
    steps = 1
#    seq = variant_tx_seqs['tumor']
    print("Tumor sequence: {}".format(seq))
    print("Sequence length: {}".format(len(seq)))
    print("InDel start: {}".format(indel_start))
    print("Non-step based indel start: {}".format(seq[indel_start]))
    print("idx\tbase")
    for i in range(indel_start-5, min(len(seq), indel_start+5)):
        print("{}\t{}".format(i, seq[i]))
    for pos_idx, pos in enumerate(list(seq)):
        if int(steps) == int(indel_start+1):
            print("pos_idx: {}".format(pos_idx))
            print("upstream offset: {}".format(upstream_offset))
            print("downstream offset: {}".format(downstream_offset))
            print("indel len: {}".format(indel_len))
            print("lb: {}".format(pos_idx-snv_pos+upstream_offset))
            print("ub: {}".format(pos_idx+snv_pos+1 + downstream_offset + indel_len))
            if indel_type == 'inframe':
                target_seq = seq[max(0, pos_idx - snv_pos + upstream_offset - 3)\
                                 :min(pos_idx+snv_pos + downstream_offset + indel_len + 3, \
                                      len(seq))]
            elif indel_type == 'frameshift':
                print("pos_idx: {}".format(pos_idx))
                print("snv_pos: {}".format(snv_pos))
                print("upstream_offset: {}".format(upstream_offset))
                target_seq = seq[max(0, pos_idx-snv_pos+upstream_offset):]
            break
        elif re.search('[A-Za-z]', pos):
            steps += 1
    print(target_seq)
    target_seq = ''.join(list(target_seq))
    print(target_seq)
    target_seq = target_seq.replace('X', '')
    print("Target seq: {}".format(target_seq))
    print("Target AA: {}".format(Seq(target_seq).translate()))
    tumor_aas = []
    tumor_nts = []
    tumor_aas.append("{}".format(Seq(target_seq).translate(to_stop=True)))
    tumor_nts.append("{}".format(Seq(target_seq)))
    print("within function tumor aas: {}".format(tumor_aas))
    print("tumor nts: {}".format(tumor_nts))

    return tumor_aas, tumor_nts


def make_indel_peptides_context(args):
    """
    """
    transcript_protein_seqs = parsers.parse_fasta(args.pep_ref)

    three_prime_seqs = parsers.parse_fasta_full_header(args.three_prime_seqs)

    three_prime_tx_to_seqs = {}

    with open(args.three_prime_keys) as ofo:
        for line in ofo.readlines():
            line = line.strip().split('\t')
            tx = line[0].split(':')[0]
            coords = line[1]
            three_prime_tx_to_seqs[tx] = three_prime_seqs[coords]

    pprint(three_prime_tx_to_seqs.items())


    mutant_peptides = {}
    mutant_seqs = {}

    somatic_txs = []

    print("Loading all transcripts harboring somatic variants...")
    with open(args.somatic_txs) as somatic_txs_fo:
        for line in somatic_txs_fo.readlines():
            somatic_txs.append(line.strip().replace('"', ''))
    print("Loaded set of transcripts harboring somatic variants.")


    # Loading actual somatic variants of interest (all InDels)
    # We will need to retrieve the nucleotide sequence records as well.
    print("Getting indels...")
    conserv_inframe_indels = var_utils.extract_inframe_indels(args.somatic_vcf, 'conservative')
    print("Retrieved {} conservative inframe indels...".format(len(conserv_inframe_indels)))
    disrupt_inframe_indels = var_utils.extract_inframe_indels(args.somatic_vcf, 'disruptive')
    print("Retrieved {} disruptive inframe indels...".format(len(disrupt_inframe_indels)))
    frameshift_indels = var_utils.extract_frameshift_indels(args.somatic_vcf)
    print("Retrieved {} frameshift indels...".format(len(frameshift_indels)))

    inframe_indels = dict(conserv_inframe_indels, **disrupt_inframe_indels)

    variant_txs_metadata = annots.extract_variant_tx_metadata(somatic_txs, args.gtf)

    for entry in inframe_indels.keys():
        for record in inframe_indels[entry]:
            variant_tx = record['transcript']
            if variant_tx not in somatic_txs or re.search('N', record['nt_change']) or \
                               variant_tx not in transcript_protein_seqs.keys() :
                continue
            print("NEW INFRAME RECORD: {}".format(record))
            record_coords = "{}_{}".format(record['meta'].CHROM, record['meta'].POS)
            print(record_coords)
            print("Variant tx: {}".format(variant_tx))
            print("Coding seqs: {}".format(variant_txs_metadata[variant_tx]['cds']))

            exon_seqs = {}
            exon_seqs['tumor'] = {}
            if not glob(os.path.join(args.var_tx_seqs,'*{}_{}.tumor.fa'.format(record['transcript'],
                                                                               record_coords))):
                continue
            for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs,
                                                            '*{}_{}.tumor.fa'
                                                            .format(record['transcript'], \
                                                                    record_coords)))[0], "fasta"):
                exon_seqs['tumor'][seq_record.description] = seq_record.seq

            variant_tx_seqs = {}

            samples = ['tumor']
            for i in samples:
                variant_tx_seqs[i] = ''
            if variant_txs_metadata[variant_tx]['strand'] == '+':
                print('positive_strand')
                for i in samples:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
                        variant_tx_seqs[i] += exon_seqs[i][cds]
            elif variant_txs_metadata[variant_tx]['strand'] == '-':
                print('negative_strand')
                for i in samples:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
                        variant_tx_seqs[i] += exon_seqs[i][cds].reverse_complement()

            print("pre-variant_tx_seqs:{}\n{}".format(i, variant_tx_seqs[i]))

            indel_start = 0
            indel_stop = 0
            indel_len = 0

            print(record)

            # DELETIONS #
            if re.search("del", record['nt_change']):
                print("INFRAME DELETION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                del_rec_coords = record['nt_change'].partition('del')[0]
                if '_' in del_rec_coords:
                    indel_start = int(del_rec_coords.split('_')[0]) - 1
                    indel_stop = int(del_rec_coords.split('_')[1]) - 1
                    indel_len = int(indel_stop) - int(indel_start) + 1
                else:
                    indel_start = del_rec_coords - 1
                    indel_stop = del_rec_coords - 1
                    indel_len = 1

            # INSERTIONS #
            elif re.search("ins", record['nt_change']):
                print("INFRAME INSERTION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')
                left_flank = int(ins_rec_nt.split('_')[0])
                right_flank = int(ins_rec_nt.split('_')[1])
                indel_len = len(ins_nt)
                indel_start = left_flank
                indel_stop = indel_start + indel_len

            # DUPLICATIONS #
            elif re.search("dup", record['nt_change']):
                print("INFRAME DUPLICATION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                dup_rec_nt, bufr, ins_nt = record['nt_change'].partition('dup')
                if '_' in dup_rec_nt:
                    indel_start = dup_rec_nt.split('_')[1]
                    indel_stop = dup_rec_nt.split('_')[1]
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = dup_rec_nt
                    indel_stop = dup_rec_nt
                    indel_len = 1

            print("Indel start: {}".format(indel_start))
            print("Indel start % 3: {}".format(float(indel_start) % 3))
            if int(indel_start) > len(variant_tx_seqs['tumor']):
                print("InDel insertion coordinates exceed the coding sequence of the transcript.")
                print("This likely means the variant actually in a UTR. Skipping variant.")
                continue

            tumor_aas, tumor_nts = make_indel_targeted_peptides(indel_start, 'inframe', \
                                                                variant_tx_seqs, args.length, \
                                                                indel_len)
            context_tumor_aas, context_tumor_nts = make_indel_targeted_peptides(indel_start, \
                                                                                'inframe', \
                                                                                variant_tx_seqs, \
                                                                                args.length + 50, \
                                                                                indel_len)

            wt_pep_seq = str(transcript_protein_seqs[record['transcript']])
            print("outer tumor aas: {}".format(tumor_aas[0]))
            if variant_tx_seqs['tumor'][:3] != 'ATG':
                # If the transcripts coding sequence isn't properly annotated
                # then we can't handle it for now.
                continue
            norm_upstream_offset = 0
            norm_downstream_offset = 0
            if (float(indel_start)) % 3 == 0:
                print("norm Direct translate")
                norm_upstream_offset += 0
                norm_downstream_offset += 0
            elif (float(indel_start)) % 3 == 1:
                print("norm Minus one")
                norm_upstream_offset += -1
                norm_downstream_offset += -1
            elif (float(indel_start)) % 3 == 2:
                print("norm Minus two")
                norm_upstream_offset += -2
                norm_downstream_offset += -2
            indel_start = float(indel_start)
            indel_stop = float(indel_stop)
            print("norm lb: {}".format((indel_start+norm_upstream_offset)/3+1 - int(args.length)))
            print("norm ub: {}".format((indel_start+norm_upstream_offset)/3+2 + int(args.length)))
            print("wt pep seq: {}".format(wt_pep_seq))
            naas_lb = max(0, int((indel_start+norm_upstream_offset)/3+1 - int(args.length) -1))
            naas_ub = int((indel_stop+norm_downstream_offset)/3+2 + int(args.length))
            norm_aas = wt_pep_seq[naas_lb:naas_ub]
            print("norm aas: {}".format(norm_aas))

            print("Aligning normal and mutant amino acids to find variant sequences.")


            # This is an attempt to correctly align the mutant and wildtype
            # sequences. There are some cases in which the mutations starting
            # in the third base of a codon may be effects beyond the expected
            # amino acids.
            correct_align = False
            not_sufficient = False
            prev_mut_aa_locs = ''
            # Adding the norm_aas and tumor_aas requirements since sometimes
            # tumor_aas canont be translated due to upstream frameshift
            # variants.
            while not correct_align and not not_sufficient and norm_aas and tumor_aas[0]:
                mut_aa_locs = []
                aligner = Align.PairwiseAligner()
                aligner.open_gap_score = -5
                aligner.mismatch_score = -5
                targeted_ref_vs_mt_align = aligner.align(norm_aas, tumor_aas[0])[0]
                ref_vs_mt_align = aligner.align(wt_pep_seq, tumor_aas[0])[0]
                tumor_aa_len = len(str(targeted_ref_vs_mt_align).split()[-1].rstrip('-'))
                print(targeted_ref_vs_mt_align)
                for idx, i in enumerate(str(targeted_ref_vs_mt_align).split()[1][:tumor_aa_len]):
                    if i in ['-']:
                        mut_aa_locs.append('{}'.format((int(idx+1))))
                print(mut_aa_locs)
                if not mut_aa_locs:
                    break
                if int(min(mut_aa_locs)) > (int(args.length)):
                    norm_aas = norm_aas[1:]
                    tumor_aas[0] = tumor_aas[0][1:]
                if len(tumor_aas[0]) - int(max(mut_aa_locs)) + 1 > int(args.length):
                    norm_aas = norm_aas[:-1]
                    tumor_aas[0] = tumor_aas[0][:-1]
                elif int(min(mut_aa_locs)) == (int(args.length)):
                    correct_align = True
                elif int(min(mut_aa_locs)) == prev_mut_aa_locs:
                    not_sufficient = True
                prev_mut_aa_locs = int(min(mut_aa_locs))

            if correct_align:
                for i in str(targeted_ref_vs_mt_align).split():
                    print("{}".format(i[:tumor_aa_len]))
                print("Mutant locations: {}".format(mut_aa_locs))

                valid_ref_orf = "F"
                print("")
                print("Checking reference ORF alignment...")
                print("External ref vs. mutant alignment:\n{}".format(ref_vs_mt_align))
                check_range = ''
                # Have to deal with cases where there may be no matched bases.
                if '|' in str(ref_vs_mt_align).split('\n')[1]:
                    pipe_start_index = str(ref_vs_mt_align).split('\n')[1].index('|')
                    print("First | index: {}".format(pipe_start_index))
                    check_range = str(ref_vs_mt_align).split('\n')[1][pipe_start_index:\
                                                                      pipe_start_index +\
                                                                          args.length - 1]
                print("Check range: {} ".format(check_range))
                if check_range == '|' * (int(args.length) - 1):
                    print("Valid ORF compared to reference protein!")
                    valid_ref_orf = "T"

                mut_aa_ranges = numbers_to_ranges(mut_aa_locs)
                print("Mutant position ranges: {}".format(mut_aa_ranges))
                if args.debug_output:
                    with open(args.debug_output, 'a') as deofo:
                        deofo.write("-------------------\n")
                        for i,j in record.items():
                            deofo.write("{}: {}\n".format(i, j))
                        deofo.write("Focused ref. vs. mt. alignment:\n")
                        deofo.write("{}".format(str(targeted_ref_vs_mt_align)))
                        deofo.write("Full ref. vs. mt. alignment:\n")
                        deofo.write("{}".format(str(ref_vs_mt_align)))

                protein_acceptable = 1
                # Making checksum identifier

                checksum_components = [str(record['meta'].CHROM),
                                       str(record['meta'].POS),
                                       str(record['transcript']),
                                       str(record['meta'].REF),
                                       str(record['meta'].ALT)]

                var_md5 = hashlib.md5("{}".format(':'.join(checksum_components)).encode('utf-8'))\
                                                                                .hexdigest()[:15]

                if protein_acceptable:

                    meta_line_comp1 = "{}\n;antigen_source:INDEL variant_coords:{}:{} "\
                                      .format(var_md5,
                                              record['meta'].CHROM,
                                              record['meta'].POS)

                    meta_line_comp2 = "transcript_id:{} indel_ref_allele:{} indel_alt_allele:{} "\
                                      .format(record['transcript'],
                                              record['meta'].REF,
                                              record['meta'].ALT)

                    meta_line_comp3 = "indel_type:{} mut_aa_range:{} valid_ref_orf:{} "\
                                      .format('inframe',
                                              ','.join(mut_aa_ranges),
                                              valid_ref_orf)

                    meta_line_comp4 = "pep_context:{} nt_context:{}"\
                                      .format(context_tumor_aas[0],
                                              context_tumor_nts[0])
                    meta_line = meta_line_comp1 + meta_line_comp2 + meta_line_comp3 + \
                                meta_line_comp4

                    mutant_peptides[meta_line] = tumor_aas[0]


                    nt_meta_line_comp_1 = "{} SRC:INDEL VARIANT_POS:{}:{} CDNA_POS:{} "\
                                          .format(var_md5,
                                                  record['meta'].CHROM,
                                                  record['meta'].POS,
                                                  record['cdna_pos'])

                    nt_meta_line_comp_2 = "TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} "\
                                          .format(record['transcript'],
                                                  record['meta'].REF,
                                                  record['meta'].ALT,
                                                  'inframe')

                    nt_meta_line = nt_meta_line_comp_1 + nt_meta_line_comp_2 + \
                                   "PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA"

                    mutant_seqs[nt_meta_line] =  tumor_nts[0]

    for entry in frameshift_indels.keys():
        print(entry)
        for record in frameshift_indels[entry]:
            print(record)
            variant_tx = record['transcript']
            if variant_tx not in somatic_txs or variant_tx not in transcript_protein_seqs.keys():
                continue
            print("NEW FRAMESHIFT RECORD: {}".format(record))
            record_coords = "{}_{}".format(record['meta'].CHROM, record['meta'].POS)
            print(record_coords)
            print(variant_tx)

            print(glob(os.path.join(args.var_tx_seqs,
                                    '*{}_{}.tumor.fa'
                                    .format(record['transcript'], record_coords))))
            exon_seqs = {}
            exon_seqs['tumor'] = {}
            if not glob(os.path.join(args.var_tx_seqs,'*{}_{}.tumor.fa'\
                                                       .format(record['transcript'],\
                                                               record_coords))):
                continue
            for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs,
                                                            '*{}_{}.tumor.fa'
                                                            .format(record['transcript'],
                                                                    record_coords)))[0], "fasta"):
                exon_seqs['tumor'][seq_record.description] = seq_record.seq


            variant_tx_seqs = {}

            samples = ['tumor']
            for i in samples:
                variant_tx_seqs[i] = ''
            if variant_txs_metadata[variant_tx]['strand'] == '+':
                print('positive_strand')
                for i in samples:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
                        print(cds)
                        variant_tx_seqs[i] += exon_seqs[i][cds]
                    variant_tx_seqs[i] += three_prime_tx_to_seqs[variant_tx]
            elif variant_txs_metadata[variant_tx]['strand'] == '-':
                print('negative_strand')
                for i in samples:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
                        print(cds)
                        variant_tx_seqs[i] += exon_seqs[i][cds].reverse_complement()
                    variant_tx_seqs[i] += Seq(three_prime_tx_to_seqs[variant_tx]).reverse_complement()

            print(variant_tx_seqs)
            print("pre-variant_tx_seqs:{}\n{}".format(i, variant_tx_seqs[i]))

            indel_start = 0
            indel_stop = 0
            indel_len = 0

            # DELETIONS #
            if re.search("del", record['nt_change']):
                print("FRAMESHIFT DELETION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                del_rec_coords = record['nt_change'].partition('del')[0]

                if '_' in del_rec_coords:
                    indel_start = int(del_rec_coords.split('_')[0]) - 1
                    indel_stop = int(del_rec_coords.split('_')[1]) - 1
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = int(del_rec_coords) - 1
                    indel_stop = int(del_rec_coords) - 1
                    indel_len = 1

            # INSERTIONS #
            elif re.search("ins", record['nt_change']):
                print("INFRAME INSERTION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')

                indel_start = int(ins_rec_nt.split('_')[0])
                indel_stop = int(ins_rec_nt.split('_')[1])
                indel_len = int(indel_stop) - int(indel_start) + 2

            # DUPLICATIONS #
            elif re.search("dup", record['nt_change']):
                dup_rec_coords, bufr, ins_nt = record['nt_change'].partition('dup')
                if '_' in dup_rec_coords:
                    indel_start = dup_rec_coords.split('_')[0]
                    indel_stop = dup_rec_coords.split('_')[1]
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = int(dup_rec_coords)
                    indel_stop = int(dup_rec_coords)
                    indel_len = 1
            if int(indel_start) > len(variant_tx_seqs['tumor']):
                print("InDel insertion coordinates exceed the coding sequence of the transcript.")
                print("This likely means the variant actually in a UTR. Skipping variant.")
                continue

            tumor_aas, tumor_nts = make_indel_targeted_peptides(indel_start, 'frameshift',
                                                                variant_tx_seqs, args.length,
                                                                indel_len)
            context_tumor_aas, context_tumor_nts = make_indel_targeted_peptides(indel_start,
                                                                                'frameshift',
                                                                                variant_tx_seqs,
                                                                                args.length + 50,
                                                                                indel_len)

            print("indel start: {}".format(indel_start))
            indel_start = float(indel_start)
            wt_pep_seq = str(transcript_protein_seqs[record['transcript']])
            print("tumor aas: {}".format(tumor_aas[0]))
            if variant_tx_seqs['tumor'][:3] != 'ATG':
                # If the transcripts coding sequence isn't properly annotated
                # then we can't handle it for now.
                continue
            norm_upstream_offset = 0
            norm_downstream_offset = 0
            if (float(indel_start)) % 3 == 0:
                print("norm Direct translate")
                norm_upstream_offset += 0
                norm_downstream_offset += 0
            elif (float(indel_start)) % 3 == 1:
                print("norm Minus one")
                norm_upstream_offset += -1
                norm_downstream_offset += -1
            elif (float(indel_start)) % 3 == 2:
                print("norm Minus two")
                norm_upstream_offset += -2
                norm_downstream_offset += -2
            indel_start = float(indel_start)
            indel_stop = float(indel_stop)
            print("norm lb: {}".format((indel_start-1+norm_upstream_offset)/3 - int(args.length)))
            print("norm ub: {}".format((indel_start-1+norm_upstream_offset)/3 + int(args.length)))
            norm_aas = wt_pep_seq[int((indel_start+norm_upstream_offset)/3+1 - int(args.length)):\
                                  int((indel_stop+norm_downstream_offset)/3+2 + int(args.length))]
            print("norm aas: {}".format(norm_aas))

            print("Aligning normal and mutant amino acids to find variant sequences.")

            aligner = Align.PairwiseAligner()
            aligner.open_gap_score = -5

            if norm_aas and tumor_aas[0]:

                targeted_ref_vs_mt_align = aligner.align(norm_aas, tumor_aas[0])[0]
                ref_vs_mt_align = aligner.align(wt_pep_seq, tumor_aas[0])[0]

                print(targeted_ref_vs_mt_align)

                if args.debug_output:
                    with open(args.debug_output, 'a') as deofo:
                        deofo.write("-------------------\n")
                        for i,j in record.items():
                            deofo.write("{}: {}\n".format(i, j))
                        deofo.write("Focused ref. vs. mt. alignment:\n")
                        deofo.write("{}".format(str(targeted_ref_vs_mt_align)))
                        deofo.write("Full ref. vs. mt. alignment:\n")
                        deofo.write("{}".format(str(ref_vs_mt_align)))

                valid_ref_orf = "F"
                if re.search('\|', str(ref_vs_mt_align).split('\n')[1]):
                    print("")
                    print("Checking reference ORF alignment...")
                    print("External ref vs. mutant alignment:\n{}".format(ref_vs_mt_align))
                    pipe_start_index = str(ref_vs_mt_align).split('\n')[1].index('|')
                    print("First | index: {}".format(pipe_start_index))
                    check_range = str(ref_vs_mt_align).split('\n')[1][pipe_start_index:\
                                                                      pipe_start_index +\
                                                                          args.length - 1]
                    print("Check range: {} ".format(check_range))
                    if check_range == '|' * (int(args.length) - 1):
                        print("Valid ORF compared to reference protein!")
                        valid_ref_orf = "T"
                else:
                    print("Alignment does not appear valid. Please investigate.")

                # Need to deal with cases where frameshift results in stop codon
                # immediately after reference codon.

                checksum_components = [str(record['meta'].CHROM),
                                       str(record['meta'].POS),
                                       str(record['transcript']),
                                       str(record['meta'].REF),
                                       str(record['meta'].ALT)]

                var_md5 = hashlib.md5("{}".format(':'.join(checksum_components))\
                                                     .encode('utf-8')).hexdigest()[:15]

                if valid_ref_orf:
                    meta_line_comp1 = "{}\n;antigen_source:INDEL variant_coords:{}:{} "\
                                      .format(var_md5,
                                              record['meta'].CHROM,
                                              record['meta'].POS)

                    meta_line_comp2 = "transcript_id:{} indel_ref_allele:{} indel_alt_allele:{} "\
                                      .format(record['transcript'],
                                              record['meta'].REF,
                                              record['meta'].ALT)

                    meta_line_comp3 = "indel_type:{} valid_ref_orf:{} "\
                                      .format('frameshift',
                                              valid_ref_orf)

                    meta_line_comp4 = "pep_context:{} nt_context:{}"\
                                      .format(context_tumor_aas[0],
                                              context_tumor_nts[0])
#                    meta_line_comp4 = "pep_context:{} nt_context:{}"\
#                                      .format('NA',
#                                              'NA')
                    meta_line = meta_line_comp1 + meta_line_comp2 + meta_line_comp3 + \
                                meta_line_comp4

                    mutant_peptides[meta_line] = tumor_aas[0]

                    nt_meta_line_comp_1 = "{} SRC:INDEL VARIANT_POS:{}:{} CDNA_POS:{} "\
                                          .format(var_md5,
                                                  record['meta'].CHROM,
                                                  record['meta'].POS,
                                                  record['cdna_pos'])

                    nt_meta_line_comp_2 = "TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} "\
                                          .format(record['transcript'],
                                                  record['meta'].REF,
                                                  record['meta'].ALT,
                                                  'frameshift')

                    nt_meta_line = nt_meta_line_comp_1 + nt_meta_line_comp_2 + \
                                   "PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA"

                    mutant_seqs[nt_meta_line] =  tumor_nts[0]

    with open(args.output, 'w') as ofo:
        for mutant_peps, mutant_pep_seq in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(mutant_peps, mutant_pep_seq))
    with open(args.nt_output, 'w') as ofo:
        for mutant_nts, mutant_nts_seq in mutant_seqs.items():
            ofo.write('>{}\n{}\n'.format(mutant_nts, mutant_nts_seq))


def make_erv_peptides(args):
    """
    """
    herv_id_to_utr = {}
    post_met_to_full = {}

    with open(args.geve_reference) as geve_fo:
        for line_idx, line in enumerate(geve_fo.readlines()):
            if line_idx != 0:
                line = line.split('\t')
                post_met = re.sub('.M$', '', line[8])
                if len(post_met.split('.')) > 1 :
                    post_met_s = post_met.split('.')
                    chrom, first, second, strand = (post_met_s[1], post_met_s[2], post_met_s[3],
                                                  post_met_s[-1])
                    post_met_to_full["{}:{}-{}".format(chrom, first, second)] = line[0]

                    if strand == '+':
                        herv_id_to_utr[line[0]] = int(post_met.split('.')[2]) - \
                                                  int(line[0].split('.')[2])
                    elif strand == '-':
                        herv_id_to_utr[line[0]] = int(line[0].split('.')[3]) - \
                                                  int(post_met.split('.')[3])


    expressed_hervs = []
    ervs_metad = {}
    col_map = {}

    with open(args.expressed_ervs) as exp_ervs_fo:
        for line_idx, line in enumerate(exp_ervs_fo.readlines()):
            line = line.rstrip().split(',')
            if line_idx == 0:
                for idx, i in enumerate(line):
                    col_map[i] = idx
            else:
                herv_id = line[col_map['Name']]
                expressed_hervs.append(herv_id)
                ervs_metad[herv_id] = {}
                ervs_metad[herv_id]['tumor_cpm'] = line[col_map['Tumor_CPM']]
                ervs_metad[herv_id]['norm_cpm'] = line[col_map['Norm_CPM']]
                ervs_metad[herv_id]['delta'] = line[col_map['log2(Tumor_CPM+1)-log2(Norm_CPM+1)']]
                ervs_metad[herv_id]['tpm'] = line[col_map['TPM']]
                ervs_metad[herv_id]['numreads'] = line[col_map['NumReads']]

    expressed_hervs_seqs = {}

    for seq_record in SeqIO.parse(args.patient_ervs_fasta, "fasta"):
        print(seq_record)
        chrom = seq_record.description.split(':')[0]
        start = seq_record.description.split(':')[1].split('-')[0]
        stop = seq_record.description.split(':')[1].split('-')[1]
        if args.species != 'mm':
            new_seq_record_description = 'Hsap38.{}.{}.{}'.format(chrom, start, stop)
        else:
            new_seq_record_description = 'Mmus38.{}.{}.{}'.format(chrom, start, stop)
        regexed_herv_id = post_met_to_full["{}:{}-{}".format(chrom, start, stop)]
        expressed_hervs_seqs[regexed_herv_id] = seq_record.seq

    expressed_hervs_aas = {}

    for id, seq in expressed_hervs_seqs.items():
        if id.endswith('-'):
            aa_seq = seq.reverse_complement().translate(to_stop=True)
            expressed_hervs_seqs[id] = seq.reverse_complement()
        else:
            aa_seq = seq.translate(to_stop=True)
        expressed_hervs_aas["{}".format(id)] = aa_seq

    with open(args.output, 'w') as ofo:
        for k, v in expressed_hervs_aas.items():
            header_content = (k,
                              ervs_metad[k]['tumor_cpm'],
                              ervs_metad[k]['norm_cpm'],
                              ervs_metad[k]['delta'],
                              ervs_metad[k]['tpm'],
                              ervs_metad[k]['numreads'],
                              v)
            hash = hashlib.md5(str(k).encode('utf-8')).hexdigest()[:15]
            header = f"{hash}"
            comment_1 = "antigen_source:ERV erv_orf_id:{} erv_tumor_cpm:{} erv_norm_cpm:{} "\
                        .format(*header_content[:3])
            comment_2 = "erv_tumor_cpm_to_norm_cpm_delta:{} erv_orf_tpm:{} "\
                        .format(*header_content[3:5])
            comment_3 = "erv_orf_raw_read_count:{} pep_context:{}"\
                        .format(*header_content[5:])

            comment = comment_1 + comment_2 + comment_3
            ofo.write(">{}\n;{}\n{}\n".format(header, comment, v))

    with open(args.nt_output, 'w') as ofo:
        for k, v in expressed_hervs_seqs.items():
            header = "{} NAME:{}".format(hashlib.md5(str(k).encode('utf-8')).hexdigest()[:15], k)
            ofo.write(">{}\n{}\n".format(header, v))


def make_self_antigen_peptides(args):
    """
    """
    print(args)
    expressed_selfs_exon_seqs = {}
    for seq_record in SeqIO.parse(args.selfs_seqs_fasta, "fasta"):
        expressed_selfs_exon_seqs[seq_record.description] = seq_record.seq

    print(expressed_selfs_exon_seqs.keys())

    expressed_selfs = []
    with open(args.expressed_selfs) as exp_selfs_fo:
        for line in exp_selfs_fo.readlines():
            expressed_selfs.append(line.split(':')[1].rstrip('\n'))

    expressed_selfs_tx_metadata = {}
    for expressed_self in expressed_selfs:
        print(expressed_self)
        expressed_selfs_tx_metadata[expressed_self] = {}
        expressed_selfs_tx_metadata[expressed_self]['cds'] = []
        expressed_selfs_tx_metadata[expressed_self]['UTR'] = []
        expressed_selfs_tx_metadata[expressed_self]['strand'] = ''


    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            for expressed_self in expressed_selfs:
                if re.search(expressed_self, line):
                    if re.search('\tCDS\t', line):
                        splitl = line.split('\t')
                        chr, start, stop, strand = (splitl[0], splitl[3], splitl[4], splitl[6])
                        expressed_selfs_tx_metadata[expressed_self]['cds'].append("{}:{}-{}"
                                                                                  .format(chr,
                                                                                          start,
                                                                                          stop))
                        expressed_selfs_tx_metadata[expressed_self]['strand'] = strand
                    #Assumes GENCODE GTF here...
                    elif re.search('\tUTR|five_prime_utr\t', line):
                        splitl = line.split('\t')
                        chr, start, stop, strand = (splitl[0], splitl[3], splitl[4], splitl[6])
                        expressed_selfs_tx_metadata[expressed_self]['UTR'].append("{}:{}-{}"
                                                                                  .format(chr,
                                                                                          start,
                                                                                          stop))

    expressed_selfs_nts = {}
    expressed_selfs_utr_buffers = {}

    # Is there a way to make this more compact while still readable?
    for expressed_self in expressed_selfs_tx_metadata.keys():
        expressed_self_tx_seq = ''
        utr_buffer = 0
        if expressed_selfs_tx_metadata[expressed_self]['strand'] == '+':
            if 'UTR' in expressed_selfs_tx_metadata[expressed_self].keys():
                for utr in expressed_selfs_tx_metadata[expressed_self]['UTR']:
                    if (int(utr.split(':')[1].split('-')[1]) <
                        int(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'])[0]\
                            .split(':')[1].split('-')[0])):
                        utr_buffer += (int(utr.split(':')[1].split('-')[1]) + 1 \
                                       - int(utr.split(':')[1].split('-')[0]))
            for cds in sorted(expressed_selfs_tx_metadata[expressed_self]['cds']):
                expressed_self_tx_seq += expressed_selfs_exon_seqs[cds]
        elif expressed_selfs_tx_metadata[expressed_self]['strand'] == '-':
            if 'UTR' in expressed_selfs_tx_metadata[expressed_self].keys():
                for utr in expressed_selfs_tx_metadata[expressed_self]['UTR']:
                    if (int(utr.split(':')[1].split('-')[0]) >\
                        int(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'],
                                   reverse=True)[0].split(':')[1].split('-')[1])):
                        utr_buffer += (int(utr.split(':')[1].split('-')[1]) + 1  \
                                       - int(utr.split(':')[1].split('-')[0]))
            for cds in sorted(expressed_selfs_tx_metadata[expressed_self]['cds']):
                expressed_self_tx_seq += expressed_selfs_exon_seqs[cds]
            expressed_self_tx_seq = str(expressed_self_tx_seq.reverse_complement())

        expressed_selfs_nts[expressed_self] = expressed_self_tx_seq
        expressed_selfs_utr_buffers[expressed_self] = utr_buffer

    selfs_peps = {}
    for self, self_seq in expressed_selfs_nts.items():
        seq = Seq(str(self_seq))
        pep_seq = seq.translate(to_stop=True)
        if pep_seq.startswith('M'):
            selfs_peps[self] = pep_seq

    with open(args.output, 'w') as ofo:
        for expressed_self, expressed_self_seq in sorted(selfs_peps.items()):
            ofo.write(">{}\n;antigen_source:CTA/SELF transcript_id:{} pep_context:{}\n{}\n"\
                      .format(hashlib.md5("{}".format(expressed_self).encode('utf-8'))\
                                                                     .hexdigest()[:15],
                                                      expressed_self,
                                                      expressed_self_seq,
                                                      expressed_self_seq))

    with open(args.nt_output, 'w') as ofo:
        for self, self_seq in expressed_selfs_nts.items():
            seq = Seq(str(self_seq))
            ofo.write(">{} SRC:CTA/SELF NAME:{} UTR_BUFFER:{}\n{}\n"\
                      .format(hashlib.md5("{}".format(self).encode('utf-8')).hexdigest()[:15],
                              self,
                              expressed_selfs_utr_buffers[self],
                              seq))


def filter_mutant_peptides(args):
    """
    Filter "mutant" sequences that do not actually have mutant amino acid
    included. Peptides are produced in such a way that a sliding window
    approach can be used to calculcate binding affinity and binding stability
    across all possible 8mers, 9mers, 10mers, etc.
    """

    # Loading this so we can tell whether an InDel is inframe (in which case we
    # only want peptides with variant sequence) or frameshift (in which case we
    # want to allow all 3' peptides).
    indel_pep_type = {}
    pep_len = {}
    if args.mt_fasta:
        for seq_record in SeqIO.parse(args.mt_fasta, "fasta"):
            id = seq_record.id
            md5 = seq_record.description.split(' ')[0].replace(':', '_')[:15]
            pep_type = seq_record.description.split(' ')[5].replace('INDEL_TYPE:', '')
            indel_pep_type[md5] = pep_type
            pep_len[md5] = len(seq_record.seq)

    threshold = int(args.max_peptide_length)
    print(threshold)
    with open(args.mt_netmhcpan) as mt_netmhcpan_fo:
        header = ''
        with open(args.mt_output, 'w') as output_fo:
            for line in mt_netmhcpan_fo.readlines():
                line = line.rstrip().split()
                if len(line) > 0 and line[0] == 'Pos' and not header :
                    header = line
                    print(header)
                    header.remove('BindLevel')
                    output_fo.write("{}\n".format('\t'.join(header)))
                elif len(line) >= 16 and line[0] not in ['Protein', 'Pos']:
                    line = parsers.sanitize_netmhc_record(line)
                    internal_id = line[10]
                    total_len = pep_len[internal_id]
                    if (int(line[0]) + len(line[2])) > threshold:
                        if (int(line[0]) < (total_len - threshold + 2) or
                             indel_pep_type[internal_id] == 'frameshift'):
                            output_fo.write("{}\n".format('\t'.join(line)))
                    elif int(line[0]) == threshold:
                        output_fo.write("{}\n".format('\t'.join(line)))

    if args.wt_netmhcpan:
        with open(args.wt_netmhcpan) as wt_netmhcpan_fo:
            header = ''
            with open(args.wt_output, 'w') as output_fo:
                for line in wt_netmhcpan_fo.readlines():
                    line = line.rstrip().split()
                    if len(line) > 0 and line[0] == 'Pos' and not header :
                        header = line
                        header.remove('BindLevel')
                        output_fo.write("{}\n".format('\t'.join(header)))
                    elif len(line) >= 16 and line[0] not in ['Protein', 'Pos']:
                        line = parsers.sanitize_netmhc_record(line)
                        internal_id = line[10]
                        total_len = pep_len[internal_id]
                        if int(int(line[0]) + len(line[2])) > threshold:
                            if int(line[0]) < int(total_len - threshold + 2):
                                output_fo.write("{}\n".format('\t'.join(line)))
                        elif int(line[0]) == threshold:
                            output_fo.write("{}\n".format('\t'.join(line)))


def make_viral_peptides(args):
    """
    """
    patient_nts = {}
    for seq_record in SeqIO.parse(args.patient_viral_fasta, "fasta"):
        patient_nts[seq_record.description.split()[0].split('_cds_')[1][:-2]] = seq_record.seq

    patient_aas = {}

    for id, seq in patient_nts.items():
        aa_seq = seq.translate(to_stop=True)
        patient_aas[id] = aa_seq

    with open(args.output, 'w') as ofo:
        for id, seq in patient_aas.items():
            vir_md5 = hashlib.md5("{}".format(id).encode('utf-8')).hexdigest()[:15]
            ofo.write(">{}\n;antigen_source:VIRUS virus_id:{} pep_context:{}\n{}\n"\
                      .format(vir_md5, id.partition('_1:1')[0], seq, seq))


def make_fusion_peptides_context(args):
    """
    """

    # Load set of relevant fusion transcripts.
    fusion_txs = []
    with open(args.fusion_txs) as fusion_txs_fo:
        for line in fusion_txs_fo.readlines():
            fusion_txs.append(line.rstrip())

    # Load relevant exonic sequences from fusion transcripts.
    exon_seqs = {}
    if args.exons_fasta:
        for seq_record in SeqIO.parse(args.exons_fasta, "fasta"):
            exon_seqs[seq_record.description] = seq_record.seq

    # Loading transcript metadata to be paired with exonic sequences.
    fusion_txs_metadata = {}
    for fusion_tx in fusion_txs:
        fusion_txs_metadata[fusion_tx] = {}
        fusion_txs_metadata[fusion_tx]['cds'] = []
        fusion_txs_metadata[fusion_tx]['strand'] = ''
    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            for fusion_tx in fusion_txs:
                if re.search(fusion_tx, line):
                    if re.search('\tCDS\t', line):
                        print(fusion_tx)
                        chr = line.split('\t')[0]
                        start = line.split('\t')[3]
                        stop = line.split('\t')[4]
                        strand = line.split('\t')[6]
                        coords = "{}:{}-{}".format(chr, start, stop)
                        fusion_txs_metadata[fusion_tx]['strand'] = strand
                        print("Actual dict: {}".format(fusion_txs_metadata[fusion_tx]['strand']))
                        if coords not in fusion_txs_metadata[fusion_tx]['cds']:
                            print("Adding CDS")
                            fusion_txs_metadata[fusion_tx]['cds'].append("{}:{}-{}".format(chr,
                                                                                           start,
                                                                                           stop))
                            print("Actual dict: {}".format(fusion_txs_metadata[fusion_tx]['cds']))

    # Combining the exonic sequences to construct each component transcripts
    # sequences. This includes sequence modifications like germline variants
    # that may modify the protein sequence (beyond the fusion event).
    fusion_tx_seqs = {}
    if exon_seqs:
        for fusion_tx in fusion_txs:
            fusion_tx_seqs[fusion_tx] = ''
            print(fusion_tx)
            if fusion_txs_metadata[fusion_tx]['strand'] == '+':
                for cds in sorted(fusion_txs_metadata[fusion_tx]['cds']):
                    fusion_tx_seqs[fusion_tx] += exon_seqs[cds]
            elif fusion_txs_metadata[fusion_tx]['strand'] == '-':
                for cds in sorted(fusion_txs_metadata[fusion_tx]['cds'], reverse=True):
                    fusion_tx_seqs[fusion_tx] += exon_seqs[cds].reverse_complement()
        print("Loaded variant transcripts metadata.")

    valid_fusions = {}
    with open(args.fusions) as fusions_fo:
        header = []
        col_map = {}
        for line_idx, line in enumerate(fusions_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx  == 0:
                col_map = {j:i for i, j in enumerate(line)}
            else:
                if line[col_map['FUSION_TRANSL']] != '.':
                    hash_components = ':'.join([line[col_map['CDS_LEFT_ID']],
                                                line[col_map['CDS_LEFT_RANGE']],
                                                line[col_map['CDS_RIGHT_ID']],
                                                line[col_map['CDS_RIGHT_RANGE']]])
                    fus_md5 = hashlib.md5("{}".format(hash_components).encode('utf-8'))\
                                                                      .hexdigest()[:15]
                    fused_tx = ''
                    left_tx = line[col_map['CDS_LEFT_ID']]
                    right_tx = line[col_map['CDS_RIGHT_ID']]
                    if left_tx in fusion_tx_seqs.keys() and right_tx in fusion_tx_seqs.keys():
                        left_tx_end_pos = int(line[col_map['CDS_LEFT_RANGE']].split('-')[1])
                        right_tx_start_pos = int(line[col_map['CDS_RIGHT_RANGE']].split('-')[0])\
                                             - 1
                        print("Both transcripts observed.")
                        print("{}\t{}".format(left_tx, left_tx_end_pos))
                        print("{}\t{}".format(right_tx, right_tx_start_pos))
                        print(line[col_map['#FusionName']])
                        try:
                            fused_tx = fusion_tx_seqs[left_tx][left_tx_end_pos].lower() +\
                                       fusion_tx_seqs[right_tx][right_tx_start_pos:].upper()\
                                       .remove('X')
                            print(fused_tx)
                        except:
                            print("Indexing issues with fused transcript.")
                    start_prot = 0
                    nuc_sub = 48
                    print(line[col_map['#FusionName']])
                    print(line[col_map['CDS_LEFT_RANGE']])
                    print(line[col_map['CDS_LEFT_RANGE']].split('-')[1])
#                    valid_fusions[line[col_map['#FusionName']]] = {}
                    valid_fusions[fus_md5] = {}
                    valid_fusions[fus_md5]['fusion_id'] = line[col_map['#FusionName']]
                    valid_fusions[fus_md5]['left_breakpoint'] = line[col_map['LeftBreakpoint']]
                    valid_fusions[fus_md5]['right_breakpoint'] = line[col_map['RightBreakpoint']]
                    valid_fusions[fus_md5]['left_tx'] = line[col_map['CDS_LEFT_ID']]
                    valid_fusions[fus_md5]['right_tx'] = line[col_map['CDS_RIGHT_ID']]
                    valid_fusions[fus_md5]['left_gn'] = line[col_map['LeftGene']]
                    valid_fusions[fus_md5]['right_gn'] = line[col_map['RightGene']]
                    valid_fusions[fus_md5]['annot'] = line[col_map['annots']]
                    if float(line[col_map['CDS_LEFT_RANGE']].split('-')[1]) % 3 == 0:
                        print("Direct translate")
                        start_prot = int(int(int(line[col_map['CDS_LEFT_RANGE']].split('-')[1]))/3)
                    elif (float(line[col_map['CDS_LEFT_RANGE']].split('-')[1]) - 1) % 3 == 0:
                        print("Minus one")
                        start_prot = int(int(int(line[col_map['CDS_LEFT_RANGE']]\
                                                              .split('-')[1]) - 1)/3)
                        nuc_sub = nuc_sub + 1
                    elif (float(line[col_map['CDS_LEFT_RANGE']].split('-')[1]) - 2) % 3 == 0:
                        print("Minus two")
                        start_prot = int(int(int(line[col_map['CDS_LEFT_RANGE']]\
                                                              .split('-')[1]) - 2)/3)
                        nuc_sub = nuc_sub + 2
                    else:
                        print("No luck")
                    print(start_prot)

                    full_nuc_seq = ''
                    full_pep_seq = ''
                    # At this point, if the fused transcript was able to be
                    # generated using the exonic sequences, then use that
                    # sequence as it is most relevant to the patient.
                    # Otherwise, use the sequence generated by STARFusion.
                    if fused_tx:
                        print("Using patient-derived sequence.")
                        full_nuc_seq = fused_tx
                        full_pep_seq = str(full_nuc_seq.translate())
                    else:
                        print("Using STARFusion-derived sequence.")
                        full_nuc_seq = line[col_map['FUSION_CDS']]
                        full_pep_seq = line[col_map['FUSION_TRANSL']]

                    if line[col_map['PROT_FUSION_TYPE']] == 'INFRAME':
                        valid_fusions[fus_md5]['type'] = 'INFRAME'
                        print('inframe')
                        print(full_pep_seq)
                        print(full_nuc_seq)
                        peptide = full_pep_seq[max(0, start_prot - 8):\
                                               min(len(full_pep_seq), start_prot + 8)]
                        peptide_context = full_pep_seq[max(0, start_prot - 16):\
                                                       min(len(full_pep_seq), start_prot + 16)]
                        nuc_context = full_nuc_seq[max(0, int(line[col_map['CDS_LEFT_RANGE']]\
                                                                  .split('-')[1]) - nuc_sub, 0):\
                                                   min(int(line[col_map['CDS_LEFT_RANGE']]\
                                                     .split('-')[1]) + nuc_sub, len(full_nuc_seq))]
                        nuc_seq = full_nuc_seq[max(0, int(line[col_map['CDS_LEFT_RANGE']]\
                                                                .split('-')[1]) - (nuc_sub - 24)):\
                                               min(int(line[col_map['CDS_LEFT_RANGE']]\
                                              .split('-')[1]) + (nuc_sub - 24 - 1), len(full_nuc_seq))]
                        print("peptide: {}".format(peptide))
                        print("peptide_context: {}".format(peptide_context))
                        print("nuc_context: {}".format(nuc_context))
                        print("nuc_seq: {}".format(nuc_seq))
                        if not re.search('\*', peptide):
                            valid_fusions[fus_md5]['peptide'] = peptide
                            valid_fusions[fus_md5]['peptide_context'] = peptide_context
                            valid_fusions[fus_md5]['nuc_context'] = nuc_context
                            valid_fusions[fus_md5]['nuc_seq'] = nuc_seq
                        else:
                            valid_fusions[fus_md5]['peptide'] = peptide
                            valid_fusions[fus_md5]['peptide_context'] = peptide_context
                            valid_fusions[fus_md5]['nuc_context'] = nuc_context
                            valid_fusions[fus_md5]['nuc_seq'] = nuc_seq

                    elif line[col_map['PROT_FUSION_TYPE']] == 'FRAMESHIFT':
                        valid_fusions[fus_md5]['type'] = 'FRAMESHIFT'
                        print('frameshift')
                        print(full_pep_seq)
                        print(start_prot)
                        init_peptide = full_pep_seq[max(0, start_prot - 8):]
                        init_peptide_context = full_pep_seq[max(0, start_prot - 16):]
                        init_nuc_context = full_nuc_seq[max(0, int(line[col_map['CDS_LEFT_RANGE']]\
                                                                 .split('-')[1]) - nuc_sub):]
                        print("init_nuc_context: {}".format(init_nuc_context))
                        
                        init_nuc_seq = full_nuc_seq[int(line[col_map['CDS_LEFT_RANGE']]\
                                                                 .split('-')[1]) - 24 - 1:]
                        print(init_peptide)
                        if re.search('\*', init_peptide):
                            print("Observed stop codon in frameshifted sequence")
                            stop_codon = init_peptide.index('*')
                            peptide = init_peptide[:stop_codon]
                            peptide_context = init_peptide_context[:stop_codon+8]
                            print("peptide: {}".format(peptide))
                            print("peptide_context: {}".format(peptide_context))
                            #valid_fusions[line[col_map['#FusionName']]]['peptide'] = peptide
                            valid_fusions[fus_md5]['peptide'] = peptide
                            valid_fusions[fus_md5]['peptide_context'] = peptide_context
                            valid_fusions[fus_md5]['nuc_context'] = init_nuc_context
                            valid_fusions[fus_md5]['nuc_seq'] = init_nuc_seq
                        else:
                            print("Did NOT observe stop codon in frameshifted sequence")
                            valid_fusions[fus_md5]['peptide'] = init_peptide
                            valid_fusions[fus_md5]['peptide_context'] = init_peptide_context
                            valid_fusions[fus_md5]['nuc_context'] = init_nuc_context
                            valid_fusions[fus_md5]['nuc_seq'] = init_nuc_seq

    with open(args.output, 'w') as ofo:
        for valid_fusion_id in valid_fusions.keys():
            print(valid_fusion_id)
            meta = valid_fusions[valid_fusion_id]
            print(meta)
            # Context sequences are too long.
            if meta['peptide']:
                ofo.write(">{}\n;antigen_source:FUSION fusion_id:{} fusion_type:{} nt_context:{} "\
                          .format(valid_fusion_id,
                                  valid_fusions[valid_fusion_id]['fusion_id'],
                                  valid_fusions[valid_fusion_id]['type'],
                                  #Dropping nuc_context now and replacing with meta['nuc_seq'].
                                  valid_fusions[valid_fusion_id]['nuc_context']))
                                  #meta['nuc_seq']))
                ofo.write("pep_context:{} fusion_left_breakpoint:{} fusion_right_breakpoint:{} "\
                          .format(valid_fusions[valid_fusion_id]['peptide_context'],
                                  valid_fusions[valid_fusion_id]['left_breakpoint'],
                                  valid_fusions[valid_fusion_id]['right_breakpoint']))
                ofo.write("fusion_left_gene:{} fusion_right_gene:{} fusion_left_transcript:{} "\
                          .format(valid_fusions[valid_fusion_id]['left_gn'],
                                  valid_fusions[valid_fusion_id]['right_gn'],
                                  valid_fusions[valid_fusion_id]['left_tx']))
                ofo.write("fusion_right_transcript:{} fusion_annotation:{}\n{}\n"\
                          .format(valid_fusions[valid_fusion_id]['right_tx'],
                                  valid_fusions[valid_fusion_id]['annot'],
                                  meta['peptide']))

    with open(args.nt_output, 'w') as ofo:
        for valid_fusion_id in valid_fusions.keys():
            print(valid_fusion_id)
            meta = valid_fusions[valid_fusion_id]
            print(meta)
            # Context sequences are too long.
            if meta['nuc_seq']:
                ofo.write(">{} FUSION_ID:{} NUCLEOTIDE_CONTEXT:{} PEPTIDE_CONTEXT:{}\n{}\n"\
                          .format(valid_fusion_id,
                                  valid_fusions[valid_fusion_id]['fusion_id'],
                                  valid_fusions[valid_fusion_id]['nuc_context'],
                                  valid_fusions[valid_fusion_id]['peptide_context'],
                                  meta['nuc_seq']))


def get_expressed_selfs_bed(args):
    """
    """
    expressed_self_transcripts = []
    with open(args.expressed_selfs) as exp_selfs_fo:
        for line in exp_selfs_fo.readlines():
            expressed_self_transcripts.append(line.rstrip().split(':')[1])

    with open(args.output, 'w') as output_fo:
        with open(args.gff) as gff_fo:
            for line in gff_fo.readlines():
                for expressed_self_transcript in expressed_self_transcripts:
                    if re.search(expressed_self_transcript, line):
                        if re.search('\tCDS\t', line):
                            parted_list = line.split('\t')
                            chrom, start, stop = (parted_list[0], parted_list[3], parted_list[4])
                            output_fo.write("{}\t{}\t{}\n".format(chrom, start, stop))


def get_expressed_transcripts_bed(args):
    """
    """
    expressed_transcripts = []
    with open(args.expressed_transcripts) as exp_tx_fo:
        for line in exp_tx_fo.readlines():
            expressed_transcripts.append(line.rstrip())

    with open(args.output, 'w') as output_fo:
        with open(args.gff) as gff_fo:
            for line in gff_fo.readlines():
                line = line.split('\t')
                if len(line) > 1 and line[2] == 'CDS':
                    for expressed_transcript in expressed_transcripts:
                        if re.search(expressed_transcript, line[8]):
                            chrom, start, stop = (line[0], line[3], line[4])
                            output_fo.write("{}\t{}\t{}\n".format(chrom, start, stop))


def get_expressed_ervs_bed(args):
    """
    """
    full_to_post_met = {}

    with open(args.geve_reference) as geve_fo:
        for line_idx, line in enumerate(geve_fo.readlines()):
            if line_idx != 0:
                line = line.split('\t')
                if line[8] != '-':
                    full_to_post_met[line[0]] = line[8].replace('.M', '')
    print(full_to_post_met)

    with open(args.expressed_ervs) as exp_ervs_fo:
        with open(args.output, 'w') as output_fo:
            for line_idx, line in enumerate(exp_ervs_fo.readlines()):
                if line_idx != 0:
                    print("line: {}".format(line))
                    if line.split(',')[0] in full_to_post_met.keys():
                        post_met = full_to_post_met[line.split(',')[0]]
                        print("post_met: {}".format(post_met))
                        chrom = post_met.split('.')[1]
                        first = post_met.split('.')[2]
                        second = post_met.split('.')[3]
                        output_fo.write("{}\t{}\t{}\n".format(chrom, first, second))


def get_expressed_viral_bed(args):
    """
    """

    viral_stops = {}

    for seq_record in SeqIO.parse(args.viral_cds_ref, "fasta"):
        viral_stops[seq_record.description.split()[0]] = len(seq_record.seq)

    with open(args.expressed_viruses) as exp_vir_fo:
        with open(args.output, 'w') as output_fo:
            for line in exp_vir_fo.readlines():
                output_fo.write("{}\t{}\t{}\n"\
                                .format(line.split()[0], 1, viral_stops[line.split()[0]]))


def numbers_to_ranges(data):
    s = pd.Series(data)
    v = s.astype(int).diff().bfill().ne(1).cumsum()
    ranges = (
        s.groupby(v).apply(
            lambda x: '-'.join(x.values[[0, -1]]) if len(x) > 1 else x.item()).tolist())
    return ranges


def filter_peptides(args):
    """
    """
    print("Loading CTAs...")
    ctas_gene_symbols = []
    with open(args.ctas, 'r') as ctafo:
        for line in ctafo.readlines():
            line = line.strip()
            ctas_gene_symbols.append(line)
    print("CTA List loaded.")

    print("Loading proteosome FASTA...")
#    pep_ref_seqs = [str(x) for x in parsers.parse_fasta(args.pep_ref).values()]
    pep_ref_seqs = parsers.parse_fasta_full_header(args.pep_ref)
    print("Original ref count: {}".format(len(pep_ref_seqs.keys())))
    print("Loaded proteosome FASTA...")

    print("Excluding CTAs from user-provided list...")

    no_cta_pep_ref_seqs = {}

    for k,v in pep_ref_seqs.items():
        if k.split('|')[-2] not in ctas_gene_symbols:
            no_cta_pep_ref_seqs[k] = v
    print("Excluded CTAs from user-provided list...")
    print("CTA-less ref count: {}".format(len(no_cta_pep_ref_seqs.keys())))

    pep_col = ''

    output_lines = []


    with open(args.input_report, 'r') as ifo:
        for line_idx, line in enumerate(ifo.readlines()):
            line = line.strip()
            line = line.split('\t')
            if line_idx == 0: 
                line.extend(['all_transcript_ids_encoding_peptide',
                             'all_gene_ids_encoding_peptide',
                             'all_gene_names_encoding_peptide'])
                output_lines.append(line)
                pep_col = line.index('peptide')
                src_col = line.index('antigen_source')
            else:
                if any([re.search(line[pep_col], v) for v in no_cta_pep_ref_seqs.values()]):
                    pass
                else:
                    if line[src_col] in ["CTA/SELF", "Self-Antigen"]:
                        print("Found CTA line.")
                        peptide_tx_ids = []
                        peptide_gene_ids = []
                        peptide_gene_symbols = []
                        hits = [k for k,v in pep_ref_seqs.items() if re.search(line[pep_col], v)]
                        if hits:
                            peptide_tx_ids = ','.join(list(set([x.split('|')[1] for x in hits])))
                            peptide_gene_ids = ','.join(list(set([x.split('|')[2] for x in hits])))
                            peptide_gene_symbols = ','.join(list(set([x.split('|')[-2] for x in hits])))
                            line.extend([peptide_tx_ids, peptide_gene_ids, peptide_gene_symbols]) 
                            output_lines.append(line)
                        else:
                            line.extend(['NA', 'NA', 'NA']) 
                            output_lines.append(line)
                        
                    else:
                        line.extend(['NA', 'NA', 'NA']) 
                        output_lines.append(line)

    with open(args.output, 'w') as ofo:
        for output_line in output_lines:
            print(output_line)
            ofo.write("{}\n".format('\t'.join(output_line)))


def prioritize_peptides(args):
    """
    """
    output_lines = []

    max_reads = 0.0
    prim_alns_max_reads = 0.0

    with open(args.input_report, 'r') as ifo:
        for line_idx, line in enumerate(ifo.readlines()):
            line = line.split('\t')
            if line_idx == 0:
                read_col = line.index('rna_reads_covering_genomic_origin_with_peptide_cds')
                prim_alns_read_col = line.index('primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds')
            else:
                print(line[read_col])
                print(line[prim_alns_read_col])

                transformed = log2(float(line[read_col]) + 1)
                if transformed > max_reads:
                    max_reads = log2(float(line[read_col]) + 1)

                ra_transformed = log2(float(line[prim_alns_read_col]) + 1)
                if ra_transformed > prim_alns_max_reads:
                    prim_alns_max_reads = log2(float(line[prim_alns_read_col]) + 1)

    with open(args.input_report, 'r') as ifo2:
        mhcflurry_data = False
        netmhcpan_data = False
        for line_idx, line in enumerate(ifo2.readlines()):
            line = line.strip()
            line = line.split('\t')
            if line_idx == 0:
#                line.append('priority_score')
#                line.append('priority_score_no_ccf')
                read_col = line.index('rna_reads_covering_genomic_origin_with_peptide_cds')
                prim_alns_read_col = line.index('primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds')
                if 'ccf' in line:
                    ccf_col = line.index('ccf')
                elif 'cellular_prevalence' in line:
                    ccf_col = line.index('cellular_prevalence')
                else:
                    ccf_col = ''
                if [x for x in line if re.search('netmhcpan.*aff_nm', x)]:
                    netmhcpan_ba_col = line.index([x for x in line if re.search('netmhcpan.*aff_nm', x)][0])
                    netmhcpan_data = True
                if [x for x in line if re.search('mhcflurry.*aff$', x)]:
                    mhcflurry_ba_col = line.index([x for x in line if re.search('mhcflurry.*aff$', x)][0])
                    mhcflurry_data = True

                if mhcflurry_data:
                    line.append('priority_score_mhcflurry')
                    line.append('priority_score_mhcflurry_no_ccf')
                    line.append('priority_score_mhcflurry_prim_alns')
                    line.append('priority_score_mhcflurry_no_ccf_prim_alns')
                if netmhcpan_data:
                    line.append('priority_score_netmhcpan')
                    line.append('priority_score_netmhcpan_no_ccf')
                    line.append('priority_score_netmhcpan_prim_alns')
                    line.append('priority_score_netmhcpan_no_ccf_prim_alns')
                line.append('priority_score_maximum')
            

            else:
                all_scores = []
                if mhcflurry_data:
                    if ccf_col and line[ccf_col] != 'NA' and line[read_col] != 'NA' and line[mhcflurry_ba_col] != 'NA':
                        print("Non-zero CCF")
                        print("priority_score_mhcflurry:{} {} {}".format(float(line[ccf_col]), \
                                                                  float(log2(float(line[read_col]) + 1)/max_reads), \
                                                                  abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)))
                        pr_scr = float(line[ccf_col]) * \
                                 float(log2(float(line[read_col]) + 1)/max_reads) * \
                                 abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)
    
                        pr_scr_non_ccf = float(log2(float(line[read_col]) + 1)/max_reads) * \
                                               abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)
                        line.append(pr_scr)
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr)
                        all_scores.append(pr_scr_non_ccf)
                    elif (ccf_col == '' or line[ccf_col] == 'NA') and line[read_col] != 'NA' and line[mhcflurry_ba_col] != 'NA':
                        print("N/A CCF")
                        print("{} {}".format(float(log2(float(line[read_col]) + 1)/max_reads), \
                                         abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)))
                        pr_scr_non_ccf = float(log2(float(line[read_col]) + 1)/max_reads) * \
                                         abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)
                        line.append("NA")
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr_non_ccf)
                    else:
                      print("Insufficient metrics")
                      line.append(["NA", "NA"])
                    if ccf_col and line[ccf_col] != 'NA' and line[prim_alns_read_col] != 'NA' and line[mhcflurry_ba_col] != 'NA':
                        print("Non-zero CCF")
                        print("priority_score_mhcflurry:{} {} {}".format(float(line[ccf_col]), \
                                                                  float(log2(float(line[prim_alns_read_col]) + 1)/max_reads), \
                                                                  abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)))
                        pr_scr = float(line[ccf_col]) * \
                                 float(log2(float(line[prim_alns_read_col]) + 1)/max_reads) * \
                                 abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)
    
                        pr_scr_non_ccf = float(log2(float(line[prim_alns_read_col]) + 1)/max_reads) * \
                                               abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)
                        line.append(pr_scr)
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr)
                        all_scores.append(pr_scr_non_ccf)
                    elif (ccf_col == '' or line[ccf_col] == 'NA') and line[prim_alns_read_col] != 'NA' and line[mhcflurry_ba_col] != 'NA':
                        print("N/A CCF")
                        print("{} {}".format(float(log2(float(line[prim_alns_read_col]) + 1)/max_reads), \
                                         abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)))
                        pr_scr_non_ccf = float(log2(float(line[prim_alns_read_col]) + 1)/max_reads) * \
                                         abs((float(line[mhcflurry_ba_col]) - 1000.0)/1000.0)
                        line.append("NA")
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr_non_ccf)
                    else:
                        print("Insufficient metrics")
                        line.append(["NA", "NA"])
                if netmhcpan_data:
                    if ccf_col and line[ccf_col] != 'NA' and line[read_col] != 'NA' and line[netmhcpan_ba_col] != 'NA':
                        print("Non-zero CCF")
                        print("priority_score_netmhcpan:{} {} {}".format(float(line[ccf_col]), \
                                                                  float(log2(float(line[read_col]) + 1)/max_reads), \
                                                                  abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)))
                        pr_scr = float(line[ccf_col]) * \
                                 float(log2(float(line[read_col]) + 1)/max_reads) * \
                                 abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)
    
                        pr_scr_non_ccf = float(log2(float(line[read_col]) + 1)/max_reads) * \
                                               abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)
                        line.append(pr_scr)
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr)
                        all_scores.append(pr_scr_non_ccf)
                    elif (ccf_col == '' or line[ccf_col] == 'NA') and line[read_col] != 'NA' and line[netmhcpan_ba_col] != 'NA':
                        print("N/A CCF")
                        print("{} {}".format(float(log2(float(line[read_col]) + 1)/max_reads), \
                                         abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)))
                        pr_scr_non_ccf = float(log2(float(line[read_col]) + 1)/max_reads) * \
                                         abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)
                        line.append("NA")
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr_non_ccf)
                    else:
                      print("Insufficient metrics")
                      line.append(["NA", "NA"])
                    if ccf_col and line[ccf_col] != 'NA' and line[prim_alns_read_col] != 'NA' and line[netmhcpan_ba_col] != 'NA':
                        print("Non-zero CCF")
                        print("priority_score_netmhcpan:{} {} {}".format(float(line[ccf_col]), \
                                                                  float(log2(float(line[prim_alns_read_col]) + 1)/max_reads), \
                                                                  abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)))
                        pr_scr = float(line[ccf_col]) * \
                                 float(log2(float(line[prim_alns_read_col]) + 1)/max_reads) * \
                                 abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)
    
                        pr_scr_non_ccf = float(log2(float(line[prim_alns_read_col]) + 1)/max_reads) * \
                                               abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)
                        line.append(pr_scr)
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr)
                        all_scores.append(pr_scr_non_ccf)
                    elif (ccf_col == '' or line[ccf_col] == 'NA') and line[prim_alns_read_col] != 'NA' and line[netmhcpan_ba_col] != 'NA':
                        print("N/A CCF")
                        print("{} {}".format(float(log2(float(line[prim_alns_read_col]) + 1)/max_reads), \
                                         abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)))
                        pr_scr_non_ccf = float(log2(float(line[prim_alns_read_col]) + 1)/max_reads) * \
                                         abs((float(line[netmhcpan_ba_col]) - 1000.0)/1000.0)
                        line.append("NA")
                        line.append(pr_scr_non_ccf)
                        all_scores.append(pr_scr_non_ccf)
                    else:
                        print("Insufficient metrics")
                        line.append(["NA", "NA"])
                line.append(max(all_scores))

            output_lines.append('\t'.join([str(x) for x in line]))

    with open(args.output, 'w') as ofo:
        for output_line in output_lines:
            ofo.write("{}\n".format(output_line))
