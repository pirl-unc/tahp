import re
from Bio.SeqUtils import seq1
import vcf
import numpy as np

from src import parsers
from src import rna_utils



def extract_missense_snvs(input_vcf):
    """
    Extracts missense SNVs from annotated somatic VCF file.
    Assumes annotated VCF follows ANN standard:
    http://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf

    Args:
        input_vcf: Annotated VCF.

    Returns:
        missense_records (dict): keys are <chr>:<pos> and values are
                                 dictionaries with variant metadata.
    """
    missense_records = {}

    vcf_reader = parsers.load_vcf(input_vcf)

    for record in vcf_reader:
        for annotation in record.INFO['ANN']:
            effects = annotation.split('|')
            variant_effect = effects[1]
            if variant_effect == 'missense_variant' and \
                    effects[13] and \
                    effects[-1] != 'WARNING_TRANSCRIPT_NO_START_CODON':
                transcript = effects[6]

                nt_change = effects[9].lstrip('c.')
                nt_change = nt_change[-3:].partition('>')

                ref_nt = nt_change[0]
                alt_nt = nt_change[2]

                nt_pos = effects[12].split('/')[0]
                nt_len = effects[12].split('/')[1]

                cdna_pos = effects[11].split('/')[0]

                aa3_change = effects[10].lstrip('p.')
                aa3_change = re.split('(\d+)', aa3_change)

                ref_aa3 = aa3_change[0]
                alt_aa3 = aa3_change[2]

                aa_pos = aa3_change[1]

                ref_aa = seq1(ref_aa3)
                alt_aa = seq1(alt_aa3)
                aa_len = effects[13].split('/')[1]

                codon_pos = int((int(effects[13].split('/')[0])*3) - 2)


                #Variant-specific metadata directory
                var_meta = {'transcript': transcript,
                            'cdna_pos': cdna_pos,
                            'nt_pos': nt_pos,
                            'nt_len': nt_len,
                            'ref_nt': ref_nt,
                            'alt_nt': alt_nt,
                            'aa_pos': aa_pos,
                            'aa_len': aa_len,
                            'ref_aa': ref_aa,
                            'alt_aa': alt_aa,
                            'codon_pos': codon_pos,
                            'meta': record}

                if var_meta:
                    if "{}:{}".format(record.CHROM, record.POS) in missense_records.keys():
                        missense_records["{}:{}".format(record.CHROM, record.POS)].append(var_meta)
                    else:
                        missense_records["{}:{}".format(record.CHROM, record.POS)] = [var_meta]

    return missense_records


def extract_inframe_indels(input_vcf, indel_type):
    """
    Extracts conservative inframe InDels from annotated somatic VCF file.

    Args:
        input_vcf: Annotated VCF.
        indel_type(str): Type of in-frame InDels to capture (conservative or disruptive)

    Returns:
        inframe_indels(dict): keys are <chr>:<pos> and values are
                              dictionaries with variant metadata.
    """
    inframe_indels = {}

    vcf_reader = parsers.load_vcf(input_vcf)

    capture_regex = ''
    if indel_type == 'conservative':
        capture_regex = '^conservative_inframe_[a-z]+$'
    elif indel_type == 'disruptive':
        capture_regex = '^disruptive_inframe_[a-z]+$'

    for record in vcf_reader:
        inframe_indels["{}:{}".format(record.CHROM, record.POS)] = []
        for annotation in record.INFO['ANN']:
            effects = annotation.split('|')
            if re.search(capture_regex, effects[1]) and \
                    effects[13] and \
                    effects[-1] != 'WARNING_TRANSCRIPT_NO_START_CODON':
                transcript = effects[6]
                aa3_change = effects[10].lstrip('p.')
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                cdna_pos = effects[11].split('/')[0]

                #Variant-specific metadata
                var_meta = {'transcript': transcript,
                            'aa3_change': aa3_change,
                            'nt_change': nt_change,
                            'aa_len': aa_len,
                            'nt_len': nt_len,
                            'cdna_pos': cdna_pos,
                            'meta': record}

                inframe_indels["{}:{}".format(record.CHROM, record.POS)].append(var_meta)

    # Removing InDels without annotation metadata.
    ks_to_del = []
    for inframe_indel, annot in inframe_indels.items():
        if not annot:
            ks_to_del.append(inframe_indel)
    for inframe_indel in ks_to_del:
        del inframe_indels[inframe_indel]

    return inframe_indels


def extract_frameshift_indels(input_vcf):
    """
    Focusing on deletions now, need to find a good example of a
    conservative_inframe_insertion
    """
    frameshift_indels = {}
    vcf_reader = parsers.load_vcf(input_vcf)

    for record in vcf_reader:
        frameshift_indels["{}:{}".format(record.CHROM, record.POS)] = []
        for annotation in record.INFO['ANN']:
            effects = annotation.split('|')
            if (re.search('^frameshift_variant$', effects[1]) and
                effects[13] and effects[7] == 'protein_coding' and
                effects[-1] != 'WARNING_TRANSCRIPT_NO_START_CODON'):
                transcript = effects[6]
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                cdna_pos = effects[11].split('/')[0]

                #Variant-specific metadata
                var_meta = {'transcript': transcript,
                            'nt_change': nt_change,
                            'aa_len': aa_len,
                            'nt_len': nt_len,
                            'cdna_pos': cdna_pos,
                            'meta': record}
                frameshift_indels["{}:{}".format(record.CHROM, record.POS)].append(var_meta)

    # Removing InDels without annotation metadata.
    ks_to_del = []
    for frameshift_indel, annot in frameshift_indels.items():
        if not annot:
            ks_to_del.append(frameshift_indel)
    for frameshift_indel in ks_to_del:
        del frameshift_indels[frameshift_indel]
    return frameshift_indels

def expressed_variants(args):
    """
    Given an annotated VCF, quants.sf file, and expression percentile (or
    threshold), emit SNV/InDel-harboring transcripts that exceed expression
    percentile and variants associated with those transcripts..

    Args:

    Returns:

    """
    tx_abundances = rna_utils.load_tx_abundances(args.quants, args.metric, args.exclude_zeros)

    tx_threshold = 0
    if not args.abundance_threshold:
        tx_threshold = rna_utils.get_tx_threshold(args.percentile, tx_abundances)
        print("Expression Threshold: {}".format(tx_threshold))
    else:
        tx_threshold = float(args.abundance_threshold)
        print("Expression Threshold: {}".format(tx_threshold))
    expressed_txids = rna_utils.get_expressed_txs(args.quants, args.metric, tx_threshold)
    print("# of expressed transcripts: {}".format(len(expressed_txids)))
    filtered_records, somatic_transcripts = filter_vcf_by_expression(args.vcf, expressed_txids)
    print("# of filtered records: {}".format(len(filtered_records)))

    write_expressed_vcf(args.vcf, args.output, filtered_records)
    with open(args.somatic_txs, 'w') as output_fo:
        for transcript in sorted(somatic_transcripts):
            output_fo.write('"{}"\n'.format(transcript))

def filter_vcf_by_expression(input_vcf, expressed_txids):
    """
    Given an annotated VCF and lis tof expressed transcript identifiers, return
    variant records associated with those transcripts and subset of transcripts
    harboring SNV/InDels.
    """
    somatic_transcripts = []
    filtered_records = []
    vcf_reader = parsers.load_vcf(input_vcf)
    for record in vcf_reader:
        annotations = record.INFO['ANN']
        for annotation in annotations:
            effects = annotation.split('|')
            transcript = effects[6]
            if transcript in expressed_txids and record not in filtered_records:
                filtered_records.append(record)
                if transcript not in somatic_transcripts:
                    somatic_transcripts.append(transcript)

    return filtered_records, list(set(somatic_transcripts))


def write_expressed_vcf(input_vcf, output_vcf, filtered_records):
    """
    """
    vcf_reader = parsers.load_vcf(input_vcf)
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)
    for filtered_record in filtered_records:
        vcf_writer.write_record(filtered_record)


def write_isolated_vcf(args, filtered_records):
    """
    """
    vcf_reader = vcf.Reader(filename=args.somatic_vcf, compressed=True)
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    for filtered_record in filtered_records:
        vcf_writer.write_record(filtered_record)


def isolated_variants(args):
    """
    """
    cand_vars = get_candidate_variants(args.somatic_vcf)
    print("# of candidate variants: {}".format(len(cand_vars)))
    isolated_vars = get_all_isolated_variants(args)
    print("# of isolated variants: {}".format(len(isolated_vars)))
    isolated_cand_vars = []
    crowded_cand_vars = []
    for record in cand_vars:
        if record in isolated_vars:
            isolated_cand_vars.append(record)
        else:
            crowded_cand_vars.append(record)
    print("# of isolated candidate variants: {}"
          .format(len(isolated_cand_vars)))
    write_isolated_vcf(args, isolated_cand_vars)


def get_candidate_variants(somatic_vcf):
    """
    """
    cand_vars = []
    vcf_reader = vcf.Reader(filename=somatic_vcf, compressed=True)
    for record in vcf_reader:
        cand_vars.append(record)
    return cand_vars


# Can likely remove this function.
def get_all_isolated_variants(args):
    """Assumes variants are sorted for now.
    """
    isolated_vars = []
    germline_vcf_reader = vcf.Reader(filename=args.germline_vcf, compressed=True)
    somatic_vcf_reader = vcf.Reader(filename=args.somatic_vcf, compressed=True)
    all_records = []
    for record_idx, record in enumerate(germline_vcf_reader):
        all_records.append(record)
    for record_idx, record in enumerate(somatic_vcf_reader):
        all_records.append(record)
    print("# of total variants: {}".format(len(all_records)))

    #A few thoughts here - should really only be checking for germline variants
    #around somatic variants. It may be worth checking for proximal somatic
    #variants, but should probably result in a warning, if anything. This is also
    #where the het phased and hom/het synonymous checks can occur. The incorporation
    #of het phased variants will have to be in the make_*_peptides step.
    for record_idx, record in enumerate(all_records[:-1]):
        upstream_record = all_records[record_idx - 1]
        downstream_record = all_records[record_idx + 1]

        upstream_isolated = False
        downstream_isolated = False

        if upstream_record.CHROM == record.CHROM:
            if np.absolute(upstream_record.POS - record.POS) > int(args.proximity):
                upstream_isolated = True
        else:
            upstream_isolated = True

        if downstream_record.CHROM == record.CHROM:
            if np.absolute(downstream_record.POS - record.POS) > int(args.proximity):
                downstream_isolated = True
        else:
            downstream_isolated = True


        if upstream_isolated and downstream_isolated:
            isolated_vars.append(record)

    return isolated_vars
