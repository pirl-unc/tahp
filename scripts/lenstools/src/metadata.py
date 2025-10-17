"""
A module for handling pMHC metadata in LENSTools.
"""

from glob import glob
import gzip
import json
import re
import os
from pprint import pprint
import pandas as pd
from Bio import SeqIO

from src import parsers

def add_tcga_data(args):
    """
    """

    abbrev_to_type = parsers.load_tcga_dicts()

    tumor_summs = parsers.load_tcga_tx_summ(args.tumor_type, args.tcga_transcript_summary)
    tumor_summs_nov = {}

    all_tumor_spec_means = [v['mean'] for k,v in tumor_summs.items()]

    # The transcript version from the TCGA reference may differ from the
    # version in the GTF, so want to match transcripts without worrying about
    # version conflicts.
    initial_txs = tumor_summs.keys()
    for transcript in initial_txs:
        tx_nov = transcript.split('.')[0]
        tumor_summs_nov[tx_nov] = tumor_summs[transcript]

    output_lines = []

    # Loading LENS report.
    col_idx_map = {}
    with open(args.report) as report_fo:
        for line_idx, line in enumerate(report_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    col_idx_map[col] = col_idx
                line.extend(['tcga_tpm_mean', 'tcga_tpm_median', 'tcga_tpm_max',
                             'tcga_tpm_iqr', 'tcga_tpm_percentile'])
                output_lines.append(line)
            else:
                if line[col_idx_map['antigen_source']] in ['SNV', 'InDel']:
                    tx_id = line[col_idx_map['transcript_identifier']].split('.')[0]
                    if tx_id in tumor_summs_nov.keys():
                        tumor_specs_gt_mean = sum([1 for x in all_tumor_spec_means if \
                                              float(x) < float(tumor_summs_nov[tx_id]['mean'])])
                        percentile = (tumor_specs_gt_mean / float(len(all_tumor_spec_means)))
                        line.extend([tumor_summs_nov[tx_id]['mean'],
                                     tumor_summs_nov[tx_id]['median'],
                                     tumor_summs_nov[tx_id]['max'],
                                     tumor_summs_nov[tx_id]['iqr'],
                                     percentile])
                else:
                    line.extend(['NA', 'NA', 'NA', 'NA', 'NA'])
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        for line in output_lines:
            ofo.write("{}\n".format('\t'.join([str(x) for x in line])))


def make_lens_report(args):
    """
    """
    reports = glob(os.path.join(args.metadata_dir, "*metadata*"))

    #Initiating dataframe with snvs...

    report_df = pd.DataFrame()
    for report in reports:
        tmp_report = ''
        tmp_report = pd.read_csv(report, sep='\t')
        if re.search('self_antigen', report):
            tmp_report = tmp_report[tmp_report["binding_affinity"] <= 50.0]
        if not(tmp_report.empty) and report_df.empty:
            report_df = tmp_report
        elif not tmp_report.empty:
            report_df = pd.concat([report_df, tmp_report])

    report_df.astype({'binding_affinity': 'float64'})

    filtered_df = report_df[report_df["binding_affinity"] <= 1000.0]
    filtered_df.to_csv(args.output, sep='\t', index=False, na_rep='NA')


def aggregate_lens_reports(args):
    """
    """
    reports = glob(os.path.join(args.reports_dir, "*report*"))

    report_df = pd.DataFrame()
    for report in reports:
        dataset = ''.join(report.split('/')[-1].split('.')[0].split('-')[0])
        pat_id = ''.join(report.split('/')[-1].split('.')[0].split('-')[1])
        samples = '-'.join(report.split('/')[-1].split('.')[0].split('-')[2:])
        print("{}\t{}\t{}\t".format(dataset, pat_id, samples))
        tmp_report = ''
        tmp_report = pd.read_csv(report, sep='\t')
        tmp_report['dataset'] = dataset
        tmp_report['patient_identifier'] = pat_id
        tmp_report['sample_identifiers'] = samples
        report_df = pd.concat([report_df, tmp_report], join='outer')
    report_df.to_csv(args.output, sep='\t', index=False, na_rep='NA')


def aggregate_pmhc_outputs(args):
    """
    """
    pmhc_summs = {}

    output_header = []

    # Getting combinations of pmhcs from netmhcpan, mhcflurry, and
    # netmhcstabpan, in that order.

    if args.netmhcpan:
        print("Detected NetMHCpan input. Using NetMHCpan as anchor point for pMHCs.")
        hdr_map = {}
        with open(args.netmhcpan, 'r') as nfo:
            for line_idx, line in enumerate(nfo.readlines()):
                if line_idx == 0:
                    line = line.rstrip().split('\t')
                    hdr_map = {col:idx for idx, col in enumerate(line)}
                else:
                    line = line.rstrip().split('\t')
                    id = '-'.join([line[hdr_map['identity']], \
                                   line[hdr_map['peptide']], \
                                   line[hdr_map['allele']], \
                                   line[hdr_map['pos']]])
                    pmhc_summs[id] = {}
                    pmhc_summs[id]['identity'] = line[hdr_map['identity']]
                    pmhc_summs[id]['peptide'] = line[hdr_map['peptide']]
                    pmhc_summs[id]['allele'] = line[hdr_map['allele']]
                    pmhc_summs[id]['pos'] = line[hdr_map['pos']]

    if args.netmhcpan:
        print("Detected NetMHCpan input. Including NetMHCpan pMHC data.")
        netmhcpan_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.netmhcpan, \
                                                                                 output_header)
        for key in keys_to_add:
            for id in pmhc_summs.keys():
                netmhcpan_id = "-".join([pmhc_summs[id]['identity'], \
                                         pmhc_summs[id]['peptide'], \
                                         pmhc_summs[id]['allele']])
                pmhc_summs[id][key] = netmhcpan_summs[netmhcpan_id][key]


    if args.netctlpan:
        print("Detected NetCTLpan input. Including NetCTLpan pMHC data.")
        netctlpan_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.netctlpan, \
                                                                                 output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            for id in pmhc_summs.keys():
                netctlpan_id = "-".join([pmhc_summs[id]['identity'], \
                                         pmhc_summs[id]['peptide'], \
                                         pmhc_summs[id]['allele']])
                if netctlpan_id in netctlpan_summs.keys():
                    pmhc_summs[id][key] = netctlpan_summs[netctlpan_id][key]
                else:
                    pmhc_summs[id][key] = 'NA'

    if args.netmhcstabpan:
        print("Detected NetMHCstabpan input. Including NetMHCstabpan pMHC data.")
        netmhcstabpan_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.netmhcstabpan, \
                                                                                     output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            for id in pmhc_summs.keys():
                netmhcstabpan_id = "-".join([pmhc_summs[id]['identity'], \
                                             pmhc_summs[id]['peptide'], \
                                             pmhc_summs[id]['allele']])
                if netmhcstabpan_id in netmhcstabpan_summs.keys():
                    pmhc_summs[id][key] = netmhcstabpan_summs[netmhcstabpan_id][key]
                else:
                    pmhc_summs[id][key] = 'NA'

    if args.mhcflurry and args.netmhcpan:
        print("Detected MHCFlurry input. Including MHCFlurry pMHC data.")
        mhcflurry_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.mhcflurry, \
                                                                                 output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            for id in pmhc_summs.keys():
                mhcflurry_id = "-".join([pmhc_summs[id]['identity'], \
                                         pmhc_summs[id]['peptide'], \
                                         pmhc_summs[id]['allele']])
                if mhcflurry_id in mhcflurry_summs.keys():
                    pmhc_summs[id][key] = mhcflurry_summs[mhcflurry_id][key]
                else:
                    pmhc_summs[id][key] = 'NA'
    elif args.mhcflurry and not args.netmhcpan:
        print("Detected MHCFlurry input, but not NetMHCpan input. Including MHCFlurry pMHC data.")
        mhcflurry_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.mhcflurry, \
                                                                                 output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            print(list(mhcflurry_summs.keys())[:10])
            for mhcflurry_id in mhcflurry_summs.keys():
                if mhcflurry_id not in pmhc_summs.keys():
                    pmhc_summs[mhcflurry_id] = {}
                pmhc_summs[mhcflurry_id][key] = mhcflurry_summs[mhcflurry_id][key]
        with open(args.mhcflurry, 'r') as nfo:
            for line_idx, line in enumerate(nfo.readlines()):
                if line_idx == 0:
                    line = line.rstrip().split('\t')
                    hdr_map = {col:idx for idx, col in enumerate(line)}
                else:
                    line = line.rstrip().split('\t')
                    id = '-'.join([line[hdr_map['identity']], \
                                   line[hdr_map['peptide']], \
                                   line[hdr_map['allele']]])
                    pmhc_summs[id]['identity'] = line[hdr_map['identity']]
                    pmhc_summs[id]['peptide'] = line[hdr_map['peptide']]
                    pmhc_summs[id]['allele'] = line[hdr_map['allele']]
                    pmhc_summs[id]['pos'] = 'NA'
        example_set = list(pmhc_summs.keys())[:10]
        for k in example_set:
            print("{}: {}".format(k, pmhc_summs[k]))

    if args.deephlapan:
        print("Detected DeepHLAPan input. Including DeepHLAPan pMHC data.")
        deephlapan_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.deephlapan, \
                                                                                  output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            for id in pmhc_summs.keys():
                deephlapan_id = "-".join([pmhc_summs[id]['identity'], \
                                          pmhc_summs[id]['peptide'], \
                                          pmhc_summs[id]['allele']])
                if deephlapan_id in deephlapan_summs.keys():
                    pmhc_summs[id][key] = deephlapan_summs[deephlapan_id][key]
                else:
                    pmhc_summs[id][key] = 'NA'
    if args.hlapollo:
        print("Detected HLApollo input.")
        hlapollo_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.hlapollo, \
                                                                                 output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            print(list(hlapollo_summs.keys())[:10])
            for hlapollo_id in hlapollo_summs.keys():
                if hlapollo_id not in pmhc_summs.keys():
                    print("hlapollo id: {}".format(hlapollo_id))
                    pmhc_summs[hlapollo_id] = {}
                pmhc_summs[hlapollo_id][key] = hlapollo_summs[hlapollo_id][key]
        with open(args.hlapollo, 'r') as nfo:
            for line_idx, line in enumerate(nfo.readlines()):
                if line_idx == 0:
                    line = line.rstrip().split('\t')
                    hdr_map = {col:idx for idx, col in enumerate(line)}
                else:
                    line = line.rstrip().split('\t')
                    id = '-'.join([line[hdr_map['identity']], \
                                   line[hdr_map['peptide']],
                                   line[hdr_map['allele']]])
                    pmhc_summs[id]['identity'] = line[hdr_map['identity']]
                    pmhc_summs[id]['peptide'] = line[hdr_map['peptide']]
                    pmhc_summs[id]['allele'] = line[hdr_map['allele']]
        example_set = list(pmhc_summs.keys())[:10]
        for k in example_set:
            print("{}: {}".format(k, pmhc_summs[k]))

    if args.pepsickle:
        print("Detected pepsickle input.")
        pepsickle_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.pepsickle, \
                                                                                 output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            print(list(pepsickle_summs.keys())[:10])
            for pepsickle_id in pepsickle_summs.keys():
                pepsickle_ident = pepsickle_id.split('-')[0]
                pepsickle_pep = pepsickle_id.split('-')[1]
                for pmhc_summ_id in pmhc_summs.keys():
                    pmhc_ident = pmhc_summ_id.split('-')[0]
                    pmhc_pep = pmhc_summ_id.split('-')[1]
                    if pepsickle_ident == pmhc_ident and pepsickle_pep == pmhc_pep:
#                    if pepsickle_id not in pmhc_summs.keys():
#                        pmhc_summs[pepsickle_id] = {}
#                    pmhc_summs[pepsickle_id][key] = pepsickle_summs[pepsickle_id][key]
                        pmhc_summs[pmhc_summ_id][key] = pepsickle_summs[pepsickle_id][key]
#        with open(args.pepsickle, 'r') as nfo:
#            for line_idx, line in enumerate(nfo.readlines()):
#                if line_idx == 0:
#                    line = line.rstrip().split('\t')
#                    hdr_map = {col:idx for idx, col in enumerate(line)}
#                else:
#                    line = line.rstrip().split('\t')
#                    id = '-'.join([line[hdr_map['identity']], \
#                                   line[hdr_map['peptide']]])
#                    pmhc_summs[id]['identity'] = line[hdr_map['identity']]
#                    pmhc_summs[id]['peptide'] = line[hdr_map['peptide']]
#        example_set = list(pmhc_summs.keys())[:10]
#        for k in example_set:
#            print("{}: {}".format(k, pmhc_summs[k]))

    if args.ag_dissim:
        print("Detected antigen.garnish dissimilarity input.")
        print("Including antigen.garnish dissimilarity pMHC data.")
        ag_dissim_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.ag_dissim, \
                                                                                 output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            for id in pmhc_summs.keys():
                ag_dissim_id = pmhc_summs[id]['peptide']
                if ag_dissim_id in ag_dissim_summs.keys():
                    pmhc_summs[id][key] = ag_dissim_summs[ag_dissim_id][key]
                else:
                    pmhc_summs[id][key] = 'NA'

    if args.ag_foreign:
        print("Detected antigen.garnish foreignness input.")
        print("Including antigen.garnish foreignness pMHC data.")
        ag_foreign_summs, output_header, keys_to_add = parse_formatted_pmhc_summs(args.ag_foreign, \
                                                                                  output_header)
        for key in keys_to_add:
            print(f"Adding {key}")
            for id in pmhc_summs.keys():
                ag_foreign_id = pmhc_summs[id]['peptide']
                if ag_foreign_id in ag_foreign_summs.keys():
                    pmhc_summs[id][key] = ag_foreign_summs[ag_foreign_id][key]
                else:
                    pmhc_summs[id][key] = 'NA'

    pprint(pmhc_summs)


    write_final_pmhc_summs(pmhc_summs, output_header, args.output)


def parse_formatted_pmhc_summs(pmhc_summ_inf, output_header):
    keys_of_interest = []
    file_pmhc_summs = {}
    hdr_map = {}
    with open(pmhc_summ_inf, 'r') as nfo:
        for line_idx, line in enumerate(nfo.readlines()):
            if line_idx == 0:
                line = line.rstrip().split('\t')
                hdr_map = {col:idx for idx, col in enumerate(line)}
                output_header.extend([x for x in line if x not in output_header])
                keys_of_interest = [x for x in hdr_map.keys() if x not in ['allele', 'peptide', \
                                                                           'pos', 'identity']]
            else:
                line = line.rstrip().split('\t')
                if 'allele' in hdr_map.keys() and 'identity' in hdr_map.keys() and \
                    'peptide' in hdr_map.keys():
                    ident = '-'.join([line[hdr_map['identity']], \
                                      line[hdr_map['peptide']], \
                                      line[hdr_map['allele']]])
                elif 'peptide' and 'identity' in hdr_map.keys():
                    ident = '-'.join([line[hdr_map['identity']], \
                                      line[hdr_map['peptide']]])
                elif 'peptide' in hdr_map.keys():
                    ident = line[hdr_map['peptide']]
                print(ident)
                file_pmhc_summs[ident] = {}
                for key in hdr_map.keys():
                    file_pmhc_summs[ident][key] = line[hdr_map[key]]
    return file_pmhc_summs, output_header, keys_of_interest


def write_final_pmhc_summs(pmhc_summs, output_header, output_file):
    """
    """
    with open(output_file, 'w') as outf:
        outf.write('\t'.join(output_header))
        outf.write('\n')
        for pmhc_id in pmhc_summs.keys():
            outf.write('\t'.join([pmhc_summs[pmhc_id][i] for i in output_header]))
            outf.write('\n')

def annot_pmhcs(args):
    """
    """
    pmhcs_meta = {}
    orig_hdr = ''
    hdr_map = {}
    hdr_to_add = []
    with open(args.pmhc_fasta) as ffo:
        for line in ffo:
            if re.search('^>', line):
                ident = line.rstrip().lstrip('>')
                nline = next(ffo)
                nline = nline.lstrip(';').split(' ')
                print(nline)
                metad = {x.partition(':')[0]:x.partition(':')[2].strip() for x in nline}
                hdr_to_add.extend([x.split(':')[0].strip() for x in nline if x not in hdr_to_add])
                pmhcs_meta[ident] = metad

    annotated_lines = []

    hdr_to_add = sorted(list(set(hdr_to_add)))
    with open(args.pmhc_summaries) as pfo:
        for line_idx, line in enumerate(pfo.readlines()):
            if line_idx == 0:
                orig_hdr = line.strip().split('\t')
                line = line.rstrip().split('\t')
                hdr_map = {col:idx for idx, col in enumerate(line)}
            else:
                line = line.rstrip().split('\t')
                ident = line[hdr_map['identity']]
                for metric in hdr_to_add:
                    if metric in pmhcs_meta[ident].keys():
                        line.append(pmhcs_meta[ident][metric])
                    else:
                        line.append('NA')
                annotated_lines.append(line)

    orig_hdr.extend(hdr_to_add)

    with open(args.output, 'w') as ofo:
        ofo.write('\t'.join(orig_hdr))
        ofo.write('\n')
        for annot_line in annotated_lines:
            print(annot_line)
            ofo.write('\t'.join(annot_line))
            ofo.write('\n')

def add_agreto_metadata(args):
    """
    """
    mt_to_wt_map = {}
    with open(args.mt_wt_map, 'r') as mtwtfo:
        for line in mtwtfo:
            if re.search('^>', line):
                mtl = line.rstrip().lstrip('>')
                wtl = next(mtwtfo).strip()
                mt_to_wt_map[mtl] = wtl


    wt_metadata = {}
    agreto_hdr_map = {}
    with open(args.agreto, 'r') as afo:
        for line_idx, line in enumerate(afo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                agreto_hdr_map = {col:idx for idx, col in enumerate(line)}
            if line[agreto_hdr_map['peptide']] not in wt_metadata.keys():
                wt_metadata[line[agreto_hdr_map['peptide']]] = {}
            wt_metadata[line[agreto_hdr_map['peptide']]][line[agreto_hdr_map['allele']]] = {}
            for col in agreto_hdr_map:
                wt_metadata[line[agreto_hdr_map['peptide']]][line[agreto_hdr_map['allele']]][col] = line[agreto_hdr_map[col]]

    pprint(wt_metadata.items())

    out_lines = []
    out_header = []
    pmhcs_hdr_map = {}
    netmhcpan_col_name = ''
    mhcflurry_col_name = ''
    with open(args.pmhcs, 'r') as pfo:
        for line_idx, line in enumerate(pfo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                out_header = line
                print(out_header)
                pmhcs_hdr_map = {col:idx for idx, col in enumerate(line)}
                netmhcpan_col_name = [x for x in line if re.search('netmhcpan.*aff_nm', x)]
                mhcflurry_col_name = [x for x in line if re.search('mhcflurry.*aff', x)]
                if netmhcpan_col_name:
                    out_header.extend(['netmhcpan_agretopicity'])
                    netmhcpan_col_name = netmhcpan_col_name[0]
                if mhcflurry_col_name:
                    out_header.extend(['mhcflurry_agretopicity'])
                    mhcflurry_col_name = mhcflurry_col_name[0]

            elif line[pmhcs_hdr_map['peptide']] in mt_to_wt_map.keys() and mt_to_wt_map[line[pmhcs_hdr_map['peptide']]] in wt_metadata.keys():
                print("matched wt: {}".format(mt_to_wt_map[line[pmhcs_hdr_map['peptide']]]))
                print(line[pmhcs_hdr_map['peptide']])
                print(netmhcpan_col_name)
                print(mhcflurry_col_name)

                if netmhcpan_col_name:
                    if line[pmhcs_hdr_map[netmhcpan_col_name]] != 'NA' and wt_metadata[mt_to_wt_map[line[pmhcs_hdr_map['peptide']]]][line[pmhcs_hdr_map['allele']]][netmhcpan_col_name] != 'NA':
                        netmhcpan_mt = float(line[pmhcs_hdr_map[netmhcpan_col_name]])
                        netmhcpan_wt = float(wt_metadata[mt_to_wt_map[line[pmhcs_hdr_map['peptide']]]][line[pmhcs_hdr_map['allele']]][netmhcpan_col_name])
                        netmhcpan_agreto = str(netmhcpan_mt/netmhcpan_wt)
                        line.append(netmhcpan_agreto)
                    else:
                        line.append("NA")

                if mhcflurry_col_name:
                    if line[pmhcs_hdr_map[mhcflurry_col_name]] != 'NA' and wt_metadata[mt_to_wt_map[line[pmhcs_hdr_map['peptide']]]][line[pmhcs_hdr_map['allele']]][mhcflurry_col_name] != 'NA':
                        mhcflurry_mt = float(line[pmhcs_hdr_map[mhcflurry_col_name]])
                        mhcflurry_wt = float(wt_metadata[mt_to_wt_map[line[pmhcs_hdr_map['peptide']]]][line[pmhcs_hdr_map['allele']]][mhcflurry_col_name])
                        mhcflurry_agreto = str(mhcflurry_mt/mhcflurry_wt)
                        line.append(mhcflurry_agreto)
                    else:
                        line.append("NA")

                out_lines.append('\t'.join(line))
            else:
                default_out = []
                if netmhcpan_col_name:
                    default_out.append('NA')
                if mhcflurry_col_name:
                    default_out.append('NA')
                line.append("{}".format('\t'.join(default_out)))
                out_lines.append('\t'.join(line))


    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(out_header)))
        for out_line in out_lines:
            ofo.write(f"{out_line}\n")

def generic_annotation(args):
    """
    Adds annotations to to a metadata file based on a user-specific column for
    identification.
    """
    out_header = []

    columns_to_add = args.columns_to_add.split(',')
    tba_col_to_idx = {}
    tba = []
    with open(args.to_be_annotated, 'r') as tbao:
        for line_idx, line in enumerate(tbao.readlines()):
            line = line.rstrip().split('\t')
            if line_idx ==  0:
                for col_idx, col in enumerate(line):
                    tba_col_to_idx[col] = col_idx
                out_header = line
                out_header.extend(columns_to_add)
                print("new header: {}".format(out_header))
            else:
                tba.append(line)


    annot_data = {}
    with open(args.annotation_data, 'r') as ado:
        col_to_idx = {}
        for line_idx, line in enumerate(ado.readlines()):
            line = line.replace("b'", '').replace("'", '')
            line = line.rstrip().split('\t')
            if line_idx ==  0:
                for col_idx, col in enumerate(line):
                    col_to_idx[str(col)] = col_idx
                pprint(col_to_idx.items())
            else:
                annot_col = line[col_to_idx[str(args.annot_column)]]
                annot_data[annot_col] = {}
                for col, idx in col_to_idx.items():
                    annot_data[annot_col][str(col)] = line[idx]

    out_lines = []
    for line in tba:
        if args.orig_column in tba_col_to_idx.keys() and line[tba_col_to_idx[args.orig_column]] in annot_data.keys():
            mut = line[tba_col_to_idx[args.orig_column]]
            for column in columns_to_add:
                print("{} {}".format(column, annot_data[mut][column]))
                line.append(annot_data[mut][column])
            out_lines.append(line)
        else:
            for column in columns_to_add:
                line.append("NA")
            out_lines.append(line)


    with open(args.output, 'w') as ofo:
        ofo.write("{}".format('\t'.join(out_header)))
        for out_line in out_lines:
            ofo.write("\n{}".format('\t'.join(out_line)))

def generic_annotation_multikey(args):
    """
    Adds annotations to to a metadata file based on a user-specific columns for
    identification. This should be rolled up into generic_annotation(), but
    want to test it before integration.
    """
    out_header = []

    columns_to_add = args.columns_to_add.split(',')
    tba_col_to_idx = {}
    tba = []
    with open(args.to_be_annotated, 'r') as tbao:
        for line_idx, line in enumerate(tbao.readlines()):
            line = line.rstrip().split('\t')
            if line_idx ==  0:
                for col_idx, col in enumerate(line):
                    tba_col_to_idx[col] = col_idx
                out_header = line
                out_header.extend(columns_to_add)
                print("new header: {}".format(out_header))
            else:
                tba.append(line)


    annot_data = {}
    with open(args.annotation_data, 'r') as ado:
        col_to_idx = {}
        for line_idx, line in enumerate(ado.readlines()):
            line = line.replace("b'", '').replace("'", '')
            line = line.rstrip().split('\t')
            if line_idx ==  0:
                for col_idx, col in enumerate(line):
                    col_to_idx[str(col)] = col_idx
                pprint(col_to_idx.items())
            else:
                #annot_col = line[col_to_idx[str(args.annot_column)]]
                annot_cols = ','.join([line[col_to_idx[str(x)]] for x in args.annot_column.split(',')])
                annot_data[annot_cols] = {}
                for col, idx in col_to_idx.items():
                    annot_data[annot_cols][str(col)] = line[idx]
            print("annot data:")
#            pprint(annot_data)

    out_lines = []
    for line in tba:
        if all([x in tba_col_to_idx.keys() for x in args.orig_column.split(',')]) and ','.join([line[tba_col_to_idx[x]] for x in args.orig_column.split(',')]) in annot_data.keys():
#            mut = line[tba_col_to_idx[args.orig_column]]
            mut = ','.join([line[tba_col_to_idx[x]] for x in args.orig_column.split(',')])
            for column in columns_to_add:
                print("{} {}".format(column, annot_data[mut][column]))
                line.append(annot_data[mut][column])
            out_lines.append(line)
        else:
            for column in columns_to_add:
                line.append("NA")
            out_lines.append(line)


    with open(args.output, 'w') as ofo:
        ofo.write("{}".format('\t'.join(out_header)))
        for out_line in out_lines:
            ofo.write("\n{}".format('\t'.join(out_line)))


def get_fda_stats(args):
    """
    Get FDA-specific statistics for form NGS_Metadata_v0.3.xlsx.
    """
    fda_columns = ['submission_type',
                   'submission_number',
                   'sample_name',
                   'library_ID',
                   'sequencing_type',
                   'read_length',
                   'total_read_count',
                   'undetermined_bases_pct',
                   'Q30_bases_pct',
                   'filtered_read_count',
                   'mapped_reads_pct',
                   'unique_mapped_reads_pct',
                   'multi_mapping_reads_pct',
                   'mapped_read_duplication_pct',
                   'mapped_exonic_reads_pct',
                   'insert_size_in_bases',
                   'mean_coverage_depth',
                   'targeted_coverage_depth',
                   'GC_content_pct',
                   'reference_genome_version',
                   'reference_genome_file',
                   'genome_annotation_version',
                   'genome_annotation_file',
                   'filetype',
                   'filename1',
                   'filename2',
                   'filename3',
                   'md5_checksum_filename1',
                   'md5_checksum_filename2',
                   'md5_checksum_filename3',
                   'normalization_method']

    fda_stats = {'read_length':'',
                 'total_read_count':'',
                 'undetermined_bases_pct':'',
                 'Q30_bases_pct':'',
                 'filtered_read_count':'',
                 'mapped_reads_pct':'',
                 'unique_mapped_reads_pct':'',
                 'multi_mapping_reads_pct':'',
                 'mapped_read_duplication_pct':'',
                 'mapped_exonic_reads_pct':'',
                 'insert_size_in_bases':'',
                 'mean_coverage_depth':'',
                 'targeted_coverage_depth':'',
                 'GC_content_pct':''}

    fastp_stats = ''
    with open(args.fastq_stats) as fastp_fo:
        fastp_stats = json.load(fastp_fo)

    exome_stats = {}
    with open(args.exome_stats) as exome_fo:
        for line in exome_fo.readlines():
            if line.startswith('SN'):
                line = line.strip().split('\t')
                exome_stats[line[1].replace(':', '')] = line[2]
    
    bam_stats = {}
    with open(args.bam_stats) as bam_fo:
        for line in bam_fo.readlines():
           if line.startswith('SN'):
               line = line.strip().split('\t')
               bam_stats[line[1].replace(':', '')] = line[2]

    cov_stats = {'exon_read_count': 0,
                 'depth_sum': 0.0,
                 'depth_n': 0}

    with open(args.cov_stats) as cov_fo:
        for line in cov_fo.readlines():
            line = line.strip().split()
            cov_stats['exon_read_count'] += int(line[3])
            cov_stats['depth_sum'] += float(line[6])
            cov_stats['depth_n'] += 1

    fq_n_stats = {'tot': 0,
                  'n': 0}

    with open(args.n_stats) as n_fo:
        for line in n_fo.readlines():
            line = line.strip().partition(':')
            print(line)
            if line[0] == 'Total bases':
                fq_n_stats['tot'] = line[2].strip()
            if line[0] == 'Total Ns':
                fq_n_stats['n'] = line[2].strip()

    fda_stats['submission_type'] = 'IND'
    fda_stats['submission_number'] = args.submission_number

    print(fastp_stats['summary']['sequencing'].split(' '))
    fda_stats['read_length'] = "{} + {}".format(fastp_stats['summary']['sequencing'].split(' ')[-2].replace('(', ''), fastp_stats['summary']['sequencing'].split(' ')[5])
    fda_stats['total_read_count'] = fastp_stats['summary']['before_filtering']['total_reads']
    fda_stats['Q30_bases_pct'] = fastp_stats['summary']['before_filtering']['q30_rate'] * 100
    fda_stats['filtered_read_count'] = fastp_stats['summary']['after_filtering']['total_reads']
    fda_stats['GC_content_pct'] = fastp_stats['summary']['after_filtering']['gc_content'] * 100

    fda_stats['mapped_reads_pct'] = float(bam_stats['reads mapped'])/float(fda_stats['filtered_read_count']) * 100
    fda_stats['unique_mapped_reads_pct'] = (float(bam_stats['reads mapped']) - float(bam_stats['non-primary alignments']))/float(fda_stats['filtered_read_count']) * 100
    fda_stats['multi_mapping_reads_pct'] = float(bam_stats['non-primary alignments'])/float(fda_stats['filtered_read_count']) * 100
    fda_stats['insert_size_in_bases'] = bam_stats['insert size average']
    fda_stats['mapped_read_duplication_pct'] = float(bam_stats['reads duplicated'])/float(bam_stats['reads mapped']) * 100

    fda_stats['mapped_exonic_reads_pct'] = float(exome_stats['reads mapped'])/float(bam_stats['reads mapped']) * 100
    fda_stats['mean_coverage_depth'] = float(cov_stats['depth_sum'])/float(cov_stats['depth_n'])

    fda_stats['sample_name'] = args.patient_name
    fda_stats['library_ID'] = "{}-{}".format(args.patient_name, args.run)

    if re.search('ad-', args.run):
        fda_stats['sequencing_type'] = 'WES/WXS'
        fda_stats['targeted_coverage_depth'] = args.tumor_dna_depth
    elif re.search('nd-', args.run):
        fda_stats['sequencing_type'] = 'WES/WXS'
        fda_stats['targeted_coverage_depth'] = args.norm_dna_depth
    elif re.search('ar-', args.run):
        fda_stats['sequencing_type'] = 'RNA-seq'
        fda_stats['targeted_coverage_depth'] = "N/A"
        fda_stats['normalization_method'] = 'TPM'

    fda_stats['filetype'] = 'fastq'
    fda_stats['filename1'] = args.fastq1_path
    fda_stats['md5_checksum_filename1'] = args.fastq1_md5sum
    fda_stats['filename2'] = args.fastq2_path
    fda_stats['md5_checksum_filename2'] = args.fastq2_md5sum

    fda_stats['reference_genome_version'] = args.ref_version
    fda_stats['reference_genome_file'] = args.ref_path

    fda_stats['genome_annotation_version'] = args.gtf_version
    fda_stats['genome_annotation_file'] = args.gtf_path
    fda_stats['undetermined_bases_pct'] = float(fq_n_stats['n'])/float(fq_n_stats['tot']) * 100

    pprint(fda_stats)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(fda_columns)))
        for fda_column in fda_columns:
            if fda_column in fda_stats.keys():
                ofo.write("{}\t".format(fda_stats[fda_column]))
            else:
                ofo.write("N/A\t")
    

def get_n_count_and_freq(args):
    """
    Get the number of uncalled bases ("N"s) and the frequency of uncalled bases
    from two FASTQs.

    TODO: Check to see if there's a faster way of acheiving this as Python is
    frankly too slow.
    """
    tot = 0
    n_count = 0
    with gzip.open(args.fastq, 'rt') as ifo:
        print("Reading FASTQ 1...")
        for seq_record in SeqIO.parse(ifo, "fastq"):
            tot += len(str(seq_record.seq))
            n_count += str(seq_record.seq).count('N')
    if args.fastq2:
        print("Reading FASTQ 2...")
        with gzip.open(args.fastq2, 'rt') as ifo2:
            for seq_record in SeqIO.parse(ifo2, "fastq"):
                tot += len(str(seq_record.seq))
                n_count += str(seq_record.seq).count('N')

    with open(args.output, 'w') as ofo:
        ofo.write("Total bases: {}\n".format(tot))
        ofo.write("Total Ns: {}\n".format(n_count))
        ofo.write("Percent Ns: {}\n".format(float(n_count/tot + 0.0)))
