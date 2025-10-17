"""
Misc. functions for LENSTools.
"""

import os
import re


def consolidate_multiqc_stats(args):
    """
    """

    sample_level_stats = {}

    stats_files = ["multiqc_fastqc.txt",
                   "multiqc_picard_RnaSeqMetrics.txt",
                   "multiqc_samtools_stats.txt"]

    for stats_file in stats_files:
        with open(os.path.join(os.getcwd(), 'multiqc_data', stats_file)) as sfo:
            metric_to_idx = {}
            for line_idx, line in enumerate(sfo.readlines()):
                line = line.rstrip().split('\t')
                if line_idx == 0:
                    for metric_idx, metric in enumerate(line[1:]):
                        metric_to_idx[metric] = metric_idx + 1
                else:
                    sample = line[0]
                    suffix = ''
                    if re.search("_[12]\.", sample) or re.search("_[12]$", sample):
                        print("Reformating sample name")
                        sample = '_'.join(line[0].split('_')[:-1])
                        suffix = ' (file {})'.format(line[metric_to_idx['Filename']].split('_')[-1])
                    print("New Sample: {}".format(sample))
                    print("Suffix Sample: {}".format(suffix))

                    if sample not in sample_level_stats.keys():
                        sample_level_stats[sample] = {}
                    for metric, metric_idx in metric_to_idx.items():
                        sample_level_stats[sample]["{}{}".format(metric, suffix)] = line[metric_idx]

    with open(args.output, 'w') as ofo:
        header = []
        for sample in sample_level_stats.keys():
            for metric in sample_level_stats[sample]:
                if metric not in header:
                    header.append(metric)

        header = sorted(header)

        header.insert(0, 'sample')

        ofo.write("{}\n".format('\t'.join(header)))

        for sample in sample_level_stats.keys():
            print(sample_level_stats[sample])
            sample_line = []
            sample_line.append(sample)
            for metric in header[1:]:
                if metric in sample_level_stats[sample].keys():
                    sample_line.append(sample_level_stats[sample][metric])
                else:
                    sample_line.append('NA')
            ofo.write("{}\n".format('\t'.join(sample_line)))


def complement(i):
    comp = {'a': 't',
            't': 'a',
            'g': 'c',
            'c': 'g'}
    return comp[i.lower()]


def add_rna_norms(args):
    """
    """
    hmap = {}
    manifest_entries = []
    with open(args.manifest, encoding='utf-8-sig') as mani_fo:
        for line_idx, line in enumerate(mani_fo.readlines()):
            if line_idx == 0:
                line = line.strip('\n').strip('\r').split('\t')
                print(line)
                for col_idx, col in enumerate(line):
                    hmap[col] = col_idx
            else:
                line = line.strip('\n')
                line = line.strip('\r') # Dealing with manifests made in Windows.
                line = line.split('\t')
                print(line)
                trunc_line = [line[hmap['Patient_Name']], line[hmap['Run_Name']],
                              line[hmap['Dataset']], line[hmap['File_Prefix']],
                              line[hmap['Sequencing_Method']], line[hmap['Normal']]]
                manifest_entries.append('\t'.join(trunc_line))
        new_line = [args.pat_name, "nr-{}".format(args.pat_name),
                    args.dataset, args.prefix, 'RNA-Seq', 'TRUE']

        manifest_entries.append('\t'.join(new_line))


    with open(args.output, 'w') as ofo:
        ofo.write('{}\n'.format('\t'.join(['Patient_Name', 'Run_Name',
                                           'Dataset', 'File_Prefix',
                                           'Sequencing_Method', 'Normal'])))
        for entry in manifest_entries:
            ofo.write("{}\n".format(entry))
    return manifest_entries


def iupac_conversion(base):
    emission = []
    conversion = {'R': ['A', 'G'],
                  'Y': ['C', 'T'],
                  'S': ['G', 'C'],
                  'W': ['A', 'T'],
                  'K': ['G', 'T'],
                  'M': ['A', 'C']}
    if base in conversion.keys():
        emission = conversion[base]
    return emission


def make_lens_bed(report, gtf, samp_id, outdir):
    """
    """
    col_map = {}

    output_lines = []

    with open(report) as report_fo:
        for line_idx, line in enumerate(report_fo.readlines()):
            line = line.rstrip().split('\t')
            print(line)
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    col_map[col] = col_idx
            else:
                if line[col_map['antigen_source']] in ['SNV', 'INDEL']:
                    iid = line[col_map['identity']]
                    chrom = line[col_map['variant_coords']].split(':')[0]
                    pos = int(line[col_map['variant_coords']].split(':')[1])
                    snapshot_lbound = pos - 25
                    snapshot_ubound = pos + 25
                    output_lines.append("{}\t{}\t{}\t{}-{}-{}.png\n".format(chrom, snapshot_lbound, \
                                                                        snapshot_ubound, samp_id, \
                                                                        iid, \
                                                                        line[col_map['antigen_source']]))
                elif line[col_map['antigen_source']] in ['FUSION']:
                    iid = line[col_map['identity']]
                    lchr = line[col_map['fusion_left_breakpoint']].split(':')[0]
                    lpos = int(line[col_map['fusion_left_breakpoint']].split(':')[1])
                    llbound = lpos - 25
                    lubound = lpos + 25
                    rchr = line[col_map['fusion_right_breakpoint']].split(':')[0]
                    rpos = int(line[col_map['fusion_right_breakpoint']].split(':')[1])
                    rlbound = rpos - 25
                    rubound = rpos + 25
                    output_lines.append("{}\t{}\t{}\t{}-{}-{}_L.png\n".format(lchr, llbound, \
                                                                              lubound, samp_id, \
                                                                              iid, 'FUSION'))
                    output_lines.append("{}\t{}\t{}\t{}-{}-{}_R.png\n".format(rchr, rlbound, \
                                                                              rubound, samp_id, \
                                                                              iid, 'FUSION'))
                elif line[col_map['antigen_source']] in ['ERV']:
                    iid = line[col_map['identity']]
                    chrom = line[col_map['erv_orf_id']].split('.')[1]
                    lbound = int(line[col_map['erv_orf_id']].split('.')[2]) - 25
                    ubound = int(line[col_map['erv_orf_id']].split('.')[3]) + 25
                    output_lines.append("{}\t{}\t{}\t{}-{}-{}.png\n".format(chrom, lbound, ubound, \
                                                                            samp_id, iid, 'ERV'))
                elif line[col_map['antigen_source']] in ['SPLICE']:
                    iid = line[col_map['identity']]
                    chrom = line[col_map['splice_coords']].split(':')[0]
                    lbound = int(line[col_map['splice_coords']].split(':')[1].split('-')[0].replace("(",'')) - 25
                    ubound = int(line[col_map['splice_coords']].split(':')[1].split('(')[0].split('-')[1]) + 25
                    if lbound > ubound:
                        lbound, ubound = ubound, lbound
                    output_lines.append("{}\t{}\t{}\t{}-{}-{}.png\n".format(chrom, lbound, ubound, \
                                                                            samp_id, iid, 'SPLICE'))

#    if outdir:
#        with open("{}/{}.lens.igv_snapshot.bed".format(outdir, samp_id), 'w') as output_fo:
#            for output_line in sorted(list(set(output_lines))):
#                output_fo.write(output_line)
#    else:
    with open("{}.lens.tumor_antigens.bed".format(samp_id), 'w') as output_fo:
        for output_line in sorted(list(set(output_lines))):
            output_fo.write(output_line)
