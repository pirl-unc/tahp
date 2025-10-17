import os
from glob import glob
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def make_antigens_barplot(args):
    """
    """
    patients = []
    counts = {}

    reports = glob(os.path.join(args.reports_dir, "*lens_report*.txt"))
    print(reports)

    antigen_sources = ['SNV', 'InDel', 'Virus', 'ERV', 'SpliceVariant',
                       'FusionEvent', 'Self-Antigen']
    for antigen_source in antigen_sources:
        counts[antigen_source] = []

    for patient_report in reports:
        print(patient_report)
        patient_id = patient_report.split('/')[-1].replace('.lens_report.txt', '')\
                                   .partition('-')[2]
        print(patient_id)
        patients.append(patient_id)

        tmp_report = pd.read_csv(patient_report, sep='\t')
        count_df = tmp_report['antigen_source'].value_counts()
        print(count_df)
        for antigen_source in counts.keys():
            if antigen_source in count_df.index.values:
                counts[antigen_source].append(count_df[antigen_source])
            else:
                counts[antigen_source].append(0)
    counts['CTA/Self-Antigen'] = counts['Self-Antigen']
    del counts['Self-Antigen']


    df = pd.DataFrame(counts, index=patients)
    print(df)

    a_series = (df != 0).any(axis=1)
    df = df.loc[a_series]

    df_sum = df.sum(axis=1)
    print(df_sum)
    df_sum.sort_values(ascending=True, inplace=True)
    df = df.reindex(df_sum.index)

    print(df)

    ax = df.plot.bar(stacked=True, colormap='viridis', log=True)
    ax.set_xlabel("Patients", fontweight='bold')
    ax.set_ylabel("Number of Predicted Peptides\n(<500 nM binding affinity)", fontweight='bold')
    ax.set_title("Predicted Peptides among Tumor Antigen Classes\nin {}".format(args.dataset), \
                 fontweight='bold')
    ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    ax.legend(loc=2, prop={'size': 5.5})

    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(args.output, format='pdf')

def make_antigens_heatmap(args):
    """
    """
    patients = []
    counts = {}

    reports = glob(os.path.join(args.reports_dir, "*lens_report*.txt"))
    print(reports)

    antigen_sources = ['SNV', 'InDel', 'Virus', 'ERV', 'SpliceVariant',
                       'FusionEvent', 'Self-Antigen']
    for antigen_source in antigen_sources:
        counts[antigen_source] = []

    for patient_report in reports:
        print(patient_report)
        patient_id = patient_report.split('/')[-1].replace('.lens_report.txt', '')\
                                   .partition('-')[2]
        print(patient_id)
        patients.append(patient_id)

        tmp_report = pd.read_csv(patient_report, sep='\t')
        count_df = tmp_report['antigen_source'].value_counts()
        print(count_df)
        for antigen_source in counts.keys():
            if antigen_source in count_df.index.values:
                counts[antigen_source].append(count_df[antigen_source])
            else:
                counts[antigen_source].append(0)
    counts['CTA/Self-Antigen'] = counts['Self-Antigen']
    del counts['Self-Antigen']


    df = pd.DataFrame(counts, index=patients)
    print(df)

    a_series = (df != 0).any(axis=1)
    df = df.loc[a_series]

    df_sum = df.sum(axis=1)
    print(df_sum)
    df_sum.sort_values(ascending=True, inplace=True)
    df = df.reindex(df_sum.index)

    print(df)

    plt.figure(figsize=(10, 5))

    hm = sns.heatmap(df.transpose(), norm=LogNorm(), cmap='viridis', xticklabels=False, \
                     linewidths=0.01, linecolor='black', square=True, \
                     cbar_kws={"shrink": 0.4, "orientation": "horizontal"}, \
                     annot_kws={"fontsize":5})

    hm.set_facecolor("lightgrey")
    hm.set_xlabel("Patients", fontweight='bold')
    hm.set_ylabel("Tumor\nAntigen\nSources", fontweight='bold')
    hm.set_title("Predicted Peptides among Tumor\nAntigen Sources in {}\n(<500 nM binding affinity)"\
                 .format(args.dataset), fontweight='bold')
    hm.set_yticklabels(hm.get_ymajorticklabels(), fontsize=4)
    hm.tick_params(labelsize=5)

    fig = hm.get_figure()
    fig.tight_layout()
    fig.savefig(args.output, format='pdf', bbox_inches='tight')


def make_antigens_circos_plot(args):
    """
    """

    tx_meta = {}

    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            if line.startswith('#'):
                continue
            if re.search('transcript_id ', line):
                tx_id_idx = ''
                meta_entries = line.split('\t')[8].split('; ')
                for meta_entry_idx, meta_entry in enumerate(meta_entries):
                    if re.search('transcript_id ', meta_entry):
                        tx_id_idx = meta_entry_idx
                        break
                tx = line.split('\t')[8].split('; ')[tx_id_idx]\
                                 .replace('"', '')\
                                 .replace('transcript_id ', '')

                tx_meta[tx] = {}
                tx_meta[tx]['chrom'] = line.split('\t')[0]
                tx_meta[tx]['start'] = line.split('\t')[3]
                tx_meta[tx]['stop'] = line.split('\t')[4]

    gn_meta = {}

    with open(args.gtf) as gtf_fo:
        for line in gtf_fo.readlines():
            if line.startswith('#'):
                continue
            if re.search('gene_name ', line):
                gn_name_idx = ''
                meta_entries = line.split('\t')[8].split('; ')
                for meta_entry_idx, meta_entry in enumerate(meta_entries):
                    if re.search('gene_name ', meta_entry):
                        gn_name_idx = meta_entry_idx
                        break
                gn = line.split('\t')[8].split('; ')[gn_name_idx]\
                                 .replace('"', '')\
                                 .replace('gene_name ', '')

                gn_meta[gn] = {}
                gn_meta[gn]['chrom'] = line.split('\t')[0]
                gn_meta[gn]['start'] = line.split('\t')[3]
                gn_meta[gn]['stop'] = line.split('\t')[4]

    cti = {}


    fusions = {}
    snvs = {}
    indels = {}
    splice = {}
    ervs = {}
    ctas = {}

    with open(args.report, 'r') as ifo:
        for line_idx, line in enumerate(ifo.readlines()):
            line = line.split('\t')
            if line_idx == 0:
                for i, j in enumerate(line):
                    cti[j] = i
            else:
                if line[cti['antigen_source']] == 'FUSION':
                    print("Adding fusion")
                    fusions[line[cti['identity']]] = {}
                    fusions[line[cti['identity']]]['fusion_id'] = line[cti['fusion_id']]
                    left_gene = line[cti['fusion_id']].partition('--')[0]
                    right_gene = line[cti['fusion_id']].partition('--')[2]
                    fusions[line[cti['identity']]]['left_breakpoint'] = line[cti['fusion_left_breakpoint']]
                    fusions[line[cti['identity']]]['left_start'] = gn_meta[left_gene]['start']
                    fusions[line[cti['identity']]]['left_chrom'] = gn_meta[left_gene]['chrom'].replace('chr', 'hs')
                    fusions[line[cti['identity']]]['right_breakpoint'] = line[cti['fusion_right_breakpoint']]
                    fusions[line[cti['identity']]]['right_end'] = gn_meta[right_gene]['stop']
                    fusions[line[cti['identity']]]['right_chrom'] = gn_meta[right_gene]['chrom'].replace('chr', 'hs')
                    fusions[line[cti['identity']]]['fusion_type'] = line[cti['fusion_type']]
                elif line[cti['antigen_source']] == 'SNV':
                    print("Adding SNV")
                    snvs[line[cti['identity']]] = {}
                    snvs[line[cti['identity']]]['chrom'] = line[cti['variant_pos']].split(':')[0].replace('chr', 'hs')
                    snvs[line[cti['identity']]]['variant_pos'] = line[cti['variant_pos']].split(':')[1]
                elif line[cti['antigen_source']] == 'INDEL':
                    print("Adding INDEL")
                    indels[line[cti['identity']]] = {}
                    indels[line[cti['identity']]]['chrom'] = line[cti['variant_pos']].split(':')[0].replace('chr', 'hs')
                    indels[line[cti['identity']]]['variant_pos'] = line[cti['variant_pos']].split(':')[1]
                    indels[line[cti['identity']]]['indel_type'] = line[cti['indel_type']]
                elif line[cti['antigen_source']] == 'SPLICE':
                    print("Adding splice")
                    splice[line[cti['identity']]] = {}
                    splice[line[cti['identity']]]['chrom'] = line[cti['chromosome']].replace('chr', 'hs')
                    splice[line[cti['identity']]]['splice_start'] = line[cti['tumor_splice']].replace('(','').split(',')[0]
                    splice[line[cti['identity']]]['splice_stop'] = line[cti['tumor_splice']].replace('(','').split(',')[1]
                elif line[cti['antigen_source']] == 'ERV':
                    print("Adding ERV")
                    ervs[line[cti['identity']]] = {}
                    ervs[line[cti['identity']]]['chrom'] = line[cti['erv_orf_id']].split('.')[1].replace('chr', 'hs')
                    ervs[line[cti['identity']]]['erv_start'] = line[cti['erv_orf_id']].split('.')[2]
                    ervs[line[cti['identity']]]['erv_stop'] = line[cti['erv_orf_id']].split('.')[3]
                elif line[cti['antigen_source']] == 'CTA/SELF':
                    print("Adding CTA")
                    ctas[line[cti['identity']]] = {}
                    tx = line[cti['transcript_id']]
                    ctas[line[cti['identity']]]['chrom'] = tx_meta[tx]['chrom']
                    ctas[line[cti['identity']]]['start'] = tx_meta[tx]['start']
                    ctas[line[cti['identity']]]['stop'] = tx_meta[tx]['stop']

        with open(os.path.join(args.output_dir, "circos.fusions"), 'w') as ofo:
            for fusion, fmetadata in fusions.items():
                ofo.write("{}\t{}\t{}\t{}\t{}\t{}\n"\
                          .format(fmetadata['left_chrom'], \
                                  fmetadata['left_start'], \
                                  fmetadata['left_breakpoint'], \
                                  fmetadata['right_chrom'], \
                                  fmetadata['right_breakpoint'], \
                                  fmetadata['right_end']))

        with open(os.path.join(args.output_dir, "circos.snvs"), 'w') as ofo:
            for snv in snvs.keys():
                ofo.write("{}\t{}\t{}\t0\n"\
                          .format(snvs[snv]['chrom'], snvs[snv]['variant_pos'], \
                                  snvs[snv]['variant_pos']))

        with open(os.path.join(args.output_dir, "circos.indels"), 'w') as ofo:
            for indel in indels.keys():
                ofo.write("{}\t{}\t{}\t0\n"\
                          .format(indels[indel]['chrom'], indels[indel]['variant_pos'], \
                                  indels[indel]['variant_pos']))


        with open(os.path.join(args.output_dir, "circos.splice"), 'w') as ofo:
            for ssplice in splice.keys():
                ofo.write("{}\t{}\t{}\t0\n"\
                          .format(splice[ssplice]['chrom'], \
                                  splice[ssplice]['splice_start'], \
                                  splice[ssplice]['splice_stop']))

        with open(os.path.join(args.output_dir, "circos.ervs"), 'w') as ofo:
            for erv in ervs.keys():
                ofo.write("{}\t{}\t{}\t0\n"\
                          .format(ervs[erv]['chrom'], \
                                  ervs[erv]['erv_start'], \
                                  ervs[erv]['erv_stop']))

        with open(os.path.join(args.output_dir, "circos.ctas"), 'w') as ofo:
            for cta in ctas.keys():
                ofo.write("{}\t{}\t{}\t0\n"\
                          .format(ctas[cta]['chrom'], \
                                  ctas[cta]['start'], \
                                  ctas[cta]['stop']))
