
import math
import re
import vcf


def make_pyclone_vi_inputs(args):
    """
    This requires an isec vcf, a proper vcf (with depth info), and sequenza results info.
    """
    # This dictionary will be populated with normal and tumor depth information
    # from the proper VCFs. It'll initially have keys populated by the isec VCF.

    vars = {}

    cellularity = 0
    with open(args.sequenza_solutions) as sso:
        for line_idx, line in enumerate(sso.readlines()):
            if line_idx == 1:
                line = line.split('\t')
                cellularity = line[0]

    with open(args.candidate_vcf) as cvo:
        for line in cvo.readlines():
            if line.startswith('#'):
                pass
            else:
                line = line.rstrip().split('\t')
                line = [i for i in line if i != '.']
                vars['{}'.format('_'.join(line[:2]))] = {}
                vars['{}'.format('_'.join(line[:2]))]['ref_depth'] = 0
                vars['{}'.format('_'.join(line[:2]))]['alt_depth'] = 0
                vars['{}'.format('_'.join(line[:2]))]['major_cn'] = 0
                vars['{}'.format('_'.join(line[:2]))]['minor_cn'] = 0
    print("Candidate variant count: {}".format(len(vars)))
    #Note: not all variants may not be in depth VCF. This can be fixed later.

    mutect_vcf_reader = vcf.Reader(open(args.mutect_vcf), compressed=False)
    depths = {}
    for record in mutect_vcf_reader:
        call = record.genotype([i.sample for i in record.samples if re.search('-ad[-_]', i.sample)][0])
        var_key = "{}_{}".format(record.CHROM, record.POS)
        depths[var_key] = {}
        depths[var_key]["ref_depth"] = call.data.AD[0]
        depths[var_key]["alt_depth"] = call.data.AD[1]

    copy_no = {}
    with open(args.sequenza_segments) as sfo:
        for line_idx, line in enumerate(sfo.readlines()):
            if line_idx == 0:
                continue
            line = line.split('\t')
            chrom, start, stop, maj_count, min_count = (line[0], line[1], line[2],
                                                        line[10], line[11])
            if maj_count != 'NA' and min_count != 'NA':
                print("{}".format(','.join([chrom, start, stop, maj_count, min_count])))
                if chrom not in copy_no.keys():
                    copy_no[chrom] = {}
                if "{}-{}".format(start, stop) not in copy_no[chrom].keys():
                    copy_no[chrom]["{}-{}".format(start, stop)] = {}
                copy_no[chrom]["{}-{}".format(start, stop)]["maj_count"] = maj_count
                copy_no[chrom]["{}-{}".format(start, stop)]["min_count"] = min_count

    easily_parsable = []
    for var in vars.keys():
        print(var)
        chrom, pos = var.split('_')
        easily_parsable.append([chrom, pos])

    pyclone_inp = []

    for chrom in copy_no.keys():
        chr_vars = sorted([x for x in easily_parsable if x[0] == chrom], key=lambda x: x[1])
        for segment in sorted(copy_no[chrom].keys(), key=lambda segment: int(segment.split('-')[0])):
            print("SEGMENT: {}".format(segment))
            start, stop = segment.split('-')
            segment_vars = sorted([x for x in chr_vars if (int(x[1]) > int(start) and
                                                           int(x[1]) < int(stop))],
                                  key=lambda x: int(x[1]))
            print(segment_vars)
            for segment_var in segment_vars:
                var = "{}_{}".format(segment_var[0], segment_var[1])
                if var in depths.keys():
                    ref_depth = str(depths[var]['ref_depth'])
                    alt_depth = str(depths[var]['alt_depth'])
                    ref_cn = str(copy_no[chrom][segment]["maj_count"])
                    alt_cn = str(copy_no[chrom][segment]["min_count"])
                    var = var.replace('_', ':')
                    pyclone_inp.append("{}\n".format('\t'.join([var,
                                                                args.samp_id,
                                                                ref_depth,
                                                                alt_depth,
                                                                ref_cn,
                                                                alt_cn,
                                                                '2',
                                                                '0.001',
                                                                cellularity])))

    for var in easily_parsable:
        var = "{}:{}".format(var[0], var[1])
        if var not in [i.split('\t')[0] for i in pyclone_inp]:
            pyc_line_to_append = '\t'.join([var,
                                           args.samp_id,
                                           '0',
                                           '0',
                                           '2',
                                           '0',
                                           '2',
                                           '0.001',
                                           cellularity])
            pyclone_inp.append("{}\n".format(pyc_line_to_append))

    with open(args.output, 'a') as ofo:
        for i in pyclone_inp:
            ofo.write(i)


def calculate_ccf(args):
    """
    """
    ccfs = {}
    gene_cns = {}

    with open(args.scna) as scnao:
        for line_idx,line in enumerate(scnao.readlines()):
            if line_idx > 0:
                line = line.split()
                genes = line[3].split(',')
                for gene in genes:
                    gene_cns[gene.partition('.')[0]] = line[6]

    print(gene_cns.keys())


    vcf_reader = vcf.Reader(open(args.vcf))
    for record in vcf_reader:
        nt = 'NA' # Assume diploid in case missing from CNVKit
        gene_id = record.INFO['ANN'][0].split('|')[4]
        call = record.genotype([i.sample for i in record.samples if re.search('-ad-', i.sample)][0])
        var_key = "{}:{}".format(record.CHROM, record.POS)
        print("variant: {}".format(var_key))
        freq = 'NA'
        if hasattr(call.data, 'AF'):
            freq = call.data.AF
        elif hasattr(call.data, 'VAF'):
            freq = call.data.VAF
        elif hasattr(call.data, 'FREQ'):
            freq = "0.{}".format(call.data.FREQ.replace('.', '').replace('%', ''))
        elif hasattr(call.data, 'TIR'):
            freq = call.data.TIR[0] / (call.data.TIR[0] + call.data.TAR[0])
        else:
            if record.alleles[0] == 'A':
                ref_count = call.data.AU[0]
            elif record.alleles[0] == 'C':
                ref_count = call.data.CU[0]
            elif record.alleles[0] == 'G':
                ref_count = call.data.GU[0]
            elif record.alleles[0] == 'T':
                ref_count = call.data.TU[0]

            if record.alleles[1] == 'A':
                alt_count = call.data.AU[0]
            elif record.alleles[1] == 'C':
                alt_count = call.data.CU[0]
            elif record.alleles[1] == 'G':
                alt_count = call.data.GU[0]
            elif record.alleles[1] == 'T':
                alt_count = call.data.TU[0]
            if alt_count > 0 and ref_count > 0:
                freq = alt_count / (alt_count + ref_count)
            else:
                freq = 0.0
        if gene_id in gene_cns.keys():
            nt = float(gene_cns[gene_id])
        elif gene_id.partition('.')[0] in gene_cns.keys():
            nt = float(gene_cns[gene_id.partition('.')[0]])
        print("nt: {}".format(nt))
        print("freq: {}".format(freq))
        if nt != 'NA' and float(freq) > 0.0:
            var_mult = _multiplicity_equation(float(freq), float(nt), float(args.purity), 1.0)
            var_ccf = _ccf_equation(float(freq), float(nt), float(args.purity), math.ceil(var_mult))
            print("freq:{}\tnt{}\tvar_mult{}\tvar_ccf{}".format(freq, nt, round(var_mult, 3), round(var_ccf, 3)))
            ccfs[var_key] = {}
            ccfs[var_key]['freq'] = freq
            ccfs[var_key]['nt'] = nt
            ccfs[var_key]['var_mult'] = var_mult
            ccfs[var_key]['var_ccf'] = var_ccf

    with open(args.output, 'w') as ofo:
        ofo.write("{}\t{}\t{}\t{}\t{}\n".format('variant', 'vaf', 'totcopynum', 'multiplicity', 'ccf'))
        for i in ccfs.keys():
            print("i: {}".format(i))
            ofo.write("{}\t{}\t{}\t{}\t{}\n".format(i,
                                                    ccfs[i]['freq'],
                                                    ccfs[i]['nt'], 
                                                    ccfs[i]['var_mult'],
                                                    ccfs[i]['var_ccf']))



def _multiplicity_equation(vaf, nt, purity, mult):
    """
    """
    print("vaf: {}\npurity: {}\nmult: {}\n".format(vaf, purity, mult))
    if mult == 0.0:
        mult = 0.00001
    leftside = float(vaf/(purity*mult))
    #print("leftside: {}".format(leftside))
    rightside = float((purity * nt) + (2 * (1-purity)))
    #print("rightside: {}".format(rightside))
    return leftside * rightside


def _ccf_equation(vaf, nt, purity, mult):
    """
    """
    return _multiplicity_equation(float(vaf), float(nt), float(purity), float(mult))
