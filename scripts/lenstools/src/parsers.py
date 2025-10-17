import re
import vcf

from Bio import SeqIO


def parse_ag_foreignness(ag_foreignness_file):
    peptide_to_foreign = {}
    with open(ag_foreignness_file) as ffo:
        foreign_col_idx = ''
        pep_col_idx = ''
        for line_idx, line in enumerate(ffo.readlines()):
            line = line.strip()
            if line_idx == 0:
                print(line)
                foreign_col_idx = line.split(',').index('foreignness_score')
                pep_col_idx = line.split(',').index('nmer')
            else:
                line = line.split(',')
                peptide_to_foreign[line[pep_col_idx]] = line[foreign_col_idx]
    return peptide_to_foreign


def parse_ag_dissimilarity(ag_dissimilarity_file):
    peptide_to_dissim = {}
    with open(ag_dissimilarity_file) as ffo:
        dissim_col_idx = ''
        pep_col_idx = ''
        for line_idx, line in enumerate(ffo.readlines()):
            line = line.strip()
            if line_idx == 0:
                dissim_col_idx = line.split(',').index('dissimilarity')
                pep_col_idx = line.split(',').index('nmer')
            else:
                line = line.split(',')
                peptide_to_dissim[line[pep_col_idx]] = line[dissim_col_idx]
    return peptide_to_dissim

def parse_netmhcstabpan(netmhcstabpan_file):
    peptide_allele_to_stability = {}
    with open(netmhcstabpan_file) as ffo:
        for line in ffo.readlines():
            if line.startswith(" ") and (re.search("HLA-", line) or re.search ("H-2-", line)):
                print(line)
                line = line.split()
                print(line)
                peptide_allele_to_stability["{}-{}".format(line[1], line[2])] = line[5]
    return peptide_allele_to_stability

def sanitize_netmhc_record(record):
    """
    """
    if 'SB' in record:
        record.remove('<=')
        record.remove('SB')
    if 'WB' in record:
        record.remove('<=')
        record.remove('WB')
    return record


def load_vcf(input_vcf):
    """
    Loads VCF file.

    Args:
        input_vcf: Input VCF file.

    Returns:
        vcf_reader (vcf Reader object): loaded VCF.
    """
    vcf_reader = ''
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')
    return vcf_reader

def load_tx_to_gene(gtf):
    """
    """
    tx_to_gene = {}
    with open(gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript' and line[1] != 'geve':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0]
                tx_to_gene[tx_id] = gene_name
                tx_to_gene[tx_id.split('.')[0]] = gene_name
    return tx_to_gene

def load_tcga_tx_summ(tumor_type, summ_file):
    """
    """
    col_idx_map = {}
    tumor_summ = {}

    # Loading TCGA expression summarization file.
    with open(summ_file) as summ_fo:
        for line_idx, line in enumerate(summ_fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    col_idx_map[col] = col_idx
            else:
                if line[col_idx_map['tumor_type']] == tumor_type:
                    tumor_summ[line[col_idx_map['identifier']]] = {}
                    tumor_summ[line[col_idx_map['identifier']]]['mean'] = line[col_idx_map['mean']]
                    tumor_summ[line[col_idx_map['identifier']]]['median'] = line[col_idx_map['median']]
                    tumor_summ[line[col_idx_map['identifier']]]['max'] = line[col_idx_map['max']]
                    tumor_summ[line[col_idx_map['identifier']]]['iqr'] = line[col_idx_map['iqr']]

    return tumor_summ
def load_tcga_dicts():
    """
    Convert TCGA tumor type abbreviations to complete name.
    """
    abbrev_to_tumor_type = {
    "LAML": "Acute Myeloid Leukemia",
    "ACC":  "Adrenocortical carcinoma",
    "BLCA": "Bladder Urothelial Carcinoma",
    "LGG":  "Brain Lower Grade Glioma",
    "BRCA": "Breast invasive carcinoma",
    "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
    "CHOL": "Cholangiocarcinoma",
    "LCML": "Chronic Myelogenous Leukemia",
    "COAD": "Colon adenocarcinoma",
    "CNTL": "Controls",
    "ESCA": "Esophageal carcinoma",
    "FPPP": "FFPE Pilot Phase II",
    "GBM":  "Glioblastoma multiforme",
    "HNSC": "Head and Neck squamous cell carcinoma",
    "KICH": "Kidney Chromophobe",
    "KIRC": "Kidney renal clear cell carcinoma",
    "KIRP": "Kidney renal papillary cell carcinoma",
    "LIHC": "Liver hepatocellular carcinoma",
    "LUAD": "Lung adenocarcinoma",
    "LUSC": "Lung squamous cell carcinoma",
    "DLBC": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
    "MESO": "Mesothelioma",
    "MISC": "Miscellaneous",
    "OV":   "Ovarian serous cystadenocarcinoma",
    "PAAD": "Pancreatic adenocarcinoma",
    "PCPG": "Pheochromocytoma and Paraganglioma",
    "PRAD": "Prostate adenocarcinoma",
    "READ": "Rectum adenocarcinoma",
    "SARC": "Sarcoma",
    "SKCM": "Skin Cutaneous Melanoma",
    "STAD": "Stomach adenocarcinoma",
    "TGCT": "Testicular Germ Cell Tumors",
    "THYM": "Thymoma",
    "THCA": "Thyroid carcinoma",
    "UCS":  "Uterine Carcinosarcoma",
    "UCEC": "Uterine Corpus Endometrial Carcinoma",
    "UVM":  "Uveal Melanoma"}

    return abbrev_to_tumor_type


def parse_fasta(fasta):
    """"
    """
    seqs = {}
    for seq_record in SeqIO.parse(fasta, 'fasta'):
        seqs[seq_record.description.split('|')[1]] = seq_record.seq

    return seqs

def parse_fasta_full_header(fasta):
    """"
    """
    seqs = {}
    for seq_record in SeqIO.parse(fasta, 'fasta'):
        seqs[seq_record.description] = str(seq_record.seq)

    return seqs
