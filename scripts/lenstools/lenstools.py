#!/usr/bin/env python3
"""
Collection of functions for running LENS workflow.
"""

import argparse
import os

from src import metadata
from src import rna_utils
from src import pep_utils
from src import misc
from src import plots
from src import onco
from src import var_utils

os.environ['OPENBLAS_NUM_THREADS'] = '1'


def get_args():
    """
    Get user-provided arguments for LENSTools.
    """
    parser = argparse.ArgumentParser(prog="lenstools",
                                     description="lenstools")
    subparsers = parser.add_subparsers(dest='command')


    # Subparser for generating peptides by coding InDels
    parser_mk_indel_peps_cntxt = subparsers.add_parser('make-indel-peptides-context',
                                                       help="Make InDel peptides FASTA.")
    parser_mk_indel_peps_cntxt.add_argument('-vts', '--var-tx-seqs',
                                            help="Variant-specific transcript sequences.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-st', '--somatic-txs',
                                            help="Expressed transcripts with somatic InDels.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-sv', '--somatic-vcf',
                                            help="Somatic variants of interest.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-g', '--gtf',
                                            help="GTF with gene annotations.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-l', '--length',
                                            help="Emitted peptide length. (default: 11)",
                                            default=11)
    parser_mk_indel_peps_cntxt.add_argument('--three-prime-keys',
                                            help="")
    parser_mk_indel_peps_cntxt.add_argument('--three-prime-seqs',
                                            help="")
    parser_mk_indel_peps_cntxt.add_argument('-pr', '--pep-ref',
                                          help="protein reference",
                                          required=True)
    parser_mk_indel_peps_cntxt.add_argument('-n', '--context-nt-length',
                                            help="Emitted sequence context length (default: 39)",
                                            default=39)
    parser_mk_indel_peps_cntxt.add_argument('-o', '--output',
                                            help="Output mutant peptides FASTA.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('--nt-output',
                                            help="Output mutant nucleotides FASTA.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-do', '--debug-output',
                                            help="Alignment debug output.",
                                            default = '')


    # Subparser for generating peptides by missense SNVs
    parser_mk_snv_peps_cntxt = subparsers.add_parser('make-snv-peptides-context',
                                                     help="Make SNV peptides FASTA.")
    parser_mk_snv_peps_cntxt.add_argument('-vts', '--var-tx-seqs',
                                          help="Directory with var-specific transcript sequences.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-st', '--somatic-txs',
                                          help="Expressed transcripts harboring somatic SNVs.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-sv', '--somatic-vcf',
                                          help="Somatic variants of interest.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-g', '--gtf',
                                          help="GTF with gene annotations.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-pr', '--pep-ref',
                                          help="protein reference",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-l', '--length',
                                          help="Emitted peptide length. (default: 11)",
                                          default=11)
    parser_mk_snv_peps_cntxt.add_argument('-n', '--context-nt-length',
                                          help="Emitted sequence context length (default: 39)",
                                          default=39)
    parser_mk_snv_peps_cntxt.add_argument('-o', '--mt-output',
                                          help="Output mutant peptides FASTA.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('--nt-output',
                                          help="Output mutant nucleotides FASTA.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-w', '--wt-output',
                                          help="Output wildtype peptides FASTA.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-do', '--debug-output',
                                          help="Alignment debug output.",
                                          default = '')


    # Subparser for making ERV peptides
    parser_mk_erv_peps = subparsers.add_parser('make-erv-peptides',
                                               help="make ERV peptides FASTA.")
    parser_mk_erv_peps.add_argument('-e', '--expressed-ervs',
                                    help="File containg list of Expressed ERVs",
                                    required=True)
    parser_mk_erv_peps.add_argument('-r', '--patient-ervs-fasta',
                                    help="Patient ERV FASTA (with homozygous germline vars).",
                                    required=True)
    parser_mk_erv_peps.add_argument('-s', '--species',
                                    help="Species (e.g. mm). Default: hs",
                                    default='hs')
    parser_mk_erv_peps.add_argument('-g', '--geve-reference',
                                    help="gEVE reference file",
                                    required=True)
    parser_mk_erv_peps.add_argument('-o', '--output',
                                    help="Output ERV peptides FASTA.",
                                    required=True)
    parser_mk_erv_peps.add_argument('-n', '--nt-output',
                                    help="Output ERV neoleotide FASTA.",
                                    required=True)


    # Subparser for making viral peptides
    parser_mk_viral_peps = subparsers.add_parser('make-viral-peptides',
                                                  help="Make viral peptides FASTA.")
    parser_mk_viral_peps.add_argument('-f', '--patient-viral-fasta',
                                      help="Patient viral fasta (with homozygous germline vars).",
                                      required=True)
    parser_mk_viral_peps.add_argument('-o', '--output',
                                      help="Output viral peptides FASTA.",
                                      required=True)


    # Subparser for making self-antigen peptides
    parser_mk_self_peps = subparsers.add_parser('make-self-antigen-peptides',
                                                help="Make CTA/self-antigen peptides FASTA.")
    parser_mk_self_peps.add_argument('-s', '--selfs-seqs-fasta',
                                     help="CTA/Self-antigen FASTA (with homozygous germline vars).",
                                     required=True)
    parser_mk_self_peps.add_argument('-g', '--gtf',
                                     help="GTF with gene annotations.",
                                     required=True)
    parser_mk_self_peps.add_argument('-e', '--expressed-selfs',
                                     help="File of expressed CTA/self-antigen genes.",
                                     required=True)
    parser_mk_self_peps.add_argument('-o', '--output',
                                     help="Output CTA/self-antigens peptide FASTA.",
                                     required=True)
    parser_mk_self_peps.add_argument('-n', '--nt-output',
                                     help="Output CTA/self-antigens nucleotide FASTA.",
                                     required=True)


    # Subparser for making fusion peptides
    parser_mk_fusion_peps_cntxt = subparsers.add_parser('make-fusion-peptides-context',
                                                        help="Make fusion peptides FASTA.")
    parser_mk_fusion_peps_cntxt.add_argument('-f', '--fusions',
                                             help="Predicted fusions (STARFusion format).",
                                             required=True)
    parser_mk_fusion_peps_cntxt.add_argument('-g', '--gtf',
                                             help="GTF with gene annotations.",
                                             required=True)
    parser_mk_fusion_peps_cntxt.add_argument('-e', '--exons-fasta',
                                             help="Fusion transcripts exons FASTA")
    parser_mk_fusion_peps_cntxt.add_argument('-t', '--fusion-txs',
                                             help="Expressed fusion transcripts",
                                             required=True)
    parser_mk_fusion_peps_cntxt.add_argument('-o', '--output',
                                             help='Output fusion peptides FASTA.',
                                             required=True)
    parser_mk_fusion_peps_cntxt.add_argument('-n', '--nt-output',
                                             help='Output fusion nucleotide FASTA.',
                                             required=True)


    # Subparser for adding SNV metadata
    parser_add_snv_metadata = subparsers.add_parser('add-snv-metadata',
                                                    help="Add SNV peptides metadata.")
    parser_add_snv_metadata.add_argument('-m', '--mutant-peptides',
                                         help="Mutant SNV peptides FASTA.",
                                         required=True)
    parser_add_snv_metadata.add_argument('-q', '--quants',
                                         help="Tumor ranscript abundance file (Salmon format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-c', '--cancer-cell-fraction',
                                         help="Cancer cell fraction file (PyClone-VI format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-g', '--gtf',
                                         help="GTF with gene annotations",
                                         required=True)
    parser_add_snv_metadata.add_argument('-b', '--binding-affinities',
                                         help="Binding affinities file (NetMHCpan-4.1b format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-s', '--binding-stabilities',
                                         help="Binding stabilities file (netMHCstabpan format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-f', '--ag-foreignness',
                                         help="Foreignness file (antigen.garnish format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-d', '--ag-dissimilarity',
                                         help="Dissimilarity file (antigen.garnish format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-o', '--output',
                                         help="Output SNV report.",
                                         required=True)


    # Subparser for adding InDel metadata
    parser_add_indel_metadata = subparsers.add_parser('add-indel-metadata',
                                                      help="Add InDel peptides metadata.")
    parser_add_indel_metadata.add_argument('-m', '--mutant-peptides',
                                           help="Mutant InDel peptides FASTA.",
                                           required=True)
    parser_add_indel_metadata.add_argument('-q', '--quants',
                                           help="Tumor ranscript abundance file (Salmon format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-c', '--cancer-cell-fraction',
                                           help="Cancer cell fraction file (PyClone-VI format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-g', '--gtf',
                                           help="GTF file with gene annotations.",
                                           required=True)
    parser_add_indel_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinities file (netMHCpan format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-s', '--binding-stabilities',
                                           help="Binding stabilities file (netMHCstabpan format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-f', '--ag-foreignness',
                                           help="Foreignness file (antigen.garnish format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-d', '--ag-dissimilarity',
                                           help="Dissimilarity file (antigen.garnish format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-o', '--output',
                                           help="Output InDel report.",
                                           required=True)


    # Subparser for adding ERV metadata
    parser_add_erv_metadata = subparsers.add_parser('add-erv-metadata',
                                                     help="Add ERV peptides metadata.")
    parser_add_erv_metadata.add_argument('-p', '--peptides',
                                         help="Annotated peptide file.",
                                         required=True)
    parser_add_erv_metadata.add_argument('-b', '--binding-affinities',
                                         help="Binding affinities file (NetMHCpan-4.1b format).",
                                         required=True)
    parser_add_erv_metadata.add_argument('-q', '--quants',
                                         help="Tumor transcript abundance file (Salmon format).",
                                         required=True)
    parser_add_erv_metadata.add_argument('-v', '--patient-vcf',
                                         help="Patient ERV VCF",
                                         required=True)
    parser_add_erv_metadata.add_argument('-n', '--nt',
                                         help="Patient-specific ERV FASTA",
                                         required=True)
    parser_add_erv_metadata.add_argument('-d', '--geve-data',
                                         help="External gEVE data.",
                                         required=True)
    parser_add_erv_metadata.add_argument('-s', '--binding-stabilities',
                                         help="Binding stabilities file (netMHCstabpan format).",
                                         required=True)
    parser_add_erv_metadata.add_argument('-f', '--ag-foreignness',
                                         help="Foreignness file (antigen.garnish format).",
                                         required=True)
    parser_add_erv_metadata.add_argument('-i', '--ag-dissimilarity',
                                         help="Dissimilarity file (antigen.garnish format).",
                                         required=True)
    parser_add_erv_metadata.add_argument('-t', '--trim-chr-prefix',
                                         help="Trim the chr prefix from gEVE references.",
                                         action='store_true')
    parser_add_erv_metadata.add_argument('-o', '--output',
                                         help="Output ERV report",
                                         required=True)


    # Subparser for getting viral metadata
    parser_add_viral_metadata = subparsers.add_parser('add-viral-metadata',
                                                      help="Add viral peptides metadata.")
    parser_add_viral_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinities data (netMHCpan-4.1b format).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-r', '--viral-cds-ref',
                                           help="Viral CDS FASTA (used for VirDetect).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-p', '--viral-pep-ref',
                                           help="Viral peptide FASTA (used for VirDetect).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-s', '--binding-stabilities',
                                           help="Binding stabilities file (netMHCstabpan format).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-f', '--ag-foreignness',
                                           help="Foreignness file (antigen.garnish format).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-d', '--ag-dissimilarity',
                                           help="Dissimilarity file (antigen.garnish format).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-o', '--output',
                                           help="Output viral report.",
                                           required=True)


    # Subparser for getting self-antigen metadata
    parser_add_self_metadata = subparsers.add_parser('add-self-antigen-metadata',
                                                     help="Add self-antigen peptides metadata.")
    parser_add_self_metadata.add_argument('-q', '--quants',
                                          help="Tumor transcript abundance file (Salmon format).",
                                          required=True)
    parser_add_self_metadata.add_argument('-b', '--binding-affinities',
                                          help="Binding affinitities data (netMHCpan-4.1b format)",
                                          required=True)
    parser_add_self_metadata.add_argument('-f', '--fasta',
                                          help="CTA/Self-antigen Peptide FASTA",
                                          required=True)
    parser_add_self_metadata.add_argument('-g', '--gtf',
                                          help="GTF with gene annotation.s",
                                          required=True)
    parser_add_self_metadata.add_argument('-l', '--gene-list',
                                          help="File with CTA/self-antigen gene list.",
                                          required=True)
    parser_add_self_metadata.add_argument('-s', '--binding-stabilities',
                                         help="Binding stabilities file (netMHCstabpan format).",
                                         required=True)
    parser_add_self_metadata.add_argument('-r', '--ag-foreignness',
                                         help="Foreignness file (antigen.garnish format).",
                                         required=True)
    parser_add_self_metadata.add_argument('-d', '--ag-dissimilarity',
                                         help="Dissimilarity file (antigen.garnish format).",
                                         required=True)
    parser_add_self_metadata.add_argument('-o', '--output',
                                          help="Output CTA/self-antigen report.",
                                          required=True)


    # Subparser for adding fusion metadata
    parser_add_fusion_metadata = subparsers.add_parser('add-fusion-metadata',
                                                       help="Add metadata for fusion peptides.")
    parser_add_fusion_metadata.add_argument('-b', '--binding-affinities',
                                            help="Binding affinities data (netMHCpan-4.1b format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-f', '--fusions',
                                            help="Predicted fusions (STARFusion format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-a', '--fasta',
                                            help="FASTA with fusion peptides.",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-s', '--binding-stabilities',
                                            help="Binding stabilities file (netMHCstabpan format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-r', '--ag-foreignness',
                                            help="Foreignness file (antigen.garnish format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-d', '--ag-dissimilarity',
                                            help="Dissimilarity file (antigen.garnish format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)


    # Subparser for filtering variants for expression
    parser_expressed_variants = subparsers.add_parser('filter-expressed-variants',
                                                 help="Filter variants for expression.")
    parser_expressed_variants.add_argument('-q', '--quants',
                                           help="Transcript abundance file (Salmon format)).",
                                           required=True)
    parser_expressed_variants.add_argument('-m', '--metric',
                                           help="Expression metric (default: TPM).",
                                           default='TPM')
    parser_expressed_variants.add_argument('-r', '--exclude-zeros',
                                           help="Exclude zeros for percentile. (default: True)",
                                           action='store_true')
    parser_expressed_variants.add_argument('-p', '--percentile',
                                           help="Expression filtering percentile (default: 90).",
                                           default=90)
    parser_expressed_variants.add_argument('-t', '--abundance-threshold',
                                           help="Expression filtering threshold (default: 0).",
                                           default=0)
    parser_expressed_variants.add_argument('-v', '--vcf',
                                           help="Annotated (snpEff) VCF.",
                                           required=True)
    parser_expressed_variants.add_argument('-o', '--output',
                                           help="Output file.",
                                           required=True)
    parser_expressed_variants.add_argument('-s', '--somatic-txs',
                                           help="Somatic transcripst output file.",
                                           required=True)


    # Subparser for filtering isolated variants (e.g. no proximal germline or somatic variants).
    parser_isolated_variants = subparsers.add_parser('filter-isolated-variants',
                                                     help="Filter isolated variants.")
    parser_isolated_variants.add_argument('-g', '--germline-vcf',
                                          help="Germline VCF.",
                                          required=True)
    parser_isolated_variants.add_argument('-s', '--somatic-vcf',
                                          help="Somatic VCF.",
                                          required=True)
    parser_isolated_variants.add_argument('-p', '--proximity',
                                          help="Clearance (in bp) around variants (default: 30).",
                                          default=30)
    parser_isolated_variants.add_argument('-a', '--allow-silent-and-homozygous',
                                          help="Allow silent/homozygous germline variants.",
                                          action="store_true")
    parser_isolated_variants.add_argument('-o', '--output',
                                          help="Output file.",
                                          required=True)


    # Subparser for filtering expressed ERVs through TPM.
    parser_expressed_ervs = subparsers.add_parser('filter-expressed-ervs',
                                                   help="Filter expressed ERVs.")
    parser_expressed_ervs.add_argument('-q', '--quants',
                                       help="Transcript abundance file (Salmon format).",
                                       required=True)
    parser_expressed_ervs.add_argument('-m', '--metric',
                                       help="Expression column. (default: TPM)",
                                       default='TPM')
    parser_expressed_ervs.add_argument('-r', '--exclude-zeros',
                                       help="Exclude zeros for expression percentile.",
                                       action='store_true')
    parser_expressed_ervs.add_argument('-z', '--abundance-threshold',
                                       help="Expression threshold for filtering. (default: 10)",
                                       default=10)
    parser_expressed_ervs.add_argument('-t', '--trim-chr-prefix',
                                       help="Trim the chr prefix from gEVE references.",
                                       action='store_true')
    parser_expressed_ervs.add_argument('-o', '--output',
                                       help="Output file.",
                                       required=True)


    # Subparser for filtering expressed ERVs through RNA coverage.
    parser_erv_rna_cov = subparsers.add_parser('check-erv-rna-coverage',
                                                help="Filter ERVs with sufficient coverage.")
    parser_erv_rna_cov.add_argument('-c', '--coverage-file',
                                    help="Coverage file created by samtools cov.",
                                    required=True)
    parser_erv_rna_cov.add_argument('-e', '--expressed-ervs',
                                    help="File containing list of expressed ERVs.",
                                    default="TPM")
    parser_erv_rna_cov.add_argument('-m', '--mean-depth',
                                    help="Required depth for expressed ERV (default: 0).",
                                    default=0)
    parser_erv_rna_cov.add_argument('-p', '--coverage',
                                    help="Required coverage for expressed ERV (default: 0).",
                                    default=0)
    parser_erv_rna_cov.add_argument('-o', '--output',
                                    help="Output file.",
                                    required=True)


    # Subparser for filtering expressed viruses through RNA coverage.
    parser_virus_rna_cov = subparsers.add_parser('check-virus-rna-coverage',
                                                      help=("Filter viruses with "
                                                            "sufficient coverage."))
    parser_virus_rna_cov.add_argument('-c', '--coverage-file',
                                      help="Coverage file created by samtools cov.",
                                      required=True)
    parser_virus_rna_cov.add_argument('-e', '--expressed-viruses',
                                      help="File containing list of expressed viruses.",
                                      required=True)
    parser_virus_rna_cov.add_argument('-m', '--mean-depth',
                                      help="Required depth for expressed virus (default: 0).",
                                      default=0)
    parser_virus_rna_cov.add_argument('-p', '--coverage',
                                      help="Required coverage for expressed virus (default: 0).",
                                      default=0)
    parser_virus_rna_cov.add_argument('-o', '--output',
                                      help="Output file.",
                                      required=True)


    # Subparser for creating expressed ERV bed file.
    parser_expressed_ervs_bed = subparsers.add_parser('get-expressed-ervs-bed',
                                                      help="Create BED file for expressed ERVs.")
    parser_expressed_ervs_bed.add_argument('-e', '--expressed-ervs',
                                           help="File containing list of expressed ERVs.",
                                           required=True)
    parser_expressed_ervs_bed.add_argument('-r', '--geve-reference',
                                           help="gEVE reference file",
                                           required=True)
    parser_expressed_ervs_bed.add_argument('-o', '--output',
                                           help="Output file.",
                                           required=True)


    # Subparser for creating expressed self-antigen/CTA BED file
    parser_expressed_selfs_bed = subparsers.add_parser('get-expressed-selfs-bed',
                                                       help=("Make BED for expressed"
                                                             " self-antigen/CTAs."))
    parser_expressed_selfs_bed.add_argument('-e', '--expressed-selfs',
                                            help="Expressed CTA/self-antigen gene list file.",
                                            required=True)
    parser_expressed_selfs_bed.add_argument('-g', '--gff',
                                            help="GFF File.",
                                            required=True)
    parser_expressed_selfs_bed.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)


    # Subparser for creating expressed self-antigen/CTA BED file
    parser_expressed_txs_bed = subparsers.add_parser('get-expressed-transcripts-bed',
                                                       help="Create BED for expressed transripts.")
    parser_expressed_txs_bed.add_argument('-t', '--expressed-transcripts',
                                          help="Expressed transcripts file.",
                                          required=True)
    parser_expressed_txs_bed.add_argument('-g', '--gff',
                                          help="GFF with gene annotations.",
                                          required=True)
    parser_expressed_txs_bed.add_argument('-o', '--output',
                                          help="Output file.",
                                          required=True)


    # Subparser for filtering expressed viruses
    parser_expressed_viral_bed = subparsers.add_parser('get-expressed-viral-bed',
                                                       help="Filter expressed hERVs.")
    parser_expressed_viral_bed.add_argument('-e', '--expressed-viruses',
                                            help="Expressed viruses file.",
                                            required=True)
    parser_expressed_viral_bed.add_argument('-r', '--viral-cds-ref',
                                            help="Viral coding sequence (CDS) FASTA.",
                                            required=True)
    parser_expressed_viral_bed.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)


    # Subparser for filtering virdetect outputs for expressed viruses
    parser_filter_virdetect_by_counts = subparsers.add_parser('filter-expressed-viruses',
                                                   help="Generate list of expressed viruses")
    parser_filter_virdetect_by_counts.add_argument('-q', '--viral-quants',
                                                   help="Viral counts file from VirDetect.",
                                                   required=True)
    parser_filter_virdetect_by_counts.add_argument('-r', '--viral-ref',
                                                   help="Viral reference FASTA (for VirDetect).",
                                                   required=True)
    parser_filter_virdetect_by_counts.add_argument('-m', '--min-threshold',
                                                   help="Minimal filtering count (default: 1).")
    parser_filter_virdetect_by_counts.add_argument('-o', '--output',
                                                   help="Output file.",
                                                   required=True)


    # Subparser for filtering virdetect outputs for expressed viruses
    parser_filter_viral_cds = subparsers.add_parser('filter-expressed-viral-cds',
                                                   help="Make expressed viral coding seqs list.")
    parser_filter_viral_cds.add_argument('-q', '--viral-cds-quants',
                                         help="Viral counts file from VirDetect.",
                                         required=True)
    parser_filter_viral_cds.add_argument('-r', '--viral-cds-ref',
                                         help="Viral reference FASTA (for VirDetect).",
                                         required=True)
    parser_filter_viral_cds.add_argument('-e', '--expressed-viruses',
                                         help="List of expressed viruses.",
                                         required=True)
    parser_filter_viral_cds.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)


    # Subparser for filtering self-antigens for expression
    parser_expressed_self_genes = subparsers.add_parser('filter-expressed-self-genes',
                                                        help=("Filter expressed "
                                                              "CTA/self-antigen genes."))
    parser_expressed_self_genes.add_argument('-q', '--quants',
                                             help="Transcript abundance file (Salmon format).",
                                             required=True)
    parser_expressed_self_genes.add_argument('-m', '--metric',
                                             help="Expression metric. (default: TPM)",
                                             default='TPM')
    parser_expressed_self_genes.add_argument('-r', '--exclude-zeros',
                                             help="Exclude zeros for percentile.",
                                             action='store_true')
    parser_expressed_self_genes.add_argument('-p', '--percentile',
                                             help="Expression filtering percentile (default: 90).",
                                             default=90)
    parser_expressed_self_genes.add_argument('-t', '--abundance-threshold',
                                             help="Expression filtering abundance (default: 0)",
                                             default=0)
    parser_expressed_self_genes.add_argument('-g', '--gene-list',
                                             help="File containing CTA/Self-antigen genes.",
                                             required=True)
    parser_expressed_self_genes.add_argument('-f', '--gtf',
                                             help="GTF file.",
                                             required=True)
    parser_expressed_self_genes.add_argument('-o', '--output',
                                             help="Output file.",
                                             required=True)


    # Subparser for calculalting agretopicity (mut BA/wt BA)
    parser_calc_agreto = subparsers.add_parser('calculate-agretopicity',
                                               help="Calcuate agreotopicity (mut BA/wt BA).")
    parser_calc_agreto.add_argument('-w', '--wt-fasta',
                                    help="Wildtype peptide FASTA (LENSTools format).",
                                    required=True)
    parser_calc_agreto.add_argument('-m', '--mt-fasta',
                                    help="Mutant peptide FASTA. (LENSTools format)",
                                    required=True)
    parser_calc_agreto.add_argument('-o', '--output',
                                    help="Output file.",
                                   required=True)

    # Subparser for adding agretopicity data to annotated pMHCs
    parser_add_agreto = subparsers.add_parser('add-agreto-metadata',
                                               help="Add agretopicity metadata.")
    parser_add_agreto.add_argument('-p', '--pmhcs',
                                    help="Annotated pMHCs.",
                                    required=True)
    parser_add_agreto.add_argument('-m', '--mt-wt-map',
                                    help="FASTA with mutant and wildtype peptides.",
                                    required=True)
    parser_add_agreto.add_argument('-a', '--agreto',
                                    help="Agretopicity metadata file",
                                    required=True)
    parser_add_agreto.add_argument('-o', '--output',
                                    help="Updated pMHCs file with agretopicity",
                                    required=True)



    # Subparser for creating PyClone-VI inputs
    parser_make_pvi_inputs = subparsers.add_parser('make-pyclone-vi-inputs',
                                                   help="Make input file required for PyClone-VI")
    parser_make_pvi_inputs.add_argument('-c', '--candidate-vcf',
                                        help="VCF containing candidate variants.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-m', '--mutect-vcf',
                                        help="Mutect VCF for variant depth information.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-s', '--sequenza-segments',
                                        help="Sequenza segments file.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('--sequenza-solutions',
                                        help="Sequenza solutions file.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('--samp-id',
                                        help="Sample identifier.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-o', '--output',
                                        help="Output file.",
                                        required=True)

    parser_make_pvi_inputs = subparsers.add_parser('calculate-ccf',
                                                   help="Calculate cancer cell fraction (CCF).")
    parser_make_pvi_inputs.add_argument('-v', '--vcf',
                                        help="Annotated somatic VCF.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-p', '--purity',
                                        help="Sample purity value.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-s', '--scna',
                                        help="CNVKit output.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-o', '--output',
                                        help="Output file.",
                                        required=True)


    # Subparser for consolidating multiqc statistics files
    parser_consol_mqc_stats = subparsers.add_parser('consolidate-multiqc-stats',
                                                   help="Consolidate multiqc stats files.")
    parser_consol_mqc_stats.add_argument('-d', '--multiqc-data',
                                         help="Directory containing multiqc stats files.",
                                         required=True)
    parser_consol_mqc_stats.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)


    # Subparser for adding RNA normals
    parser_add_rna_norms = subparsers.add_parser('add-rna-normals',
                                                 help="Add RNA normal runs to patient data.")
    parser_add_rna_norms.add_argument('-m', '--manifest',
                                      help="Manifest file with original patient data.",
                                      required=True)
    parser_add_rna_norms.add_argument('-p', '--prefix',
                                      help="RNA normal sample Run Name (e.g. nr-CTRL).",
                                      required=True)
    parser_add_rna_norms.add_argument('-d', '--dataset',
                                      help="RNA normal sample dataset.",
                                      required=True)
    parser_add_rna_norms.add_argument('-n', '--pat_name',
                                      help="RNA normal sample Patient Name (e.g. CTRL).",
                                      required=True)
    parser_add_rna_norms.add_argument('-o', '--output',
                                      help="Output manifest file.",
                                      required=True)


    # Subparser for adding splice metadata
    parser_add_splice_meta = subparsers.add_parser('add-splice-metadata',
                                                 help="Add splice metadata.")
    parser_add_splice_meta.add_argument('-s', '--splice-summary',
                                        help="Splice summary file from NeoSplice.",
                                        required=True)
    parser_add_splice_meta.add_argument('-f', '--ag-foreignness',
                                        help="Foreignness file (antigen.garnish format).",
                                        required=True)
    parser_add_splice_meta.add_argument('-d', '--ag-dissimilarity',
                                        help="Dissimilarity file (antigen.garnish format).",
                                        required=True)
    parser_add_splice_meta.add_argument('-b', '--binding-stabilities',
                                        help="Binding stabilities file (netMHCstabpan format).",
                                        required=True)
    parser_add_splice_meta.add_argument('-o', '--output',
                                        help="Output file.",
                                        required=True)


    # subparser for making lens report
    parser_mk_lens_report = subparsers.add_parser('make-lens-report')
    parser_mk_lens_report.add_argument('--metadata-dir', '-d',
                                       help="path to metadata directory.",
                                       required=True)
    parser_mk_lens_report.add_argument('--output', '-o',
                                       help="output file.",
                                       required=True)


    # Subparser for making pan-patient antigen source barplot
    parser_mk_antigen_plot = subparsers.add_parser('make-antigens-barplot')
    parser_mk_antigen_plot.add_argument('--reports-dir', '-d',
                                        help="Path containing patient reports.",
                                        required=True)
    parser_mk_antigen_plot.add_argument('--output', '-o',
                                        help="Output file.",
                                        required=True)
    parser_mk_antigen_plot.add_argument('--dataset',
                                        help="Dataset (for plot title).",
                                        default='')

    # Subparser for making pan-patient antigen source heatmap
    parser_mk_antigen_hm = subparsers.add_parser('make-antigens-heatmap')
    parser_mk_antigen_hm.add_argument('--reports-dir', '-d',
                                        help="Path containing patient reports.",
                                        required=True)
    parser_mk_antigen_hm.add_argument('--output', '-o',
                                        help="Output file.",
                                        required=True)
    parser_mk_antigen_hm.add_argument('--dataset',
                                        help="Dataset (for plot title).",
                                        default='')

    # Subparser for making patient antigen source circos
    parser_mk_antigen_circos = subparsers.add_parser('make-antigens-circos-plot')
    parser_mk_antigen_circos.add_argument('--report', '-r',
                                        help="LENS report",
                                        required=True)
    parser_mk_antigen_circos.add_argument('--output-dir', '-o',
                                        help="Output file.",
                                        required=True)
    parser_mk_antigen_circos.add_argument('--gtf',
                                        help="GTF.",
                                        required=True)


    # Subparser for adding TCGA data to SNV or InDel LENS report
    parser_add_tcga_data = subparsers.add_parser('add-tcga-data')
    parser_add_tcga_data.add_argument('--report', '-r',
                                      help="SNV or InDel LENS Report",
                                      required=True)
    parser_add_tcga_data.add_argument('--tumor-type', '-t',
                                      help="Tumor type (either full name or abbreviation).",
                                      required=True)
    parser_add_tcga_data.add_argument('--tcga-transcript-summary', '-s',
                                      help="TCGA transcript summary file.",
                                      required=True)
    parser_add_tcga_data.add_argument('--output', '-o',
                                      help="Output file.",
                                      required=True)


    # Subparser for determining read count support for Viral, ERV, and self peptides.
    parser_get_pep_rd_cnt = subparsers.add_parser('get-peptide-read-count')
    parser_get_pep_rd_cnt.add_argument('--pmhcs', '-p',
                                       help="File with pMHCs of interest.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--bam', '-b',
                                       help="RNA-Sequencing BAM.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--antigen-source', '-s',
                                       help="Antigen source of interest.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--cds-fasta', '-c',
                                       help="FASTA with peptide nucleotide coding sequence.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--output', '-o',
                                       help="Output file.",
                                       required=True)


    # Subparser for determining read count support for SNV-derived peptides.
    parser_get_snv_pep_rd_cnt = subparsers.add_parser('get-snv-peptide-read-count')
    parser_get_snv_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                           help="NetMHCpan file with peptides of interest.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--bam', '-b',
                                           help="RNA-Sequencing BAM.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--nt-fasta', '-c',
                                           help="FASTA file peptide nucleotide coding sequence.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--gtf', '-g',
                                           help="GTF with gene annotations.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--output', '-o',
                                           help="Output file.",
                                           required=True)


    # Subparser for determining read count support for InDel peptides.
    parser_get_indel_pep_rd_cnt = subparsers.add_parser('get-indel-peptide-read-count')
    parser_get_indel_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                             help="NetMHCpan file with peptides of interest.",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--bam', '-b',
                                             help="RNA-Sequencing BAM",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--nt-fasta', '-c',
                                             help="Peptide nucleotide coding sequence FASTA",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--gtf', '-g',
                                             help="GTF with gene annotations",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--output', '-o',
                                             help="Output file.",
                                             required=True)


    # Subparser for determining read count support for fusion peptides.
    parser_get_fus_pep_rd_cnt = subparsers.add_parser('get-fusion-peptide-read-count')
    parser_get_fus_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                           help="NetMHCpan file with peptides of interest.",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--fusions', '-f',
                                           help="Fusions (STARFusion format)",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--nt-fasta', '-t',
                                           help="Peptide nucleotide coding sequences FASTA",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--fusion-reads', '-r',
                                           help="Junction or discordant reads supporting fusions",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--output', '-o',
                                           help="Output file.",
                                           required=True)


    # Subparser for determining read count support for splice peptides.
    parser_get_splc_pep_rd_cnt = subparsers.add_parser('get-splice-peptide-read-count')
    parser_get_splc_pep_rd_cnt.add_argument('--neosplice-summary', '-n',
                                            help="NeoSplice summary file",
                                            required=True)
    parser_get_splc_pep_rd_cnt.add_argument('--bam', '-b',
                                            help="RNA-Sequencing BAM file",
                                            required=True)
    parser_get_splc_pep_rd_cnt.add_argument('--output', '-o',
                                            help="Output file.",
                                            required=True)


    # Subparser for filtering peptides that exist in wildtype sample
    parser_fltr_mut_peps = subparsers.add_parser('filter-mutant-peptides')
    parser_fltr_mut_peps.add_argument('--wt-netmhcpan', '-w',
                                      help="Wildtype NetMHCpan")
    parser_fltr_mut_peps.add_argument('--mt-netmhcpan', '-m',
                                      help="Mutant NetMHCpan",
                                      required=True)
    parser_fltr_mut_peps.add_argument('--wt-output', '-wo',
                                      help="Wildtype filtered output")
    parser_fltr_mut_peps.add_argument('--mt-output', '-mo',
                                      help="Mutant filtered output",
                                      required=True)
    parser_fltr_mut_peps.add_argument('--mt-fasta', '-mf',
                                      help="Mutant fasta (for InDels).")
    parser_fltr_mut_peps.add_argument('--max-peptide-length', '-l',
                                      help="Maximum peptide length (default: 11)",
                                      default=11)

    # Subparser for filtering expressed variants
    parser_expressed_txs = subparsers.add_parser('filter-expressed-transcripts',
                                                 help="Filter transcripts for expression.")
    parser_expressed_txs.add_argument('-q', '--quants',
                                      help="Transcript abundance file (Salmon format)).",
                                      required=True)
    parser_expressed_txs.add_argument('-m', '--metric',
                                      help="Expression metric (default: TPM).",
                                      default='TPM')
    parser_expressed_txs.add_argument('-r', '--exclude-zeros',
                                      help="Exclude zeros for percentile. (default: True)",
                                      action='store_true')
    parser_expressed_txs.add_argument('-p', '--percentile',
                                      help="Expression filtering percentile (default: 90).",
                                      default=90)
    parser_expressed_txs.add_argument('-t', '--abundance-threshold',
                                      help="Expression filtering threshold (default: 0).",
                                      default=0)
    parser_expressed_txs.add_argument('-o', '--output',
                                      help="Output file.",
                                      required=True)


    parser_agg_lens_reports = subparsers.add_parser('aggregate-lens-reports',
                                                    help="Combine multiple patient LENS reports.")
    parser_agg_lens_reports.add_argument('-d', '--reports-dir', required=True)
    parser_agg_lens_reports.add_argument('--no-pat-id',
                                         action='store_true')
    parser_agg_lens_reports.add_argument('-o', '--output', required=True)


    parser_mk_lens_bed = subparsers.add_parser('make-lens-bed',
                                                    help="Make BED file for LENS tumor antigens.")
    parser_mk_lens_bed.add_argument('-d', '--report', required=True)
    parser_mk_lens_bed.add_argument('-g', '--gtf')
    parser_mk_lens_bed.add_argument('-s', '--samp-id', required=True)
    parser_mk_lens_bed.add_argument('-o', '--output', default='')


    parser_agg_pmhc_outs = subparsers.add_parser('aggregate-pmhc-outputs',
                                                    help="Combine multiple tool pMHC outputs.")
    parser_agg_pmhc_outs.add_argument('--netmhcpan', default='')
    parser_agg_pmhc_outs.add_argument('--netctlpan', default='')
    parser_agg_pmhc_outs.add_argument('--netmhcstabpan', default='')
    parser_agg_pmhc_outs.add_argument('--mhcflurry', default='')
    parser_agg_pmhc_outs.add_argument('--deephlapan', default='')
    parser_agg_pmhc_outs.add_argument('--ag-dissim', default='')
    parser_agg_pmhc_outs.add_argument('--ag-foreign', default='')
    parser_agg_pmhc_outs.add_argument('--pepsickle', default='')
    parser_agg_pmhc_outs.add_argument('--hlapollo', default='')
    parser_agg_pmhc_outs.add_argument('-o', '--output', required=True)

    parser_annot_pmhcs = subparsers.add_parser('annotate-pmhcs',
                                               help="Annotate pMHC summary with relevant metadata.")
    parser_annot_pmhcs.add_argument('--pmhc-summaries', '-s', required=True)
    parser_annot_pmhcs.add_argument('--pmhc-fasta', '-f', required=True)
    parser_annot_pmhcs.add_argument('--output', '-o', required=True)


    parser_gener_annot = subparsers.add_parser('generic-annotation',
                                               help="Annotate pMHC summary with relevant metadata.")
    parser_gener_annot.add_argument('--to-be-annotated', '-t', required=True)
    parser_gener_annot.add_argument('--annotation-data', '-d', required=True)
    parser_gener_annot.add_argument('--orig-column', '-a', required=True)
    parser_gener_annot.add_argument('--annot_column', '-b', required=True)
    parser_gener_annot.add_argument('--columns-to-add', '-c', required=True)
    parser_gener_annot.add_argument('--output', '-o', required=True)

    parser_gener_annot_mk = subparsers.add_parser('generic-annotation-multikey',
                                               help="Annotate pMHC summary with relevant metadata.")
    parser_gener_annot_mk.add_argument('--to-be-annotated', '-t', required=True)
    parser_gener_annot_mk.add_argument('--annotation-data', '-d', required=True)
    parser_gener_annot_mk.add_argument('--orig-column', '-a', required=True)
    parser_gener_annot_mk.add_argument('--annot_column', '-b', required=True)
    parser_gener_annot_mk.add_argument('--columns-to-add', '-c', required=True)
    parser_gener_annot_mk.add_argument('--output', '-o', required=True)

    parser_priortze_peps = subparsers.add_parser('prioritize-peptides',
                                               help="Prioritize peptides.")
    parser_priortze_peps.add_argument('--input-report', '-i', required=True)
    parser_priortze_peps.add_argument('--output', '-o', required=True)

    parser_filtr_peps = subparsers.add_parser('filter-peptides',
                                               help="Filter peptides against peptidome.")
    parser_filtr_peps.add_argument('--input-report', '-i', required=True)
    parser_filtr_peps.add_argument('--pep-ref', '-p', required=True)
    parser_filtr_peps.add_argument('--ctas', '-c', required=True)
    parser_filtr_peps.add_argument('--output', '-o', required=True)

    parser_get_ns = subparsers.add_parser('get-n-count-and-freq',
                                           help="Get N count and frequency from FASTQs.")
    parser_get_ns.add_argument('--fastq', '-f', required=True)
    parser_get_ns.add_argument('--fastq2', '-f2')
    parser_get_ns.add_argument('--output', '-o', required=True)

    parser_get_fda_stats = subparsers.add_parser('get-fda-stats',
                                           help="Get FDA stats.")
    parser_get_fda_stats.add_argument('--fastq-stats', '-f', required=True)
    parser_get_fda_stats.add_argument('--bam-stats', '-b', required=True)
    parser_get_fda_stats.add_argument('--cov-stats', '-c', required=True)
    parser_get_fda_stats.add_argument('--exome-stats', '-e', required=True)
    parser_get_fda_stats.add_argument('--n-stats', '-n', required=True)
    parser_get_fda_stats.add_argument('--submission-number', required=True)
    parser_get_fda_stats.add_argument('--fastq1-path', required=True)
    parser_get_fda_stats.add_argument('--fastq1-md5sum', required=True)
    parser_get_fda_stats.add_argument('--fastq2-path', required=True)
    parser_get_fda_stats.add_argument('--fastq2-md5sum', required=True)
    parser_get_fda_stats.add_argument('--ref-path', required=True)
    parser_get_fda_stats.add_argument('--ref-version', required=True)
    parser_get_fda_stats.add_argument('--ref-md5sum', required=True)
    parser_get_fda_stats.add_argument('--gtf-path', required=True)
    parser_get_fda_stats.add_argument('--gtf-version', required=True)
    parser_get_fda_stats.add_argument('--gtf-md5sum', required=True)
    parser_get_fda_stats.add_argument('--run', required=True)
    parser_get_fda_stats.add_argument('--patient-name', required=True)
    parser_get_fda_stats.add_argument('--tumor-dna-depth', required=True)
    parser_get_fda_stats.add_argument('--norm-dna-depth', required=True)
    parser_get_fda_stats.add_argument('--output', '-o', required=True)

    return parser.parse_args()


def main():
    args = get_args()

    # Generating peptide sequences
    if args.command == 'make-snv-peptides-context':
        pep_utils.make_snv_peptides_context(args)
    if args.command == 'make-indel-peptides-context':
        pep_utils.make_indel_peptides_context(args)
    if args.command == 'make-fusion-peptides-context':
        pep_utils.make_fusion_peptides_context(args)
    if args.command == 'make-erv-peptides':
        pep_utils.make_erv_peptides(args)
    if args.command == 'make-viral-peptides':
        pep_utils.make_viral_peptides(args)
    if args.command == 'make-self-antigen-peptides':
        pep_utils.make_self_antigen_peptides(args)

    # Peptide quantification
    if args.command == 'get-peptide-read-count':
        rna_utils.get_peptide_read_count(args)
    if args.command == 'get-fusion-peptide-read-count':
        rna_utils.get_fusion_peptide_read_count(args)
    if args.command == 'get-splice-peptide-read-count':
        rna_utils.get_splice_peptide_read_count(args)
    if args.command == 'get-snv-peptide-read-count':
        rna_utils.get_snv_peptide_read_count(args)
    if args.command == 'get-indel-peptide-read-count':
        rna_utils.get_indel_peptide_read_count(args)

    # Peptide-related
    if args.command == 'filter-mutant-peptides':
        pep_utils.filter_mutant_peptides(args)
    if args.command == 'calculate-agretopicity':
        pep_utils.calculate_agretopicity(args)
    if args.command == 'add-agreto-metadata':
        metadata.add_agreto_metadata(args)

    # Metadata generation
    if args.command == 'add-snv-metadata':
        metadata.add_snv_metadata(args)
    if args.command == 'add-indel-metadata':
        metadata.add_indel_metadata(args)
    if args.command == 'add-erv-metadata':
        metadata.add_erv_metadata(args)
    if args.command == 'add-self-antigen-metadata':
        metadata.add_self_antigen_metadata(args)
    if args.command == 'add-viral-metadata':
        metadata.add_viral_metadata(args)
    if args.command == 'add-fusion-metadata':
        metadata.add_fusion_metadata(args)
    if args.command == 'add-splice-metadata':
        metadata.add_splice_metadata(args)

    # Additional metadata
    if args.command == 'add-tcga-data':
        metadata.add_tcga_data(args)

    # Reporting
    if args.command == 'make-lens-report':
        metadata.make_lens_report(args)
    if args.command == 'make-fda-stats':
        metadata.get_fda_stats(args)

    # Combining reports
    if args.command == 'aggregate-lens-reports':
        metadata.aggregate_lens_reports(args)

    # Variant filtering
    if args.command == 'filter-expressed-variants':
        var_utils.expressed_variants(args)
    if args.command == 'filter-isolated-variants':
        var_utils.isolated_variants(args)

    # BED generation
    if args.command == 'get-expressed-transcripts-bed':
        pep_utils.get_expressed_transcripts_bed(args)
    if args.command == 'get-expressed-ervs-bed':
        pep_utils.get_expressed_ervs_bed(args)
    if args.command == 'get-expressed-viral-bed':
        pep_utils.get_expressed_viral_bed(args)
    if args.command == 'get-expressed-selfs-bed':
        pep_utils.get_expressed_selfs_bed(args)

    if args.command == 'make-lens-bed':
        misc.make_lens_bed(args.report, args.gtf, args.samp_id, args.output)

    # Expression filtering
    if args.command == 'filter-expressed-self-genes':
        rna_utils.expressed_self_genes(args)
    if args.command == 'filter-expressed-viruses':
        rna_utils.filter_virdetect_by_counts(args)
    if args.command == 'check-erv-rna-coverage':
        rna_utils.check_erv_rna_coverage(args)
    if args.command == 'check-virus-rna-coverage':
        rna_utils.check_virus_rna_coverage(args)
    if args.command == 'filter-expressed-ervs':
        rna_utils.filter_expressed_ervs(args)
    if args.command == 'filter-expressed-transcripts':
        rna_utils.filter_expressed_txs(args.quants, args.metric, args.output, args.exclude_zeros, \
                                       args.percentile)

    # Onco-related
    if args.command == 'make-pyclone-vi-inputs':
        onco.make_pyclone_vi_inputs(args)
    if args.command == 'calculate-ccf':
        onco.calculate_ccf(args)

    # Misc
    if args.command == 'add-rna-normals':
        misc.add_rna_norms(args)
    if args.command == 'consolidate-multiqc-stats':
        misc.consolidate_multiqc_stats(args)
    if args.command == 'get-n-count-and-freq':
        metadata.get_n_count_and_freq(args)
    if args.command == 'get-fda-stats':
        metadata.get_fda_stats(args)

    # Plots
    if args.command == 'make-antigens-barplot':
        plots.make_antigens_barplot(args)
    if args.command == 'make-antigens-heatmap':
        plots.make_antigens_heatmap(args)
    if args.command == 'make-antigens-circos-plot':
        plots.make_antigens_circos_plot(args)

    # PMHC-related
    if args.command == 'aggregate-pmhc-outputs':
        metadata.aggregate_pmhc_outputs(args)
    if args.command == 'annotate-pmhcs':
        metadata.annot_pmhcs(args)
    if args.command == 'generic-annotation':
        metadata.generic_annotation(args)
    if args.command == 'generic-annotation-multikey':
        metadata.generic_annotation_multikey(args)
    if args.command == 'filter-peptides':
        pep_utils.filter_peptides(args)
    if args.command == 'prioritize-peptides':
        pep_utils.prioritize_peptides(args)


if __name__=='__main__':
    main()
