# lenstools


# make-snv-peptides-context
Makes peptide sequences based on somatic SNVs observed in patient.

Inputs:
-vts/--var-tx-seqs: Directory containing FASTA files. FASTA files are named by
the following convention:
<SAMPLE_NAME>.<TRANSCRIPT>_<VARIANT_POS>.<"tumor"|"normal">.fa. For example,

```
<SAMPLE_NAME>.ENST00000248846_chr22_50233350.normal.fa
<SAMPLE_NAME>.ENST00000248846_chr22_50233350.tumor.fa
<SAMPLE_NAME>.ENST00000264709_chr2_25234373.normal.fa
<SAMPLE_NAME>.ENST00000264709_chr2_25234373.tumor.fa
<SAMPLE_NAME>.ENST00000328090_chr10_5730622.normal.fa
<SAMPLE_NAME>.ENST00000328090_chr10_5730622.tumor.fa
<SAMPLE_NAME>.ENST00000453712_chr4_98418725.normal.fa
<SAMPLE_NAME>.ENST00000453712_chr4_98418725.tumor.fa
<SAMPLE_NAME>.ENST00000597229_chr19_7520446.normal.fa
<SAMPLE_NAME>.ENST00000597229_chr19_7520446.tumor.fa
```

These files are generated using the `lenstools_make_variant_specific_tx_seqs`
process in the `lenstools.nf` Nextflow module. Briefly, transcripts are created
by creating variant-specific VCFs consisting of the variant of interest and any
proximal phased or homozygous germline variants. These variants are integrated
into the transcript sequencing using `bcftools consensus`.

Entries in the FASTA are exon sequences. These are reconstructed into the
transcript sequence within the function.


-st/-somatic-txs: List of transcripts that are expressed and harbor somatic variants. For example,

```
ENST00000264709
ENST00000453712
ENST00000328090
ENST00000597229
ENST00000248846
ENST00000287295
```

-sv/-somatic-vcf: Somatic VCF containing somatic SNVs.

-g/--gtf: GTF file containing gene annotations.  Used for reconstructing
transcript sequences from the exonic sequences contained within
`-vts/--var-tx-seqs`.

-l/--length: Length of peptides to emit (default: 11 amino acids).

-n/--context-nt-length: Length of surrounding nucleotide sequences around
variant of interest. Relevant for peptide manufacturability and discovering
codind sequence within RNA-Sequencing reads.

-o/--mt-output: Mutant peptide output FASTA.

--nt-output: Mutant nucleotide output FASTA.

-w/--wt-output: Wildtype peptide output FASTA.



# make-indel-peptides-context
Makes peptide sequences based on somatic InDels observed in patient.

Inputs:
-vts/--var-tx-seqs: Directory containing FASTA files. FASTA files are named by the following convention:
<SAMPLE_NAME>.<TRANSCRIPT>_<VARIANT_POS>.<"tumor"|"normal">.fa. For example,

```
<SAMPLE_NAME>.ENST00000248846_chr22_50233350.normal.fa
<SAMPLE_NAME>.ENST00000248846_chr22_50233350.tumor.fa
<SAMPLE_NAME>.ENST00000264709_chr2_25234373.normal.fa
<SAMPLE_NAME>.ENST00000264709_chr2_25234373.tumor.fa
<SAMPLE_NAME>.ENST00000328090_chr10_5730622.normal.fa
<SAMPLE_NAME>.ENST00000328090_chr10_5730622.tumor.fa
<SAMPLE_NAME>.ENST00000453712_chr4_98418725.normal.fa
<SAMPLE_NAME>.ENST00000453712_chr4_98418725.tumor.fa
<SAMPLE_NAME>.ENST00000597229_chr19_7520446.normal.fa
<SAMPLE_NAME>.ENST00000597229_chr19_7520446.tumor.fa
```

These files are generated using the `lenstools_make_variant_specific_tx_seqs`
process in the `lenstools.nf` Nextflow module. Briefly, transcripts are created
by creating variant-specific VCFs consisting of the variant of interest and any
proximal phased or homozygous germline variants. These variants are integrated
into the transcript sequencing using `bcftools consensus`.

Entries in the FASTA are exon sequences. These are reconstructed into the
transcript sequence within the function.

Note: make-indel-peptides-context creates `normal.fa` FASTAs, but does not use
them for peptide generation.


-st/-somatic-txs: List of transcripts that are expressed and harbor somatic
variants. For example,

``` ENST00000264709 ENST00000453712 ENST00000328090 ENST00000597229
ENST00000248846 ENST00000287295 ```

-sv/-somatic-vcf: Somatic VCF containing somatic SNVs.

-g/--gtf: GTF file containing gene annotations.  Used for reconstructing
transcript sequences from the exonic sequences contained within
`-vts/--var-tx-seqs`.

-l/--length: Length of peptides to emit (default: 11 amino acids).

-n/--context-nt-length: Length of surrounding nucleotide sequences around
variant of interest. Relevant for peptide manufacturability and discovering
codind sequence within RNA-Sequencing reads.

-o/--output: Mutant peptide output FASTA.

--nt-output: Mutant nucleotide output FASTA.

