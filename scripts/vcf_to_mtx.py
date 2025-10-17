#!/usr/bin/env python3

import subprocess
import pandas as pd
import sys
import tempfile
import os

def vcf_to_sample_variant_matrix(input_vcf, output_tsv):
    """
    Convert VCF to TSV matrix where:
    - Rows = samples
    - Columns = variants  
    - Values = 0 (hom ref), 1 (het/hom alt), NA (missing)
    """
    
    print(f"Processing VCF: {input_vcf}")
    
    # Extract sample names using bcftools
    print("Extracting sample names...")
    samples_result = subprocess.run(
        ['bcftools', 'query', '-l', input_vcf], 
        capture_output=True, text=True, check=True
    )
    samples = samples_result.stdout.strip().split('\n')
    print(f"Found {len(samples)} samples")
    
    # Extract genotype matrix using bcftools
    print("Extracting genotype data...")
    genotype_result = subprocess.run([
        'bcftools', 'query', 
        '-f', '%CHROM:%POS:%REF:%ALT[\\t%GT]\\n',
        input_vcf
    ], capture_output=True, text=True, check=True)
    
    # Process genotype data
    print("Converting genotypes and creating matrix...")
    variants = []
    genotype_matrix = []
    
    for line in genotype_result.stdout.strip().split('\n'):
        if not line:
            continue
            
        fields = line.split('\t')
        variant_id = fields[0]
        genotypes = fields[1:]
        
        variants.append(variant_id)
        
        # Convert genotypes: 0/0|0|0 -> 0, het/hom_alt -> 1, missing -> NA
        converted_gts = []
        for gt in genotypes:
            if gt in ['0/0', '0|0']:
                converted_gts.append('FALSE')
            elif gt in ['./.', '.|.', '.']:
                converted_gts.append('FALSE')  # pandas will convert to NaN
            elif '/' in gt or '|' in gt:
                # Any other phased/unphased genotype (het or hom alt)
                converted_gts.append('TRUE')
            else:
                converted_gts.append('FALSE')
        
        genotype_matrix.append(converted_gts)
    
    print(f"Found {len(variants)} variants")
    
    # Create DataFrame with variants as rows, samples as columns
    print("Creating pandas DataFrame...")
    df = pd.DataFrame(genotype_matrix, index=variants, columns=samples)
    
    # Transpose so samples are rows and variants are columns
    print("Transposing matrix...")
    df_transposed = df.T
    
    # Reset index to make sample names a column
    df_transposed.reset_index(inplace=True)
    df_transposed.rename(columns={'index': 'barcode'}, inplace=True)

    df_transposed = df_transposed.sort_values(by="barcode")
    
    # Save to TSV
    print(f"Saving to {output_tsv}...")
    df_transposed.to_csv(output_tsv, sep='\t', index=False, na_rep='NA')
    
    print("Matrix conversion complete!")
    print(f"Final dimensions: {len(samples)} samples x {len(variants)} variants")
    
    # Show preview
    print(f"\nPreview (first 3 samples, first 5 variants):")
    preview = df_transposed.head(3).iloc[:, :6]  # First 3 rows, first 6 columns (including barcode)
    print(preview.to_string(index=False))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 vcf_to_matrix.py input.vcf.gz output.tsv")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_tsv = sys.argv[2]
    
    try:
        vcf_to_sample_variant_matrix(input_vcf, output_tsv)
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
