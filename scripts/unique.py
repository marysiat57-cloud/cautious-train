#!/usr/bin/env python
'''
Extract unique MLST sequence profiles and export them as FASTA.

Usage:
    unique.py --input FILE --output FILE

Arguments:
    --input   TSV file containing merged MLST results.
    --output  FASTA file with unique merged MLST locus sequences.
'''

import argparse
import polars as pl

def parse_args() -> argparse.Namespace:
    '''
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace with attributes:
            input  (str) - TSV file with merged MLST results
            output (str) - FASTA file with unique sequences
    '''
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input', type=str, required=True,
        help='TSV file with merged MLST results')
    parser.add_argument('--output', type=str, required=True,
        help='FASTA file with unique merged squences of MLST loci')

    args = parser.parse_args()
    return args

def unique() -> None:
    '''
    Extract unique MLST sequence profiles.
    
    The function:
      1. Loads MLST assignment results.
      2. Removes entries without a valid sequence type.
      3. Keeps one representative for each unique merged sequence.
      4. Exports unique sequences in FASTA format.
    '''
    
    args = parse_args()
    
    # Load MLST results.
    mlst_df = pl.read_csv(args.input, separator='\t')
    
    # Remove genomes for which no sequence type was assigned.
    mlst_df = mlst_df.filter( pl.col('ST') != -1 )
    
    # Keep only one entry for each unique merged sequence.
    mlst_df = mlst_df.unique('mergedseq')
    
    print(mlst_df)
    
    # Write unique MLST profiles to FASTA format.
    f_out = open(args.output, 'w')
    for ST, mergedseq in mlst_df.select('ST', 'mergedseq').rows():
        lines = '\n'.join(
            [ mergedseq[i:i+60] for i in range(0, len(mergedseq), 60) ]
        )
        fasta = f'>ST-{ST}\n{lines}\n'
        f_out.write(fasta)
    f_out.close()

if __name__ == '__main__':
    unique()

