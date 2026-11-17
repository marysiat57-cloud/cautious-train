#!/usr/bin/env python
'''
Filter tabular BLAST output by query coverage per subject (qcovs).

Usage:
    filter_blast.py [ --asmid ASMID ] --qcovs INT --cols "COL1 COL2 ..." \
                      --input FILE --output FILE

Arguments:
    --asmid   Genome assembly accession number. Not required. If provided the
              value will be included in the output table in a separate column
              named 'asmid'.
    --qcovs   Minimum query coverage per subject (%).
    --eqlen   Alignment length must be equal to the query sequence length.
              Presence flags True.
    --cols    Whitespace-separated column names in -outfmt order,
              must include 'qcovs'.
    --input   BLAST tabular results (tab-separated, no header).
    --output  Filtered results (tab-separated, with header).
'''

import argparse
import polars as pl
from pathlib import Path

def parse_args() -> argparse.Namespace:
    '''
    Parse command-line arguments.

    Returns:
        argparse - Namespace with attributes:
            asmid  (str)  - genome assembly accession number
            qcovs  (int)  - query coverage threshold
            eqlen  (bool) - alignment length must be equal
                            to the query sequence length
            cols   (str)  - whitespace-separated column names
            input  (str)  - path to BLAST results
            output (str)  - path for filtered results
    '''
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--asmid', type=str, required=False,
        help='Assembly accession number (assembly ID)')
    parser.add_argument('--qcovs', type=int, required=True,
        help='Query coverage by subject threshold')
    parser.add_argument('--eqlen', action='store_true',
        help='Alignment length must be equal to the query sequence length ' +
             '(if combined with --qcovs 100, ensures gapless alignment)')
    parser.add_argument('--cols', type=str, required=True,
        help='Column names for BLAST results')
    parser.add_argument('--input', type=str, required=True,
        help='BLAST results')
    parser.add_argument('--output', type=str, required=True,
        help='Filtered results')
    
    args = parser.parse_args()
    
    return args


def filter_hits() -> None:
    '''
    Load BLAST hits to a Polars DataFrame. Keep only those with
    qcovs >= threshold and (optionally) with the query length equal
    to the alignment lenght. Optionally add the assembly accession number
    to the output DataFrame. Save filtered results to a TSV file.
    '''
    
    args = parse_args()
    
    # Load BLAST results to a Polars DataFrame. BLAST tabular output has
    # no header. Column names come from --cols (args.cols). If the file size
    # is 0, create an empty DataFrame with args.cols columns.
    p = Path(args.input)
    if p.stat().st_size != 0:
        df = pl.read_csv(
            args.input,
            has_header     = False,
            new_columns    = args.cols.split(),
            separator      = '\t'
        )
    else:
        df = pl.DataFrame( schema=args.cols.split() )
    
    df = df.filter( pl.col('qcovs') >= args.qcovs )
    
    if args.eqlen:
        df = df.filter( pl.col('qlen')  == pl.col('length') )
    
    # If the genome assembly accession is provided, add its value in a new
    # column named 'asmid'.
    if args.asmid is not None:
        df.insert_column(0, pl.lit( args.asmid ).alias('asmid'))
    
    # Write the filtered results to a TSV file, write_csv() includes
    # a header row by default.
    df.write_csv(args.output, separator = '\t')

if __name__ == '__main__':
    filter_hits()

