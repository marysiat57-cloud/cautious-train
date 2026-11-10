#!/usr/bin/env python
'''
Assign MLST sequence types (STs) to a genome assembly based on detected
allele sequences.

Usage:
    mlst.py --asmid ASMID --allelesdb FILE --mlstdb FILE \
            --alleles FILE --output FILE

Arguments:
    --asmid      Genome assembly accession number.
    --allelesdb  FASTA file containing allele reference sequences.
    --mlstdb     TSV file containing MLST profiles and sequence types.
    --alleles    TSV file containing allele sequences detected in
                 the analysed genome.
    --output     Output TSV file with the assigned sequence type.
'''

import argparse
import polars as pl


def parse_args() -> argparse.Namespace:
    '''
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace with attributes:
            asmid      (str) - genome assembly accession number
            allelesdb  (str) - allele database FASTA file
            mlstdb     (str) - MLST profiles table
            alleles    (str) - detected allele sequences
            output     (str) - output file path
    '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--asmid', type=str, required=True,
        help='Assembly accession number (assembly ID)')
    parser.add_argument('--allelesdb', type=str, required=True,
        help='FASTA file with alleles database')
    parser.add_argument('--mlstdb', type=str, required=True,
        help='TSV file with MLST profiles')
    parser.add_argument('--alleles', type=str, required=True,
        help='TSV file with alleles sequences for one genome')
    parser.add_argument('--output', type=str, required=True,
        help='TSV file with STs for analysed genomes')

    args = parser.parse_args()
    return args


def load_allelesdb(fpath: str) -> dict:
    '''
    Load allele reference sequences from a FASTA file.
    
    Sequence identifiers are expected to have the format:
        GENE_ALLELE
    
    Returns:
        dict:
            {
                gene1: {sequence1: allele_no, ...},
                gene2: {sequence2: allele_no, ...},
                ...
            }
    '''
    
    seqs  = {}
    seqid = None
    f_in  = open(fpath)
    
    for line in f_in:
        if line[0] == '>':
            seqid = line[1:line.find(' ')]
            seqs[seqid] = ''
        else:
            seqs[seqid] += line[:-1]
    
    f_in.close()
    
    allelesdb = {}
    for seqid, seq in seqs.items():
        gene, allele = seqid.rsplit('_', 1)
        gene_alleles = allelesdb.get(gene, {})
        gene_alleles[seq] = int(allele)
        allelesdb[gene]   = gene_alleles
        
    return allelesdb


def mlst() -> None:
    '''
    Determine the MLST sequence type (ST) for a genome.
    
    The function:
      1. Loads the allele sequence database.
      2. Loads MLST profiles and detected allele sequences.
      3. Verifies that all seven MLST loci are present exactly once.
      4. Converts allele sequences to allele numbers.
      5. Searches the MLST profile database for a matching profile.
      6. Writes the assigned ST or -1 if no unambiguous match exists.
    '''
    
    args = parse_args()
    allelesdb  = load_allelesdb(args.allelesdb)
    mlstdb_df  = pl.read_csv(args.mlstdb, separator='\t')
    alleles_df = pl.read_csv(args.alleles, separator='\t')
    
    mlstout = 'asmid\tST\tmergedseq\n'
    seqout  = ''
    
    # A valid MLST profile must contain exactly seven loci and each
    # locus must appear only once.
    if alleles_df.shape[0] != 7 or alleles_df['qseqid'].unique().shape[0] != 7:
        mlstout +=  args.asmid + '\t-1\t\n'
    else:
        # Convert detected allele sequences into allele numbers using
        # the reference allele database.
        alleles = {}
        for gene, seq in alleles_df.select('qseqid', 'sseq').iter_rows():
            alleles[gene] = allelesdb[gene].get(seq, -1)
        
        # Find MLST profiles matching all allele assignments.
        hits_db = mlstdb_df.filter(**alleles)
        
        # Exactly one profile match is required for ST assignment.
        if hits_db.shape[0] != 1:
            mlstout +=  args.asmid + '\t-1\t\n'
        else:
            ST = hits_db.item(0, 'ST')
            
            # Concatenate allele sequences in sorted gene order to
            # create a deterministic sequence signature.
            mergedseq = ''.join(alleles_df.sort('qseqid')['sseq'])
            mlstout  +=  f'{args.asmid}\t{ST}\t{mergedseq}\n'
    
    # Write MLST result to the output file.
    with open(args.output, 'w') as f:
        f.write(mlstout)

if __name__ == '__main__':
    mlst()

