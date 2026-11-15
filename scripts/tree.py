#!/usr/bin/env python
'''
Construct and visualize a phylogenetic tree from MLST sequence
alignments.

Usage:
    tree.py --mlst FILE --seqs FILE      \
            --outtree FILE --outvis FILE

Arguments:
    --mlst     TSV file containing MLST results.
    --seqs     FASTA multiple sequence alignment (MSA).
    --outtree  Output tree in Newick format.
    --outvis   Tree visualization image file.
'''

import argparse
import polars as pl
import toytree, toyplot
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator

def parse_args() -> argparse.Namespace:
    '''
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace with attributes:
            mlst     (str) - MLST results table
            seqs     (str) - aligned sequences in FASTA format
            outtree  (str) - output Newick tree file
            outvis   (str) - output tree visualization file
    '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--mlst', type=str, required=True,
        help='TSV file with merged MLST results')
    parser.add_argument('--seqs', type=str, required=True,
        help='FASTA multiple sequnece alignment (MSA) file')
    parser.add_argument('--outtree', type=str, required=True,
        help='Computed tree in Newick format')
    parser.add_argument('--outvis', type=str, required=True,
        help='Tree visualisation')

    args = parser.parse_args()
    return args

def construct_tree(seqs_fpath: str) -> Phylo.BaseTree.Tree:
    '''
    Construct a neighbour-joining phylogenetic tree from a
    multiple sequence alignment.
    
    Args:
        seqs_fpath (str):
            FASTA alignment file.
    
    Returns:
        Bio.Phylo tree object.
    '''
    
    aln = AlignIO.read( open(seqs_fpath) , 'fasta')
    constructor = DistanceTreeConstructor()
    calculator  = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    tree = constructor.nj(dm)
    return tree

def draw_tree(tree: Phylo.BaseTree.Tree, counts: pl.DataFrame) -> toyplot.canvas.Canvas:
    '''
    Create a graphical representation of the phylogenetic tree.
    
    The figure contains:
      - a phylogenetic tree with branch lengths
      - a bar chart showing the number of isolates per ST
    
    Args:
        tree:
            Toytree tree object.
        counts:
            DataFrame containing counts of isolates per ST.
    
    Returns:
        toyplot.Canvas containing the visualisation.
    '''
    
    # Scale figure height according to the number of tips.
    height = tree.ntips * 50
    canvas = toyplot.Canvas(width=800, height=height)
    canvas.style['background-color'] = 'white'
    
    # Left panel: phylogenetic tree.
    # Right panel: isolate counts per sequence type.
    ax1 = canvas.cartesian(bounds=(10,  600, 20, height-40))
    ax2 = canvas.cartesian(bounds=(610, 790, 10, height-30))
    
    ax1.show = False
    
    ax2.y.show = False
    ax2.x.show = True
    ax2.x.ticks.show = True
    
    # Draw horizontal bars indicating the number of genomes
    # associated with each sequence type.
    ax2.bars(
        range(tree.ntips),
        [counts.filter( pl.col('ST') == label ).item(0, 'count') for label in tree.get_tip_labels()],
        color = 'lightgray',
        along = 'y'
    )
    
    # Draw the phylogenetic tree.
    tree.draw(
        axes = ax1, tip_labels_align = True,
        edge_align_style = {
            'stroke'           : 'lightgray',
            'stroke-width'     : 0.8,
            'stroke-dasharray' : '1,1'
        },
        tip_labels_style = {
            'font-family'  : 'Monospace',
            'font-weight'  : 'bold',
            'font-size'    : 20,
            'anchor-shift' : 10
        },
    )
    
    # Add branch length labels to tree edges.
    dists = tree.get_node_data('dist')
    tree.annotate.add_edge_labels(
        ax1, labels=dists.round(3), color='green', font_size=18, yshift=-10,
        mask=False, style={'font-family' : 'Monospace'}
    )
    
    # Add distance scale bar.
    tree.annotate.add_axes_scale_bar(ax1)
    
    return canvas

def tree() -> None:
    '''
    Generate a phylogenetic tree and its visualization.
    
    The function:
      1. Builds a neighbour-joining tree from aligned sequences.
      2. Saves the tree in Newick format.
      3. Counts genomes assigned to each sequence type.
      4. Generates a combined tree and abundance visualization.
      5. Saves the resulting figure.
    '''
    
    args = parse_args()
    
    # Construct phylogenetic tree from aligned sequences.
    tree = construct_tree(args.seqs)
    
    # Save the generated tree in Newick format.
    Phylo.write(tree, args.outtree, 'newick')
    
    # Load MLST assignments.
    mlst_df = pl.read_csv(args.mlst, separator='\t')
    
    # Format sequence types as ST-<number>.
    mlst_df = mlst_df.with_columns(
        (pl.lit('ST-') + pl.col('ST').cast(pl.String) ).alias('ST')
    )
    
    # Count the number of isolates per sequence type.
    counts = mlst_df['ST'].value_counts()
    
    # Load the saved tree into Toytree for visualization.
    tree = toytree.tree(args.outtree)
    
    # Create and save the final figure.
    canvas = draw_tree(tree, counts)
    toytree.save(canvas, args.outvis)

if __name__ == '__main__':
    tree()

