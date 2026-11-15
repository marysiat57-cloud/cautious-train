# MLST pipeline for Staphylococcus aureus genome assemblies
# blastn -> filter -> merge_flt -> mlst -> merge_mlst -> unique -> tree

# Indicate the config file (it defines 'genomesdir')
configfile:
    'config.yml'

# Extract query sequence names (sequence IDs) from file paths
# to use them as IDs for query files (the 7 MLST genes)
seqids, = glob_wildcards('input/queries/{seqid}.fna')

# Extract genome assembly accession numbers (assembly IDs) from file paths
# to use them as IDs for genome files
asmids, = glob_wildcards(config['genomesdir'] + '/{asmid}.fna.gz')

# Rule all: collect the final output of the pipeline:
# (the tree in the Newick format and its visualisation)
rule all:
    input:
        newick = 'output/tree.newick',
        png    = 'output/tree.png'

# Rule blastn: for each genome perform a nucleotide BLAST search using:
# - a single-sequence FASTA query file (allele no. 1 of one MLST gene)
# - a multiple-sequence FASTA file (one staphylococcal genome assembly contigs)
# Uncompress the gzipped genome assembly file and redirect the output
# to blastn's stdin while asking blastn to read the subject FASTA sequence
# from that stream ('-' instead of a subject file path)
rule blastn:
    conda:
        'envs/powb-blast.yml'
    params:
        evalue = 0.01,
        cols   = 'qseqid qlen sseqid length qcovs pident evalue sseq'
    input:
        query = 'input/queries/{seqid}.fna',
        subj  = config['genomesdir'] + '/{asmid}.fna.gz'
    output:
        'output/blastn/blastn_{asmid}_{seqid}.tsv'
    log:
        'log/blastn/blastn_{asmid}_{seqid}.log'
    shell:
        '''gunzip -c {input.subj} 2> {log} |      \
           blastn -evalue         {params.evalue} \
                  -outfmt      "6 {params.cols}"  \
                  -query          {input.query}   \
                  -subject        -               \
                  -out            {output}        \
                   >> {log} 2>&1
        '''

# Rule filter: for each search retain 100% query coverage hits (alignments)
# of the same length as the query sequence (gapless alignment)
rule filter:
    conda:
        'envs/powb-python.yml'
    params:
        qcovs = 100,
        cols  = rules.blastn.params.cols
    input:
        rules.blastn.output
    output:
        'output/filter/filtered_{asmid}_{seqid}.tsv'
    log:
        'log/filter/filter_{asmid}_{seqid}.log'
    shell:
        '''scripts/filter.py --eqlen                    \
                             --asmid  {wildcards.asmid} \
                             --qcovs  {params.qcovs}    \
                             --cols  "{params.cols}"    \
                             --input  {input}           \
                             --output {output}          \
                               > {log} 2>&1
        '''

# Rule merge_flt: for each genome join the filtered BLAST results of all
# genes into one file using a short Bash script
# ({seqid} is expanded, {asmid} stays a wildcard: one output file per genome)
rule merge_flt:
    params:
        mask = rules.filter.output[0].replace('{seqid}', '*')
    input:
        expand(rules.filter.output, seqid=seqids, allow_missing=True)
    output:
        'output/merge_flt/merged_{asmid}.tsv'
    log:
        'log/merge_flt/merge_{asmid}.log'
    shell:
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail -n +2 ${{file}}
           done >> {output} 2>> {log}
        '''

#DONE:
# Rule mlst: perform MLST analysis, one genome at a time.
# The allele sequences found in the genome are looked up in the PubMLST
# database files to get the allele numbers and then the sequence type (ST)
rule mlst:
    conda:
        'envs/powb-python.yml'
    input:
        alleles   = rules.merge_flt.output,
        allelesdb = 'input/alleles.fna',
        mlstdb    = 'input/bigsdb.tsv'
    output:
        'output/mlst/mlst_{asmid}.tsv'
    log:
        'log/mlst/mlst_{asmid}.log'
    shell:
        '''scripts/mlst.py --asmid     {wildcards.asmid} \
                           --allelesdb {input.allelesdb} \
                           --mlstdb    {input.mlstdb}    \
                           --alleles   {input.alleles}   \
                           --output    {output}          \
                             > {log} 2>&1
        '''

#DONE:
# Rule merge_mlst: join the MLST results of all genomes into one file
# using a short Bash script (no wildcards left in the output: one file)
rule merge_mlst:
    params:
        mask = rules.mlst.output[0].replace('{asmid}', '*')
    input:
        expand(rules.mlst.output, asmid=asmids)
    output:
        'output/merge_mlst/merged.tsv'
    log:
        'log/merge_mlst.log'
    shell:
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail -n +2 ${{file}}
           done >> {output} 2>> {log}
        '''

#DONE:
# Rule unique: from the merged MLST results extract unique merged MLST loci
# sequences (i.e., one sequence per ST) to a FASTA file (it serves as an MSA)
rule unique:
    conda:
        'envs/powb-python.yml'
    input:
        rules.merge_mlst.output
    output:
        'output/unique.fna'
    log:
        'log/unique.log'
    shell:
        '''scripts/unique.py --input  {input}  \
                             --output {output} \
                               > {log} 2>&1
        '''

#DONE:
# Rule tree: take the merged MLST results and the extracted merged MLST loci
# sequences, and then compute and visualise a phylogenetic tree (NJ tree)
# together with the number of genomes per ST

rule tree:
    conda:
        'envs/powb-python.yml'
    input:
        merged = rules.merge_mlst.output,
        unique = rules.unique.output
    output:
        newick = 'output/tree.newick',
        png    = 'output/tree.png'
    log:
        'log/tree.log'
    shell:
        '''scripts/tree.py --mlst    {input.merged}  \
                           --seqs    {input.unique}  \
                           --outtree {output.newick} \
                           --outvis  {output.png}    \
                             > {log} 2>&1
        '''
