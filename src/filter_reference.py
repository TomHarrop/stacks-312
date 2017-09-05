#!/usr/bin/env python3

from Bio import SeqIO

# process snakemake object
reference = snakemake.input['reference']
filtered_reference = snakemake.output['filtered_reference']
min_len = snakemake.params['min_len']

# filter scaffolds into memory
seqrecords = SeqIO.parse(reference, 'fasta')
filtered_scaffolds = list(x for x in seqrecords
                          if len(x) >= min_len)

# write output
SeqIO.write(sequences=filtered_scaffolds,
            handle=filtered_reference,
            format='fasta')
