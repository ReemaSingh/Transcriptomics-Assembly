#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = "Trinity.fasta" # Input fasta file
wanted_file = "Spikeins" # Input interesting sequence IDs, one per line
result_file = "result_file.fasta" # Output fasta file

wanted = set()
with open(wanted_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id not in  wanted:
            SeqIO.write([seq], f, "fasta")
