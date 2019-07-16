from Bio import SeqIO
import argparse
import sys

parser = argparse.ArgumentParser(description='Run flu builds')
parser.add_argument('-l', '--list', type = str, help='list of strains to select')
parser.add_argument('-f', '--fasta', type = str, help='fasta file to select from')
parser.add_argument('-o', '--output', type = str, help='name of output file')
args = parser.parse_args()

with open(args.list, 'r') as f:
    seqs = [ line.strip() for line in f.readlines() ]

records = []
for record in SeqIO.parse(args.fasta, "fasta"):
    if record.id in seqs:
        records.append(record)

with open(args.output, "w") as handle:
    SeqIO.write(records, handle, "fasta")
