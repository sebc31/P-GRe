import sys
import Bio.SeqIO

########################################################################################################################
#### FUNCTIONS

def insertNewlines(string, every=60):
    string=str(string)
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

########################################################################################################################
#### MAIN

#ARGV
#1. P-GRe annotation
#2. P-GRe sequences

long_pg = set()

file_two = sys.argv[2]

new_fasta = open(file_two.replace(".fasta", "_high_confidence.fasta"), "w")

for sequence in list(Bio.SeqIO.parse(file_two, 'fasta')):
    if len(sequence.seq) > 20:
        new_fasta.write(">" + sequence.id + "\n" + insertNewlines(sequence.seq) + "\n")
        long_pg.add(sequence.id)

file_one = sys.argv[1]
new_gff = open(file_one.replace(".gff", "_high_confidence.gff"), "w")

with open(file_one) as gff:
    for line in gff:
        toks = line.replace("\n", "").split("\t")
        pgid = toks[8].split(";")[0].split("=")[1]
        if pgid in long_pg:
            new_gff.write(line)
