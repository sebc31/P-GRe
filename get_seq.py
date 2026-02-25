import sys
from re import sub

########################################################################################################################
#### FUNCTIONS

def insertNewlines(string, every=60):
    string=str(string)
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

########################################################################################################################
#### MAIN

#ARGV
#1. P-GRe annotation
#2. Miniprot raw results

#### MAIN
######## READING P-GRe OUPUT

#Loading GFF to keep pseudogenes id

pseudogene = set()
with open(sys.argv[1]) as results:
    for line in results:
        id = line.split("=")[1].split(";")[0]
        pseudogene.add(id)

#### MAIN
######## READING RAW MINIPROT GFF

sequence_dic = {}
j = 0

with open(sys.argv[2]) as raw_mp:
    for line in raw_mp:
        if line[0:5] == "##ATA":
            sequence = line
        elif "mRNA" in line:
            id = line.split("=")[1].split(";")[0]
            if id in pseudogene:
                j += 1
                pseudogene.remove(id)
                sequence_dic[id] = sequence
                sequence = sub(r"[^a-zA-Z*]", "", sequence)
                print(">" + id + "\n" + insertNewlines(sequence[3:]))
            else:
                sequence = ""



