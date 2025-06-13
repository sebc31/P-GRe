import sys

def insertNewlines(string, every=60):
    string=str(string)
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

ids = set()
write = False
with open(sys.argv[1]) as matching_query:
    for line in matching_query:
        ids.add(line.replace("\n",""))

with open(sys.argv[2]) as large_fasta:
    for line in large_fasta:
        if line[0] == ">":
            if line.replace("\n","").replace(">","").split()[0] in ids:
                write = True
            else:
                write = False
        if write:
            print(insertNewlines(line), end="")
