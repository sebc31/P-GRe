import sys

#ARGV
# 1 = PGRe.gff
# 2 = poly(A) file

polya_found = set()
with open(sys.argv[2]) as polya_tail:
    for line in polya_tail:
        polya_found.add(line.split("\t")[0])

with open(sys.argv[1]) as gff:
    for line in gff:
        line = line.split("\t")
        if line[2] == "pseudogene":
            pg_id = line[8].split("=")[1].split(";")[0]
            if pg_id in polya_found:
                line[2] = "processed_pseudogene"
        print("\t".join(line), end = "")