import sys

#ARGV
# 1 = Diamond results
# 2 = PGRe's GFF file
# 3 = Outdir
# 4 = PGRe's pseudogenes' proteins file

if sys.argv[3][-1] != "/":
    sys.argv[3] += "/"

hqpg = set()
with open(sys.argv[1]) as diamond_results:
    for line in diamond_results:
        hqpg.add(line.replace("\n", ""))

filtered_gff = open(sys.argv[3] + "PGRe_filtered.gff", "w")
filtered_pg = set()
with open(sys.argv[2]) as pgre_gff:
    for line in pgre_gff:
        pgid = line.split("=")[1].split(";")[0]
        if pgid in hqpg:
            filtered_gff.write(line)
        else:
            filtered_pg.add(pgid)

filtered_pg = len(filtered_pg)
total_pg = len(hqpg) + filtered_pg
print("Removed", filtered_pg, "pseudogenes (", round(filtered_pg / total_pg * 100, 2), "%) with no homology with the given protein sequences.")

filtered_fasta = open(sys.argv[3] + "pseudogene_protein_filtered.fasta", "w")
with open(sys.argv[4]) as pgre_fasta:
    for line in pgre_fasta:
        if line[0] == ">":
            pgid = line.replace(">", "").replace("\n", "")
            if pgid in hqpg:
                to_save = True
            else:
                to_save = False
        if to_save:
            filtered_fasta.write(line)