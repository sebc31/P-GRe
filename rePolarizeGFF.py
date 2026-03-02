import sys

#ARGV
#1 = GFF
#2 = A TSV file indicating the scaffolds name in the first column and their length in the second column

upto = 500
scaffold_upper_limit = {}
with open(sys.argv[2]) as scaffold_length:
    for line in scaffold_length:
        scaffold, length = line.replace("\n", "").split("\t")
        scaffold_upper_limit[scaffold] = int(length) - 1  # -1 is just in case their is a +/- 1 error in the inputted GFF file.

with open(sys.argv[1]) as gff:
    for line in gff:
        line = line.split("\t")
        line[3] = int(line[3]); line[4] = int(line[4])
        if line[2] == "pseudogene":
            pg_id = line[8].split("=")[1].split(";")[0] + "_three_prime_" + str(upto) + "_prior"

            if line[6] == "+":
                right = min(line[4] + upto, scaffold_upper_limit[line[0]])
                if right - line[4] > 15:
                    print("\t".join([line[0], str(line[4]), str(right), pg_id, "0", "+"]), sep="\t")

            elif line[6] == "-":
                if line[3] > 15:
                    left = max(line[3] - 500, 1)
                    print("\t".join([line[0], str(left), str(line[3]), pg_id, "0", "-"]), sep="\t")
