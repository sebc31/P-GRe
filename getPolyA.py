import sys

#ARGV
# 1 = 3' sequences

poly_a_found = 0
with open(sys.argv[1]) as three_prime_file:
    for line in three_prime_file:
        if line[0] == ">":
            current_pseudogene = line.split("_")[0].replace(">", "")
        else:
            line = line.replace("\n", "")
            for i in range(0, len(line)):
                if line[i:i+20].upper().count("A") >= 15:
                    poly_a_found += 1
                    print(current_pseudogene, line[i:i+20], sep = "\t")
                    break
