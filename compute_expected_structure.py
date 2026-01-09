import sys
import math

#ARGV
#1 = GFF

id_has_structure_name = True
transcript_structure = {}
transcript_strand = {}
tot = {}
count = 1
with open(sys.argv[1]) as gff:
    for line in gff:
        if line[0] != "#":
            line = line.replace("\n","").split("\t")
            if line[2] == "CDS":
                value = ""
                left = int(line[3])
                right = int(line[4])
                attributes = line[8].split(";")
                for attribute in attributes:
                    key, value = attribute.split("=")
                    if key == "Parent":
                        if id_has_structure_name:
                            value = value.split(":")[1]
                        if value not in transcript_structure:
                            transcript_structure[value] = []
                            transcript_strand[value] = line[6]
                            tot[value] = 0
                        break
                if value == "":
                    print("ERROR in GFF file. Line " + str(count) + " has a CDS with no parent structure.")
                    exit(1)

                transcript_structure[value].append([min(left, right), max(left, right)])
                tot[value] += transcript_structure[value][-1][1]-transcript_structure[value][-1][0]+1
        count += 1

for transcript in transcript_structure:
    tot[transcript] /= 3
    previous_end = 0
    previous_hang = 3
    if transcript_strand[transcript] == "+":
        transcript_structure[transcript] = sorted(transcript_structure[transcript], key=lambda x: x[0])
    else:
        transcript_structure[transcript] = sorted(transcript_structure[transcript], key=lambda x: x[0], reverse=True)
    for i in range(0, len(transcript_structure[transcript])):
        true_size = transcript_structure[transcript][i][1] - transcript_structure[transcript][i][0] + 1
        non_five_multriframe_codon_size = true_size - (3 - previous_hang)
        non_ft_multriframe_coded_size = math.floor(non_five_multriframe_codon_size/3)

        if previous_hang == 3:
            transcript_structure[transcript][i][0] = previous_end + 1
            transcript_structure[transcript][i][1] = non_ft_multriframe_coded_size + previous_end
        else:
            transcript_structure[transcript][i][0] = previous_end + 2
            transcript_structure[transcript][i][1] = non_ft_multriframe_coded_size + previous_end + 1

        previous_hang = non_five_multriframe_codon_size%3
        if previous_hang == 0:
            previous_hang = 3
        previous_end = transcript_structure[transcript][i][1]

for transcript in transcript_structure:
    print(transcript, end="\t")
    for prot_part in transcript_structure[transcript]:
        print(prot_part[0],":",prot_part[1], end=";", sep="")
    print()

