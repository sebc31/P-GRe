import sys
import re

# ARGV
# 1 = results from miniprot
# 2 = mrna_coordinates.bed
# 3 = outdir
# 4 = in_organism_prot.id

if sys.argv[3][-1] != "/":
    sys.argv[3] += "/"

########################################################################################################################
#### FUNCTIONS

def addAnchors(currentProtFiltered, file):
    with open(file) as gene:
        for line in gene:
            chr, left, right = line.replace("\n", "").split("\t")
            id = chr + "_" + left + "_" + right
            currentProtFiltered[id] = {}
            currentProtFiltered[id]['chr'] = chr
            currentProtFiltered[id]['spawn'] = [int(left), int(right)]
            currentProtFiltered[id]['size'] = currentProtFiltered[id]['spawn'][1] - currentProtFiltered[id]['spawn'][0] + 1
            currentProtFiltered[id]['type'] = "gene"

def filterOverlap(currentProtFiltered, pseudogeneDict, scoring_param='size'):
    to_remove = set()
    for i in range(1, len(currentProtFiltered)):
    
        overlap_factor = 0
        overlap_size = currentProtFiltered[i - 1][1]['spawn'][1] - currentProtFiltered[i][1]['spawn'][0] + 1
        overlap_treshold = min(currentProtFiltered[i - 1][1]['size'] * overlap_factor, currentProtFiltered[i][1]['size'] * overlap_factor)
        
        if overlap_size > overlap_treshold and currentProtFiltered[i - 1][1]['chr'] == currentProtFiltered[i][1]['chr']:

            if currentProtFiltered[i - 1][1]['type'] == "gene":
                to_remove.add(currentProtFiltered[i][0])
                currentProtFiltered[i] = currentProtFiltered[i - 1]

            elif currentProtFiltered[i][1]['type'] == "gene":
                to_remove.add(currentProtFiltered[i - 1][0])
                currentProtFiltered[i - 1] = currentProtFiltered[i]

            elif currentProtFiltered[i - 1][1]['in_organism_prot'] == True and currentProtFiltered[i][1]['in_organism_prot'] == False:
                to_remove.add(currentProtFiltered[i][0])
                currentProtFiltered[i] = currentProtFiltered[i - 1]

            elif currentProtFiltered[i][1]['in_organism_prot'] == True and currentProtFiltered[i - 1][1]['in_organism_prot'] == False:
                to_remove.add(currentProtFiltered[i - 1][0])
                currentProtFiltered[i - 1] = currentProtFiltered[i]

            elif currentProtFiltered[i - 1][1][scoring_param] > currentProtFiltered[i][1][scoring_param]:
                to_remove.add(currentProtFiltered[i][0])
                currentProtFiltered[i] = currentProtFiltered[i - 1]
            else:
                to_remove.add(currentProtFiltered[i - 1][0])
                currentProtFiltered[i - 1] = currentProtFiltered[i]

    for i in range(0, len(currentProtFiltered)):
        if currentProtFiltered[i][0] not in to_remove and currentProtFiltered[i][1]['type'] != "gene":
            pseudogeneDict[currentProtFiltered[i][0]] = currentProtFiltered[i][1]

    return pseudogeneDict

def readGFF(file):
    currentProt = {}
    with open(file) as blast:
        for line in blast:

            if line[0] != "#":
                line = line.replace("\n", "").split('\t')

                # Get informations from mRNA lines
                if line[2] == "mRNA" or line[2] == "pseudogene":
                    # Get id and uses it to generate a unique pseudogene identifier
                    id = line[8].split(";")[0].split("=")[-1]

                    # Get parent name
                    parent = line[8].split("=")[-1].split(" ")[0]

                    # The next lines add necessary infomartions to the dictionnary
                    currentProt[id] = {}
                    currentProt[id]['chr'] = line[0]
                    currentProt[id]['score'] = int(line[5])
                    currentProt[id]['cds_unf'] = {}
                    currentProt[id]['total_align_length'] = []
                    currentProt[id]['identity'] = float(re.search(r"Identity=([0-1]\.[0-9]+)", line[8]).group(1))
                    currentProt[id]['type'] = "prediction"
                    currentProt[id]['spawn'] = [int(line[3]), int(line[4])]
                    currentProt[id]['size'] = int(line[4]) - int(line[3]) + 1
                    currentProt[id]['other'] = "\t".join([line[6], line[7], line[8]])

                    if parent in parent_list:
                        currentProt[id]['in_organism_prot'] = True
                    else:
                        currentProt[id]['in_organism_prot'] = False

                # Get informations from CDS lines
                if line[2] == "CDS":

                    cds_id = line[0] + "_" + line[3] + "_" + line[4]
                    currentProt[id]['cds_unf'][cds_id] = {}
                    currentProt[id]['cds_unf'][cds_id]['left'] = int(line[3])
                    currentProt[id]['cds_unf'][cds_id]['right'] = int(line[4])
                    currentProt[id]['cds_unf'][cds_id]['score'] = int(line[5])
                    currentProt[id]['cds_unf'][cds_id]['other'] = "\t".join([line[6], line[7], line[8]])

    return currentProt

########################################################################################################################
####################################### IN-ORGANISM FILTERING ##########################################################
########################################################################################################################

########################################################################################################################
#### IN-ORGANISM PROT. LIST

parent_list = set()
with open(sys.argv[4]) as parent_file:
    for line in parent_file:
        parent_list.add(line.replace("\n",""))

########################################################################################################################
#### GFF PARSING

currentProt = readGFF(sys.argv[1])

#Next block add the gene coordinates to remove eventual predictions that overlap them
addAnchors(currentProt, sys.argv[2])
currentProt = sorted(currentProt.items(),key=lambda x: (x[1]['chr'], x[1]['spawn'][0], x[1]['spawn'][1]))
# Note that the sorted() function transforms the dictionnary structure into a list of tuple like:
# [(id1,currentProt[id1]),(id2,currentProt[id2]),...]

########################################################################################################################
#### OVERLAPS FILTERING

pseudogeneDict = {}
pseudogeneDict = filterOverlap(currentProt, pseudogeneDict)
del currentProt

########################################################################################################################
################################################# OUTPUTS ##############################################################
########################################################################################################################

good_lines = open(sys.argv[3] + "miniprot_res_filtered.gff", "w")
for pseudog in pseudogeneDict:

    chr = pseudogeneDict[pseudog]['chr']
    left = pseudogeneDict[pseudog]['spawn'][0]
    right = pseudogeneDict[pseudog]['spawn'][1]
    score = pseudogeneDict[pseudog]['score']
    attributes = [chr, "miniprot", "pseudogene", str(left), str(right), str(score), pseudogeneDict[pseudog]['other']]
    line = "\t".join(attributes)
    good_lines.write(line + "\n")

    for cds in pseudogeneDict[pseudog]['cds_unf']:
        left = pseudogeneDict[pseudog]['cds_unf'][cds]['left']
        right = pseudogeneDict[pseudog]['cds_unf'][cds]['right']
        score = pseudogeneDict[pseudog]['cds_unf'][cds]['score']
        attributes = [chr, "miniprot", "pseudogenic_CDS", str(left), str(right), str(score), pseudogeneDict[pseudog]['cds_unf'][cds]['other']]
        line = "\t".join(attributes)
        good_lines.write(line + "\n")
