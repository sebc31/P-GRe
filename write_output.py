import sys

########################################################################################################################
#### FUNCTIONS

def initAtt(attributes_string):
    attributes = attributes_string.split(";")
    att_dic = {}
    for attribute in attributes:
        key, value = attribute.split("=")
        att_dic[key] = value
    return att_dic

def initDic(dictionarry, chr, type, left, right, score, strand, frame, in_org_state, att_dic):
    dictionarry = {}
    dictionarry["chr"] = chr
    dictionarry["method"] = "PGRe"
    dictionarry["type"] = type
    dictionarry["left"] = int(left)
    dictionarry["right"] = int(right)
    dictionarry["score"] = int(score)
    dictionarry["strand"] = strand
    dictionarry["frame"] = frame
    dictionarry["in_org"] = in_org_state
    dictionarry["attributes"] = attributes
    dictionarry["Frameshift"] = "0"
    for element in att_dic:
        dictionarry[element] = att_dic[element]
    return dictionarry

def insertNewlines(string, every=60):
    string=str(string)
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

def printGFFLine(type, dictionarry, current_structure):
    gff_line = [dictionarry[type][current_structure]["chr"], \
                dictionarry[type][current_structure]["method"], \
                dictionarry[type][current_structure]["type"], \
                str(dictionarry[type][current_structure]["left"]), \
                str(dictionarry[type][current_structure]["right"]), \
                str(dictionarry[type][current_structure]["score"]), \
                dictionarry[type][current_structure]["strand"], \
                dictionarry[type][current_structure]["frame"]]
    print("\t".join(gff_line))

########################################################################################################################
#### MAIN

#ARGV
#1. In-organism protein identifiers list (with Diamond homology)
#2. Protein expected alignments
#3. Miniprot filtered gff
#4. Outdir

if sys.argv[4][-1] != "/":
    sys.argv[4] += "/"

#### MAIN
######## LOADING FILES

#Loading the list of parent genes
parent_dic = {}
all_in_organism_prot = set()
with open(sys.argv[1]) as dmd:
    for line in dmd:
        key, value = line.replace("\n","").split("\t")
        parent_dic[key] = value
        all_in_organism_prot.add(value)

#Loading the expeceted proteins structures (/alignements)
protein_structure = {}
with open(sys.argv[2]) as gs:
    for line in gs:
        transcript, structure = line.replace("\n","").strip(";").split("\t")
        protein_structure[transcript] = []
        structure = structure.split(";")
        for protein_part in structure:
            protein_structure[transcript].append(set())
            begin, end = protein_part.split(":")
            for i in range(int(begin), int(end) + 1):
                protein_structure[transcript][-1].add(i)

#### MAIN
######## READING MINIPROT RESULT

fasta_file = open(sys.argv[4] + "pseudogene_protein.fasta","w")
result = {}
sequence = ""
with open(sys.argv[3]) as mini:
    for line in mini:
        #ATA lines are used to retrieve the (pseudo)protein sequences
        if "##ATA" in line:
            sequence = line.replace("#ATA","").replace(".","").replace("$","").replace("!","").replace(" ","").replace("\t","").replace("\n","")

        #Other non-commentary lines are tretaed as classic GFF lines
        if line[0] != "#":

            #Next block recovers all the informations for the line.
            chr, method, type, left, right, score, strand, frame, attributes = line.replace("\n", "").split("\t")
            att_dic = initAtt(attributes)

            #For each "mRNA-type" lines, as subdictionnary is created to contains informations about CDS and eventual stop_codon
            if type == "pseudogene":
                #Next line check if the prediction is a "good one"
                fasta_file.write(">" + att_dic["ID"] + "\n" + sequence + "\n")
                att_dic["Target"] = att_dic["Target"].split(" ")[0]

                #Next lines check if predictions were made from in-organism protein. This is made to eventually filter more strictly the other pseudogenes
                if att_dic["Target"] in all_in_organism_prot:
                    in_org = True
                else:
                    in_org = False

                #Next lines change the parent of pseudogenes when their parent didn't originate from the organism,
                #yet homology was found between the parent and an organism protein
                if att_dic["Target"] in parent_dic.keys():
                    att_dic["Target"] = parent_dic[att_dic["Target"]]
                else:
                    type = "unitary_pseudogene"
                result[att_dic["ID"]] = {}
                result[att_dic["ID"]] = initDic(result[att_dic["ID"]], chr, type, left, right, score, strand, frame, in_org, att_dic)
                result[att_dic["ID"]]["CDS"] = {}

            #"CDS-type" lines
            if type == "pseudogenic_CDS":
                in_org = "CDS"
                cds_id = chr + "_" + left + "_" + right
                covered_protein_part = att_dic["Target"].split(" ")
                att_dic["Target"] = [int(covered_protein_part[1]), int(covered_protein_part[2])]
                result[att_dic["Parent"]]["CDS"][cds_id] = {}
                result[att_dic["Parent"]]["CDS"][cds_id] = initDic(result[att_dic["Parent"]], chr, type, left, right, score, strand, frame, in_org, att_dic)
                result[att_dic["Parent"]]["CDS"][cds_id]["covered_pos"] = set()
                #Next lines recover all positions covered by the CDS. This is made to infer the pseudogene type later on
                for i in range(att_dic["Target"][0], att_dic["Target"][1] + 1):
                    result[att_dic["Parent"]]["CDS"][cds_id]["covered_pos"].add(i)

#### MAIN
######## PREDICTING PSEUDOGENE TYPE

check_type = True
if check_type:
    for pseudogene in result:
        if result[pseudogene]["type"] != "unitary_pseudogene":
            #Intron gain or loss can not be infered for monoexonic gene, thus pseudogenes predicted from monoexonic parent
            #are simply annotated as pseudogene
            if result[pseudogene]["Target"] in protein_structure:
                expected_aln = protein_structure[result[pseudogene]["Target"]]

                #Next (large) block get the structure of the pseudogene. It checks wether or not each expected protein part
                # (i.e. part of protein coded by each parent gene's CDSs) is indeed associated with a predicted CDS.
                if len(expected_aln) != 1:
                    retained_introns = 0
                    lost_introns = 0
                    result[pseudogene]["structure"] = []
                    for cds in result[pseudogene]["CDS"]:
                        is_overlapping_structure = False
                        covered_pos = result[pseudogene]["CDS"][cds]["covered_pos"]
                        #Next block check, for each expected protein part, if it is covered by a predicted CDS by at least 45%
                        #of its length. Protein part that are thus covered are considered highly conserved, and intron gains
                        #or losses will only be checked between these protein parts. Note that this threshold is < 50%,
                        #meaning that a protein part can be considered encoded by two predicted CDS. This allow to find potential
                        #predicted intron that are actually large insertion events. Also note that CDS with score < 0 are not
                        #taken into account
                        structure_index = 1
                        for structure in expected_aln:
                            threshold = len(structure)*0.45
                            intersection = structure.intersection(covered_pos)
                            if len(intersection) > threshold:
                                result[pseudogene]["structure"].append(str(structure_index))
                                is_overlapping_structure = True
                            structure_index += 1
                        if not is_overlapping_structure:
                            result[pseudogene]["structure"].append("None")
                        result[pseudogene]["structure"].append("predicted_intron")

                    #Next block remove the last "predicted intron"
                    if result[pseudogene]["structure"][-1] == "predicted_intron":
                         del result[pseudogene]["structure"][-1]

                    #Next block looks for intron retention or loss pattern
                    for i in range(0, len(result[pseudogene]["structure"]) - 2):

                        #If a single correspondance was found between a protein part and a CDS, it means that only one CDs from
                        #the (multi-exonic) parent gene is conserved, and it is therefore impossible to infer its type
                        #Next block checks for patterns
                        if len(result[pseudogene]["structure"]) != 1:
                            first_struct = result[pseudogene]["structure"][i]
                            second_struct = result[pseudogene]["structure"][i + 1]
                            if len(result[pseudogene]["structure"]) > 2:
                                third_struct = result[pseudogene]["structure"][i+2]
                                #If a predicted intron is between two conserved protein parts...
                                if second_struct == "predicted_intron" and first_struct != "None" and third_struct != "None":
                                    # ... and they are not the same, an intron is retained
                                    if first_struct != third_struct:
                                        retained_introns += 1
                                    # ... and they are the same, a potential large insertion is found
                                    else:
                                        result[pseudogene]["structure"][i + 1] = "potential_large_insertion"
                            #If two consecutive different conserved protein parts are found, an intron is loss
                            if (first_struct != second_struct and first_struct != "predicted_intron" and second_struct != "predicted_intron" and \
                                first_struct != "None" and second_struct != "None") or \
                            (second_struct != third_struct and second_struct != "predicted_intron" and third_struct != "predicted_intron" and \
                                    second_struct != "None" and third_struct != "None"):
                                lost_introns += 1

                    #If no intron loss or gain could be found (only low-conserved protein part), is it impossible to infer
                    #the pseudogene type
                    if retained_introns + lost_introns != 0:
                        #Porportion of retained intron is used to infer the pseudogene type.
                        retained_proportion = retained_introns/(retained_introns + lost_introns)
                        #If 60% or more introns are conserved, it is considered a duplicated pseudogene. This means that:
                        #1. P-GRe considers that intron loss are more indicative that intorn retention, because retropseudogenes coming
                        #from alternative transcripts can retain some introns.
                        #2. P-GRe considers that some duplicated pseudogene can have "some" lost introns. It is actually a rare events
                        #, but I think that splicing signal can easily be lost and induce errors in miniprot.
                        if retained_proportion > 0.6:
                            result[pseudogene]["type"] = "duplicated_pseudogene"
                        else:
                            result[pseudogene]["type"] = "processed_pseudogene"

#### MAIN
######## WRITING GFF OUTPUT

for pseudogene in result:

    ignore = False

    ############ ONE LAST FILTER
    # This is made to filter predictions that were made from the second set of proteins. As this set can be extremly
    # large, it can produced many short results. I don't think they are wrong, but its a truth the reviwers are not
    # ready to accept! (just kidding)
    if not result[pseudogene]["in_org"]:
        size = result[pseudogene]["right"] - result[pseudogene]["left"] + 1
        if size < 150:
            ignore = True

    if not ignore:
        gff_line = [result[pseudogene]["chr"], result[pseudogene]["method"], result[pseudogene]["type"], \
                    str(result[pseudogene]["left"]), str(result[pseudogene]["right"]), str(result[pseudogene]["score"]), \
                    result[pseudogene]["strand"], result[pseudogene]["frame"]]
        attributes_string = "ID=" + result[pseudogene]["ID"] + ";"
        if result[pseudogene]["type"] != "unitary_pseudogene":
            attributes_string += "parent_gene=" + result[pseudogene]["Target"] + ";"
        attributes_string += "identity=" + result[pseudogene]["Identity"] + ";"
        attributes_string += "positive=" + result[pseudogene]["Positive"] + ";"
        attributes_string += "frameshift=" + result[pseudogene]["Frameshift"]
        gff_line.append(attributes_string)
        print("\t".join(gff_line))

        for cds in result[pseudogene]["CDS"]:
            gff_line = [result[pseudogene]["chr"], result[pseudogene]["CDS"][cds]["method"], result[pseudogene]["CDS"][cds]["type"], \
                        str(result[pseudogene]["CDS"][cds]["left"]), str(result[pseudogene]["CDS"][cds]["right"]), str(result[pseudogene]["CDS"][cds]["score"]), \
                        result[pseudogene]["CDS"][cds]["strand"], result[pseudogene]["CDS"][cds]["frame"]]
            attributes_string = "Parent=" + result[pseudogene]["CDS"][cds]["Parent"] + ";"
            if result[pseudogene]["type"] != "unitary_pseudogene":
                attributes_string += "parent_gene=" + result[pseudogene]["Target"] + ";"
            attributes_string += "identity=" + result[pseudogene]["CDS"][cds]["Identity"] + ";"
            attributes_string += "frameshift=" + result[pseudogene]["CDS"][cds]["Frameshift"]
            gff_line.append(attributes_string)
            print("\t".join(gff_line))

