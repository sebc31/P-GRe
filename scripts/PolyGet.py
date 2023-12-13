# -*- coding: utf-8 -*-

'''Quick tools to retrieve Poly-A tail and categorize pseudogene.'''

# Find in which category each kind of PG belongs depending on their structure ('structure_cons' attribute)
# The attribute 'structure_cons' is changed into a 3-elements list, where the first element is the type of
# copy ("full" copy, noted "Copy", or partial, noted "Fragment") and the second element is the type of pseudogene
# ("Pseudogene" or "Retropseudogene"). To classify the pseudogenes, only the "high quality evidences" are taken
# into account (see next blocks of code). When no such evidences are found, and depending on the cases, other
# categories can exists. The third element is the structure from which the pseudogene originates.
#### Case A. Pseudogenes consisting of only one CDS

def catergorizePG(pseudogeneDic, pseudogene, parentDic, interval_list, parent_prot):
    #### Case A. Pseudogenes consisting of only one CDS
    if len(pseudogeneDic[pseudogene]['structure_cons']) == 0:
        ######## Case A - A. Parent also consists of only one CDS
        if parentDic[pseudogeneDic[pseudogene]['parent']]['CDS_nb'] == 1:
            ############ Case A - A - A. Length of the reconstructed pseudogene sequence is close to the parent gene
            ############ sequence (i.e. >= 80% of the length of the parent sequence)
            if (interval_list[-1][-1] - interval_list[0][0]) / 3 >= len(parent_prot) * 0.80:
                categories = ["Copy", "Unknown", pseudogeneDic[pseudogene]['parent']]
            ############ Case A - A - B. Length of the reconstructed pseudogene sequence is shorter than the parent
            ############ gene (i.e. < 80% of the length of the parent sequence)
            else:
                categories = ["Fragment or degraded copy", "Unknown", pseudogeneDic[pseudogene]['parent']]
        ######## Case A - B. Parent have more than one CDS
        else:
            categories = ["Fragment", "Unknown", pseudogeneDic[pseudogene]['parent']]
    #### Case B. Pseudogenes consisting of more than one CDS
    else:
        ######## Case B - A. Parent has only one CDS, => potential insertion into the CDS
        if parentDic[pseudogeneDic[pseudogene]['parent']]['CDS_nb'] == 1:
            ############ Case B - A - A. % of aligned amino acids between the two sequences is >= than 80%
            if pseudogeneDic[pseudogene]['structure_cons'][0] >= 0.80:
                categories = ["Copy", "Unknown", pseudogeneDic[pseudogene]['parent']]
            ############ Case B - A - B. % of aligned amino acides between the two sequences is < 80%
            else:
                categories = ["Fragment or degraded copy", "Unknown", pseudogeneDic[pseudogene]['parent']]
        ######## Case B - B. Parent has more than one CDS
        else:
            categories = ["", "", pseudogeneDic[pseudogene]['parent']]
            intron = False
            intron_loss = False
            high_quality_cds_found = False
            ######## P-GRe will look for intron or intron loss that are surrounded by "high quality evidence", i.e.
            ######## part of the pseudogenes where sequence aligns with at least 80% of corresponding parent sequence
            ######## CDS
            for i in range(1, len(pseudogeneDic[pseudogene]['structure_cons']) - 1):
                border_left = pseudogeneDic[pseudogene]['structure_cons'][i - 1]
                struct = pseudogeneDic[pseudogene]['structure_cons'][i]
                border_right = pseudogeneDic[pseudogene]['structure_cons'][i + 1]
                if str(border_left).replace(".", "").isnumeric() and \
                        str(border_right).replace(".", "").isnumeric() and \
                        border_left >= 0.80 and border_right >= 0.80:
                    if struct == "INTRON":
                        intron = True
                    elif struct == "NO INTRON":
                        intron_loss = True
                ######## Also checks if a single high-quality CDS is found
                if (str(border_left).replace(".", "").isnumeric() and \
                    border_left >= 0.80) or \
                        (str(border_right).replace(".", "").isnumeric() and \
                         border_right >= 0.80):
                    high_quality_cds_found = True
            ############ Case B - B - 1A. First and last CDS of the parent sequence have >= 80% of its a.a. aligned
            ############ with the pseudogene sequence
            if str(pseudogeneDic[pseudogene]['structure_cons'][0]).replace(".", "").isnumeric() and \
               str(pseudogeneDic[pseudogene]['structure_cons'][-1]).replace(".","").isnumeric():
                if pseudogeneDic[pseudogene]['structure_cons'][0] >= 0.80 \
                   and pseudogeneDic[pseudogene]['structure_cons'][-1] >= 0.80:
                    categories[0] = 'Copy'
            ############ Case B - B - 1B. First or last CDS of the parent sequence have <80% of its a.a. aligned
            ############ with the pseudogene sequence
            if categories[0] != "Copy":
                ################# Case B - B - 1B - A. At least one of the CDS of the parent sequence have >= 80% of
                ################# its a.a. aligned with pseudogene sequence
                if high_quality_cds_found:
                    categories[0] = 'Fragment'
                ################ Case B - B - B - 1B. No CDS of the parent sequence have its CDS where its a.a. align
                ################ with at least 80% of the pseudogene sequence
                else:
                    categories[0] = 'Fragment or degraded copy'
            ############ Case B - B - 2A. Intron AND intron loss found between two "high quality" CDS
            if intron and intron_loss:
                categories[1] = '(Iso)retropseudogene'
            ############ Case B - B - 2B. Intron AND NOT intron loss found between two "high quality" CDS
            elif intron and not intron_loss:
                categories[1] = 'Duplicated pseudogene'
            ############ Case B - B - 2C. NO intron AND intron loss found between two "high quality" CDS
            elif not intron and intron_loss:
                categories[1] = 'Retropseudogene'
            ############ Case B - B - 2D. No intron AND NO intron loss found between two "high quality" CDS
            else:
                categories[1] = 'Unknown'
    return categories


def containsPolyA(seq, length = 20, expected = 0.75, unexpected = 0.35, max_reach=500):
    # If a seq is shorter than 20 bp, it is considered insuficient to check for polyA site.
    if len(seq) < length:
        return "Unknown"
    # PolyGet will look for a (by default) 20 bp region downstream of the pseudogene that contains at least 75%
    # of A. If such region is found, a polyA site is considered found. If no 20-bp region containing at least 30%
    # of A is found, it is considered that no polyA site is present, even in a degenerated form, and pseudogene can
    # thus be considered a (non-retro)pseudogene. This categorization is based on the absence of PolyA site and is
    # very approximate, the "unexpected" threshold must therefore be very strict.
    else:
        absent = True
        limit = min(max_reach, len(seq)-length-1)
        for i in range(0, limit):
            adenylation_rate = seq[i:i+length].upper().count('A') / length
            if adenylation_rate >= expected:
                return "Found"
            if adenylation_rate >= unexpected:
                absent = False
    if absent: return 'Not found'
    else : return "Unkown"
