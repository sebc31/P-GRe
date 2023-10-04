# -*- coding: utf-8 -*-

'''Tools for Analyses from Gff Length for Intron Alignment (TAGLIA) retrieves the length of every peptides coded by
every CDS in a GFF file, or for one specific protein. Can also performs an Lindley-like process to find introns
positions on alignment.'''

import re
import math
import collections

#pepLength[Full/Fast] uses a GFF file to retrieve the length of every amino acid sequences coded by every CDS of a
# gene in a dictionary format like: pep_len[protein1][pep1]=length_of_the_aa_seq, pep_len[protein1][pep2], ...
def pepLengthFull(gff):
    pep_len=collections.defaultdict(dict)
    with open(gff) as gff:
        for line in gff:
            CDS=re.search("\tCDS\t([0-9]*)\t([0-9]*).*Parent=([^;]*);",line)
            if CDS:
                if CDS.group(3) not in pep_len:
                    i=1
                    pep_len[CDS.group(3)]["pep"+str(i)]=int(CDS.group(2))-int(CDS.group(1))+1
                else:
                    i+=1
                    pep_len[CDS.group(3)]["pep"+str(i)]=int(CDS.group(2))-int(CDS.group(1))+1

    # After retrieving the length of every CDS, %3 is used to see if the length is divisible by 3. If not, the
    # remainder is added to the next CDS, as it will be translated with the next CDS after splicing.
    for prot in pep_len:
        next_pep_follows_intron=False
        for pep in pep_len[prot]:
            index=list(pep_len[prot]).index(pep)+1
            length=pep_len[prot][pep]
            if next_pep_follows_intron:
                length+=divisible
            divisible = length % 3
            if divisible==0:
                length=length/3
                next_pep_follows_intron=False
            else :
                length=math.floor(length/3)
                next_pep_follows_intron=True
            if index==len(pep_len[prot]):  # Remove stop codon
                length-=1
            pep_len[prot][pep]=int(length)
    return pep_len

#The max_error_expected allows for the Lindley process to start looking for introns a few amino acids before the
# expected positions.
def pepLengthFastRun(gff,prot_to_look,max_error_expected=0):
    attributes_to_look = ['Parent', 'Name']
    pep_len=[]
    pep_len_exact=[]
    strand=""
    with open(gff) as gff:
        for line in gff:
            if line[0]!="#":
                CDS=line.split("\t")
                if CDS[2]=="CDS":
                    attributs_dic={}
                    #Next XX lines are used to convert attributes into a dictionnary
                    attributes=CDS[8].split(";")
                    for attribute in attributes:
                        if "=" in attribute:
                            attribute=attribute.split("=")
                            attributs_dic[attribute[0]]=attribute[1]
                    for attribute_to_look in attributes_to_look:
                        if attribute_to_look in attributs_dic:
                            if attributs_dic[attribute_to_look]==prot_to_look:
                                pep_len.append(int(CDS[4])-int(CDS[3])+1)
                                if strand=="":strand=CDS[6]
                                break

    next_pep_follows_intron=False
    for i in range(0,len(pep_len)):
        length=pep_len[i]
        if next_pep_follows_intron:
            length+=divisible
        divisible = length % 3
        if divisible==0:
            length=length/3
            next_pep_follows_intron=False
        else :
            length=math.floor(length/3)
            next_pep_follows_intron=True
        pep_len[i]=int(length)
    #Accounts for strands
    if strand=="-":pep_len.reverse()
    #Remove stop codon
    pep_len[-1]-=1
    #Sums every length by the sum of its previous length to get the "intron" positions and remove the "tolerance margin"
    for i in range(1,len(pep_len)):
        pep_len[i]=pep_len[i-1]+pep_len[i]
    pep_len_exact = pep_len.copy()
    pep_len[0]-=max_error_expected
    for i in range(1,len(pep_len)):
        pep_len[i]-=max_error_expected
    return [pep_len,pep_len_exact]

#linldeyAlign perfoms a simplified Lindley process (Lindley, 1952) to remove introns from the pseudo-protein. This
# allows for a soft filtering of introns. This (modified) Lindley process uses CDS length informations to only cumulate
# local score when an intron is expected. A gap ("-") in the parent protein add 1 to local score, -1 if not.
#Local score is also caped to ten. More examples and explanations at the end of this script.
def lindleyAlign(alignement,pep_len,error, log):
    ######DEV-TOOLS
    maxLS=10  #LS cap
    ######
    ##############################
    #STEP1. Computing Local Score#
    ##############################
    pep_len_no_last=pep_len[0][0:-1]
    state_string=""
    local_score_list=[]
    splice_pos_list=[]
    pg_aligned=alignement[0][0]  # Sequence of the aligned pseudogene
    ref_aligned=alignement[0][1]  # Sequence of the aligned parent gene
    len_align=max(len(pg_aligned),len(ref_aligned))  # Alignement length
    log.write(" [PSEUDO-PROTEIN] "+pg_aligned+"\n")
    log.write(" [PARENT PROTEIN] "+ref_aligned+"\n")
    aa_count=0
    aa_inprocess_count=0
    inProcess=False
    local_score=0
    intron=False
    for i in range(0,len_align):
        splice_pos_list.append(0)
        # Count the current positions
        if ref_aligned[i]!="-":
            aa_count+=1
        # If the current position is expected to be gap, starts the process
        if aa_count in pep_len_no_last:
            splice_pos = pep_len[1][pep_len[0].index(aa_count)]
            inProcess=True
        if inProcess:
            splice_pos_list[i] = splice_pos  #For every positions, keep the position of the closet splicing site
            state_string += "I"
            #If a gap is indeed found, add 1 to the local score
            if ref_aligned[i]=="-":
                local_score=min(maxLS,local_score+1)
                intron=True
                aa_inprocess_count=0
            else:
                #Otherwise, substract 1 to the local score
                local_score=max(0,local_score-1)
                if not intron:
                    aa_inprocess_count+=1
            #When in-Lindley Process, the first amino acid with a score of 0 stops the process
            if intron and local_score==0:
                inProcess=False
                intron=False
                aa_inprocess_count=0
            #The next blocks allows for permature stops of the process in case of retropseudogene
            if aa_inprocess_count>=10+error*2:
                inProcess=False
                aa_inprocess_count=0
        else:state_string+="O"

        local_score_list.append(local_score)
        print_score=""
        for scores in local_score_list:
            if scores<10: print_score+=str(scores)
            else: print_score+="+"
    log.write("  [LINDLEY STATE] "+state_string+"\n")
    log.write("[LINDLEY-LIKE LS] "+print_score+"\n")
    ###############################
    #STEP2. Finding relevant zones#
    ###############################
    start_of_peak=[]
    intron=False
    all_intron_pos=[]
    for i in range(0,len_align):
        #If a local score is >0, it is stocked until the first in-sequence local score=maxLS is found
        if local_score_list[i]>0 and not intron:
            start_of_peak.append(i)
        else:  # If no local score = maxLS is found, considered as a missing intron
            intron=False
            if local_score_list[i]==0:
                start_of_peak=[]
            else:
                start_of_peak=[i]
        if local_score_list[i]>maxLS-1:
            if not intron:
                intron=True
                all_intron_pos.extend(start_of_peak)
            all_intron_pos.append(i)
    ##################################
    #STEP3. Converts to binary result#
    ##################################
    #Every local score is converted to 1 if associated with an intronic positions, 0 otherwise
    bin_res=[]
    for i in range(0,len_align):
        if i in all_intron_pos:
            bin_res.append(1)
        else:
            bin_res.append(0)
    ########################################################
    #STEP4. Realign amino acids that align inside an intron#
    ########################################################
    aa_count=0
    current_bin_list=[]
    misaligned_to_left=0
    misaligned_to_right=0
    for i in range(0,len_align):
        if ref_aligned[i]!="-": aa_count+=1
        if bin_res[i]:
            current_bin_list.append(i)
            if ref_aligned[i]!="-":  #Misaligned amino acid
                if aa_count<splice_pos_list[i]+1: #Check if the misaligned aa if located before or after the splicing site
                    misaligned_to_left+=1
                else:
                    misaligned_to_right+=1
        elif len(current_bin_list)>0:
            for j in range(0,misaligned_to_left):
                bin_res[current_bin_list[j]]=0
            for j in range(1,misaligned_to_right+1):
                bin_res[current_bin_list[-j]]=0
            current_bin_list=[]
            misaligned_to_left=0
            misaligned_to_right=0
    ############################################
    #STEP5. Construction of the pseudo-proteins#
    ############################################
    pseudo_prot=""
    log.write("    [EXON/INTRON] ")
    for i in range(0,len(bin_res)):
        if bin_res[i]:
            if pg_aligned[i]!="-":
                log.write("-")
                pseudo_prot+="-"
            else: log.write("_")
        else:
            if pg_aligned[i]!="-":
                log.write(pg_aligned[i])
                pseudo_prot+=pg_aligned[i]
            else: log.write("_")
    log.write('\n\n')
    return pseudo_prot

#Check[First|Last]Exon check if the pseudoprotein has the first and last exons
def checkFirstExon(alignement,pep_len):
    if len(pep_len)==1:return True
    missingFirstExon=False
    pg_aligned=alignement[0][0]  # Sequence of the aligned pseudogene
    ref_aligned=alignement[0][1]  # Sequence of the aligned parent gene
    len_align=max(len(pg_aligned),len(ref_aligned))  # Alignement length
    aa_count=1  #pseudocount
    gap_count=1  #pseudocount
    for i in range(0,len_align):
        if aa_count==pep_len[0]:break
        if ref_aligned[i]!="-":aa_count+=1
        if pg_aligned[i]=="-":gap_count+=1
    if (gap_count/aa_count)>0.95:
        missingFirstExon=True
    return missingFirstExon

def checkLastExon(alignement,pep_len):
    if len(pep_len)==1:return True
    missingLastExon=False
    pg_aligned=alignement[0][0]  # Sequence of the aligned pseudogene
    ref_aligned=alignement[0][1]  # Sequence of the aligned parent gene
    len_align=max(len(pg_aligned),len(ref_aligned))  # Alignement length
    aa_count=1  #pseudocount
    j=0
    gap_count=1  #pseudocount
    while aa_count<pep_len[-2]: # This compute the position of the last exon
        if ref_aligned[j]!="-":aa_count+=1
        j+=1
    aa_count=0
    for i in range(j,len_align):
        if ref_aligned[i]!="-":aa_count+=1
        if pg_aligned[i]=="-":gap_count+=1
    if (gap_count/aa_count)>0.95:
        missingLastExon=True
    return missingLastExon

#alignToInterv transform an aligned protein sequence and a dictionary in the format
#dict[first_amino_acid][start_of_codon_position,stop_of_codon_position] to intervals.
def alignToInterv(pseudo_prot,aa_dic):
    new_start = 0
    intervals = []
    for i in range(0, len(pseudo_prot)):
        if pseudo_prot[i] != "-":
            if new_start == 0 or (aa_dic[str(i)]['codon_pos'][0] != aa_dic[str(i - 1)]['codon_pos'][1] \
                                  and i > 0):
                if new_start > 0:
                    new_end = aa_dic[str(i - 1)]['codon_pos'][1] - 1
                    intervals.append([new_start, new_end])
                new_start = aa_dic[str(i)]['codon_pos'][0]
        elif new_start > 0:
            new_end = aa_dic[str(i - 1)]['codon_pos'][1] - 1
            intervals.append([new_start, new_end])
            new_start = 0
        if i == len(pseudo_prot) - 1 and new_start > 0:
            new_end = aa_dic[str(i)]['codon_pos'][1] - 1
            intervals.append([new_start, new_end])
    return intervals

#######################################################################################################################
##############################MODIFIED LINDLEY PROCESS FOR INTRON REMOVAL##############################################
#######################################################################################################################

'''The modified Lindley-like process implemented in this script is used to find pseudo-introns of the reconstructed
pseudo-proteins in a pseudogene from pairewise alignement of the pseudo-proteins and the protein coded by the parent
gene of the pseudogene. To achieve this, the expected positions of the gaps are pre-computed with the pepLengthFastRun.
In short, a gap is expected after every last amino acid coded by each CDS of the parent protein. For more flexibility,
a margin of error is added, meaning that the gap can occurs two (default for P-Gre, zero for TAGLIA) amino acids before
the expected gap position. When one of the expected gaps position if found, the modified Lindley Process starts. The
examples below shows the start of the Lindley Process when a gap is expected at the 10th position (and finish at the
20th position, but this doesn't affect the start of the Lindley Process):

 1           5              10             15             20             25             30
 M  A  F  S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  A  F  S  S  Y  S  L  V  -  -  -  -  -  -  -  -  -  -  L  K  K  G  S  D  A  S  L  L  K
                      8
                      ^
              Lindley Process start
         
           
 1           5              10             15             20             25             30
 M  A  F  S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  A  F  S  S  Y  S  -  -  -  -  -  -  -  -  -  -  -  -  -  K  K  G  S  D  A  S  L  L  K
                      8
                      ^
              Lindley Process start


Note in the last example that the Lindley Process starts at the 8th amino acid, meaning that gap aren't counted
as valid positions.

 1   |   |   4              9              14             19             24             29
 M  A| F |S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  A| - |S  S  Y  S  L  V  -  -  -  -  -  -  -  -  -  -  L  K  K  G  S  D  A  S  L  L  K
     |   |               8
                         ^
                 Lindley Process start

Once the Lindley Process has started, if no gap is found in the next ten positions, it is prematurely stops, no
intron is found which can be a marker of a retropseudogene. If a gap is found, the Lindley Process give a local score
to an amino acid which is equal to the score of the previous amino acid + 1 if a gap is found, and the maximum between
0 and the score of the previous amino acid - 1 if not. Note that in this modified Lindley Process, the process can
starts and stops many times (if the score reach 0, it is stopped until the next expected gap is found), meaning that
a negative local score is not expected. The local score is also caped to 10. The example belows show the local score
of every amino acids when a gap is expected at the 10th positions.

 1           5              10             15             20             25             30  <- positions
 M  A  F  S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  A  F  S  S  Y  S  L  V  -  -  -  -  -  -  -  -  -  -  L  K  K  G  S  D  A  S  L  L  K
 0  0  0  0  0  0  0  0  0  1  2  3  4  5  5  5  5  5  5  4  3  2  1  0  0  0  0  0  0  0  <- local scores
                      ^                                               ^
              Lindley Process start                           Lindley Process stop
 
 1           5              10             15             20             25             30  <- positions
 M  A  F  S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  -  -  -  -  -  S  L  V  -  -  -  -  -  H  L  -  -  -  L  -  K  G  S  D  -  S  L  -  L
 0  0  0  0  0  0  0  0  0  1  2  3  4  5  4  3  4  5  5  4  5  4  3  2  1  2  1  0  0  0  <- local scores
                      ^                                                           ^
              Lindley Process start                                       Lindley Process stop
              
Finally, start of intron is defined as the first position with a local score >0 in a sequence of local scores
without 0 and a maximum local score of 5. End of the intron is defined as the last positions with a local score = 10
in the same sequence, as shown in the examples below :

 1           5              10             15             20             25             30  <- positions
 M  A  F  S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  A  F  S  S  Y  S  L  V  -  -  -  -  -  -  -  -  -  -  L  K  K  G  S  D  A  S  L  L  K
 0  0  0  0  0  0  0  0  0  1  2  3  4  5  5  5  5  5  5  4  3  2  1  0  0  0  0  0  0  0  <- local scores
                      ^     ^                          ^              ^
                      | Intron start               Intron stop        |
              Lindley Process start                           Lindley Process stop


 1           5              10             15             20             25             30  <- positions
 M  A  F  S  S  L  I  L  V  F  A  S  S  *  H  L  A  Y  I  L  *  K  G  C  D  A  S  L  L  L
 M  -  -  -  -  -  S  L  V  -  -  -  -  -  H  L  -  -  -  L  -  K  G  S  D  -  S  L  -  L
 0  0  0  0  0  0  0  0  0  1  2  3  4  5  4  3  4  5  5  4  5  4  3  2  1  2  1  0  0  0  <- local scores
                      ^     ^                                ^                    ^
                      | Intron start                     Intron stop              |
              Lindley Process start                                       Lindley Process stop
'''
