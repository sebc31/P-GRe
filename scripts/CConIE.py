# -*- coding: utf-8 -*-

'''Chimer Construction from Inputed Exons (CConIE) is a sub-programm invocated by P-GRe to construct all the possible
chimeric sequences between two predicted overlapping exons in different ORF.'''

import collections
import math

def chimer_construct(interval1,interval2,sequences,log):
    # From the exons coordinates, the two possibles "half"-proteins are constructed and their lengths are computed.
    chim_dic=collections.defaultdict(dict)
    seqframe1=sequences[interval1[0]:interval1[1]+1].translate()
    len_frame1=len(seqframe1)
    seqframe2=sequences[interval2[0]:interval2[1]+1].translate()
    len_frame2=len(seqframe2)
    # The (theorical) protein that should be constructed if no frame-shift were found has its length is computed. The 
    # difference between this theorical protein and the sum of the lengths of the two "half"-proteins
    # equals the length of the overlap between these two "half"-protein, and the frame shift is expected to take place
    # in this overlap, which is later called "conflict". This "conflict" subsequence is therefore the part of the
    # protein sequence where it's hard to determine when one ORF stop and another starts.
    expected_size=len(sequences[interval1[0]:interval2[1]+1].translate())
    overlap_size=len_frame1+len_frame2-expected_size
    dsh_str = ""
    for i in range(0,expected_size-len_frame1):
        dsh_str+="-"
    log.write("[FRAME ALPHA]   "+str(seqframe1)+dsh_str+"\n")
    dsh_str = ""
    for i in range(0,expected_size-len_frame2):
        dsh_str+="-"
    log.write(" [FRAME BETA]   "+dsh_str+str(seqframe2)+"\n")
    dsh_str = ""
    for i in range(0,overlap_size):
        dsh_str+="?"
    log.write("   [CONFLICT]   "+str(seqframe1[:-overlap_size])+dsh_str+str(seqframe2[overlap_size:])+"\nSolving conflict...")
    # All the possible sequences with a frame-shift in the conflict zone ("chimers") are computed here
    for i in range(0,overlap_size+1):
        chimer_name="chimer"+str(i+1)
        chim_dic[chimer_name]['sequence']=str(seqframe1[:len_frame1-i])+str(seqframe2[overlap_size-i:])
        chim_dic[chimer_name]['interval1']=[interval1[0],interval1[1]-i*3]
        chim_dic[chimer_name]['interval2']=[interval2[0]+(overlap_size-i)*3,interval2[1]]
        chim_dic[chimer_name]['conflict']=[len_frame1-overlap_size,len_frame1+overlap_size]
    return chim_dic

def chimdic_to_fasta(dic,wd):
    # This function is for writing a FASTA file, which is used for BLASTp and find the best chimeric protein
    file=open(wd+"/tmp/chim_temp.fasta", "w")
    for chimer in dic:
        file.write('>'+chimer+'\n'+dic[chimer]['sequence']+'\n')

def bestChim(file):
    # This function retrieve the best chimer from a blast file. The best chimer is the one with the lowest E-value
    lowest_eval=math.inf
    best_chim=""
    with open(file) as chim_res:
        for line in chim_res:
            line=line.split("\t")
            if float(line[10])<lowest_eval:
                lowest_eval=float(line[10])
                best_chim=line[0]
    return best_chim




