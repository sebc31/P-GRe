# -*- coding: utf-8 -*-

'''BLAST And Separate Next Inputs (BASNI) is used to...'''

import collections
import re
import argparse
import Bio.SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

########################################################################################################################
#################################################ARUGMENT PARSING#######################################################
########################################################################################################################

parser = argparse.ArgumentParser(description="Filter BLAST results")
parser.add_argument("-f","--fasta", help="FASTA files of the genome to annotate", required=False, type=str)
parser.add_argument("-r","--results", help="BLAST results obtained by blasting the proteome against the masked genome", required=True, type=str)
parser.add_argument("-w","--wd", help="Working directory for pipeline integration", required=False, type=str)
parser.add_argument("-v","--verbose", help="Verbose ON/OFF", required=True, type=bool)
args = parser.parse_args()

########################################################################################################################
####################################################FUNCTIONS###########################################################
########################################################################################################################

def Overlap(intervals1,intervals2):
    if (intervals1[0]<intervals2[0] and intervals1[1]>intervals2[0]) or \
       (intervals1[0]>intervals2[0] and intervals1[0]<intervals2[1]) or \
       (intervals1[0]<intervals2[0] and intervals1[1]>intervals2[1]) or \
       (intervals1[0]>intervals2[0] and intervals1[1]<intervals2[1]):
        return True
    else:return False

def extendToMinMax(intervals1,intervals2):
    intervals1[0]=min(intervals1[0],intervals2[0])
    intervals1[1]=max(intervals1[1],intervals2[1])
    return intervals1

def updateEval(currentProtFiltered,index):
    #best_eval = float(currentProtFiltered[index][1]['eval'])
    best_eval = currentProtFiltered[index][1]['spawn'][1]-currentProtFiltered[index][1]['spawn'][0]+1
    best_hit = currentProtFiltered[index][0]
    current_conflict = currentProtFiltered[index][1]['spawn'][1]
    return best_eval, best_hit, current_conflict

########################################################################################################################
########################################################TBLASTN#########################################################
########################################################################################################################

#After masking the protein-coding gene loci (i.e. replacing nucletodies with "N" in the sequence), the protein sequences
#are blasted against the genome. This steps "forces" the proteins to align on part of the genomes that aren't protein
#coding gene. Max_intron_length is overly exagerated to ensure that many hits from the same protein get the same
#E-value on a single loci, making it easier to find the left and right position of the pseudogene. Gap extension penalty
#are also set to the max, as some short gaps are expected in pseudogenes (indel events), extended gaps corresponding
#to pseudo-intron. This also allows to avoid (some) frame-shift as it can separates a large hit overlapping a frame-shit
#into two smaller hits, precilsely separated at the frame-shift location. Undetected frame-shift are later corrected
#by P-GRe.

''' Old code for blasting via Python:
cmd_tblastn = NcbitblastnCommandline(cmd=args.binary, query=args.fasta, out="./blast_res.blast", outfmt=6, db=db, evalue=0.01, word_size=3, gapextend=2, max_intron_length=100000)
cmd_tblastn()'''

########################################################################################################################
#################################################RESULTS PARSING########################################################
########################################################################################################################

#currentProt will be used to stock the results in a dictionnary with format:
# churrentProt[id1]
# [chr]=chromosome
# [strand]=+|-
# [eval]=E-value of the hits
# [hits]=positions of the hits on the scaffold
# [spawn]=[start of first hit, end of last hit]
# [missing]=[amino acid upstream of the first hit, last amino acid that matched]
# [total_align_length]=total length (in amino acids) of the hits
# [identity]=[% of identity of hit1, % of identity of hit2, % of identity of hit3,...]
# For a protein "A" that aligns with two loci in the genome, keys (id) of this dictionnary are generated like:
# A.1 for the first loci, A.2 for the second loci.
currentProt=collections.defaultdict(dict)
index_dic={}  # To generate unique identifier
current_query=""
current_eval=None

with open(args.results) as blast:
    for line in blast:
        line=line.split('\t')
        left = min(int(line[8]), int(line[9]))
        right = max(int(line[8]), int(line[9]))
        #If a protein found on a new line is the same as the last protein, on the same chromosome, with the same E-value
        #, its caracteristics are added to the last protein sub-dictionnary
        if line[0]==current_query and current_chrom==line[1] and \
                (Overlap(currentProt[id]['spawn'],[left-200,right+200])
                or (line[10]==current_eval and Overlap(currentProt[id]['spawn'],[left-10000,right+10000]))):
                #The same E-value should be given to every hit of a protein in a loci due to BLAST parameters.
                #However, frame-shift can cause small overlap that prevent this. To avoid this, hits from the same
                #proteins that overlaps are merged together. If E-value is equal to 0, every hits from the same protein
                #in a chromosome will merge together. To avoid this, hits from the same protein having an E-value=0.0
                #are simply merged together if they are separated by a distance < 10 kb
            currentProt[id]['identity'].append([float(line[2]),int(line[3])])
            currentProt[id]['hits'].append([left,right])
            currentProt[id]['total_align_length'] += int(line[3])
            currentProt[id]['spawn']=extendToMinMax(currentProt[id]['spawn'],[left,right])
            currentProt[id]['missing']=extendToMinMax(currentProt[id]['missing'],[int(line[6]),int(line[7])])
        #Otherwise, a new sub-dictonnary is created
        else:
            #The next 3 lines generate a unique identifier
            if line[0] not in index_dic:index_dic[line[0]] = 1
            else:index_dic[line[0]] += 1
            id = line[0] + "." + str(index_dic[line[0]])
            if int(line[8])>int(line[9]):strand="-"
            else:strand="+"
            currentProt[id]['chr'] = line[1]
            currentProt[id]['strand'] = strand
            currentProt[id]['hits'] = [[left,right]]
            currentProt[id]['eval']=line[10]
            currentProt[id]['spawn']=[left,right]
            currentProt[id]['missing']=[int(line[6]),int(line[7])]
            currentProt[id]['total_align_length'] = int(line[3])
            currentProt[id]['identity'] = [[float(line[2]),int(line[3])]]
            current_query=line[0]
            current_eval=line[10]
            current_chrom=line[1]
            current_left=left
            current_right=right
if args.verbose: print("               >> Parsing:       OK <<")

########################################################################################################################
###############################################RESULTS FILTERING########################################################
########################################################################################################################

#The next block retrieve the multi-hits having a total length > threshold and at least one hit with identity >
# threshold. the results are then sorted by chromsome, left position of the first hit and right position of the
# last hit ("spawn")
currentProtFiltered=collections.defaultdict(dict)
#Next lines are used to quickly modify the threshold equation, if needed (see next commentary block)
left_id=75
right_id=20
left_align=20
right_align=100

for multi_hits in currentProt:
    for identity in currentProt[multi_hits]['identity']:
        #Rather than having strict threshold for alignment length and identity, P-GRe uses a fonction so that long
        #alignments need less identity than short alignments. A minimum identity is still expected. The equation of
        # the identity threshold is:
        # id_thr = (left_align-align_len)/((left_id-right_id)/(right_align-left_align))+left_id
        #where left_id is the minimal identity expected for short alignment of size left_align
        #and right_id id the minimal identity expected for long alignment of size right_align
        #The equation is represented here:
        #   left_align(ex:15 aa)           ----------------------> right_align(ex:100 aa)
        #   has an expected identity of...                         has an expected identity of...
        #   left_id(ex:75.0)               ----------------------> right_id(ex:20.0)
        id_thresh=max((left_align-identity[1])/((left_id-right_id)/(right_align-left_align))\
                      +left_id,right_id)
        #A minimum total length is still expected
        if identity[0]>=id_thresh:
            currentProtFiltered[multi_hits]=currentProt[multi_hits].copy()
            break
del currentProt

currentProtFiltered=sorted(currentProtFiltered.items(),key=lambda x: (x[1]['chr'],x[1]['spawn'][0],x[1]['spawn'][1]))
#Note that the sorted() function transforms the dictionnary structure into a list of tuple like:
#[(id1,currentProt[id1]),(id2,currentProt[id2]),...]. I chose not to reconstruct the dictionary structure.
if args.verbose: print("               >> Filtering:     OK <<")

#The next block iterates over the results to find overlaping results and keep the one with the lowest E-value.
#Contrary to other pseudogene prediction softwares, P-Gre does not condier results that are separated by less than
#XX bp to be overlapping. Considering all the hits of a protein on a loci permits to avoid this, moreover this might
#make it impossible to detect chimeric retropseudogenes, as they have more than one parent gene and their hits needs
#to be treated as separate results
parents=set()
best_hit=""
#best_eval=1.0
best_eval=0
current_conflict=0
for i in range(1,len(currentProtFiltered)):
    #While a results overlap the previous one, compare the E-values and stocks the best result in the 'best_hit'
    #variable.
    if (currentProtFiltered[i-1][1]['spawn'][1]>=currentProtFiltered[i][1]['spawn'][0] or\
            current_conflict>=currentProtFiltered[i][1]['spawn'][0]) and currentProtFiltered[i-1][1]['chr']\
            ==currentProtFiltered[i][1]['chr']:
        #if float(currentProtFiltered[i-1][1]['eval'])<best_eval:
        if currentProtFiltered[i-1][1]['spawn'][1]-currentProtFiltered[i-1][1]['spawn'][0]+1>best_eval:
            best_eval, best_hit, current_conflict=updateEval(currentProtFiltered,i-1)
        #if float(currentProtFiltered[i][1]['eval'])<best_eval:
        if currentProtFiltered[i][1]['spawn'][1] - currentProtFiltered[i][1]['spawn'][0] + 1 > best_eval:
            best_eval, best_hit, current_conflict=updateEval(currentProtFiltered,i)
    else:
        #If the current result doesn't overlap the previous result and the previous result was being "compared" with
        #other overlapping results, stocks the best_hit from this previous comparison as a potential parent gene.
        if best_hit!="":
            parents.add(best_hit)
            best_hit=""
            #best_eval=1.0
            best_eval=0
            current_conflict=0
        #Otherwise, if the current result doesn't overlap the previous result and the previous result didn't olverap
        #any previous results, the previous hit is stocked as a potential parent gene, without being compared to any
        # other hits.
        else:
            parents.add(currentProtFiltered[i-1][0])
    #The next two line stocks the result at the last iteration
    if i==len(currentProtFiltered)-1:
        if best_hit!="":parents.add(best_hit)
if len(currentProtFiltered)==1: parents.add(currentProtFiltered[0][0])

if args.verbose: print("               >> Overlaping:    OK <<")
print(parents)
########################################################################################################################
###########################################RESULTS (~BED) WRITING#######################################################
########################################################################################################################

#For every results in the filtered currentProt dictionnary (transformed into a list of tuples), retrieve only the ones
#which have been noted as parent in the parents list.
pseudogeneDict=collections.defaultdict(dict)
chrDic=collections.defaultdict(dict)  # To quickly retrieve sequence from scaffold
pseudogene_count=1
for tuple in currentProtFiltered:
    if tuple[0] in parents:
        tuple[1]['hits']=sorted(tuple[1]['hits'],key=lambda x: (x[0],x[1]))
        zero=""
        for i in range(0,6-len(str(pseudogene_count))):
            zero+="0"
        id_pg="pseudogene"+zero+str(pseudogene_count)
        pseudogeneDict[id_pg]['parent']=re.sub("\.[0-9]*$","",tuple[0])
        pseudogeneDict[id_pg]['chr'] = tuple[1]['chr']
        ###The next four lines are only used to accelerate the FASTA file writing
        if tuple[1]['chr'] not in chrDic:
            chrDic[tuple[1]['chr']]=[]
        chrDic[tuple[1]['chr']].append(tuple[1]['spawn'].copy())
        chrDic[tuple[1]['chr']][-1].append(id_pg)
        ######
        pseudogeneDict[id_pg]['strand'] = tuple[1]['strand']
        pseudogeneDict[id_pg]['missing_N'] = tuple[1]['missing'][0]
        pseudogeneDict[id_pg]['missing_C'] = tuple[1]['missing'][1]
        pseudogeneDict[id_pg]['spawn'] = tuple[1]['spawn']
        pseudogeneDict[id_pg]['hits'] = tuple[1]['hits'].copy()
        pseudogeneDict[id_pg]['real_position'] = tuple[1]['spawn'][0]
        ###Convert chromosome-centred coordinates to the FASTA-centred coordinates (see FASTA later)
        pseudogeneDict[id_pg]['spawn'][1]=pseudogeneDict[id_pg]['spawn'][1]-pseudogeneDict[id_pg]['spawn'][0]+151
        if pseudogeneDict[id_pg]['spawn'][0]>150:
            pseudogeneDict[id_pg]['spawn'][0]=151
        for hits in pseudogeneDict[id_pg]['hits']:
            hits[0]=hits[0]-pseudogeneDict[id_pg]['real_position']+pseudogeneDict[id_pg]['spawn'][0]
            hits[1]=hits[1]-pseudogeneDict[id_pg]['real_position']+pseudogeneDict[id_pg]['spawn'][0]
        ######
        pseudogene_count+=1
del currentProtFiltered

#Write the output in a TSV format:
#Chromosome/scaffold   pseudogene_id   5'-position on the scaffold   strand   parent_id   missing_N   missing_C   hits
file = open(args.wd+"/tmp/tblastn_results.tempfile.tsv", "w")
for pseudogene in pseudogeneDict:
    file.write(pseudogeneDict[pseudogene]['chr']+"\t"+pseudogene+"\t"+\
               str(pseudogeneDict[pseudogene]['real_position'])+"\t"+\
               pseudogeneDict[pseudogene]['strand']+"\t"+pseudogeneDict[pseudogene]['parent']+"\t"+\
               str(pseudogeneDict[pseudogene]['missing_N'])+"\t"+str(pseudogeneDict[pseudogene]['missing_C'])+"\t")
    hit_string=""
    for hits in pseudogeneDict[pseudogene]['hits']:
        hit_string+=str(hits[0])+":"+str(hits[1])+";"
    hit_string=hit_string.strip(";")
    file.write(hit_string+"\n")
if args.verbose: print("               >> Writing TSV:   OK "+args.wd+"/tmp/tblastn_results.tempfile.tsv <<")

########################################################################################################################
################################################FASTA WRITING###########################################################
########################################################################################################################

file = open(args.wd+"/tmp/first_pg_set.tempfile.fasta", "w")
for scaffold in list(Bio.SeqIO.parse(args.fasta,'fasta')):
    for hits in chrDic[scaffold.id]:
        file.write(">"+hits[2]+"\n")
        #To work with Python, where first character in a string has an index of 0, 1 is substracted to every
        #hits found. As Python also exclude the right value in a list, 1 is added to the right value of a hit.
        #in the end, only the left value of the hit is substracted by 1. For P-GRe N-ter and C-ter reconstruction
        #, sequences need to be flanked by at least 60 bp, in which P-GRe will try to find a strat and stop codon, thus
        #, it is necessary to substract and add ~60 bp (150 bp) upstream and downstream of the sequence.
        if hits[0]-151>0: left_start=hits[0]-151
        else: left_start=0
        file.write(str(scaffold.seq[left_start:hits[1]+150])+"\n")
if args.verbose: print("               >> FASTA writing: OK "+args.wd+"/tmp/first_pg_set.tempfile.fasta <<")
