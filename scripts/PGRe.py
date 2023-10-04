# -*- coding: utf-8 -*-
import math

#P-GRe scripts
import CConIE
import VITo
import argparse
import TAGLIA

#Python utilities
import collections
import warnings
import re
import sys

#Biopython
import Bio.SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Emboss.Applications import StretcherCommandline
from Bio.Align import substitution_matrices
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio import BiopythonDeprecationWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonDeprecationWarning)
    from Bio import pairwise2

########################################################################################################################
#################################################ARUGMENT PARSING#######################################################
########################################################################################################################

parser = argparse.ArgumentParser(description="P-GRe pipeline for pseudogene reconstruction")
parser.add_argument("-f","--fasta", help="Proteins in multi-FASTA format", required=True, type=str)
parser.add_argument("-g","--gff", help="Genome's corresponding GFF3 file", required=True, type=str)
parser.add_argument("-s","--scaffold", help="Genome's FASTA file", required=True, type=str)
parser.add_argument("-w","--wd", help="Working directory for pipeline integration", required=True, type=str)
parser.add_argument("-b","--binaries", help="Path to binaries", required=False, default="../bin")
parser.add_argument("-v","--verbose", help="Verbose ON/OFF", required=True, type=bool)
args = parser.parse_args()

if args.verbose: print("pg generation:")

########################################################################################################################
####################################################FUNCTIONS###########################################################
########################################################################################################################

def createFile(to_write):
    type_of_seq = re.sub("\.f.a$", "", to_write)
    file = open(args.wd+"/res/pseudogenes." + to_write, "w")
    if args.verbose: print("                                 OK " + args.wd + "/res/pseudogenes."+to_write+" <<")
    for tuple in pseudogeneDic:
        file.write(">" + tuple[1]['id'] + "\n")
        file.write(insert_newlines(tuple[1][type_of_seq]) + "\n")

def gap_function(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    #The first BLAST assure that aligning the pseudo-protein with its parent gene protein will produce "blocks", thus
    #the exagerated gap opening penality. For some reason, overly exagerating this penality tends to produce more gaps,
    #and the best results where obtained with a gap opening penalty of 5.
    else: return -5

def insert_newlines(string, every=60):
    string=str(string)
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

def invert_coord(list_of_interv,right_lim):
    for hits in list_of_interv:
        tmp=abs(hits[1]-right_lim)
        hits[1]=abs(hits[0]-right_lim)
        hits[0]=tmp
    list_of_interv = sorted(list_of_interv, key=lambda x: (x[0], x[1]))
    return list_of_interv

def title(title,endline='\n',size=40,char="-"):
    global log
    nice_dash=""
    for i in range(0,int(size-len(title)/2)):
        nice_dash+=char
    log.write(nice_dash+title+nice_dash+endline)

########################################################################################################################
###############################################RETRIEVE BLAST RESULTS###################################################
########################################################################################################################

#Retrieve the results from the BLAST analysis (filtered and parsed) and convert them into a dictionary with format:
#pseudogeneDic[pseudogene1]
#['parent']=parent gene
#['chromosome']=chromosome
#['strand']=+|-
#['intervals']=[[start of hit1, end of hit1], [start of hit2, end of hit2], ...]
#['missing_N']=number of amino acids before the start of the fist hit
#['missing_C']=number of amino acids after the last hit
#['real_position']=real position on the chromosome, not on the FASTA file
log=open(args.wd+"/log/PGRe_pseudogene_construct_step.log", "w")
pseudogeneDic=collections.defaultdict(dict)
parentDic={} #This is used to retrieve the protein coded by the parent gene. It's not integrated into the pseudogeneDic
#for optimisation purpose
with open(args.wd+"/tmp/tblastn_results.tempfile.tsv") as blast_res:
    for line in blast_res:
        if line[0]!="#":
            line=line.split("\t")
            if args.verbose:
                last_pg = int(re.search("pseudogene[0]*([0-9]+)",line[1]).group(1))
            pseudogeneDic[line[1]]['chromosome']=line[0]
            pseudogeneDic[line[1]]['strand']=line[3]
            pseudogeneDic[line[1]]['parent'] = line[4]
            parentDic[line[4]]=""
            pseudogeneDic[line[1]]['intervals'] = []
            for hit in line[7].split(";"):
                pseudogeneDic[line[1]]['intervals'].append(hit.split(":"))
                #Note that the hits positions (i.e. "intervals") are substracted by one to work with Python, as the
                #first character of a string (after, DNA or protein sequences) is at the index 0.
                pseudogeneDic[line[1]]['intervals'][-1][0]=int(pseudogeneDic[line[1]]['intervals'][-1][0])-1
                pseudogeneDic[line[1]]['intervals'][-1][1]=int(pseudogeneDic[line[1]]['intervals'][-1][1])-1
            pseudogeneDic[line[1]]['initial_shift'] = int(pseudogeneDic[line[1]]['intervals'][0][0])+1
            pseudogeneDic[line[1]]['missing_N']=int(line[5])
            pseudogeneDic[line[1]]['missing_C']=int(line[6])
            pseudogeneDic[line[1]]['one_exon_only']=False
            pseudogeneDic[line[1]]['real_position']=int(line[2])

#These two lines add the attribute:
#['sequence']=sequence of the scaffold between the positions (+-100) where the pseudogene has been found
for sequence in list(Bio.SeqIO.parse(args.wd+"/tmp/first_pg_set.tempfile.fasta", 'fasta')):
    if sequence.id in pseudogeneDic:
        pseudogeneDic[sequence.id]['sequence']=sequence.seq

#These two lines add the attribute:
#['parent_seq']=sequence of the protein coded by the parent gene
for sequence in list(Bio.SeqIO.parse(args.fasta,'fasta')):
    if sequence.id in parentDic:
        parentDic[sequence.id]=sequence.seq
if args.verbose: print("               >> BLAST parsing: OK <<")
if args.verbose: print("               >> PG generation: 0.00000% <<", end="")
blosum62=substitution_matrices.load('BLOSUM62')

########################################################################################################################
#################################################FRAME-SHIFT LOOKUP#####################################################
########################################################################################################################

for pseudogene in pseudogeneDic:
    if args.verbose:
        percent=int(re.search("pseudogene[0]*([0-9]+)",pseudogene).group(1))/last_pg*100
        percent=str('{:.5f}'.format(percent))
    print("\r               >> PG generation: "+percent+"% <<", end="")
    sys.stdout.flush()
    parent_prot=parentDic[pseudogeneDic[pseudogene]['parent']]
    open(args.wd+"/tmp/last_parent.fasta","w").write(">last_parent\n"+str(parent_prot))
    #The two next blocks converts the sequence and coordinates if the pseudogene is located on the negative strand
    if pseudogeneDic[pseudogene]['strand']=="-":
        pseudogeneDic[pseudogene]['sequence']=pseudogeneDic[pseudogene]['sequence'].reverse_complement()
        right_lim=len(pseudogeneDic[pseudogene]['sequence'])-1
        pseudogeneDic[pseudogene]['intervals']=invert_coord(pseudogeneDic[pseudogene]['intervals'],right_lim)

    #VITo is used to merge overlapping intervals. It doesn't overlaps intervals when the merged intervals has a length
    #that is not divisible by 3, as it shows frame-shift. The intervals that overlaps after this step are considered
    #frame-shift and are computed by another subprogram (CConIE)
    interval_list=VITo.mergeOverlapWOFrameShift(pseudogeneDic[pseudogene]['intervals'])
    if len(interval_list)==1 : pseudogeneDic[pseudogene]['one_exon_only']=True

    title(pseudogene,char="=")
    title("pseudo-CDS search")
    log.write("pseudo-CDS found at positions: "+str(interval_list)+"\n\n")
    title("looking for frame-shift")

    #For each interval, if it overlaps with the next interval, CConIE functions are used to retrieve the best chimer, i.e.
    #the positions at which the first frame shift stops and the second starts that allows to get the best E-value when
    #blasting against user DB.
    interval_fs_corrected=[]
    fs_found=False
    if len(interval_list)>1:
        for i in range(0,len(interval_list)-1):
            if interval_list[i][1]>=interval_list[i+1][0]:
                fs_found=True
                log.write("Frame-shift found between positions "+str(interval_list[i+1][0])+" and "+str(interval_list[i][1])+"\n")
                chimer_dic= CConIE.chimer_construct(interval_list[i], interval_list[i + 1], pseudogeneDic[pseudogene]['sequence'],log)
                CConIE.chimdic_to_fasta(chimer_dic, args.wd)
                cmd_blastp = NcbiblastpCommandline(cmd=args.binaries + "/blastp",
                    query=args.wd+"/tmp/chim_temp.fasta",out=args.wd+"/tmp/chim_temp.tab", outfmt=6, subject=args.wd+"/tmp/last_parent.fasta", \
                        evalue=0.01, word_size=3, gapextend=2, max_target_seqs=1)
                cmd_blastp()
                best_chimer= CConIE.bestChim(args.wd+"/tmp/chim_temp.tab")
                if best_chimer=="":
                    log.write(" Unable to solve conflict\n\n")
                else:
                    #The two next lines changes the coordinates of the intervals according to CConIE results to avoid frame-shift
                    interval_list[i]=chimer_dic[best_chimer]['interval1']
                    interval_list[i+1]=chimer_dic[best_chimer]['interval2']
                    log.write(" Solved.\n[BEST CHIMER]   "+chimer_dic[best_chimer]['sequence']+"\n\n")
            interval_fs_corrected.append(interval_list[i])
            if i==len(interval_list)-2:  #Iterations stop before last interval, this line avoid it
                interval_fs_corrected.append(interval_list[i+1])

        if not fs_found: log.write("No frame-shifts found, but others might be found during the intron refinement step.\n\n")
        else: log.write("pseudo-CDS frame-shift corrected: "+str(interval_fs_corrected)+"\n\n")
        interval_list=interval_fs_corrected.copy()
        del interval_fs_corrected
    log.write('\n')
    log.flush()

    ########################################################################################################################
    ###############################################INTRON REFINEMENT########################################################
    ########################################################################################################################

    title("intron refinement")
    #Intervals in interval_fs_corrected where used in the first version of P-GRe, with good results. However, part of
    #sequences with low similarity to other sequences can't be retrieved by BLAST, leading to erronous pseudo-introns.
    #To avoid this, interval are extended to the next one, meaning that the pseudo-proteins is the one ~obtained if no
    # pseudo-introns existed. Only frame-shift are conserved. Pseudo-introns are removed at the next step.
    fs_start = []
    fs_list = []
    if len(interval_list)==1:
        log.write("Only one CDS found, skipping intron and splicing steps\n\n")
    else:
        interval_extended = VITo.extendedIntervalToNext(interval_list)
        proto_prot=""
        #Presence of frame-shift and extended intervals not divisble by three can leads to a proteins that contains more or
        #less amino acids that it should (i.e. End of last interval - Start of first interval +1 / 3). The most accurate way to
        #convert amino acids position after introns removal (see later) to their genomic positions is to translate every codon
        #inside the intervals and to stock their positions.
        aa_dic={}
        aa_count=0
        for intervals in interval_extended:
            i=0
            while intervals[0]+i<intervals[1]+1:
                start_of_codon=intervals[0]+i
                stop_of_codon=intervals[0]+i+3
                if stop_of_codon<=intervals[1]+1:
                    aa=str(pseudogeneDic[pseudogene]['sequence'][start_of_codon:stop_of_codon].translate())
                    aa_dic[str(aa_count)]={}
                    aa_dic[str(aa_count)]["aa"]=aa
                    aa_dic[str(aa_count)]["codon_pos"]=[start_of_codon,stop_of_codon]
                    proto_prot += aa  # Construct protein codon by codon
                    aa_count += 1
                i+=3
        log.write("Parent gene: "+pseudogeneDic[pseudogene]['parent']+'\n')
        #TAGLIA is used to find the expected gap positions during the FASTA alignment. This will allows for the modified-
        #Lindley-process to avoid cumulating score when a gap isn't at an expected position. In other words, Lindley process
        # will only cumulate score when gap are due to intron and not to alignment artefact.
        #The "error margin", can be modified here
        error_expected_gap=5
        pep_len=TAGLIA.pepLengthFastRun(args.gff,pseudogeneDic[pseudogene]['parent'],error_expected_gap)
        log.write("Gaps expected in the FASTA alignement at approximative positions: "+str(pep_len[0][0:-1])+"\n")
        log.write("Gaps exact positions: " + str(pep_len[1][0:-1]) + "\n")
        #Alignment is performed with global alignement with exagerated gap penality
        if max(len(parent_prot),len(proto_prot))<=300:  # Gotho's algorithm (O. Gotoh, 1982) for short sequence
            alignement=pairwise2.align.globaldc(proto_prot,parent_prot,blosum62,gap_function,gap_function, one_alignment_only=True)
        else:  # EMBOSS stretcher algorithm (W. Myers & W. Miller, 1988) for long sequence
            alignement=[[]]  #List in list to get the same structure as the pairwise2 result
            open(args.wd + "/tmp/last_long_pg.fasta", "w").write(">"+pseudogene+"\n" + str(proto_prot))
            stretch = StretcherCommandline(cmd=args.binaries + "/stretcher", \
                                           asequence=args.wd + "/tmp/last_long_pg.fasta", \
                                           bsequence=args.wd + "/tmp/last_parent.fasta", gapopen=5, gapextend=0, \
                                           outfile=args.wd + "/tmp/last_stretcher.fasta", aformat='fasta', auto=True)
            stretch()
            for sequence in list(Bio.SeqIO.parse(args.wd + "/tmp/last_stretcher.fasta", 'fasta')):
                alignement[0].append(str(sequence.seq))
        pseudo_prot=TAGLIA.lindleyAlign(alignement,pep_len,error_expected_gap,log)
        #The reconstitued pseudo-protein and the amino-acid dictionnary, where each amino acid is associated with
        #its genomic coordinates, are used to retrieve the exon/intron structures.
        interval_list=TAGLIA.alignToInterv(pseudo_prot,aa_dic)

        ########################################################################################################################
        ##############################################FRAME-SHIFT REFINMENT#####################################################
        ########################################################################################################################

        fs_list=[]
        fs_start=[]
        for i in range (1,len(interval_list)):
            fs=False
            while interval_list[i][0]<=interval_list[i-1][1]:
                interval_list[i-1][1]-=3
                fs=True
            if fs:
                fs_list.append([interval_list[i-1][1]+1,interval_list[i][0]-1])
                fs_start.append(interval_list[i-1][1])
        title("pseudo-CDS intron-corrected")
        log.write("pseudo-CDS found at positions: "+str(interval_list)+"\n")
        log.write("frame-shift: "+str(fs_list)+"\n\n")

        ########################################################################################################################
        ##############################################SPLICING SITE LOOK-UP#####################################################
        ########################################################################################################################

        title("looking for splicing sites")
        previous_free_nucleotide=0
        for i in range(0, len(interval_list)-1):
            #splicing sites are stocked in the ss dictionary, with format:
            #ss[left splicing site id]...
            #...['left_splice'] = [left coordinate of splicing site, right coordinate of splicing site]
            #...['frame'] = number of overhanging nucleotide, i.e. nucleotide at the right extermity which aren't in a triplet
            #...['score'] = distance of the (left) splicing site motif to the previously predicted intron position
            #...['right_splice'][right splicing site id]...
            #... ... ['right_splice'] = [left coordinate of splicing site, right coordinate of splicing site]
            #... ... ['score'] = distance of the (right + left) splicing site motifs to the previously predicted intron position
            ss = {}
            modification = {}
            best_score=math.inf
            best_left=0
            best_right=0
            if interval_list[i][1] not in fs_start:
                #####DEV-TOOLS
                up=9  # Start looking for splicing site 9 bp before start of intron, because word size of BLASTp is 3
                #####
                for j in range(-up,up+1):
                    if pseudogeneDic[pseudogene]['sequence'][interval_list[i][1]+1+j:interval_list[i][1]+3+j]=="GT":
                        left=interval_list[i][1]+1+j
                        right=interval_list[i][1]+3+j
                        id=str(left)+":"+str(right)
                        ss[id] = {}
                        ss[id]['left_splice']=[left,right-1]
                        ss[id]['right_splice']={}
                        ss[id]['length']=left-interval_list[i][0]
                        ss[id]['frame'] = ss[id]['length']%3 #Number of free nucleotide (i.e. not in a triplet). This is used to avoid
                        #frame-shift during translation
                        ss[id]['score']=abs(interval_list[i][1]-left) #Distance to predicted intron start
                        for k in range(-up-2+ss[id]['frame'],up-2,3):
                            if pseudogeneDic[pseudogene]['sequence'][interval_list[i+1][0]+k\
                                    :interval_list[i+1][0]+2+k]=="AG":
                                left=interval_list[i+1][0]+k
                                right=interval_list[i+1][0]+2+k
                                id_r=str(left)+":"+str(right)
                                ss[id]['right_splice'][id_r]={}
                                ss[id]['right_splice'][id_r]['right_splice']=[left,right-1]
                                ss[id]['right_splice'][id_r]['score']=abs(right-interval_list[i+1][0])
                #print('\n'+pseudogene+" ["+str(interval_list[i][1]+1)+":"+str(interval_list[i+1][0]-1)+"] "+str(ss))
                for left_splice_sites in ss:
                    for right_splices_sites in ss[left_splice_sites]['right_splice']:
                        if ss[left_splice_sites]['right_splice'][right_splices_sites]['score']<best_score:
                            best_left=ss[left_splice_sites]['left_splice'][0]
                            best_right=ss[left_splice_sites]['right_splice'][right_splices_sites]['right_splice'][1]
                if best_left!=0 and best_right!=0:
                    modification[str(interval_list[i][1])]=best_left-1
                    modification[str(interval_list[i+1][0])]=best_right+1
                    log.write("Splicing sites found: GT ["+str(best_left)+", "+str(best_left+1)+"]"+\
                              ", AG ["+str(best_right-1)+", "+str(best_right)+"]\n")
        for i in range(0, len(interval_list)):
            key1=str(interval_list[i][0])
            key2=str(interval_list[i][1])
            if key1 in modification:
                interval_list[i][0]=modification[key1]
            if key2 in modification:
                interval_list[i][1]=modification[key2]
        log.write('\n')
        log.flush()
    ########################################################################################################################
    ###############################################N-TER CONSTRUCTION#######################################################
    ########################################################################################################################

    #To finely reconstruct the pseudoprotein, N-ter and C-ter reconstruction are performed. To do so, the hits with the
    #"parent protein" (i.e. the protein encoded by the parent gene) are used to estimate how many amino acids are missing
    #from the reconstituted pseudo-protein in both direction. In the N-ter direction, P-Gre will look for a start codon
    # where a start codon is expected. If none is found, P-Gre will check the next amino acid until a "M" is found.
    aa_up=pseudogeneDic[pseudogene]['missing_N']-1
    start_found=0  # 0=not found, 1=ATG found, 2=mutated start found
    title("n-ter reconstruction")
    log.write("Parent protein has "+str(aa_up)+" amino acids upstream.\n")
    #Before trying to find a start codon, TAGLIA checks if the first exon is missing from the pseudogene, which can
    #be due to many biological phenomenons (genetic recombination, template-switching during retrotranscription, etc)
    if (pseudogeneDic[pseudogene]['one_exon_only'] and aa_up<=50) or \
       (not pseudogeneDic[pseudogene]['one_exon_only'] and not TAGLIA.checkFirstExon(alignement, pep_len[0])):
        aa=""
        i=0
        coordinates=[]
        mutated_start=["ATA","ATC","ATT","AAG","ACG","AGG","TTG","CTG","GTG"]
        upstream=interval_list[0][0]-(aa_up*3)
        #Before trying to find a start codon ("M"), P-Gre allows for "mutated start codon". If any of the codons that
        #can be created by one substitution event on an ATG (see the mutated_start list) is found EXACTLY where the start
        #codon is expected, it is considered a degenerated start codon, else P-GRe will continue to try reconstructing
        #the N-ter part of the pseudo-protein by striclty searching an ATG codon.
        if pseudogeneDic[pseudogene]['sequence'][upstream:upstream + 3].upper() in mutated_start:
            log.write("Mutated start codon found at the exact expected position: ATG ->"+ \
                str(pseudogeneDic[pseudogene]['sequence'][upstream:upstream + 3].upper()) + "\n\n")
            interval_list[0][0]=upstream
            start_found=2
        else:
            #Adding +2 to the expected start position allows for more flexible searching, as P-GRe will look for
            #a start codon a bit further than expected.
            aa_up += 2
            aa = pseudogeneDic[pseudogene]['sequence'][upstream + i:upstream + i + 3].translate()
            coordinates = [upstream + i, upstream + i + 2]
            while aa!="M" and i/3<aa_up:
                i += 1
                aa= pseudogeneDic[pseudogene]['sequence'][upstream + i:upstream + i + 3].translate()
                coordinates = [upstream + i,upstream + i + 2]
            if aa=="M":
                log.write("Start codon found: "+str(coordinates)+"\n\n")
                if i%3==0:
                    interval_list[0][0]=coordinates[0]
                else:
                    unmatched=3-i%3
                    previous_start=interval_list[0][0]
                    interval_list.insert(0,[coordinates[0],previous_start-unmatched-1])
                    fs_list.insert(0,[interval_list[0][1]+1,interval_list[1][0]-1])
                    fs_start.insert(0,interval_list[0][1]+1)
                start_found=1
            else: log.write("Couldn't reconstruct N-ter\n\n")
    else: log.write("Couldn't reconstruct N-ter.\n\n")
    log.flush

    ########################################################################################################################
    ###############################################C-TER CONSTRUCTION#######################################################
    ########################################################################################################################

    #To finely reconstruct the pseudoprotein, N-ter and C-ter reconstruction are performed. To do so, the hits with the
    #"parent protein" (i.e. the protein encoded by the parent gene) are used to estimate how many amino acids are missing
    #from the reconstituted pseudo-protein in both direction. In the C-ter direction, P-Gre will look for a stop codon
    # where a stop codon is expected. If none is found, P-Gre will check the previous amino acid until a "*" is found.
    stop_found=False
    coordinates=[]
    aa_down=len(parentDic[pseudogeneDic[pseudogene]['parent']])-pseudogeneDic[pseudogene]['missing_C']+1
    title("c-ter reconstruction")
    log.write("Parent protein has "+str(aa_down)+" amino acids downstream.\n")
    #Before trying to find a stop codon, TAGLIA checks if the last exon is missing from the pseudogene, which can
    #be due to many biological phenomenons (genetic recombination, template-switching during retrotranscription, etc)
    if (pseudogeneDic[pseudogene]['one_exon_only'] and aa_down <= 50) or \
       (not pseudogeneDic[pseudogene]['one_exon_only'] and not TAGLIA.checkLastExon(alignement, pep_len[0])):
        aa=""
        i=0
        aa_down+=2
        downstream=interval_list[-1][1]+(aa_down*3)
        aa = pseudogeneDic[pseudogene]['sequence'][downstream - i - 2:downstream - i + 1].translate()
        coordinates = [downstream - i - 2, downstream - i]
        while aa!="*" and i/3<aa_down:
            i+=1
            aa= pseudogeneDic[pseudogene]['sequence'][downstream - i - 2:downstream - i + 1].translate()
            coordinates=[downstream-i-2,downstream - i]
        if aa=="*":
            if i % 3 == 0:
                interval_list[-1][1] = coordinates[1]
            else:
                unmatched = 3 - i % 3
                previous_stop = interval_list[-1][1]
                interval_list.append([previous_stop + unmatched + 1, coordinates[1]])
                fs_list.append([interval_list[-2][1] + 1, interval_list[-1][0] - 1])
                fs_start.append(interval_list[-2][1] + 1)
            log.write("Stop codon found: "+ str(coordinates) + "]\n\n")
            stop_found=True
        else: log.write("Couldn't reconstruct C-ter\n\n")
    else: log.write("Couldn't reconstruct C-ter\n\n")
    log.flush()

    ########################################################################################################################
    ###################################################SAVING RESULTS#######################################################
    ########################################################################################################################

    pseudogeneDic[pseudogene]['genomic']=pseudogeneDic[pseudogene]['sequence'][interval_list[0][0]:interval_list[-1][1] + 1]
    sequence=""
    for intervals in interval_list:
        sequence+= pseudogeneDic[pseudogene]['sequence'][intervals[0]:intervals[1] + 1]
        if intervals[1] in fs_start:
            overhang=(intervals[1]-intervals[0]+1)%3
            if overhang>0:
                for i in range(0,overhang):
                    sequence+="N"  #Add N where a frame-shift was found to avoid frame-shift during translation
    pseudogeneDic[pseudogene]['cds']=sequence
    pseudogeneDic[pseudogene]['protein']=pseudogeneDic[pseudogene]['cds'].translate()

    pseudogeneDic[pseudogene]['intervals']=interval_list.copy()
    #The next block is used to retrieve the true coordinates if the pseudogene is on the negative strand
    if pseudogeneDic[pseudogene]['strand']=="-":
        pseudogeneDic[pseudogene]['intervals']=invert_coord(pseudogeneDic[pseudogene]['intervals'],right_lim)
        fs_list=invert_coord(fs_list,right_lim)
    #The next block transforms the intervals from "Python coordinate" (i.e. first character in a string is at index 0)
    #to regular intervals. The real position on the scaffold is also used to retrieve the position of the first
    #nucleotide from the FASTA sequence, which is substracted by the initial shift, i.e. the number of nucleotide
    #that were added upstream of this FASTA sequence by the previous subprogram.
    for intervals in pseudogeneDic[pseudogene]['intervals']:
        intervals[0] += 1 + pseudogeneDic[pseudogene]['real_position'] - pseudogeneDic[pseudogene]['initial_shift']
        intervals[1] += 1 + pseudogeneDic[pseudogene]['real_position'] - pseudogeneDic[pseudogene]['initial_shift']
    for intervals in fs_list:
        intervals[0] += 1 + pseudogeneDic[pseudogene]['real_position'] - pseudogeneDic[pseudogene]['initial_shift']
        intervals[1] += 1 + pseudogeneDic[pseudogene]['real_position'] - pseudogeneDic[pseudogene]['initial_shift']
    pseudogeneDic[pseudogene]['frame_shift']=fs_list
    pseudogeneDic[pseudogene]['start_found']=start_found
    pseudogeneDic[pseudogene]['stop_found']=stop_found

    id="pseudogene_" + pseudogeneDic[pseudogene]['chromosome'] + "_" + str(pseudogeneDic[pseudogene]['intervals'][0][0]) + \
       "_" + str(pseudogeneDic[pseudogene]['intervals'][-1][1])
    id=re.sub("[^a-zA-Z0-9\.\_]","_",id)
    pseudogeneDic[pseudogene]['id']=id
    log.flush()

#The above line transforms the dictionnary in a list of tuple with format:
#tuple[0] = pseudogene name
#tuple[1]['genomic'] = genomic sequence
#tuple[1]['cds'] = CDS sequence
#tuple[1]['protein'] = protein sequence
#tuple[1]['intervals'] = CDS coordinates
#tuple[1]['strand'] = strand
#tuple[1]['frame_shift'] = frame-shift coordinates
#tuple[1]['start_found'] = start found (0, 1 or 2)
#tuple[1]['stop_found'] =  stop found (True/False)
#tuple[1]['chromosome'] = scaffold name
#tuple[1]['id'] = identifier

pseudogeneDic=sorted(pseudogeneDic.items(), key=lambda x: (x[1]['chromosome'],x[1]['intervals'][0][0]))
del parentDic
if args.verbose: print("\r               >> PG generation: OK <<             ")

########################################################################################################################
#################################################CHIMER DETECTION#######################################################
########################################################################################################################

#Next block aggregate pseudogenes that are separated by less than 2500 bases into chimeric pseudogenes
in_overlap=False
new_tuple=()
pseudogeneChimerDic = []
for i in range(0, len(pseudogeneDic)-1):
    if pseudogeneDic[i+1][1]['intervals'][0][0] < pseudogeneDic[i][1]['intervals'][-1][1] + 2500 and \
            pseudogeneDic[i+1][1]['chromosome']==pseudogeneDic[i][1]['chromosome'] and \
            pseudogeneDic[i][1]['stop_found']==False:
        if not in_overlap:
            in_overlap=True
            new_tuple=pseudogeneDic[i]
            new_tuple[1]['genomic']=[pseudogeneDic[i][1]['intervals'][0][0],pseudogeneDic[i+1][1]['intervals'][-1][1]]
        new_tuple[1]['genomic'][1]=pseudogeneDic[i+1][1]['intervals'][-1][1]
        new_tuple[1]['cds']+=pseudogeneDic[i+1][1]['cds']
        new_tuple[1]['parent']+=","+pseudogeneDic[i+1][1]['parent']
        new_tuple[1]['protein']+=pseudogeneDic[i+1][1]['protein']
        for intervals in pseudogeneDic[i+1][1]['intervals']:
            new_tuple[1]['intervals'].append(intervals)
        if new_tuple[1]['strand']!=pseudogeneDic[i+1][1]['strand']:new_tuple[1]['strand']="."
        for fs in pseudogeneDic[i+1][1]['frame_shift']:
            new_tuple[1]['frame_shift'].append(fs)
        if new_tuple[1]['strand']=="+": new_tuple[1]['stop_found']=pseudogeneDic[i+1][1]['stop_found']
        else: new_tuple[1]['start_found']=pseudogeneDic[i+1][1]['start_found']
        if i==len(pseudogeneDic)-2:
            pseudogeneChimerDic.append(new_tuple)
    else:
        if in_overlap:
            in_overlap=False
            pseudogeneChimerDic.append(new_tuple)
        else:
            in_overlap=False
            pseudogeneChimerDic.append(pseudogeneDic[i])
        if i==len(pseudogeneDic)-2:
            pseudogeneChimerDic.append(pseudogeneDic[i+1])

if len(pseudogeneDic)>1:pseudogeneDic=pseudogeneChimerDic.copy()  #If added for exemple data
del pseudogeneChimerDic

# Retreieving genomic sequence for newly merged chimeric pseudogenes
for scaffold in list(Bio.SeqIO.parse(args.scaffold,'fasta')):
    for tuple in pseudogeneDic:
        if isinstance(tuple[1]['genomic'],list) and scaffold.id==tuple[1]['chromosome']:
            tuple[1]['genomic']=scaffold.seq[tuple[1]['intervals'][0][0]+1:tuple[1]['intervals'][-1][1]+2]

########################################################################################################################
##############################################PSEUDOGENEISATION CHECK###################################################
########################################################################################################################

if args.verbose: print("               >> Files:         OK "+args.wd+"/res/ <<")
file=open(args.wd+"/res/true_genes.id", "w")
for pseudogene in pseudogeneDic:
    try:
        prot_seq_to_try=pseudogene[1]['protein']
        if prot_seq_to_try[0]=="M" and len(pseudogene[1]['frame_shift'])==0 and len(prot_seq_to_try)>=80:
            while prot_seq_to_try[-1]=="*":
                prot_seq_to_try=prot_seq_to_try[:-1]
            if prot_seq_to_try.count("*")==0:
                file.write(pseudogene[0]+'\n')
    except:
        print(pseudogene, prot_seq_to_try, pseudogene[1]['frame_shift'])

if args.verbose: print("                                 OK "+args.wd+"/res/true_genes.id <<")

########################################################################################################################
#################################################FILES GENERATION#######################################################
########################################################################################################################

###FASTA
createFile("genomic.fna")
createFile("cds.fna")
createFile("protein.faa")

###GFF
file=open(args.wd+"/res/pseudogenes.gff", "w")
for tuple in pseudogeneDic:
    structure_list=[]
    #Exon and frame-shift are put in the same list, and the type of structure (i.e. "frame-shift" or "pseudo_CDS") is added
    for intervals in tuple[1]['intervals']:
        intervals.append("pseudo_CDS")
        structure_list.append(intervals)
    for intervals in tuple[1]['frame_shift']:
        intervals.append("frame_shift")
        structure_list.append(intervals)
    # The next 4 lines are used to count the number of each feature. If a CDS is splitted by a frame-shift, for example
    # CDS1, the first part of the CDS will be annotated as CDS1.1 and the second part as CDS1.2
    structure_list[0].append("1.1")
    CDS_num=1
    CDS_frag_num=1
    fs_num=1

    structure_list=VITo.sortIntervals(structure_list)
    for i in range(1,len(structure_list)):
        if structure_list[i][2]=='pseudo_CDS':
            if structure_list[i-1][2]=="pseudo_CDS":
                CDS_frag_num=1
                CDS_num+=1
                structure_list[i].append(str(CDS_num)+"."+str(CDS_frag_num))
            else:
                CDS_frag_num+=1
                structure_list[i].append(str(CDS_num) + "." + str(CDS_frag_num))
                structure_list[i-1].append("."+str(fs_num))
                fs_num+=1
    #Generating the GFF file
    #First, write the full pseudogene positions
    file.write(tuple[1]['chromosome'] + "\tP_GRe\tpseudogene\t" + str(tuple[1]['intervals'][0][0]) + "\t" + \
               str(tuple[1]['intervals'][-1][1]) +"\t.\t"+tuple[1]['strand']+"\t.\tID=" + tuple[1]['id'] +";Parent_gene=" + \
               str(tuple[1]['parent']) +"\n")
    #Then start_codon, if found
    if tuple[1]['start_found']>0:
        if tuple[1]['start_found']==1:tuple[1]['start_found']="pseudo_start_codon"
        elif tuple[1]['start_found']==2:tuple[1]['start_found']="degenerated_start_codon"
        if tuple[1]['strand']=="+":
            left=str(structure_list[0][0])
            right=str(structure_list[0][0]+2)
        else:
            left=str(structure_list[-1][1] - 2)
            right=str(structure_list[-1][1])
        file.write(tuple[1]['chromosome'] + "\tP_GRe\t" + tuple[1]['start_found'] + "\t" + \
                   left + "\t" + right +"\t.\t"+tuple[1]['strand']+"\t.\tID=" + tuple[1]['id'] +"." \
                   + tuple[1]['start_found'] +".1;Parent=" + tuple[1]['id'] +"\n")
    #Then CDS and frame-shifts
    for intervals in structure_list:
        file.write(tuple[1]['chromosome'] + "\tP_GRe\t" + intervals[2] + "\t" + str(intervals[0]) + "\t" + str(intervals[1]) + "\t" +\
                "\t.\t"+tuple[1]['strand']+"\t.\tID=" + tuple[1]['id'] +"." + intervals[2] + str(intervals[3]) \
                   +";Parent=" + tuple[1]['id'] +"\n")
    #Then stop_codon, if found
    if tuple[1]['stop_found']:
        if tuple[1]['strand']=="-":
            left=str(structure_list[0][0])
            right=str(structure_list[0][0]+2)
        else:
            left=str(structure_list[-1][1] - 2)
            right=str(structure_list[-1][1])
        file.write(tuple[1]['chromosome'] + "\tP_GRe\tpseudo_stop_codon\t" + left + "\t" + right +"\t.\t"+\
                   tuple[1]['strand']+"\t.\tID=" + tuple[1]['id'] +".pseudo_stop_codon.1;Parent=" + tuple[1]['id']\
                   +"\n")
if args.verbose: print("                                 OK "+args.wd+"/res/pseudogenes.gff <<")
