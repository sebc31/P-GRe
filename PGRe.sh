#!/bin/bash

##############################################################################################################################
#### HELP ####################################################################################################################
##############################################################################################################################


echo -e "P-GRe pipeline v1.0\n\nP-GRe is a pipeline dedicated to the automatic detection and annotation of pseudogenes. If you are usign P-GRe, please cite:\n"

Help() {
  echo "usage: PGRe.sh -f fasta -g gff -p proteins [-o output directory] [-u other proteins] [-d path to P-GRe] [-t thread]"
  echo "MANDATORY ARGUMENTS"
  echo "       -f   path to the genome to annotate, in FASTA format"
  echo "       -g   path to the genome annotation, in GFF format"
  echo "       -p   path to the organism protein sequences, in FASTA format"
  echo "OPTIONAL ARGUMENTS"
  echo "       -o   output directory. If none is given results will be written in the directory containing the genome to annotate [fasta]"
  echo "       -u   path to protein sequences from other organisms, used for unitary pseudogene predictions"
  echo "       -d   path to the directory containing the PGRe.sh script. This might be necessary if running P-GRe on a cluster [.]"
  echo "       -t   number of threads to use for multi-threading [1]"
  echo "       -h   displays this help and exit"
  exit
}

#Set option variables
GENOME=""; GFF=""; OUTDIR=""; PROTEOME=""; PGRE_PATH=$(dirname $(readlink -f "$0")); THREAD=1; PROTEOME2=""; COMMAND="PGRe.sh"

#Get the options
while getopts "f:g:p:m:o:u:d:t:h" option; do
  case $option in
  
    h)#display help
      Help;;
      
    f)#Retrieve genome to annotate
      GENOME=$(readlink -f $OPTARG);
      if [ ! -f $GENOME ]
      	then echo "ERROR with -f option: $GENOME file does not exist."
      	exit
      fi;
      COMMAND+=" -f $GENOME";;
      
    g)#Retrieve genome annotation
      GFF=$(readlink -f $OPTARG);
      if [ ! -f $GFF ]
      	then echo "ERROR with -g option: $GFF file does not exist."
      	exit
      fi;
      COMMAND+=" -g $GFF";;
      
    o)#Retrieve output directory
      OUTDIR=$OPTARG;
      COMMAND+=" -o $OUTDIR";;
      
    p)#Retrieve proteome
      PROTEOME=$(readlink -f $OPTARG);
      if [ ! -f $PROTEOME ]
      	then echo "ERROR with -p option: $PROTEOME file does not exist."
      	exit
      fi;
      COMMAND+=" -p $PROTEOME";;
      
    u)#Retrieve second proteome
      PROTEOME2=$(readlink -f $OPTARG);
      if [ ! -f $PROTEOME2 ]
      	then echo "ERROR with -u option: $PROTEOME2 file does not exist."
      	exit
      fi;
      COMMAND+=" -u $PROTEOME2";;
      
    d)#Retrieve PGRe.sh directory
      PGRE_PATH=$(readlink -f $OPTARG);
      if [ ! -d $PGRE_PATH ]
      	then echo "ERROR with -d option: $PGRE_PATH path does not exist."
      	exit
      fi;
      COMMAND+=" -d $PGRE_PATH";;
      
    t)#Retrieve number of threads
      THREAD=$OPTARG;
      if [[ ! $THREAD =~ ^[0-9]+$ ]]
      	then echo "ERROR with option -t: Number of threads must be an integer, not $THREAD."
      	exit
      fi;
      COMMAND+=" -t $THREAD";;
      
    \?)#Invalid option
      exit;;
  esac
done

#Check if no arguments are given
if [ $OPTIND -eq 1 ]; then Help; fi

#Retrieve script path
SCRIPTPATH=$(dirname $(realpath "$0"))


##############################################################################################################################
#### WORKING FILES ###########################################################################################################
##############################################################################################################################

#Create WD if none is given
if [ -z $OUTDIR ]
	then OUTDIR=$(dirname $GENOME)/PGRe_$(date +%d_%m_%Y_%Hh_%Mm_%Ss)
	mkdir $OUTDIR
	else mkdir -p $OUTDIR
fi

echo "Running P-GRe. Command used: $COMMAND"

#Retrieve gene coordinates
echo -e "\nCreating working files..."
awk '
BEGIN{OFS="\t"} 
/\tmRNA\t/{
	if ($4<$5){print $1,$4,$5} 
	else {print $1,$5,$4}}
' $GFF > $OUTDIR/mrna_coordinates.bed

#Mask genes sequences in genome
bedtools maskfasta -fi $GENOME -bed $OUTDIR/mrna_coordinates.bed -fo $OUTDIR/masked_genome.fasta

#Get the in-organism protein id from the FASTA file
awk '/>/{print substr($1,2,length($1))}' $PROTEOME > $OUTDIR/in_organism_prot.id
echo -e "Done."

#Concat the two proteome files
PROTEOME_FULL=$PROTEOME
if [ ! -z $PROTEOME2 ]
	then cat $PROTEOME $PROTEOME2 > $OUTDIR/all_prot.fa
	PROTEOME_FULL=$OUTDIR/all_prot.fa
fi

##############################################################################################################################
#### MAIN ####################################################################################################################
##############################################################################################################################

echo -e "\nRunning miniprot... This may take a few minutes."
miniprot -L 15 -B 5 -J 25 -F 23 -I -t $THREAD -p 0.6 --outs 0.7 --outc 0.1 --gff --aln $OUTDIR/masked_genome.fasta $PROTEOME_FULL > $OUTDIR/miniprot_res.gff

echo -e "Done.\n\nFiltering results... This may take a few minutes."
python3 $SCRIPTPATH/overlap_filter.py $OUTDIR/miniprot_res.gff $OUTDIR/mrna_coordinates.bed $OUTDIR $OUTDIR/in_organism_prot.id
sed -rn '/pseudogene/s/.+Target=([^ ]+).+/\1/gp' $OUTDIR/miniprot_res_filtered.gff | grep -Fxv -f $OUTDIR/in_organism_prot.id - > $OUTDIR/check_unitarity.list

touch $OUTDIR/diamond.res
if [ ! -z $PROTEOME2 ]	
	then echo -e "Done.\n\nLooking for unitary pseudogenes with DIAMOND... This may take a few minutes."
	python3 $SCRIPTPATH/extract_unitary.py $OUTDIR/check_unitarity.list $PROTEOME2 > $OUTDIR/unitary_check.fasta
	diamond makedb --in $PROTEOME -d $OUTDIR/diamond.db -p $THREAD
	diamond blastp -p $THREAD -k 1 -f 6 qseqid sseqid -d $OUTDIR/diamond.db --query $OUTDIR/unitary_check.fasta --out $OUTDIR/diamond.res
fi
awk 'BEGIN{OFS="\t"} {print $1,$1}' $OUTDIR/in_organism_prot.id | cat - $OUTDIR/diamond.res > $OUTDIR/in_op_full.id

echo -e "Done\n\nChecking pseudogene's parent structures to determine their types..."
python3 $SCRIPTPATH/compute_expected_structure.py $GFF > $OUTDIR/exp_struc.tsv

##############################################################################################################################
#### OUTPUTS #################################################################################################################
##############################################################################################################################

echo -e "Done.\n\nInferring pseudogene types and writing output..."
python3 $SCRIPTPATH/write_output.py $OUTDIR/in_op_full.id $OUTDIR/exp_struc.tsv $OUTDIR/miniprot_res_filtered.gff $OUTDIR > $OUTDIR/PGRe.unsorted.res
sort --version-sort -k1,1 -k4,4 $OUTDIR/PGRe.unsorted.res > $OUTDIR/PGRe.gff

python3 $SCRIPTPATH/get_seq.py $OUTDIR/PGRe.gff $OUTDIR/miniprot_res.gff > $OUTDIR/pseudogene_protein.fasta

#echo -e "Done.\n\nCleaning working directory..."
mkdir -p $OUTDIR/tmp
mv $OUTDIR/* $OUTDIR/tmp 2>/dev/null
mv $OUTDIR/tmp/PGRe.gff $OUTDIR/tmp/pseudogene_protein.fasta $OUTDIR 2>/dev/null

echo -e "Done.\n\nThank you for using P-GRe."
