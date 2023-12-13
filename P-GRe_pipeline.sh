###########################################################
###########################################################
# P-GRe pipeline v1.0                                     #
###########################################################
###########################################################


#P-GRe is a pipeline dedicated to the automatic detection and annotation of pseudogenes. If you are usign P-GRe, please cite:
#XXXXXXX.

###########################################################
# Help & options                                          #
###########################################################

Help() {
  echo "usage: PGRe_pipeline.sh -f fasta -g gff [-o output directory] [-p proteome] [-b blast results] [-d path to P-GRe] [-t thread] [-v]"
  echo "MANDATORY ARGUMENTS"
  echo "       -f   path to the genome to annotate, in FASTA format"
  echo "       -g   path to the genome annotation, in GFF format"
  echo "OPTIONAL ARGUMENTS"
  echo "       -o   output directory. If none is given results will be written in the directory containing the genome to annotate"
  echo "       -p   path to the proteome corresponding to the genome to annotate. Specifying this option avoid the proteome creation step"
  echo "       -b   path to the previous BLAST results obtained by tblastn with proteome as the query and masked genome as the database"
  echo "       -d   path to the directory containing the P-GRe_pipeline.sh. This might be necessary if running P-GRe on a cluster"
  echo "       -t   number of threads to use for multithreading. Currently only affects the tBLASTn step"
  echo "       -v   enable verbose mode"
  echo "       -h   displays this help and exit"
  exit
}

#Set option variables
GENOME=""; GFF=""; OUTDIR=""; PROTEOME=""; BLAST_RES=""; PGRE_PATH=$(dirname $(readlink -f "$0")); THREAD=1; VERBOSE=false

#Get the options
while getopts "f:g:o:p:b:d:t:vh" option; do
  case $option in
    h)#display help
      Help;;
    f)#Retrieve genome to annotate
      GENOME=$(readlink -f $OPTARG);;
    g)#Retrieve genome annotation
      GFF=$(readlink -f $OPTARG);;
    o)#Retrieve output directory
      OUTDIR=$(readlink -f $OPTARG);;
    p)#Retrieve proteome
      PROTEOME=$(readlink -f $OPTARG);;
    b)#Retrieve BLAST results
      BLAST_RES=$(readlink -f $OPTARG);;
    d)#Retrieve P-Gre_pipeline.sh directory
      PGRE_PATH=$(readlink -f $OPTARG);;
    t)#Retrieve number of threads
      THREAD=$OPTARG;;
    v)#Retrieve verbose option
      VERBOSE=true;;
    \?)#Invalid option
      echo "Error: Invalid option"
      Help;;
  esac
done

#Checking options and testing if no options are specified
AUTO_CREATE_WD=false
if [ -z $OUTDIR ]; then OUTDIR=$(dirname $GENOME)/PGRe_$(date +%d_%m_%Y_%Hh_%Mm_%Ss); AUTO_CREATE_WD=true; fi
if [[ ! $THREAD =~ ^[0-9]+$ ]]; then echo 'Error: Number of threads must be an integer, not '$THREAD; Help; fi
if [ $OPTIND -eq 1 ]; then Help; fi
if [ -z $GENOME ]; then >&2 echo "Error: No FASTA file"; Help; fi
if [ -z $GFF ]; then >&2 echo "Error: No GFF file"; Help; fi

COMMAND="P-GRe_pipeline.sh -f "$GENOME" -g "$GFF" -o "$OUTDIR
if [ ! -z $PROTEOME ]; then COMMAND+=" -p "$PROTEOME; fi
if [ ! -z $BLAST_RES ]; then COMMAND+=" -b "$BLAST_RES; fi
if [ ! -z $PGRE_PATH ]; then COMMAND+=" -d "$PGRE_PATH; fi
if [ ! -z $THREAD ]; then COMMAND+=" -t "$THREAD; fi
if [ $VERBOSE = true ]; then COMMAND+=" -v"; fi


###########################################################
# Errors                                                  #
###########################################################

MissingScript(){
  >&2 echo ""
  >&2 echo "Error: At least one of the P-GRe Python script is missing"
  >&2 echo "Please, check that all P-GRe scripts are located in the following directory:"
  >&2 echo $PGRE_PATH/script/
  >&2 echo "If this path isn't the one where P-GRe has been installed, you might be using a job-scheduled cluster to run P-GRe. In this case, you should provid the P-GRe_pipeline.sh directory with the -d option"
  exit 1
}

MissingDep(){
  >&2 echo ""
  >&2 echo "Error: At least one of the dependency binary is missing"
  >&2 echo "Please, check that tblastn, blastp, bedtools and XXXX are located in the following directory:"
  >&2 echo $PGRE_PATH/bin/
  >&2 echo "If this path isn't the one where P-GRe has been installed, you might be using a job-scheduled cluster to run P-GRe. In this case, you should provid the P-GRe_pipeline.sh directory with the -d option"
  exit 1
}

MissingFiles(){
  >&2 echo ""
  >&2 echo "Error: At least one of the input files doesn't exist"
  exit 1
}

########################################################### 
###########################################################
# Main program                                            #
###########################################################
###########################################################


if [ $VERBOSE = true ]; then
  echo "*********************************";
  echo " _____         _____ _____"
  echo "|  __ \       / ____|  __ \\"
  echo "| |__) |_____| |  __| |__) |___ "
  echo "|  ___/______| | |_ |  _  // _ \\"
  echo "| |          | |__| | | \ \  __/"
  echo "|_|           \_____|_|  \_\___|"
  echo ""
  echo "*********************************"            
  echo "P-GRe v1.0 starting."
  echo "Command:"
  echo $PGRE_PATH/$COMMAND
  echo ""
  fi

###########################################################
# Checking P-GRe installation                             #
###########################################################

INSTAL_IS_OK=true
if [ $VERBOSE = true ]; then echo "Checking P-GRe installation..."; fi
if [ -f $PGRE_PATH/script/BASNI.py ]; then if [ $VERBOSE = true ]; then echo "BASNI:         OK "$PGRE_PATH/script/BASNI.py; fi; else echo "BASNI:         Error"; INSTAL_IS_OK=false; fi
if [ -f $PGRE_PATH/script/CConIE.py ]; then if [ $VERBOSE = true ]; then echo "CConIE:        OK "$PGRE_PATH/script/CConIE.py; fi; else echo "CConIE:        Error"; INSTAL_IS_OK=false; fi
if [ -f $PGRE_PATH/script/PGRe.py ]; then if [ $VERBOSE = true ]; then echo "PGRe:          OK "$PGRE_PATH/script/PGRe.py; fi; else echo "P-GRe:         Error"; INSTAL_IS_OK=false; fi
if [ -f $PGRE_PATH/script/TAGLIA.py ]; then if [ $VERBOSE = true ]; then echo "TAGLIA:        OK "$PGRE_PATH/script/TAGLIA.py; fi; else echo "TAGLIA:        Error"; INSTAL_IS_OK=false; fi
if [ -f $PGRE_PATH/script/VITo.py ]; then if [ $VERBOSE = true ]; then echo "VITo:          OK "$PGRE_PATH/script/VITo.py; fi; else echo "VITo:          Error"; INSTAL_IS_OK=false; fi
if [ -f $PGRE_PATH/script/PolyGet.py ]; then if [ $VERBOSE = true ]; then echo "PolyGet:       OK "$PGRE_PATH/script/PolyGet.py; fi; else echo "PolyGet:       Error"; INSTAL_IS_OK=false; fi
if [ $INSTAL_IS_OK = false ]; then MissingScript; fi

DEPENDENCIES_ARE_OK=true
if [ $VERBOSE = true ]; then echo ""; echo "Checking dependencies..."; fi
if [ -f $PGRE_PATH/bin/tblastn ]; then if [ $VERBOSE = true ]; then echo "tblastn:       OK "$PGRE_PATH/bin/tblastn; fi; else echo "tblastn:       Error"; DEPENDENCIES_ARE_OK=false; fi
if [ -f $PGRE_PATH/bin/blastp ]; then if [ $VERBOSE = true ]; then echo "blastp:        OK "$PGRE_PATH/bin/blastp; fi; else echo "blastp:        Error"; DEPENDENCIES_ARE_OK=false; fi
if [ -f $PGRE_PATH/bin/makeblastdb ]; then if [ $VERBOSE = true ]; then echo "makeblastdb:   OK "$PGRE_PATH/bin/makeblastdb; fi; else echo "makeblastdb:   Error"; DEPENDENCIES_ARE_OK=false; fi
if [ -f $PGRE_PATH/bin/bedtools ]; then if [ $VERBOSE = true ]; then echo "bedtools:      OK "$PGRE_PATH/bin/bedtools; fi; else echo "bedtools:      Error"; DEPENDENCIES_ARE_OK=false; fi
if [ -f $PGRE_PATH/bin/gffread ]; then if [ $VERBOSE = true ]; then echo "gffread:       OK "$PGRE_PATH/bin/gffread; fi; else echo "gffread:       Error"; DEPENDENCIES_ARE_OK=false; fi
if [ -f $PGRE_PATH/bin/stretcher ]; then if [ $VERBOSE = true ]; then echo "stretcher:     OK "$PGRE_PATH/bin/stretcher; fi; else echo "stretcher:     Error"; DEPENDENCIES_ARE_OK=false; fi
if [ $DEPENDENCIES_ARE_OK = false ]; then MissingDep; fi

FILES_ARE_OK=true
if [ $VERBOSE = true ]; then echo ""; echo "Checking input files..."; fi
#Mandatory inputs
if [ -f $GENOME ]; then if [ $VERBOSE = true ]; then echo "FASTA:         OK" $GENOME; fi; else echo "FASTA:         Error"; FILES_ARE_OK=false; fi
if [ -f $GFF ]; then if [ $VERBOSE = true ]; then echo "GFF:           OK" $GFF; fi; else echo "GFF:           Error"; FILES_ARE_OK=false; fi
#Optional inputs
if [ ! -z $PROTEOME ]; then if [ -f $PROTEOME ]; then if [ $VERBOSE = true ]; then echo "PROTEOME:      OK" $PROTEOME; fi; else echo "PROTEOME:      Error"; FILES_ARE_OK=false; fi; fi
if [ ! -z $BLAST_RES ]; then if [ -f $BLAST_RES ]; then if [ $VERBOSE = true ]; then echo "BLAST RESULTS: OK" $BLAST_RES; fi; else echo "BLAST RESULTS: Error"; FILES_ARE_OK=false; fi; fi
if [ $FILES_ARE_OK = false ]; then MissingFiles; fi

###########################################################
# Creating work files and directories                     #
###########################################################

##### Creating directories
if [ $VERBOSE = true ]; then echo ""; echo "Creating working directories and files..."; fi
if [ $VERBOSE = true ] && [ $AUTO_CREATE_WD = true ]; then echo "Warning: No directory was given for outputting the results. Writting in "$OUTDIR; fi
mkdir -p $OUTDIR
mkdir -p $OUTDIR/tmp; mkdir -p $OUTDIR/data; mkdir -p $OUTDIR/res; mkdir -p $OUTDIR/log

##### Verifications
#Check if files are already present in the data directory. When copying files to the working directory, P-GRe uses cat -v instead of cp to remove every binary from text file, as Python tends to throw an error when reading binary with "open()" function
##Mandatory files
if [ -f $OUTDIR/data/$(basename $GENOME) ]; then if [ $VERBOSE = true ]; then echo "Warning: the requested genome is already in the data/ folder. P-GRe will use the already existing genome."; fi;
else cat -v $GENOME > $OUTDIR/data/$(basename $GENOME); fi
if [ -f $OUTDIR/data/$(basename $GFF) ]; then if [ $VERBOSE = true ]; then echo "Warning: the specified GFF file is already in the data/ folder. P-GRe will use the already existing GFF file."; fi;
else cat -v $GFF > $OUTDIR/data/$(basename $GFF); fi
##Optionial files
###PROTEOME
if [ ! -z $PROTEOME ]; then
  if [ -f $OUTDIR/tmp/$(basename $PROTEOME) ]; then
    if [ $VERBOSE = true ]; then echo "Warning: the specified proteome file is already in the tmp/ folder. P-GRe will use the already existing proteome file."; fi
  else cat -v $PROTEOME > $OUTDIR/tmp/$(basename $PROTEOME); fi
  PROTEOME=$OUTDIR/tmp/$(basename $PROTEOME); fi
###BLAST RESULTS
if [ ! -z $BLAST_RES ]; then
  if [ -f $OUTDIR/tmp/$(basename $BLAST_RES) ]; then
    if [ $VERBOSE = true ]; then echo "Warning: the specified BLAST results file is already in the tmp/ folder. P-GRe will use the already existing proteome file."; fi
  else cat -v $BLAST_RES > $OUTDIR/tmp/$(basename $BLAST_RES); fi
  BLAST_RES=$OUTDIR/tmp/$(basename $BLAST_RES); fi

##### Generating working files
if [ -z $BLAST_RES ]; then
  #Generating the gene-only GFF file to mask the genome
  GO_GFF=$OUTDIR/tmp/$(basename $GFF | sed 's/\.[^\.]*$/.gene_only.tempfile.gff/g')  #Same name as the input GFF file, with extension replaced by .gene_only.tempfile.gff
  grep -P '\tgene\t' $OUTDIR/data/$(basename $GFF) > $GO_GFF
  if [ $VERBOSE = true ]; then echo "gene-only GFF: OK "$GO_GFF; fi
  if [ $VERBOSE = true ]; then echo "               >> INFO: "$(wc -l < $GO_GFF)" genes extracted from the GFF file <<"; fi
  
  #Masking the genome on the gene loci
  MASKED_GENOME=$OUTDIR/tmp/$(basename $GENOME | sed 's/\.[^\.]*$/.hard_masked.tempfile.fna/g')  #Same name as the input FASTA file, with extension replaced by .hard_masked.tempfile.fna
  rm -f $OUTDIR/log/bedtools_genome_masking_step.log
  $PGRE_PATH/bin/bedtools maskfasta -fi $OUTDIR/data/$(basename $GENOME) -bed $GO_GFF -fo $MASKED_GENOME >> $OUTDIR/log/bedtools_genome_masking_step.log 2>&1
  if [ ! $? -eq 0 ]; then echo "fasta masking: Error, check: "$OUTDIR"/log/bedtools_genome_masking_step.log"; exit 1; fi
  if [ $VERBOSE = true ]; then echo "fasta masking: OK "$MASKED_GENOME; fi
  
  #Making masked-genome db for BLAST
  rm -f $OUTDIR/log/makeblast_db_genome_masking_step.log
  $PGRE_PATH/bin/makeblastdb -in $MASKED_GENOME -dbtype nucl >> $OUTDIR/log/makeblast_db_genome_masking_step.log 2>&1
  if [ ! $? -eq 0 ]; then echo "makeblastdb:   Error, check: "$OUTDIR"/log/makeblast_db_genome_masking_step.log"; exit 1; fi
  if [ $VERBOSE = true ]; then echo "makeblastdb:   OK "$MASKED_GENOME; fi
  if [ $VERBOSE = true ]; then for file in $MASKED_GENOME.n*; do echo "               OK "$file; done; fi
  if [ $VERBOSE = true ]; then echo "               >> INFO: "$(sed -n 's/.*added \([0-9]*\) sequences.*/\1/gp' $OUTDIR/log/makeblast_db_genome_masking_step.log)" chromosomes/scaffolds added <<"; fi

elif [ $VERBOSE = true ]; then
  echo "gene-only GFF: SKIP "; echo "fasta masking: SKIP"; echo "makeblastdb:   SKIP"
  echo "               >> INFO: BLAST results provided: "$BLAST_RES", skipping all genome masking related steps <<"; fi

#Generating proteome from genome and annotation file
rm -f $OUTDIR/log/gffrad_proteome_generation_step.log
if [ -z $PROTEOME ]; then
  PROTEOME=$OUTDIR/tmp/$(basename $GENOME | sed 's/\.[^\.]*$/.prot.tempfile.faa/g')  #Same name as the input genome FASTA file, with extension replaced by .prot.tempfile.faa
  $PGRE_PATH/bin/gffread -S -y $PROTEOME -g $OUTDIR/data/$(basename $GENOME) $OUTDIR/data/$(basename $GFF) >> $OUTDIR/log/gffread_proteome_generation_step.log 2>&1
  if [ ! $? -eq 0 ]; then echo "gffread:       Error, check: "$OUTDIR"/log/gffread_proteome_generation_step.log"; exit 1; fi
  if [ $VERBOSE = true ]; then echo "gffread:       OK "$PROTEOME; fi
  if [ $VERBOSE = true ]; then echo "               >> INFO: "$(grep -c '>' $PROTEOME)" proteins generated <<"; fi
elif [ $VERBOSE = true ]; then
  echo "gffread:       SKIP "
  echo "               >> INFO: Proteome provided: "$PROTEOME" , skipping GFF reading step <<"; fi

###########################################################
# Running P-GRe                                           #
###########################################################

if [ $VERBOSE = true ]; then
  echo ""
  echo "*********************************"            
  echo "P-GRe running."
  echo ""
fi

#Blasting proteome against masked genome
if [ -z $BLAST_RES ]; then
  rm -f $OUTDIR/log/tblastn_proteome_vs_masked_genome_step.log
  BLAST_RES=$OUTDIR/tmp/proteome_vs_masked_genome_tblastn.blast
  if [ $VERBOSE = true ]; then echo "tblastn:       RUNNING..."; fi
  $PGRE_PATH/bin/tblastn -query $PROTEOME -db $MASKED_GENOME -seg 'yes' -db_soft_mask dust -soft_masking True -outfmt 6 -evalue 0.01 -word_size 3 -gapextend 2 -max_intron_length 10000 -out $BLAST_RES -num_threads $THREAD >> $OUTDIR/log/tblastn_step.log 2>&1
  if [ ! $? -eq 0 ]; then echo -e "\033[1A\033[Ktblastn:       Error, check: "$OUTDIR"/log/tblastn_step.log"; exit 1; fi
  if [ $VERBOSE = true ]; then echo -e "\033[1A\033[Ktblastn:   OK "$BLAST_RES; fi
elif [ $VERBOSE = true ]; then
  echo "tblastn:       SKIP "$BLAST_RES
  echo "               >> INFO: BLAST results provided: "$BLAST_RES", skipping tBLASTn step <<"; fi

#Filtering hits according to BASNI algorithm
rm -f $OUTDIR/log/PGRe_BASNI_step.err
if [ $VERBOSE = true ]; then echo "filtering hits:"; fi
python $PGRE_PATH/script/BASNI.py -f $OUTDIR/data/$(basename $GENOME) -r $BLAST_RES -w $OUTDIR -v $VERBOSE 2>> $OUTDIR/log/PGRe_BASNI_step.err  #Redirecting stderr only
if [ ! $? -eq 0 ]; then echo "filtering hits:Error, check: "$OUTDIR"/log/PGRe_BASNI_step.err"; exit 1; fi

#Reconstructing pseudogenes
rm -f $OUTDIR/log/PGRe_pseudogene_construct_step.err
python $PGRE_PATH/script/PGRe.py -f $PROTEOME -g $OUTDIR/data/$(basename $GFF) -s $OUTDIR/data/$(basename $GENOME) -w $OUTDIR -b $PGRE_PATH/bin -v $VERBOSE 2>> $OUTDIR/log/PGRe_pseudogene_construct_step.err  #Redirecting stderr only, PGRe.py will write other outputs in stdout/.log file
if [ ! $? -eq 0 ]; then echo ""; echo "pg generation: Error, check: "$OUTDIR"/log/PGRe_pseudogene_construct_step.err"; exit 1; fi

#Terminating
if [ $VERBOSE = true ]; then
  echo ""
  echo "*********************************"            
  echo "P-GRe finished."
  echo ""
fi
