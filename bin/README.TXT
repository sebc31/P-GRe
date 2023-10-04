#This folder must contains the binary files from the different software used by P-GRe.

#At the time of P-GRe testing, these software/dependencies versions were used:

#########################################
### GFFRead 0.12.7 ######################
#########################################

#GFFRead can be downloaded from https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz and installed with make:

wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz
tar xzf gffread-0.12.7.Linux_x86_64.tar.gz
cd gffread-0.12.7.Linux_x86_64
make

#########################################
### BEDTools 2.30.0 #####################
#########################################

#A precompiled version of BEDTools binary can be downloaded from https://github.com/arq5x/bedtools2/releases.
#Additionally, you should rename the binary file to simply "bedtools" so that P-GRe can work properly:

wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools

#Alternatively, you can compile it from source:

wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
tar -zxvf bedtools-2.30.0.tar.gz
cd bedtools2
make

#########################################
### NCBI BLAST+ 2.13.0+ #################
#########################################

# NCBI BLAST+ latest version can be installed with:

sudo apt-get install ncbi-blast+

# Alternatively, you can download the 2.13.0+ version on which P-GRe was developed from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/:

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz

#########################################
### BioPython 1.81 ######################
#########################################

#To install BioPython, you first need to install Python3 package manager with:

sudo apt-get install python3-pip

#You can then install BioPython by using the pip3 package manager:

sudo pip3 install biopython

#########################################
### EMBOSS:6.6.0.0 ######################
#########################################

#EMBOSS tool suite can be downloaded from ftp://emboss.open-bio.org/pub/EMBOSS/.
#Once downloaded you will need to compile with `make`:

wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar -xvf EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0/
./configure
make

#Binary files will be located in the emboss folder. This folder contains a lot of binary files, but only stretcher is used by P-GRe.
