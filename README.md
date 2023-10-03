# P-GRe v.0
Contacts and authors
====================
Sébastien Cabanac, PhD student, sebastien.cabanac@univ-tlse3.fr

Christophe Dunand, professor at Université Paul Sabatier Toulouse 3, christophe.dunand@univ-tlse3.fr

Catherine Mathé, associate professor at Université Paul Sabatier Toulouse 3, catherine.mathe-dehais@univ-tlse3.fr

All members of the Laboratoire de Recherche en Sciences Végétales, Université de Toulouse, CNRS, Université Paul Sabatier Toulouse 3, Toulouse INP, Auzeville-Tolosane, France.

Funding and acknowledgements
============================
The authors are thankful to the Paul Sabatier-Toulouse 3 University and to the Centre National de la Recherche Scientifique (CNRS) for granting their work. Sébastien Cabanac is the recipient of a fellowship from the “École Universitaire de Recherche (EUR)” TULIP-GS (ANR-18-EURE-0019). This study is set within the framework of the “Laboratoires d’Excellences (LABEX)” TULIP (ANR-10-LABX-41).

Contents
========

- [Contacts and authors](#contacts-and-authors)
- [Funding and acknowledgements](#funding-and-acknowledgements)

About P-GRe
===========
Pseudogenes are genomic sequences with homology to functional genes but that harbor deleterious mutations, such as loss of the start codon, loss of coding sequence, gain of stop or frame-shifts. No longer coding for a functional protein, pseudogenes are rarely transcribed and are often described as having no function. It is now known that part of the pseudogenes are transcribed, this part representing for example 15% of all the pseudogenes in mice (Sisu et al., 2020). It has also been shown that these transcripts, originating from pseudogenes, could form duplexes with the mRNAs of homologous functional genes and thus participate in post-transcriptional regulation through the RNAi pathway (Tam et al., 2008; Watanabe et al., 2008; Guo et al., 2009). In the other hand, the exhaustive prediction of pseudogenes allow a better understanding of dynamic of gene evolution in multigenenic families often subjected to duplication and pseudogenisation events.

The goal of PseudoGene REtriever (P-GRe [/pɛɡʁ/]) is to find the position of pseudogenes on a genome, as well as to infer their structures in pseudo-exons and pseudo-introns. P-GRe aims to be more user-friendly, with a limited number of dependencies, ease of use and total automaticity, while producing qualitative results and having greater sensitivity than other software with the same goal.

Requirements and input
======================

- Use a high quality genome assembly in FASTA or multi-FASTA format.

- Use an annotation file corresponding to the genome to be annotated in GFF3 format. The GTF format is currently not supported. Genes must be designated by the type "gene" in the third column of this file. CDS must be designated by the type “CDS”, and the gene to which each CDS belongs must be designated by the attribute "Parent" or "Name".

- If a proteome file is provided, it must be in FASTA or multi-FASTA format. For now, each protein present in this file **must** have its correspondence on the annotation file in GFF format. Conversely, each gene in the GFF format annotation file **must** have a correspondence in the proteome file.

How P-GRe works
===============
<p align="center"><img src="docs/figs/figure_article.jpg" width="50%" height="50%"></p>
P-GRe works in two main stages, themselves divided into several sub-stages. The first main step is the search for the position of pseudogenes on the genome, and the second main step is the refinement of the structure of the pseudogenes found.

## Finding the pseudogenes
(See the numbered sub-steps on gray background on the top figure) The search for pseudogenes is carried out in four sub-steps:

- **[1]** Gene annotation and genome files are used to generate the proteome file. Alternatively, a proteome file can be provided to P-GRe.
- **[2]** Gene positions are used to hard-mask gene loci on the genome.
- **[3]** All protein sequences are aligned locally to the genome. Because the original loci of the protein-coding genes are masked, the proteins will be aligned to areas of the genome with similarity to the original coding gene. It is assumed that the areas of similarity thus found are potential pseudogenes derived from the coding genes, called "parent genes" in the literature.
- **[4]** The results of local alignments are filtered according to several criteria, in particular the percentage of homology between a potential pseudogene and the sequence of the aligned protein. Close hits obtained from alignments of the same sequence are merged together.

## Finding the structure of pseudogenes
(See the numbered sub-steps on black background on the top figure) Pseudogene structure inference is divided into three substeps:

- **[1]** Hits that overlap, and whose overlap length is not divisible by 3, are considered frame-shift markers, because this can show that two parts of a protein sequence encoded by the same CDS are aligned in different reading frames. The presice position of the frame-shift is found by a so-called “chimera” approach. More information on this approach can be found [here](docs/figs/chimeras.md).
- **[2]** The hits obtained may have the defect of not completely covering the (pseudo-)exons of the pseudogenes. To correct this, each hit obtained for a pseudogene is extended to the next hit and translated. Thus, a peptide sequence devoid of frame-shift is obtained, but retaining the introns. By aligning this sequence with the sequence of the protein encoded by the parent gene, extended gaps are expected at the introns loci (Note that an absence of gap marks the presence of a retropseudogene). The alignments obtained are corrected by a process inspired by the Lindley process. More information on this "Lindley-inspired" process can be found [here]().
- **[3]** The ends of the pseudogenes are refined by searching for a start codon and a stop codon. For start codons, P-GRe accepts "degenerate" start codons (*i.e.* which have a single substitution) provided that they are at a precise position upstream of the pseudogene. This position is determined from the alignment between the amino acid sequences encoded by the pseudogene and its parent gene.

Installation
============

At the time of P-GRe testing, these software/dependencies versions were used:

### GFFRead 0.12.7

[GFFRead](https://github.com/gpertea/gffread) can be downloaded from [here](https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz) and installed with `make`:
```
wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz
tar xzf gffread-0.12.7.Linux_x86_64.tar.gz
cd gffread-0.12.7.Linux_x86_64
make
```
<br>

### BEDTools 2.30.0

A precompiled version of [BEDTools](https://bedtools.readthedocs.io/en/latest/) binary can be downloaded from [here](https://github.com/arq5x/bedtools2/releases). Additionally, you should rename the binary file to simply "bedtools" so that P-GRe can work properly:
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
```
Alternatively, you can compile it from source:
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
tar -zxvf bedtools-2.30.0.tar.gz
cd bedtools2
make
```
<br>

### NCBI BLAST+ 2.13.0+

[NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) latest version can be installed with:
```
sudo apt-get install ncbi-blast+
```
Alternatively, you can download the 2.13.0+ version on which P-GRe was developed from [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/):
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz
```
<br>

### BioPython 1.81
To install [BioPython](https://biopython.org/), you first need to install Python3 package manager with:
```
sudo apt-get install python3-pip
```
You can then install BioPython by using the pip3 package manager:
```
sudo pip3 install biopython
```
<br>

### EMBOSS:6.6.0.0
[EMBOSS](https://emboss.sourceforge.net/) tool suite can be downloaded from [this]( ftp://emboss.open-bio.org/pub/EMBOSS/) FTP server. Once downloaded you will need to compile with `make`:
```
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar -xvf EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0/
./configure
make
```
Binary files will be located in the `emboss` folder. This folder contains a lot of binary files, but only `stretcher` is used by P-GRe.

### Other
P-GRe requires a Linux-system and a Python3 installation. It is known to work with CentOs 7 distribution and Python 3.10.0 and 3.6.3 versions.
