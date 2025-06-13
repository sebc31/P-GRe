# P-GRe v.0.2

Contacts and authors
====================
Sébastien Cabanac, PhD student, sebastien.cabanac@univ-tlse3.fr

Christophe Dunand, professor at Université Paul Sabatier Toulouse 3, christophe.dunand@univ-tlse3.fr

Catherine Mathé, associate professor at Université Paul Sabatier Toulouse 3, catherine.mathe-dehais@univ-tlse3.fr

Contents
========

- [About P-GRe](#about-p-gre)
- [Input requirements](#requirements-and-input)
- [Installation requirements](#installation)
- [Running P-GRe](#running-p-gre)
- [Citing P-GRe and associated software](#citing-p-gre-and-associated-software)
- [Licence](#licence)

About P-GRe
===========
Pseudogenes are genomic sequences with homology to functional genes but that harbor deleterious mutations, such as loss of the start codon, loss of coding sequence, gain of stop or frame-shifts. The goal of PseudoGene REtriever (P-GRe) is to find the position of pseudogenes on a genome, as well as to infer their structures in pseudo-CDSs. P-GRe aims to be more user-friendly, with a limited number of dependencies, ease of use and total automaticity, while producing qualitative results and having greater sensitivity than other software with the same goal. P-GRe relies on miniprot to align user-provided protein sequences on the genome and then filter the overlapping results. P-GRe also categorizes all predictions into the three main categories of pseudogenes (unitary, duplicated or processed).

Input requirements
======================
P-GRe requires the genome in FASTA format and the structural annotation of the genes present on this genome, as well as the set of protein sequences encoded by these genes. Optionally, a second set of protein sequences from other organisms can be provided to increase the sensitivity and to annotate unitary pseudogenes. Input files must follow some common formatting rules: the genome and protein sequence sets must be in FASTA format, while the structural annotation must be in GFF3 format. <br/><br/>The GFF3 file must contain the "CDS" information in the type field (3rd column). More importantly, on CDS lines, the attributes field (9th column) must contain the predefined tag "Parent" and the value associated with this tag must be identical to the identifier of the corresponding sequence in the protein sequences FASTA file. This is necessary because P-GRe uses this tag to establish the correspondence between the structural annotation of each transcript and the corresponding protein sequence. Comparing the structure of the transcripts to those obtained on the pseudogenes then allows P-GRe to identify intron loss events, which is needed for the categorization of pseudogenes.

Installation requirements
============
P-GRe is designed to run on a Unix system and requires some dependencies. Before running P-GRe, make sure you have installed the following software: `bedtools` [(arq5x/bedtools2)](https://github.com/arq5x/bedtools2), `miniprot` [(lh3/miniprot)](https://github.com/lh3/miniprot) and `diamond` [(bbuchfink/diamond)](https://github.com/bbuchfink/diamond). P-GRe also requires a `python3` interpreter.

Running P-GRe
=============
To predict the structure of all pseudogenes in an organism's genome (`genome_A.fasta`) from that organism's protein sequences (`proteins_A.fasta`), a typical command line is:<br/><br/>
`PGRe.sh -f genome_A.fasta -g genome_A_annotation.gff -p proteins_A.fasta`<br/>

To add a set of protein sequences from multiple organisms (`proteins_B.fasta`) to increase sensitivity and predict unitary pseudogenes, the typical command line would look like this:<br/><br/>
`PGRe.sh -f genome_A.fasta -g genome_A_annotation.gff -p proteins_A.fasta -u proteins_B.fasta`<br/>

Other options include multithreading and output directory options. The main output files are the structural annotation of the pseudogenes (`PGRe.gff`), and the sequences of peptides encoded by the predicted pseudogenes (`pseudogene_protein.fasta`).

Citing P-GRe
============
To be added.

Licence
=======
All source code contained in the `scripts` folder and the `P-GRe_pipeline.sh` file are under the [Artistic Licence 2.0](https://opensource.org/license/artistic-2-0/).
