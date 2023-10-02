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
![pgre-main-a\[fig1\]](docs/figs/figure_article.png)
