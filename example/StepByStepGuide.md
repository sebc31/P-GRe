# Exemple data Step by Step Guide
The purpose of this guide is to allow you to understand how P-GRe works and to launch it for your future analyses. For this, the data available in the `example` subfolder will be used as an example.

Descritpion of the data
=======================
The data provided for this example contains a fraction of chromosome 1 of the model plant *Arabidospis thaliana* ([TAIR10, release 56](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/)) extending from positions 3,823,000 to 3,825,500. A simplified version of the corresponding GFF file is provided, and only describes this part of the genome. The GFF file can be opened with the command `cat AtChr1_part.gff`.

| Chromosome | Version   | Type | Strat | End  |   | Sense | Frame | Attributes                                                 |
|------------|-----------|------|-------|------|---|-------|-------|------------------------------------------------------------|
| 1          | araport11 | gene | 1510  | 2231 | . | -     | .     | ID=gene:AT1G11362;Name=AT1G11362                           |
| 1          | araport11 | CDS  | 1630  | 2184 | . | -     | 0     | ID=CDS:AT1G11362.1;Parent=AT1G11362;protein_id=AT1G11362.1 |

The GFF file contains only one gene, AT1G11362, itself made up of a single CDS. Note that the gene coordinates match the coordinates on the FASTA file `AtChr1_part.fasta`, the true coordinates of this gene on chromosome 1 are 3,824,510 to 3,825,231. The interesting thing about this part of the genome of *A. thaliana* is that **~800 bases upstream of this gene there is a pseudogene similar to AT1G11362**, and it is even possible that one of these structures is a pseudogenized copy of the other. Thus, on this small fraction of the genome, the AT1G11362 gene will be used for pseudogene prediction.

Running P-GRe
=============
After [installing P-GRe](https://github.com/sebc31/P-GRe#installation), go to the `example` folder and use the following command to use P-GRe on the provided example data:
```
../P-GRe_pipeline.sh -f AtChr1_part.fasta -g AtChr1_part.gff -v
```
Here, the options `-f`, `-g` and `-v` are respectively used to specify the path to the genome FASTA file, to specify the path to the corresponding GFF file and to activate the verbose mode. **If you are using P-GRe on a cluster, you may need to provide the absolute path to the `P-GRe` main folder via the option `-d`.** The process should only take a few seconds given the small size of the dataset. The P-GRe output is divided into several parts. The first corresponds to the verification that P-GRe makes of its correct installation, i.e. the verification of the presence of scripts and dependencies. The second gives information regarding the creation of working files and folders. Finally, the third part gives information on the progress of the different stages of the P-GRe scripts. As the `-o` option was not specified, a new folder should be created in the folder where the FASTA file is located, in our case, the current `example` folder. This folder, automatically named `PGRe_` followed by the date, is divided into four subfolders, including the `res` subfolder containing the prediction results.

Looking at the results
======================
The `res` folder contains the (pseudo)coding sequence and the genomic sequence of a unique pseudogene predicted upstream of the AT1G11362 gene. More interestingly, it also contains the sequence of the protein virtually encoded by the pseudogene:
```
>pseudogene_1_151_714
MFSNVKSCLLIHIVFFLVFVVSSSARFSTEITKSEINSICTHTEANTNASLCFEFLNSSP
EIASLDFYGLTKYILTYNSRKTSDMLKQFQALVRSTTDPNAKGSYHVCAEMFDKATDCFD
DAFTTLASKDYITLEQRVVCTFDMGDGCKEELLTFTPSPQLFKDVSIVKNLSSMVLVILD
YYLKKK*
```
P-GRe uses GFFRead to generate the proteome's sequences if no proteome file is provided. The sequence of the protein encoded by the AT1G11362 gene is thus located in the working file `AtChr1_part.prot.tempfile.faa
`, itself located in the `tmp` folder. By globally aligning these two sequences using the [Needleman–Wunsch algorithm](https://www.ebi.ac.uk/Tools/psa/emboss_needle/), the similarity between these two sequences appears obvious:
```
# Length: 188
# Identity:     140/188 (74.5%)
# Similarity:   157/188 (83.5%)
# Gaps:           4/188 ( 2.1%)
# Score: 701.5

pseudogene_1_      1 MFSNVKSCLLI-HIVFFLVFVVSSSARFSTEITKSEINSICTHTEANTNA     49
                     |..|||||.|. :||||||||||:|||||||:||||||||||||:  .:|
AT1G11362          1 MCLNVKSCFLFANIVFFLVFVVSASARFSTEVTKSEINSICTHTD--VDA     48

pseudogene_1_     50 SLCFEFLNSSPEIASLDFYGLTKYILTYNSRKTSDMLKQFQALVRSTTDP     99
                     |||||||||||:||:|||||||||::||.|||.||||||||:||.|||||
AT1G11362         49 SLCFEFLNSSPQIAALDFYGLTKYLITYESRKFSDMLKQFQSLVNSTTDP     98

pseudogene_1_    100 NAKGSYHVCAEMFDKATDCFDDAFTTLASKDYITLEQRVVCTFDMGDGCK    149
                     :||||||||...|||.|.||||||..||||||||||..|.|||||..||:
AT1G11362         99 SAKGSYHVCVGTFDKGTGCFDDAFRHLASKDYITLEWSVECTFDMAAGCE    148

pseudogene_1_    150 EELLTFTPSPQLFKDVSIVKNLSSMVLVILDYYLKKK*    187
                     :||.||.|:||||||:|||||||.|.|||:..:|||. 
AT1G11362        149 DELSTFKPNPQLFKDISIVKNLSIMGLVIVKLFLKK*-    185
```
The similarity between these two sequences, and the absence of a stop codon in the sequence encoded by the pseudogene, could lead one to believe that P-GRe predicted a functional gene. To have more information on the structure of the pseudogene, the `pseudogenes.gff` file, in the `res` folder, can be checked.
| Chromosome | Version | Type               | Strat | End |   | Sense | Frame | Attributes                                                               |
|------------|---------|--------------------|-------|-----|---|-------|-------|--------------------------------------------------------------------------|
| 1          | P_GRe   | pseudogene         | 151   | 714 | . | -     | .     | ID=pseudogene_1_151_714; Parent_gene=AT1G11362                            |
| 1          | P_GRe   | pseudo_start_codon | 712   | 714 | . | -     | .     | ID=pseudogene_1_151_714.pseudo_start_codon.1; Parent=pseudogene_1_151_714 |
| 1          | P_GRe   | pseudo_CDS         | 151   | 156 | . | -     | .     | ID=pseudogene_1_151_714.pseudo_CDS1.1; Parent=pseudogene_1_151_714        |
| 1          | P_GRe   | frame_shift        | 157   | 157 | . | -     | .     | ID=pseudogene_1_151_714.frame_shift.1; Parent=pseudogene_1_151_714        |
| 1          | P_GRe   | pseudo_CDS         | 158   | 670 | . | -     | .     | ID=pseudogene_1_151_714.pseudo_CDS1.2; Parent=pseudogene_1_151_714        |
| 1          | P_GRe   | frame_shift        | 671   | 672 | . | -     | .     | ID=pseudogene_1_151_714.frame_shift.2; Parent=pseudogene_1_151_714        |
| 1          | P_GRe   | pseudo_CDS         | 673   | 714 | . | -     | .     | ID=pseudogene_1_151_714.pseudo_CDS1.3; Parent=pseudogene_1_151_714        |
| 1          | P_GRe   | pseudo_stop_codon  | 151   | 153 | . | -     | .     | ID=pseudogene_1_151_714.pseudo_stop_codon.1; Parent=pseudogene_1_151_714  |

As you can see, the single CDS is actually "broken" by two frame-shift events and divided into three sub-CDS: CDS1.1, CDS1.2 and CDS1.3. These frame-shift events most likely caused the loss of function of the pseudogene.
