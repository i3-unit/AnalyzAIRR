# RepSeq

This package is for analysing high-thoughput sequencing repertoire data. It focus on clonotype data sets pre-processd by ClonotypeR, rTCR and MiXCR. Clonotype tables from each sample are concatenated, filtered, normalized and tested for differential expression between groups of samples. Diversity indices are also studied. Results could be visualized in different ways.

Contact: ph.phuong-a-gmail.com

# Description

RepSeq is an open source R package which aims to facilitate Immune Repertoire data analyzed by high-thougthput sepquencing technologies.


# Package installation

**RepSeq** depends on the following packages: *data.table, pbapply, pheatmap, 
DESeq2, Rcpp, vegan, ggplot2, naturalsort, scales, magick*. 
These above packages could be installed using the following scripts:
```r
list.pkgs <- c("data.table", "pbapply", "pheatmap", "DESeq2", "Rcpp", "vegan", "ggplot2", "naturalsort", "scales", "magick")
pkgs <- list.pkgs[!(list.pkgs %in% installed.packages()[,"Package"])]
if(length(pkgs)>0) install.packages(pkgs)
```

The latest developement package could be installed from GitHub using **devtools**.

```r
#install.packages("devtools")  # if necessary
#devtools::install_github("ph-pham/RepSeq")
library(RepSeq)
```

# Data format

The package can load the following clonotype tables:

MiXCR output:
~~~~~~~~~~~~~
cloneId	cloneCount	cloneFraction	clonalSequence	clonalSequenceQuality	allVHitsWithScore	allDHitsWithScore	allJHitsWithScore	allCHitsWithScore	allVAlignments	allDAlignments	allJAlignments	allCAlignments	nSeqFR1	minQualFR1	nSeqCDR1	minQualCDR1	nSeqFR2	minQualFR2	nSeqCDR2	minQualCDR2	nSeqFR3	minQualFR3	nSeqCDR3	minQualCDR3	nSeqFR4	minQualFR4	aaSeqFR1	aaSeqCDR1	aaSeqFR2	aaSeqCDR2	aaSeqFR3	aaSeqCDR3	aaSeqFR4	refPoints
0	55	9.70873786407767E-4	TGTGTTGTGAGCTATAACCAGGGAGGAAAGCTTATCTTC	IIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIII	TRAV8-2*01(893,1),TRAV8-2*02(888,1),TRAV8-4*01(834,2),TRAV8-4*04(834,2)		TRAJ23*01(289),TRAJ23*02(275)	270|281|304|0|11||55.0;270|280|318|0|10||50.0;270|281|304|0|11|SC274T|41.0;270|281|304|0|11|SC274T|41.0	25|52|83|12|39||135.0;25|52|83|12|39||135.0	TGTGTTGTGAGCTATAACCAGGGAGGAAAGCTTATCTTC	39	CVVSYNQGGKLIF	:::::::::0:-3:11:::::12:-5:39:::
1	44	7.766990291262136E-4	TGTGCTGTGGATAGCAACTATCAGTTAATCTGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	TRAV2*01(744,8),TRAV2*02(730,5)	TRAJ33*01(280,5)	249|260|283|0|11||55.0;249|259|284|0|10||50.0	20|46|77|7|33||130.0	TGTGCTGTGGATAGCAACTATCAGTTAATCTGG	40	CAVDSNYQLIW	:::::::::0:-3:11:::::7:0:33:::
2	41	7.237422771403354E-4	TGCATCCTGAGAGGTGGGGGTTACCAGAAAGTTACCTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	TRAV26-2*01(937,2)	TRAJ13*01(279,7),TRAJ13*02(265,7)	261|274|296|0|13||65.0	27|52|83|14|39||125.0;27|52|83|14|39||125.0	TGCATCCTGAGAGGTGGGGGTTACCAGAAAGTTACCTTT	40	CILRGGGYQKVTF	:::::::::0:-2:13:::::14:-7:39:::
3	41	7.237422771403354E-4	TGTGCCGTTTTCAACAATGACATGCGCTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIII	TRAV39*01(985,6)	TRAJ43*01(244,3)	264|272|297|0|8||40.0	25|43|74|12|30||90.0	TGTGCCGTTTTCAACAATGACATGCGCTTT	40	CAVFNNDMRF	:::::::::0:-5:8:::::12:-5:30:::
4	37	6.531332744924978E-4	TGTGCTGTGGAGGACCTTTATAACCAGGGAGGAAAGCTTATCTTC	IIIIIIIIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIII	TRAV2*01(907,5),TRAV2*02(872,4)		TRAJ23*01(299,6),TRAJ23*02(285,6)	249|267|283|0|17|DT263|73.0;249|259|284|0|10||50.0	23|52|83|16|45||145.0;23|52|83|16|45||145.0												TGTGCTGTGGAGGACCTTTATAACCAGGGAGGAAAGCTTATCTTC	39	CAVEDLYNQGGKLIF	:::::::::0:4:17:::::16:-3:45:::

~~~~~~~~~~~~~

rTCR output:
~~~~~~~~~~~~

Number of reads	Amino acid sequence	V gene	J gene	Junction nucleotide sequence	V gene end position	J gene start position	Frame	Number of stop codons	Minimum Phred	Quality
58	CVVSYNQGGKLIF	TRAV8-2*01	TRAJ23*01	TGTGTTGTGAGCTATAACCAGGGAGGAAAGCTTATCTTC	11	13	0	0	39	40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|39|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40
51	CAVDSNYQLIW	TRAV21*01	TRAJ33*01	TGTGCTGTGGATAGCAACTATCAGTTAATCTGG	9	10	0	0	40	40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40
43	CILRGGGYQKVTF	TRAV26-2*01	TRAJ13*01	TGCATCCTGAGAGGTGGGGGTTACCAGAAAGTTACCTTT	13	15	0	0	40	40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40
43	CAVFNNDMRF	TRAV39*01	TRAJ43*01	TGTGCCGTTTTCAACAATGACATGCGCTTT	8	13	0	0	40	40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40
41	CAVEDLYNQGGKLIF	TRAV2*01	TRAJ23*01	TGTGCTGTGGAGGACCTTTATAACCAGGGAGGAAAGCTTATCTTC	14	17	0	0	39	40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|39|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40
38	CAVPSPYSGAGSYQLTF	TRAV21*01	TRAJ28*01	TGTGCTGTTCCGTCGCCATACTCTGGGGCTGGGAGTTACCAACTCACTTTC	8	17	0	0	39	40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|39|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40|40

~~~~~~~~~~~~~

ClonotypeR output:
~~~~~~~~~~~~~~~~~~
cdna35	TRAV1-1	TRAJ38	115	185	SN7001332:388:H5YJNBCX2:1:1108:14316:30492	GCTGTGAGAGCCTTAGGTGGCAACAACTGTAAGCTGAT	IHHIIIIIIIIIIIIIIIHIIIIHIIIIIIIIIIIIHH	AVRALGGNNCKL
cdna35	TRAV1-1	TRAJ38	115	185	SN7001332:388:H5YJNBCX2:1:1113:12029:52603	GCTGTGAGAGCCTTAGGTGGCAACAACCGTAAGCTGAT	IIIIIHHHHIIIIIIHIHGIIIIIIIIIIIIIIHIIII	AVRALGGNNRKL
cdna35	TRAV1-1	TRAJ38	115	185	SN7001332:388:H5YJNBCX2:1:1202:12045:74829	GCTGTGAGAGCCTTAGGTGGCAACAACCGTAAGCTGAT	IIIIIIIIIIIIIIIHIHGIIIIIIIIIIIIIIIIIII	AVRALGGNNRKL
cdna35	TRAV1-1	TRAJ38	115	185	SN7001332:388:H5YJNBCX2:1:2204:11922:37401	GCTGTGAGAGCCTTAGGTGGCAACAACCGTAAGCTGAT	IIIIIIHIIIIIHIIIIHHIIIIIIIIIIIIIIIHIII	AVRALGGNNRKL
cdna35	TRAV1-1	TRAJ38	176	115	SN7001332:388:H5YJNBCX2:2:1216:1394:6329	GCTGTCAGAATGCTGGCAACAACCGTAAGCTGAT	IIHIHHIIIHIIIHFIIIIIIIIIIHIIIIIIII	AVRMLATTVS*
cdna35	TRAV1-1	TRAJ38	179	211	SN7001332:388:H5YJNBCX2:1:1203:1842:95471	GCTGTTCAATGCTGGCAACAACCGTAAGCTGAT	IIIIIHHGHIIIHEIIIIIIHIIIHHIIIIHFH	AVQCWQQP*AD
cdna35	TRAV1-1	TRAJ38	183	212	SN7001332:388:H5YJNBCX2:1:1208:9879:19838	GCTGTTCAATGCTGGCAACAACCGTAAGCTGAT	IIIIIHIIHIGIHGIIIIIGIIIIIGIHHEG@E	AVQCWQQP*AD
cdna35	TRAV1-1	TRAJ38	183	212	SN7001332:388:H5YJNBCX2:2:2115:2594:20108	GCTGTTCAATGCTGGCAACAACCGTAAGCTGAT	IIIIIIIIIIIIIIIIIIIIIIIIIIHIIIIIH	AVQCWQQP*AD
cdna35	TRAV1-1	TRAJ38	183	212	SN7001332:388:H5YJNBCX2:2:1108:10381:30399	GCTGCCCGATGCTGGCAACAACCGTAAGCTGAT	IIHIIIIIIIIIHGIIIIIIIIIIIIHIIIIIG	AARCWQQP*AD
cdna35	TRAV1-1	TRAJ38	183	212	SN7001332:388:H5YJNBCX2:2:2105:4465:27811	GCTGCCCGATGCTGGCAACAACCGTAAGCTGAT	HIIIHIIIIIIIIHIIIIIIIIIIIIIIIIIII	AARCWQQP*AD
~~~~~~~~~~~~~~~~~~


# The RepSeqExperiment object

The RepSeqExperiment object is a **R** **S4** class which stores the clonotype tables in long format along with experimantal data. 

## Anatomy of a RepSeqExperiment

The ```RepSeqExperiment``` class contains four slots. 
* **assayData** a **data.table** that contains clonotype tables in long format with the following columns: 
  * **lib**: sample name
  * **V**: V-gene nomenclature
  * **J**: J-gene nomenclature 
  * **CDR3aa**: peptide chain 
  * **CDR3dna**:  
  * **VpJ**: clonotype defined as the combination of V-gene, peptide chain & J-gene 
  * **VJ**: combination of V & J-genes 
  * **score**: alignment score 
  * **count**: clonotype count (count of VpJ); 
* **sampleData** a **data.frame** containing information related to samples. Sample names are stored in ```rownames``` and must be found in the column **lib** of the slot **assay** (with the same order);
* **metaData** a **list** containing all other meta data related to experiment, empty by default;
* **History** a **data.frame** containing history of treatment applied to the RepSeqExperiment object.

## Getter & setter
...

For more details 

```r
?RepSeqExperiment
```

## Example of data in **assayData** and in **sampleData** slots

assayData slot:
~~~~~~~

             lib     V       J             CDR3aa                                                CDR3dna                              VpJ            VJ score count
      1: cdna100 TRBV1 TRBJ2-5        CTSEAEETQYF                      TGCACCAGCGAAGCGGAAGAGACCCAGTACTTC        TRBV1 CTSEAEETQYF TRBJ2-5 TRBV1 TRBJ2-5    37     1
      2: cdna100 TRBV1 TRBJ2-7   CTSNWGLAGGTYEQYF       TGCACCAGCAACTGGGGGCTAGCGGGGGGGACCTACGAGCAGTACTTC   TRBV1 CTSNWGLAGGTYEQYF TRBJ2-7 TRBV1 TRBJ2-7    27     1
      3: cdna100 TRBV1 TRBJ1-1     CTSSPSGSQGNLIF             TGCACCAGCAGCCCCTCGGGAAGCCAAGGAAATCTCATCTTT     TRBV1 CTSSPSGSQGNLIF TRBJ1-1 TRBV1 TRBJ1-1    39     1
      4: cdna100 TRBV1 TRBJ2-1    CTSSQACSSYNEQFF          TGCACCAGCAGCCAAGCCTGCAGCTCCTACAATGAGCAGTTCTTC    TRBV1 CTSSQACSSYNEQFF TRBJ2-1 TRBV1 TRBJ2-1    34     1
      5: cdna100 TRBV1 TRBJ1-1 CTSSQDGRDRKGNTEAFF TGCACCAGCAGCCAAGATGGTCGGGACAGGAAAGGGAACACTGAAGCTTTCTTT TRBV1 CTSSQDGRDRKGNTEAFF TRBJ1-1 TRBV1 TRBJ1-1    36     1
     ---                                                                                                                                                           
1772579:  cdna99 TRBV9 TRBJ2-3    GASSVSGGVRDTQYF          GGTGCCAGCAGCGTGAGCGGGGGGGTCAGAGATACGCAGTATTTT    TRBV9 GASSVSGGVRDTQYF TRBJ2-3 TRBV9 TRBJ2-3    38     1
1772580:  cdna99 TRBV9 TRBJ1-4     RASSRTVTNEKLFF             CGTGCCAGCAGCCGGACGGTTACTAATGAAAAACTGTTTTTT     TRBV9 RASSRTVTNEKLFF TRBJ1-4 TRBV9 TRBJ1-4    38     1
1772581:  cdna99 TRBV9 TRBJ2-7    SASSPRDRGIHEQYF          AGTGCCAGCAGCCCCCGGGACAGGGGTATTCACGAGCAGTACTTC    TRBV9 SASSPRDRGIHEQYF TRBJ2-7 TRBV9 TRBJ2-7    36     1
1772582:  cdna99 TRBV9 TRBJ1-2    YASSGRVSVDYGYTF          TATGCCAGCAGCGGCAGGGTCTCAGTGGACTATGGCTACACCTTC    TRBV9 YASSGRVSVDYGYTF TRBJ1-2 TRBV9 TRBJ1-2    38     1
1772583:  cdna99 TRBV9 TRBJ2-3     YASSVGTYTDTQYF             TATGCCAGCAGCGTCGGGACGTACACAGATACGCAGTATTTT     TRBV9 YASSVGTYTDTQYF TRBJ2-3 TRBV9 TRBJ2-3    38     2

~~~~~~~

sampleData slot:
~~~~~~~

                Sample   Cell  Organ Patient nReads   VpJ  V  J  VJ
cdna100 cdna100.tsv.gz  nTeff    pLN     p31  56445 34606 58 13 613
cdna101 cdna101.tsv.gz amTeff    pLN     p31  35159 19694 56 13 572
cdna102 cdna102.tsv.gz amTreg    pLN     p31  43550  8160 53 13 529
cdna103 cdna103.tsv.gz  nTreg spleen     p32  72796  7906 51 13 513
cdna104 cdna104.tsv.gz  nTeff spleen     p32  83631 24795 55 13 587
cdna105 cdna105.tsv.gz amTeff spleen     p32  82018 20744 53 13 581
cdna106 cdna106.tsv.gz amTreg spleen     p32  34961  3240 50 13 453
...
cdna94   cdna94.tsv.gz amTreg spleen     p31  14162  5931 50 13 495
cdna95   cdna95.tsv.gz  nTreg    mLN     p31  51652  2671 48 13 420
cdna96   cdna96.tsv.gz  nTeff    mLN     p31  42961 17784 56 13 578
cdna97   cdna97.tsv.gz amTeff    mLN     p31  31177 12666 55 13 552
cdna98   cdna98.tsv.gz amTreg    mLN     p31  52426  7250 53 13 522
cdna99   cdna99.tsv.gz  nTreg    pLN     p31  69638  7433 53 13 518

~~~~~~~

# Getting started 

The format of aligner output file is usually tab-delimited, one tab-delimited file for each sample. 


```r
# load library in memory
library(RepSeq)
# list of aligner output files (suppose to be stored in /MiXCR_output) 
inputFolder <- list.files("/MiXCR_output/", full.name = TRUE, pattern = ".tsv")
# Create an object of class RepSeqExperiment using the wrapper function readClonotypeSet
datatab <- readClonotypeSet(inputFolder, cores=2L, aligner="MiXCR", chain="A", sampleinfo=NULL, keep.ambiguous=FALSE, keep.unproductive=FALSE, aa.th=8) 
```

If you have **assayData** data1 and **sampleData** data2:
```r
datatab <- methods::new("RepSeqExperiment", assayData=data1, sampleData=data2, metaData=list(), History=data.frame())
```

# Web-based interface

A web-based interface was implemented using the Shiny tools which source code can be found at: https://github.com/ph-pham/DiversiTR

 

