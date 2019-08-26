# SplicingCompass

SplicingCompass: differential splicing detection using RNA-Seq data.

Moritz Aschoff(1,2), Agnes Hotz-Wagenblatt(1), Karl-Heinz Glatting(1), Matthias Fischer(4), Roland Eils(2,3) and Rainer König(2,3,5,6)

(1) Bioinformatics "HUSAR", Genomics Proteomics Core Facility, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 580, 69120 Heidelberg, Germany
(2) Department of Theoretical Bioinformatics, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 280, 69120 Heidelberg, Germany
(3) Department of Bioinformatics and Functional Genomics, Institute of Pharmacy and Molecular Biotechnology, and Bioquant, University of Heidelberg, Im Neuenheimer Feld 267, 69120 Heidelberg, Germany
(4) Department of Pediatric Oncology and Hematology and Center of Molecular Medicine Cologne, University Children's Hospital, Joseph-Stelzmann-Strasse 9, D-50924 Cologne, Germany
(5) Leibniz Institute for Natural Product Research and Infection Biology, Hans Knöll Institute Jena, Beutenbergstrasse 11a, 07745 Jena, Germany
(6) Integrated Research and Treatment Center, Center for Sepsis Control and Care (CSCC), Jena University Hospital, Erlanger Allee 101, 07747 Jena, Germany

corresponding authors: m.aschoff@dkfz.de and Rainer.Koenig@uni-jena.de

Paper Abstract


Motivation: Alternative splicing is central for cellular processes and substantially increases transcriptome and proteome diversity. Aberrant splicing events often have pathological consequences and are associated with various diseases and cancer types. The emergence of next-generation RNA sequencing (RNA-seq) provides an exciting new technology to analyse alternative splicing on a large scale. However, algorithms that enable the analysis of alternative splicing from short-read sequencing are not fully established yet and there are still no standard solutions available for a variety of data analysis tasks.
Results: We present a new method and software to predict genes that are differentially spliced between two different conditions using RNA-seq data. Our method uses geometric angles between the high dimensional vectors of exon read counts. With this, differential splicing can be detected even if the splicing events are composed of higher complexity and involve previously unknown splicing patterns. We applied our approach to two case studies including neuroblastoma tumour data with favourable and unfavourable clinical courses. We show the validity of our predictions as well as the applicability of our method in the context of patient clustering. We verified our predictions by several methods including simulated experiments and complementary in silico analyses. We found a significant number of exons with specific regulatory splicing factor motifs for predicted genes and a substantial number of publications linking those genes to alternative splicing. Furthermore, we could successfully exploit splicing information to cluster tissues and patients. Finally, we found additional evidence of splicing diversity for many predicted genes in normalized read coverage plots and in reads that span exon-exon junctions.

## Software

SplicingCompass is licensed under the GNU General Public License and available as a package in GNU R.

## Downloads

    Installation guide and tutorial
    GNU R package SplicingCompass Version 2.0.1: SplicingCompass2.1_1.0.tar.gz (licensed under GNU GENERAL PUBLIC LICENSE Version 3) (VersionChangelog)
    GTF annotation of union transcripts based on the UCSC CCDS annotation (please download and see tutorial)
    CountTable object file comparing brain vs liver samples on chromosome 8 (saved R object for demonstration purposes, please download and see tutorial)
