# Brachypodium dehydrins
> This repository contains data and instructions to run JAVA programs and scripts used in different analyses of Brachypodium dehydrins (Bdhn) data and other phenotypic and climate niche data included in the paper *"Evolution and functional dynamics of dehydrins in model Brachypodium grasses"* coauthored by Maria Ángeles Decena, Sergio Galvez-Rojas, Federico Agostini, Rubén Sancho, Bruno Contreras-Moreira, David L. Des Marais, Pilar Hernández and Pilar Catalán.

## Table of Contents
* [Brachypodium_Bdhn_sequences_structure_chromosomal mapping](#brachypodium_bdhn_sequences_structure_chromosomal-mapping-)
* [Brachypodium_Bdhn_genes_CREs](#brachypodium_bdhn_genes_cres)
* [Brachypodium_Bdhn_MSAs_phylogenies](#brachypodium_bdhn_msas_phylogenies)
* [Bdistachyon_Bdhn expression_drought_phenotypes_climate niche_statistics_phylogenetic signal](#bdistachyon_bdhn-expression_drought_phenotypes_climate-niche_statistics_phylogenetic--signal)
* [Supplementary_Data](#supplementary_data)


## Brachypodium_Bdhn_sequences_structure_chromosomal mapping ![](https://img.shields.io/badge/Code-Java&nbsp;11-informational?style=plastic&logo=Java&logoColor=white&color=2bbc8a)
This folder contains short Java programs to deal with *Brachypodium distachyon* ecotypes. 
This project contains three programs to extract information from data files extracted from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). In particular:
* **es.uma.keyword.LookForKeyword**: Allows iterating through the fasta and gff3 files associated to several species and extract information related to a family of genes. Before working with this package, you have to download some files from Phytozome, configure some lines of source Java files and recompile. The program can be executed as is by means of the sample files already provided in the folder `"genomes/Phytozome/PhytozomeV13"`.
* **es.uma.html.GenerateSVG**: Takes as input the previously extracted data and generates an HTML file with the distribution of the genes. This file is specifically prepared to work with the 10 dehydrins found in *Brachypodium distachyon*. This program can be executed as is.
* **es.uma.motif.HTMLDecorator**: Iterates over a fasta file and converts it into an HTML file with highlighted motives. This program can be executed as is.
These programs are not intented to simply execute them but, instead, the code must be modified slightly if any developer wants to adapt to their needs. In any case, you need to clone this project by using the next command from any directory of your own:
```
git clone https://github.com/Bioflora/Brachypodium_dehydrins
cd Brachypodium_dehydrins/01_Brachypodium_Bdhn_sequences_structure_chromosomal_mapping
```
In the next sections are explained how to work with them.
## es.uma.keyword.LookForKeyword
This program uses as input an external structure of files that the user has to download from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). This has been tested in a case of use in which 54 ecotypes of *Brachypodium distachyon* have been downloaded. Before compilation, you have to change the line:
```Java
    private static final String DIR_NAME = "genomes/Phytozome/PhytozomeV13";
```
to adapt to your own directory that contains the files taken from Phytozome. In particular, you have to download the files with extension `.annotation_info.txt`, `.gene.gff3.gz`, `.protein_primaryTranscriptOnly.fa.gz` and `genes.fasta`. Please, note that the file `genes.fasta` is not in Phytozome and it must be generated by downloading also the `assembly` file from Phtyzome and executing additional programs and applications to extract genes from it (e.g., `gffread -x <out.fasta> -g <genome.fasta> <annotation.gff>`). In addition, the file `genomes/Order.txt` contains the order in which you want to generate the data on the extracted genes. Finally, the genes extracted are those that verify the sentences:
```Java
if (line.toUpperCase().contains("PF00257") || 
    line.toUpperCase().contains("DEHYDRIN")) { ...
```
By default this refers to the LEA2 gene family of Dehydrins whose PFAM code is PF00257. You have to modify these lines to adapt them to your needs. Genes are marked with different shapes according to the text found in it; this is done through the sentence:
```
String shape = line.toUpperCase().contains("PF00257")?"RECT":"TRIANGLE";
```
so a gene will be drawn as a square (PFAM code is found) or a triangle (PFAM code not found but, instead, a literal keyword) depending on its annotation in Phytozome.
Once these changes have been performed, compilation is carried by executing the next commands:
```
cd src
javac es/uma/keyword/LookForKeyword.java
```
and executed using 
```
cd ..
java -cp ./src es.uma.keyword.LookForKeyword
```
The result is displayed by console.
The file `datafiles/dehydrins_ecotypes.data` contains an example of the output of this program once executed over 54 ecotypes of *Brachypodium distachyon* downloaded from Phytozome. In addition, it is also given a file `datafiles/Order.txt` that refers to these 54 ecotypes.
## es.uma.html.GenerateSVG
This program uses the data extracted from the previous one (es.uma.keyword.LookForKeyword). As an example, a file `datafiles/dehydrins_ecotypes.data` is provided so you could run it from scratch. Compilation is carried by executing the next commands:
```
cd scr
javac es/uma/html/GenerateSVG.java
```
and executed using 
```
cd ..
java -cp ./src es.uma.html.GenerateSVG
```
This will display on console a trace of the genes processed and, in addition, a file `result.html` will be generated in the root directory; the order of the 54 ecotypes is also give by the file `datafiles/Order.txt`. Such a file has the next appearance:
![Dehydrins](01_Brachypodium_Bdhn_sequences_structure_chromosomal_mapping/images/DHN_Ecotypes.png?raw=true "Dehydrins")
## es.uma.motif.HTMLDecorator
This is a more general program very useful to colour motives in fasta files and to show the results in HTML. As an example, a file `datafiles/all_varieties_DHN_peptide.fasta` is provided so you could run it from scratch. It searches for 6 different motives, as can be seen in the next lines of the code:
```Java
    	StringBuilder sb =
    	process("datafiles/all_varieties_DHN_peptide.fasta",
    	                new PerfectMatchArgument[]{
    	                    new PerfectMatchArgument("EKKGIMDKIKEKLPG", "KSegment", 4),
    	                    new PerfectMatchArgument("VDEYGNP", "YSegment", 3),
    	                    new PerfectMatchArgument("EDDGQGR", "InterSegment", 2),
    	                    new PerfectMatchArgument("KKDKKKKKEKK", "NLSSegment", 2),
    	                    new PerfectMatchArgument("DRGLFDKFIGKK", "FSegment", 4), // For Brachypodium
    	                },
    	                new PatternArgument[]{
    	                    new PatternArgument(Pattern.compile("SSSS+"), "SSegment"),
    	                });
```
A search can be carried out by means of a regular Expression (`PatternArgument`) or by a perfect match of a given String (`PerfectMatchArgument`). In this last case, a third argument is added in order to specify how many mismatches we allow at most. The second argument of both `PatternArgument` and `PerfectMatchArgument` is a CSS class name used to colour the found motives. These CSS classes must appear in the file `resources/hubVarietiesDHN.html`; this HTML file is used as a skeleton to contain the result.
Compilation is carried by executing the next commands:
```
cd scr
javac es/uma/motif/HTMLDecorator.java
```
and executed using 
```
cd ..
java -cp ./src es.uma.motif.HTMLDecorator
```
This will display on console a trace of the motives found in each genes processed and, in addition, a file `DHNSegments.html` will be generated in the root directory. Such a file has the next appearance:
![Dehydrins](01_Brachypodium_Bdhn_sequences_structure_chromosomal_mapping/images/DHN_Segments.png?raw=true "Dehydrins")

## Brachypodium_Bdhn_genes_CREs
- 01_Bdhn_Brachypodium_promoters.fasta: Bdhn promoter sequences, with a window of  -500bp to +200pb around the inferred Transcription Start Site (TSS).
- 02_Brachypodium_Bdhn_Peak.tf: output of Peak-motif in TRANSFAC format.
- 03_CRE_Bdhn.tf: final clusters ID from RSAT matrix-clustering
- 04_Brachypodium_Bdhn_scan: gff3 files from RSAT matrix-scan
- README file describing how to carry this analysis

## Brachypodium_Bdhn_MSAs_phylogenies
This folder contains the alignments of the 10 Bdhn genes retrieved from the genomes of 4 Brachypodium species (*B.distachyon*, B.stacei, B.hybridum* subgenomes D and S, and *B. sylvaticum*) plus 5 outgroup species (*A. tauschii, H. vulgare, S. bicolor, O. sativa* , and *Z. mays*). Additionally, it also contains the alignments of the 10 Bdhn genes retrieved from 54 *B. distachyon* ecotypes. 
These files were used in the phylogenetic analyses of the dehydrin genes. 
- READMe file: description of multiple sequence alignment (MSA) folders.
This folder contains:
1.	Brachypodium.spp_Bdhn_MSAs_phylogenies: 
	a folder containing a global Bdhn alignment of the 10 Bdhn genes (exon+intron sequences) for the 4 Brachypodium species studied plus 5 outgropus species and separate
	MSAs for each Bdhn gene. These data sets were used to describe the Bdhn genes in Brachypodium and to perform the phylogenetic analysis. The global Bdhn alignment was
	used to compute the Brachypodium Bdhn ML tree shown in Fig. 3. 
	Additionally, there is another folder (Brachypodium_Bdhn_Annotation)  containing annotation files for each Bdhn gene.
2.	B.distachyon_MSAs:
	a folder containing alignments (exon+intron sequences) of each of the 10 Bdhn genes for the 54 *B. distachyon* ecotypes plus outgroup *B.stacei*. These files were used
	to perform phylogenetic analyses. The resulting Bdhn ML gene trees are shown in Supplementary Fig. S4. 
	Additionally, there is another folder (Bdistachyon_Bdhn_Anotations)  containing annotation files for each Bdhn gene and a B.distachyon_Bdhns_concatenated_MSA.fna, an
	aligment of 6 Bdhn genes (*Bdhn1, Bdhn2, Bdhn3, Bdhn6, Bdhn7, Bdhn8*) that showed a congruent topology for the 54 *B.distachyon* ecotypes. This file was used to compute
	the *B. distachyon* Bdhn tree shown in Supplementary Fig. 5Sa.

## Bdistachyon_Bdhn expression_drought_phenotypes_climate niche_statistics_phylogenetic  signal 
This folder contains five subfolders with input files and Rscripts used in different analyses of Bdhn expression, drought-related phenotypic traits, climate niche data and phylogenetic signals in *B. distachyon* lines. 
-	 01_Bdhn_Drough&TemperatureExpressions
This folder contains 12 input files of Bdhn1a, Bdhn2, Bdhn3 and Bdhn7 gene expressions under drought (watered-W vs dry-D) and temperature (cold-C vs hot-H) conditions and their combination (CD, CW, HD, HW) across 32 ecotypes of *B. distachyon* and 3 R scrips used to compute boxplots and pairwise Wilcoxon tests between conditions in combined and separate treatments. The 3 output files All_Bdhn_CD_CW_HD_HW_v3.png, All_Bdhn_D_W_v3.png and All_Bdhn_C_H_v3.png showing the results of these analyses correspond to Supplementary Figs S7a,b, and c, respectively. 
-	02_BdhnExpressionStatsDEtests
This folder contains two input csv files: Bdistachyon_Bdhn_TPM_W&D.csv with Bdhn1a, Bdhn2, Bdhn3 and Bdhn4 expression data (transcripts per million, TPM) for 32 B. distachyon ecotypes under watered (W) and drought (D) conditions with up to 4 replicas per accession and condition and Bdistachyon_Bdhn_Phenotypictraits_W&D.csv   containing the previous Bdhn1a, Bdhn2, Bdhn3 and Bdhn4 gene expression data plus data for 12 drought-related phenotypic traits [(leaf_rwc: relative water content in leaf; leaf_wc: water content in leaf; lma: leaf mass per area; pro: leaf proline content; abvrgd: above ground biomass; blwgrd: below ground biomass; ttlmass: total plant mass; rmr: root mass ratio; delta13c: carbon isotope, a proxy for lifetime integrated WUE; leafc: leaf carbon content; leafn: leaf nitrogen content; cn: leaf carbon/nitrogen ratio)] for the same set of B. distachyon ecotypes and up to 4 replicas per ecotype and condition.
It also includes the file Bdistachyon_Statistics.R with the R scripts used to compute basic statistics and boxplots for Bdhn gene expressions under W and D conditions per ecotype, pairwise Wilcoxon tests of DEs (W vs D) per ecotype with p-values adjusted for multiple comparison (BH method), Kruskal-Wallis rank test for each group (Bdhn1a, Bdhn2, Bdhn3 and Bdhn4, W vs D) and Tukey tests for multiple comparisons of W and D expressions across ecotypes.
-	03_Bdhn&PhenotraitsRegressions
This folder contains the input file "Bdistachyon_Bdhn_PhenotypicTraitsW&D.csv”   with data on the 4 Bdhn genes (Bdhn1a, Bdhn2, Bdhn3, Bdhn7) expression levels (TPM) and 12 drought-related phenotypic traits values under W and D conditions across the 32 B. distachyon ecotypes with up to 4 replicas per ecotype and the the BdistachyonRegressionsRscript file with commands to perform the linear model regression analyses   of pairwise Bdhn expressions and pairwise Bdhn expressions and phenotypic trait changes.  
-	04_ClimateNichePCA
This folder contains the input file BdistachyonLongLatAltClimate.csv showing geographic (longitude, latitude, altitude) and climate data (19 Bioclim temperature and precipitation variables) for the *B. distachyon* ecotypes and the Bdistachyon_Climate_PCA.R file with the commands to compute a principal component analysis (PCA) of the 19 climate variables and to extract the PCA1 values for the ecotypes. The results of these analyses are shown in Supplementary Fig. S6 and Supplementary Table 3. 
-	05_PhylogeneticSignal 
This folder contains 10 input files of averaged expression values (TPM) of the 4 Bdhn genes under W and D conditions, averaged values of the 12 drought-related phenotypic traits under W and D conditions and climate niche PCA1 values for the B. distachyon ecotypes, 2 input tree files corresponding to the B.distachyonBdhntree and the BdistachyonSpeciesTree, and the Brachypodium_Phylosig.R file with commands to compute the phyloheatmaps of the three sets of variables (Bdhn expressions, phenotypic traits changes, climate PCA1 values) on each phylogenetic tree and to perform phylogenetic signal tests using Blomberg’s K and Pagel’s lambda tests with the phylosig function of the package *phytools*. The results of the phylogenetic signal tests on the B. distachyon Bdhn tree are shown in table 5 and Fig. 6 and those on the B. distachyon species tree in Supplementary Table 12 and Supplementary Fig. S13. 

## Supplementary_Data
This folder contains the Supplementary Tables, Supplementary Figures and Supplementary Materials of the article. 
