This folder contains data and results of the cis-regulatory elements of de novo discovery analysis of 47 Bdhn gene promoters from four Brachypodium species using Rsat::plants.
The starting dataset is in file "01_Bdhn_Brachypodium_promoters.fasta" containing 47 Bdhn promoter sequences with a length of -500bp to +200 bp around the initial of the 
transcription site. 
Analyses were performed by http://rsat.eead.csic.es/plants/

1. Peak-motif tool was used to detect over-represented motif at Bdhn promoter sequences.Brachypodium genomes from the species sampled were used as background.
    The analysis was run four time, using one genome background each time (B. distachyon Bd21 v3.0.46, B. stacei ABR114 v1.1.JGI, B. hybridum ABR113 v1.1.JGI, 
    B. sylvaticum Ain1v1.1.JGI) 
    Significance of the motif discovery were evaluated using negative controls, consisting in 47 random promoter sequences from each background genome analysed with the same 
    conditions as the Bdhn promoters.

  Results from the 4 analyses are in file "02_Brachypodium_Bdhn_Peak.tf" . This file is the input in the next step.

2. Matrix-clustering tool. Significance motif selected, based on their k-mer and number of sites were clustered to avoid redundant motif. 
	The result is in file "03_CRE_Bdhn.tf"

3. Matrix-scan tool. Consensous motif from matrix-clustering result were scanned along the 47 Bdhn promoter sequences in order to locate potential CREs with a maximum threshold of
9 nucleotides on the median length of the 3 selected motif.
Results are in gff3 files from each Brachypodium at folder: 04_Brachypodium_Bdhn_scan. Frecuency results are in table 3 
