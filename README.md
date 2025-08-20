# VCF-Compare-and-Consensus
A tool to compare 2 or more single or multi-sample VCF files created with different callers, and generate a consensus VCF and a table summarizing variants.  
__________________________________________________________________________________________________________________________________________________________

Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT:

 This program takes 2 or 3 single or multi-sample VCFs and performs a comparison. A new VCF is created that contains the variants found in common to all input files.  A summary table of shared and unique variants is also produced.  This program was tested using VCFs from illumina short-read data using a Linux workstation. The program supports compressed and uncompressed VCF files.

RATIONAL: 

Different programs used to call sequence variants use different algorithms and thus have different false positive and false negative errors. One approach to prioritizing potentially true variants is to generate VCFs of the same samples with different callers and select those variants found by all callers for further evaluation. While many variant calling tools exist, this program is limited to the evaluation of no more than 3 VCFs based on previous experience (1,2). 

PREREQUISITES:

1) zcat (usually available with gzip) is required for compressed inputs 
2) bgzip is required for compressed outputs
3) C++ compiler (for example GCC/G++) 

INSTALLATION:

This program should work on all systems. Download the .cpp file and compile according to your system. For Linux Ubuntu, compile with g++ (g++ -o intersect2 vcf_compare_consensus_V1_2.cpp -std=c++11).


COMMAND LINE OPTIONS:  

<b>Basic usage (uncompressed)</b>

./vcf_compare file1.vcf file2.vcf

<b>With compression</b>

./vcf_compare --compress file1.vcf.gz file2.vcf.gz file3.vcf.gz

<b>With mixed input formats</b>

./vcf_compare --compress file1.vcf file2.vcf.gz file3.vcf


REFERENCES

1) Gupta P, Reddaiah B, Salava H, Upadhyaya P, Tyagi K, Sarma S, Datta S, Malhotra B, Thomas S, Sunkum A, Devulapalli S, Till BJ, Sreelakshmi Y, Sharma R. Next-generation sequencing (NGS)-based identification of induced mutations in a doubly mutagenized tomato (Solanum lycopersicum) population. Plant J Cell Mol Biol. 2017 Nov;92(3):495–508.

2) Hawliczek A, Bolibok L, Tofil K, Borzęcka E, Jankowicz-Cieślak J, Gawroński P, Kral A, Till BJ, Bolibok-Brągoszewska H. Deep sampling and pooled amplicon sequencing reveals hidden genic variation in heterogeneous rye accessions. BMC Genomics. 2020 Nov 30;21(1):845. 
