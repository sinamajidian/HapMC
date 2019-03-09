# HapMC


We apply matrix completion methods for haplotype assembly from NGS reads to develop the new HapSVT, HapNuc, and HapOPT algorithms.

## Cite us
```
Sina Majidian, Mohammad Hossein Kahaei, "NGS based haplotype assembly using matrix completion", Submitted to PLOS ONE, 2018
```

## Input
Our developed program works with the fragment file used in AltHap and SDhaP with the following format:
```
Number of reads
Number of columns 
"Number of parts" "Read ID" "index of the first allele of the 1st part" "Alleles of 1st part" "index of first allele of the 2nd part" "Alleles of the 2nd part" ... "Quality scores of alleles of all parts" 
...
```
An example of fragment file:
```
3
9
1 chr_1 1 010 @@@@@
2 chr_2 3 011 7 010 @@@@@@ 
1 chr_3 1 01 7 010 @@@@@@ 
```
When you have Bam and VCF files, we suggest to use ExtractHAIRS as the following command-line:
```
extractHAIRS --bam reads.sorted.bam --VCF variants.VCF --out fragment.txt
```
This works for diploids. If you want to use it for polyploids, artificially convert the genotypes (GP tags) to "0/1" using sed command.


## Output
The output of algorithm is a text file named as `Reconstructed_Haplotype.txt`. For the all-heterozygous diploid case, it is like this:

```
Block 1 3 
1 0
2 1
3 0
```

in which the first column is the index of variant in the VCF file and the second column is the allele at that SNP. It is obvious that another haplotype can be found by flipping the second column. For the case where both heterozygous and homozygous variants exist, the output file contains three columns.

# How to use
One can run the MATLAB code using the following command-line:
```
matlab -r "HapMC(fragment_file,Hap_algorithm)"
```
Hap_algorithm can be either 'O', 'S' or 'N' corresponding to 'HapOPT', 'HapSVT' or 'HapNuc', respectively,

Test use
Here, we can easily test the example described in Table 2 of our paper. The input file is provided in the data folder.
```
git clone https://github.com/smajidian/HapMC.git
cd HapMC
matlab -r "HapMC('data/fragment_sample.txt','O');exit";
cat Reconstructed_Haplotype.txt
```
## File description
The MATLAB code convert_frag_mat.m can be used to convert a fragment file to a matrix of MATLAB format with .mat extension. first_block_extractor.m is used within the MATLAB code, HapMC.m, to extract overlapped reads which are considered as a connected read block. To convert the text file of the simulated data, we used two codes in R programming language:

convert_simulated_data_to_mat-hap.R and convert_simulated_data_to_mat-read.R. 

To convert the HapOPT’s output to HapCUT format, one can use convert_hapcut.py.

