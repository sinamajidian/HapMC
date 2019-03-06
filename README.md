# HapMC



We apply matrix completion methods for haplotype assembly from NGS reads to develop the new HapSVT, HapNuc, and HapOPT algorithms.

## Cite us

```
Sina Majidian, Mohammad Hossein Kahaei, "NGS based haplotype assembly using matrix completion",
Submitted to PLOS ONE, 2018
```



## Input

Our developed program works with the fragment file used in [AltHap](https://github.com/realabolfazl/AltHap) and SDhaP  with the following format:


```
Number of reads
Number of columns 
"Number of parts" "Read ID"  "index of the first allele of the 1st part"  "Alleles of 1st part"  "index of first allele of 2nd part"  "Alleles of 2nd part" ... "Quality scores of alleles of all parts  " 
```
The last line is repeated for each (paired/mate/singel) read. It is noted that out algorithm can also work with 10x data. I this case each line corresponds to each barcode.

An example of fragment file:
```
3
9
1 chr1_1 1 010 @@@@@
2 chr1_2 3 011 7 010 @@@@@@ 
1 chr1_2 1 01 7 010 @@@@@@ 
```


When you have a Bam and VCF files, we suggest you to use [ExtractHAIRS](https://github.com/vibansal/HapCUT2) to create a fragment file.

```
./extractHAIRS  --bam reads.sorted.bam --VCF variants.VCF --out fragment.txt
```
This works for diploids, if you want to use it for polyploids, artificially convert the genotypes (GP tags) to a "0/1" using `sed`.


## Output

The output of algorithm is a text file named as `Reconstructed_Haplotype.txt`. For the all-heterzygous case 
```
Block 1 3  
1 0
2 1
3 0
````
in which first column is the index of variant in VCF file and the second column is the allele at that SNP. It is obvious that another haplotype can be found by flipping this.  For the case where heterzygous and homozygous variants exists the output file contains three column.

# How to use


```
matlab -r "HapMC(fragment_file,Hap_algorithm)"
```
Hap_algorithm can be either 'O', 'S' or 'N' correspond to 'HapOPT', 'HapSVT' or 'HapNuc'  respectively,



## Test use
Here, we can test the example of paper in Table 2. Run this lines in Terminal!

```
git clone https://github.com/smajidian/HapMC.git
cd HapMC
matlab -r "HapMC('data/fragment_sample.txt','O');exit"; #if it is not in bash /Applications/MATLAB_R2018b.app/bin/matlab 
cat Reconstructed_Haplotype.txt
```



``convert_frag_mat.m`` can be used for converting a fragment file  to a matrix of MATLAB format with .mat extension.

``first_block_extractor.m`` is used within the matlab code, HapMC, for extracting overlaped reads which are considered as a connected read block.

The `` HapMC.m `` is the core of algorithm. 
For converting the matrix file of the simulated data, one can use the R code.
For comparison we can use compare.py. To convert the output to a more standard format, on can use the convert_hapcut.py.


