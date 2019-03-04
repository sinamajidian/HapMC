# HapMC



We apply matrix completion methods for haplotype assembly from NGS reads to develop the new HapSVT, HapNuc, and HapOPT algorithms.

## Cite

```

Sina Majidian, Mohammad Hossein Kahaei, "NGS based haplotype assembly using matrix completion",
Submitted to PLOS ONE, 2018

```



## Input



Our developed program works with the fragment file in which each line is with the following format:

```
"Number of parts" "Read ID"  "index of the first allele of the 1st part"  "Alleles of 1st part"  "index of first allele of 2nd part"  "Alleles of 2nd part" ... "Quality scores of alleles of all parts  " 
```

An example of fragment file:
```
1 chr1_1 1 010 @@@@@
2 chr1_2 3 011 7 010 @@@@@@ 
1 chr1_2 1 01 7 010 @@@@@@ 

```


When you have a Bam and VCF files, we suggest you to use [ExtractHAIRS](https://github.com/vibansal/HapCUT2).
 
```
./extractHAIRS  --bam reads.sorted.bam --VCF variants.VCF --out fragment.txt
```
This works for diploids, if you want to use it for polyploids, artificially convert the genotypes (GP tags) to a "0/1" using `sed`.

You can set the argument in order to set which algorithm you want to run.

```
matlab -r 'HapMC(fragment_file,Hap_algorithm)';
```
Hap_algorithm can be either 'O', 'S' or 'N' correspond to 'HapOPT', 'HapSVT' or 'HapNuc'  respectively,


## Output


```
Block
1 0 1
2 1 0
3 0 1
````



# Test use

```
git clone https://github.com/smajidian/HapMC.git
cd HapMC
matlab -r 'HapMC('data/fragment_sample.txt','O')';
cat Reconstructed_Haplotype.txt

```



``convert_frag_mat.m`` can be used for converting a fragment file in hapcut format to a matrix of MATLAB format with .mat extension.

``first_block_extractor.m`` is used within the matlab code, HapMC, for extracting overlaped reads which are considered as a connected read block.

The `` HapMC.m `` is the core of algorithm. The input of this code is a matrix of matlab. The output will be three haplotype fles, named as ``chr1_opt.hap``, ``chr1_svt.hap``, and  ``chr1_nuc.hap``.
