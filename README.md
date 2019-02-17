# HapMC



We apply matrix completion methods for haplotype assembly from NGS reads to develop the new HapSVT, HapNuc, and HapOPT algorithms.


```

Sina Majidian, Mohammad Hossein Kahaei, "NGS based haplotype assembly using matrix completion",
Submitted to PLOS ONE, 2018

```


``convert_frag_mat.m`` can be used for converting a fragment file in hapcut format to a matrix of MATLAB format with .mat extension.

``first_block_extractor.m`` is used within the matlab code, HapMC, for extracting overlaped reads which are considered as a connected read block.
