# FGTpartitioner
Partitions genome data using the 4-gamete test into a minimal number of blocks which contain no recombinations (=FGT conflicts)

Input is a VCF file, which can represent unphased genotypes, and output is a parsimonious set of breakpoints separating non-overlapping intervals which do not show evidence of recombination, as tested using the four-gamete test. 

### Status
FGTpartitioner is currently working properly, and finds the same FGT conflicts as other programs that I have tested. However, it is currently very slow! I've sped it up slightly by Cython-izing a major bottleneck, and enabling a parallel search for FGT conflicts using the multiprocess module. 

### Dependencies
Requires Python 3 and the following modules:
- pyVCF 
- pySAM
- intervaltree
- Cython > 0.27 
- multiprocess

You will additionally need tabix installed to clock-compress and index your VCF file (which enables me to parse it more quickly).

The easiest way to install all of the dependencies is through conda:
```
conda install -c conda-forge -c bioconda pyvcf intervaltree cython multiprocess pysam tabix
```

If you don't have conda installed, go [here](https://conda.io/en/latest/miniconda.html) and choose the correct Python3 installer for your system.

### Installation

To prep FGTpartioner for running, you will first need to pre-compile the cythonized portions of the code:
```
python setup.py build_ext --inplace
```
After that, FGTpartitioner is ready to run. You can view the help menu by typing:
```
./FGTpartitioner.py -h
```

### Inputs

The input file (provided via -v) is a standard VCF file, including contig lengths in the ##contig headers. You can find an example of a VCF file in the examples/ directory. In short, a minimally conforming VCF should have the following structure:
```
##fileformat=VCFv4.2
##contig=<ID=chr1.scaffold1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">
##FORMAT=<ID=PL,Number=G,Type=Float,Description="Phred-scaled Genotype Likelihoods">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMP001 SAMP002 SAMP003 SAMP004
chr1.scaffold1  100     rs11449 G       A       .       PASS    .       GT      0/0     0/0     1/1     1/1
chr1.scaffold1  200     rs11449 T       A       .       PASS    .       GT      0/0     1/1     0/0     1/1
chr1.scaffold1  300 rs84825 A   T       .       PASS    .       GT      0/0     1/1     1/1     1/1
chr1.scaffold1  400 rs84825 A   G       .       PASS    .       GT      1/1
...
...
...
```

You will also need to block-compress and index your VCF file to speed up parsing:
```
#Run bgzip (NOT gzip) to compress your joint VCF file
bgzip file.vcf

#tabix to index it
tabix -h -f -p vcf file.vcf.gz
```
The result will be a binary compressed-VCF ".vcf.gz" file, and a ".vcf.gz.tbi" index file. The ".vcf.gz" will be the input provided to FGTpartitioner using the -v flag, and the '.vcf.gz.tbi" file should be in the same directory, and with the same prefix.

### Usage
You can view all of the possible options by calling the help menu in the command-line interface:
```
./FGTpartitioner.py -h
```













### License
Copyright © 2019 Tyler K. Chafin <tylerkchafin@gmail.com>

This work is free. You can redistribute it and/or modify it under the
terms of the Do What The Fuck You Want To Public License, Version 2,
as published by Sam Hocevar. See http://www.wtfpl.net/ for more details.
