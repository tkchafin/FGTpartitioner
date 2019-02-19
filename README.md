# FGTpartitioner
Partitions an alignment using the 4-gamete test

Input is a VCF file, which can represent unphased genotypes, and output is a parsimonious set of breakpoints separating non-overlapping intervals which do not show evidence of recombination, as tested using the four-gamete test. 

### Dependencies
Requires Python 3 and the following modules:
- pyVCF 
- pySAM
- intervaltree
- Cython > 0.27 
- multiprocess

You will additionally need tabix installed to pre-process your input file.

The easiest way to install all of the dependencies is through conda:
```
conda install -c conda-forge -c bioconda pyvcf intervaltree cython multiprocess pysam tabix
```

If you don't have conda installed, go [here](https://conda.io/en/latest/miniconda.html) and choose the correct Python3 installer for your system.

To prep FGTpartioner for running, you will first need to pre-compile the cythonized portions of the code:
```
python setup.py build_ext --inplace
```


### Inputs

Coming Soon

You will need to block-compress and index your VCF file to speed up parsing:
```
#Run bgzip (NOT gzip) to compress your joint VCF file
bgzip file.vcf

#tabix to index it
tabix -h -f -p vcf file.vcf.gz
```















### License
Copyright © 2019 Tyler K. Chafin <tylerkchafin@gmail.com>

This work is free. You can redistribute it and/or modify it under the
terms of the Do What The Fuck You Want To Public License, Version 2,
as published by Sam Hocevar. See http://www.wtfpl.net/ for more details.
