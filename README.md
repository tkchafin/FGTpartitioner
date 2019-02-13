# FGTpartitioner
Partitions an alignment using the 4-gamete test

Input is a VCF file, which can represent unphased genotypes, and output is a parsimonious set of breakpoints separating non-overlapping intervals which do not show evidence of recombination, as tested using the four-gamete test. 

### Dependencies
- Python >3
- pyVCF
- intervaltree

The easiest way to install the dependencies is through conda:
```
conda install -c conda-forge -c bioconda pyvcf intervaltree
```
You will additionally need tabix installed to pre-process your input file:
```
conda install -c bioconda tabix 
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
