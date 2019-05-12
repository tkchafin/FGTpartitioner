# FGTpartitioner
Partitions genome data using the 4-gamete test into a minimal number of blocks which contain no recombinations (=FGT conflicts)

Input is a VCF file, which can represent unphased genotypes, and output is a parsimonious set of breakpoints separating non-overlapping intervals which do not show evidence of recombination, as tested using the four-gamete test. 

Of note, there are multiple options for partioning a genome using the four-gamete test. Here are a few, sorry if I left anyone out:
- https://github.com/RILAB/rmin_cut
- https://github.com/YichaoOU/genome_partition
- http://www.csbio.unc.edu/mcmillan/pubs/BCB10_Wang.pdf

FGTpartitioner is my implementation, if you find it useful for your research, please just cite this GitHub page:
```
Chafin, TK. 2019. FGTpartitioner: https://github.com/tkchafin/FGTpartitioner
```


### Status
FGTpartitioner is currently working properly, and finds the same FGT conflicts as other programs that I have tested. I may work on continuing to optimize it a bit, but for right now the best solutions for speeding up FGTpartitioner with very large alignments is to use multiprocessing <-t> and to set a maximum physical distance for calculating FGTs <-d>. The justification for the latter is that there exists a certain physical map distance which is almost guaranteed to have spanned multiple recombinations, thus there is no point in continuing to search for FGT conflicts. Ideally you can find some published studies calculating linkage disequilibrium for your species or something similar, and inform this parameter with the larger-end of the expected linkage block size distribution.

### Dependencies
Requires Python 3 and the following modules:
- pyVCF 
- pySAM
- intervaltree
- Cython > 0.27 
- pathos
- tabix

You will additionally need tabix installed to block-compress and index your VCF file (which enables me to parse it more quickly).

The easiest way to install all of the dependencies is through conda:
```
conda install -c conda-forge -c bioconda pyvcf intervaltree cython pathos pysam tabix
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
tyler:FGTpartitioner $ ./FGTpartitioner.py -h

Must provide VCF file <-v,--vcf>

FGTpartitioner.py

Contact:Tyler K. Chafin, tylerkchafin@gmail.com

Usage:  FGTpartitioner.py -v <input.vcf> -r <1|2|3> [-c chr1]

Description: Computes parsimonious breakpoints to partition chromosomes into recombination-free blocks

	Arguments:
		-v	: VCF file for parsing
		-r	: Strategy to treat heterozygotes. Options:
			   1: Randomly haploidize heterozygous sites
			   2: Fail FGT if heterozygote might be incompatible
			   3: Pass FGT if heterozygote might be compatible
		-c	: Chromosome or contig in VCF to partition
		-o	: Output file name [default: regions.out]
		-t	: Number of threads for parallel execution
		-m	: Minimum number of individuals genotyped to keep variant [default=2]
		-a	: Maximum number of alleles allowed per locus [default=2]
		-h	: Displays help menu
```
One important option which you will need to consider is <-r>, which determines how FGTpartitioner behaves when it encounters heteroozygotes. The four-gamete test has several assumptions, the most important being: 1) That you have sampled haploid chromosomes; and 2) an [infinite-sites](https://en.wikipedia.org/wiki/Infinite_sites_model) mutation model (e.g. all mutations occur at a new site- no back mutation, or multiple mutations per site). You can find more details below in the "Four-Gamete Test" section.

In order to meet assumption #1, we need to manipulate our [unphased genotype data](https://www.biostars.org/p/7846/). FGTpartioner allows 3 ways in which this can be accomplished (passed as an integer option to the -r flag): 1) <-r 1> Randomly choose one allele and treat the sample as homozygous for that allele, at that position; 2) <-r 2> Ask if **either** allele causes a failure of the four-gamete test, and treat the comparison as failed if so (e.g. a pessimistic/safe approach); or 3) <-r 3> ask if **either** allele could possibly be consistent with the four-gametes assumption, and pass the comparison if so (e.g. an optimistic approach). This pessimistic/optimistic approach was inspired by [Wang et al (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5690570/)

### Output
The output of FGTpartitioner is a list of regions in GATK format (chromosome:start-end), in a 1-based indexing scheme:
```
CM000001.3:1-218120
CM000001.3:218121-218432
CM000001.3:218433-221068
CM000001.3:221069-224784
CM000001.3:224785-228655
...
...
...
```
You can then use these to subset your VCF using your choice of tool. I usually use [vcftools](http://vcftools.sourceforge.net/) and [bcftools](https://samtools.github.io/bcftools/bcftools.html) for manipulating VCF files. 

If you want to ultimately parse multiple sequence alignments from your VCF file and the FGTpartitioner outputs, I would recommend my tool [vcf2msa.py](https://github.com/tkchafin/vcf2msa.py).

To get an average region size, you can use the following bash command:
```
less regions.out | sed "s/.*://g" | awk 'BEGIN{FS="-"}{print $2-$1}' |awk 'BEGIN{sum=0; count=0; sumsq=0}{count+=1; sum +=$1; sumsq+=$1*$1}END{print(sum/count); print sqrt(sumsq/count - (sum/count)**2)}'
```

### The Four-Gamete Test
The general principle of the [four-gamete test (FGT)](https://en.wikipedia.org/wiki/Four-gamete_test) is this: If we've sampled two positions (SNPs) along a chromosome, **assuming that multiple mutations at a site never occur** (= 'infinite sites' assumption), the ONLY way that we could possibly sample haploid chromosomes exhibiting **all** combinations of alleles at those two positions is if recombination had occured. 

For example, if we sampled two SNPs on a chromosome, one at position 1028 which is either an A or a T, and another at position 4589 which is either a G or a C:

| Position  | Allele1 | Allele2 |
|:-----: |:-----:| :-----:|
| 1028     | A | T |
| 4589     | G | C |

It should be impossible to observe all possible combinations of A/T and G/C (i.e. A at position 1 and G at position 2, A as position 1 and C at position 2, and so on...) without one of two things occuring: 1) There have been multiple convergent substitutions at a site (see: [homoplasy](https://en.wikipedia.org/wiki/Homoplasy)); or 2) A crossover event has occured, allowing chromosomes to exchange alleles. Since the FGT assumes that #1 never occurs (the infinite sites assumption), then we must conclude that #2 has occured. 

Conducting the test is simple. If we sample at least 4 haploid chromosomes, we shouldn't see all four possible combinations of alleles. For example, the following sampled chromosomes do NOT violate the four-gamete test:

| Sample | 1028 Allele | 4589 Allele |
|:-----: |:-----:| :-----:|
| 1     | A | C |
|2     | A | C |
|3     | A | C |
|4     | T | G |

Here, only 2 combinations are seen here. One haplotype, found in samples 1-3, has ----A--------C----, while the other haplotype (sample 4) has ----T--------G----. However, if we sampled the following chromosomes:

| Sample | 1028 Allele | 4589 Allele |
|:-----: |:-----:| :-----:|
| 1     | A | G |
|2     | A | C |
|3     | T | C |
|4     | T | G |

We see that all possible combinations were sampled: 
```
----A--------C----
----A--------G----
----T--------C----
       and
----T--------G----
```
This is **only** possible (again, given an assumption that multiple substitutions per site never occurs), if at some point there had been a crossover event duyring meiosis allowing recombination between these two sites:
```
----A--------C----   -->   ----A---\/----C----   -->  ----A--------G----  
----T--------G----   -->   ----T---/\----G----   -->  ----T--------C---- 
```
Meaning, in this sample of 4 chromosomes, all possible gametes resulting from this recombination event are seen. Understanding this should demonstrate the difficulty of applying the four-gamete test to an unphased diploid sample: What do we do with heterozygotes?

If we sample the following unphased genotype:
----A/T--------C/G---- (where A/T and C/G represent heterozygous genotypes). 

How do we know the haplotypes of the two individual chromosomes? The short answer is, given this data alone, that we can't. That is why I give 3 options in FGTpartitioner:

**-r 1**
This option randomly resolves heterozygotes, yielding one of the following as a pseudo-haplotype for this sample:
```
----A--------C----
----A--------G----
----T--------C----
----T--------G----
```
**-r 2**
This option considers all possible pseudo-haplotypes for this sample, and fails the FGT if **any** of these would violate the FGT. 

**-r 3**
This option considers all possible pseudo-haplotypes for the sample, and optimistically assumes that if **any** of them could be consistent with the FGT, that recombination did not occur. 

While not ideal, in the absence of phasing information (as is often the case when dealing with non-model data) we don't have many options. 


### Algorithm
The basic idea is that we perform a minimal number of pairwise comparisons between SNPs so as to characterize any possible recombination events which might have occured, and then place a minimal number of breakpoints along the chromosome to yield a parimonious set of blocks which do not contain variants violating the FGT. 

That is to say that, none of these blocks appear to contain recombination events according to a relatively low-power test (the FGT). I also use a 'greedy' approach in that I retain all variant information. As a result, these blocks likely **do** contain recombinations, but the FGT lacks the power to identify it. 

When identifying FGT violations, I place these as intervals in an [interval tree](https://en.wikipedia.org/wiki/Interval_tree) data structure, which is an efficient way of storing and parsing this sort of data. For each SNP-pair which violate the FGT, I also assign a layer (k) representing the number of SNPs (since these are the only informative columns in the alignment) which are subsumed by the interval. For example, if the A/C and A/T SNPs in the following violated the FGT, the interval length (k), or the number of SNPs spanned, would be 5.
```
        _________________________________________________
--------A/T----------T/C----T/C----G/A---C/T------------A/C--------T/G--
```
With increasing distance between SNPs, the higher probability that multiple recombination events have occured. For this reason, if multiple "conflict intervals" overlap, I place a higher weight on the interval with the minimum-k:
```
        _________________________________________________
	              _______
--------A/T----------T/C----T/C----G/A---C/T------------A/C--------T/G--
```
In this example, two conflicts were found: one of k=1 and another of k=5. To "solve" this interval tree (i.e. place the smallest number of breakpoints so as to eliminate all FGT conflicts), I solve from min(k) to max(k). Here, FGTpartitioner would find the following solution:
```
        __________________|_______________________________
	              ____|___
--------A/T----------T/C--|--T/C----G/A---C/T------------A/C--------T/G--
```
Where "|" represents the breakpoint. The algorithm first considers the k=1 interval, and places a breakpoint in the center between the flanking SNPs. This breakpoint also resolves the k=5 interval, thus solving the alignment to create two blocks (one breakpoint), which is the most parsimonious solution. In reality, we do not know where the breakpoint actually occured: it could be anywhere within the k=1 interval. My solution was to evenly divide these monomorphic nucleotides between the two breakpoints. This might not be your desired behavior... But fortunately this would be very easy to change in the code- have at it. 

In reality, FGTpartitioner actually considers every possible breakpoint (centerpoint between SNPs, or "nodes"), so as to maximally resolve intervals. For example, with the following case:
```
        _______________________________________________________
	                        _________________________________________
	              ____________________
--------A/T----------T/C--(1)--T/C--(2)--G/A---C/T------------A/C--------T/G--
```
The algorithm would look at centerpoint #1, count the number of intervals which would be resolved (=2), and then do the same for centerpoint #2 (=3). In this case, the most parsimonious solution is to place a single breakpoint at candidate #2. 

Because of this strategy, you might have noticed that it would be a waste of time to consider any pairwise comparisons after identifying a conflict, because those intervals would always be of larger k (and thus never considered, since we solve from min(k) to max(k)). So, I don't actually consider all possible comparisons, but instead consider sequentially more distance SNPs until I find a conflict- which will be guaranteed to be the minimal-k conflict for the focal SNP. Starting from the first SNP, I would make passes until discovering a conflict:

```
Pass #1 ==============>
Pass #2 =====================>
Pass #3 ============================>
...
...
--------A/T----------T/C----T/C----G/A---C/T------------A/C--------T/G--
```

The same procedure is then performed for each SNP:
```
Pass #1               =======>
Pass #2               ==============>
Pass #3               ====================>
...
...
--------A/T----------T/C----T/C----G/A---C/T------------A/C--------T/G--
```
There are probably better ways to do it, but this accomplishes my goal so ¯\_(ツ)_/¯

### Profiling
To efficiently process very large alignments, FGTpartitioner is best used on an HPC system.

#### Memory efficiency

For example, with a very large full-chromosome alignment for a mammalian dataset, comprising >2,000,000 variants, FGTpartitioner peaked at about 18GB memory usage. 
![](https://raw.githubusercontent.com/tkchafin/FGTpartitioner/master/images/mem_profile.png)
Note that this will be effected by the number of cores, with less cores requiring less memory.

#### Runtimes

Predicting the runtimes for FGTpartitioner is difficult, as it will depend on both the number of variants, and the number of them showing FGT conflicts (which increases the time for finding the most parsimonious breakpoints). 

Runtimes in general will scale n^2 with dataset size:
![](https://raw.githubusercontent.com/tkchafin/FGTpartitioner/master/images/size_scaling.png)

This can be helped by using the <-d> makimum allowable distance, which produced a linear speed-up:
![](https://raw.githubusercontent.com/tkchafin/FGTpartitioner/master/images/distance_scaling.png)
Note that these were calculated on a fairly large alignment of ~2 million variants. See above for how dataset size will effect runtimes.

#### Multiprocess scaling

FGTpartitioner shows diminishing returns with added processors, with adding a second core generally halving the runtime, while adding a 16th core does very little:
![](https://raw.githubusercontent.com/tkchafin/FGTpartitioner/master/images/cpu_scaling.png)

As such, if you have limited cores available but very large amount of sequence to process, a better investment of your available cores would be to use a modest amount (2-4) per chromosome (selected using -c), and to run individual chromosomes in separate FGTpartitioner runs. 

One way you could do this would be with the excellent [GNU Parallel](https://www.gnu.org/software/parallel/) tool:
```
parallel `python3 FGTpartitioner.py -i input.vcf.gz -t 2 -c {} -o {}_regions.out' ::: {chr1 chr2 chr3 chr4}  
```
The above would dedicate 2 cores each to 4 separate FGDpartitioner runs (one for each of 4 chromosomes). You can read more on the use of GNU Parallel [here](https://www.gnu.org/software/parallel/parallel_tutorial.html)
### License
Copyright © 2019 Tyler K. Chafin <tylerkchafin@gmail.com>

This work is free. You can redistribute it and/or modify it under the
terms of the Do What The Fuck You Want To Public License, Version 2,
as published by Sam Hocevar. See http://www.wtfpl.net/ for more details.
