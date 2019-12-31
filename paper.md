---
title: 'FGTpartitioner: A rapid method for parsimonious delimitation of ancestry breakpoints in large genome-wide SNP datasets'
tags:
  - Python
  - phylogenomics
  - recombination
  - four-gamete test
  - SNP analysis
authors:
  - name: Tyler K. Chafin
    orcid: 0000-0001-8687-5905
    affiliation: 1
affiliations:
 - name: Tyler K. Chafin, Ph.D Candidate, University of Arkansas
   index: 1
date: 31 December 2019
bibliography: paper.bib

---


# Summary

Partitioning large (e.g. chromosomal) alignments into ancestry blocks is a common step in 
phylogenomic analyses `[@Dutheil2009]`. However, current solutions require complicated 
analytical assumptions, or are difficult to implement due to excessive runtimes.  Multiple approaches 
have been proposed for delimiting ancestry blocks in genomes (i.e. establishing recombination 
breakpoints), which generally fall into one of two categories: those which require dense or phased 
genotypic data (Liu et al., 2013); and those with complex analytical assumptions which require 
the definition of informative prior probability distributions and are computationally intensive 
(Dutheil et al., 2009). Both conditions are problematic for genome-scale studies of non-model 
species, where large-scale resequencing and phased reference data are unavailable, 
and genomes are often sequenced at low coverage. 

I here describe a solution, ``FGTpartitioner``, which is specifically designed for use with 
non-model genomic data without the need for high-quality phased reference data or dense 
population-scale sampling. ``FGTpartitioner`` delimits chromosome scale alignments using a 
fast interval-tree approach which detects pairwise variants which violate the four-gametes 
assumption (Hudson & Kaplan, 1985), and rapidly resolves a most parsimonious set of recombination 
events to yield non-overlapping intervals which are both unambiguously defined and consistent 
regardless of processing order. These sub-alignments are then suitable for separate phylogenetic 
analysis, or as a ‘first pass’ which may facilitate parallel application of finer-resolution (yet 
more computationally intensive) methods.

For ease of application, inputs are required to follow the widely used VCF format (Danecek et al., 2011). 
Users may provide parameter settings as arguments in the command-line interface which can restrict block 
delimitation to a certain chromosome (<-c> flag), with the option to additionally target a region via 
start (<-s>) and end (<-e>) coordinates. Parallel computation is also possible (<-t>) for particularly 
large alignments. After parsing user-inputs, the workflow of ``FGTpartitioner`` is as follows:

(1)	For each SNP, perform four-gamete tests sequentially for rightward neighboring records, up to a 
maximal physical distance (if defined; <-d>) and stopping when a conflict (=’interval’) is found. Intervals are
stored in a self-balancing tree. When using multiprocessing (<-t>), daughter processes are each provided 
an offset which guarantees a unique pairwise SNP comparison for each iteration 
(2)	Merge interval trees of daughter processes (if <-t 2 or greater>)
(3)	Assign rank k per-interval, defined as the number of SNP records (indexed by position) spanned by each 
interval 
(4)	Order intervals by k; starting at min(k), resolve conflicts as follows: For each candidate recombination 
site (defined as the mid-point between SNPs), compute the depth d of spanning intervals. The  most parsimonious
breakpoint is that which maximizes d

These algorithm choices have several implications: indexing SNPs by physical position guarantees that 
the same recombination sites will be chosen given any arbitrary ordering of SNPs; and defining breakpoints 
as physical centerpoints between nodes means that monomorphic sites will be evenly divided on either side 
of a recombination event. Because monomorphic sites by definition lack phylogenetic information, they 
cannot be unambiguously assigned to any particular ancestry block, thus my solution is to evenly divide them.
Heterozygous sites in diploid genomes are dealt with in multiple ways. By default, FGTpartitioner will 
randomly resolve haplotypes. The user can select an alternate resolution strategy (<-r>) which will either 
treat a SNP pair as failing if any resolution meets the four-gamete condition, or as passing if any possible 
resolution passes (i.e. the 'pessimistic' and 'optimistic' strategies of Wang et al., 2010).

In conclusion, ``FGTpartitioner`` has several advantages over similar methods: 1) algorithmic and performance enhancements 
allow it to perform orders of magnitude faster, thus extending application to larger genomes; and 2) the 
flexibility of diploid resolution strategies precludes the need for haplotype phasing a priori. Validation 
using empirical data indicated the suitability of FGTpartitioner for highly distributed work on high-performance
computing clusters, with parallelization easily facilitated by built-in options in the command-line interface. 
Additionally, runtime and memory profiling indicate its applicability on modern desktop workstations as well, 
when applied to moderately sized datasets. Thus, it provides an efficient and under-friendly solution to alignment 
pre-processing for phylogenomic studies, or as a method of breaking up large alignments in order to efficiently 
distribute computation for more rigorous recombination tests.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"


# Acknowledgements
I would like to thank vonHoldt et al. (2016) and Kukekova et al. (2018) for making their raw datasets 
available via the NCBI SRA, and the Arkansas High Performance Computing Center and my Ph. D. advisors 
Drs. Michael and Marlis Douglas for access to computational resources which I used for benchmarking this program. 

# References
