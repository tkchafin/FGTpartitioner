#!/usr/bin/python

import re
import sys
import os
import getopt
import vcf
import pysam
import collections
import intervaltree
from intervaltree import Interval, IntervalTree
from collections import OrderedDict 

def main():
	params = parseArgs()

	vfh = vcf.Reader(filename=params.vcf)

	#grab contig sizes
	contigs = dict()
	for c,s in vfh.contigs.items():
		if params.chrom and s.id != params.chrom:
			continue
		contigs[s.id] = s.length

	if len(contigs) < 1:
		print("No chromosomes found. Please check your input or verify that options provided with -c are correct")
		sys.exit(1)
		
	'''
	1. for each chromosome (or 1 if region is set)
	2. build list where value=SeqRecord
	3. next, loop through ordered dict to build IntervalTree
		Interval = start, end, index
		k_lookup[index] = k (order, or number of nodes encompassed)
	4. after examining last node, resolve tree
		from smallest k to largest k. Maybe sort k_lookup by value
	
	tree resolution: 
	For each layer from k=1 to kmax:
		for each interval in layer k
			query center point of interval to get interval depth
			if overlaps exist:
				query each SNP-node centerpoint between interval start and end
				place break at centerpoint which maximizes depth, or center of maximum depth region
				delete all intervals intersecting with breakpoint
	'''
	

	#for each chromosome
	for this_chrom in contigs: 
		#initialize data structures
		tree = IntervalTree()
		k_lookup = dict()
		nodes = list()
		start = 0
		stop =1
		
		#gather records from VCF
		records = vfh.fetch(this_chrom)
		for rec in records:
			if rec.CHROM != this_chrom:
				continue
			else:
				#if this SNP
				if rec.is_snp and not rec.is_monomorphic:
					#Check if parsimony-informative
					#if no current lastNode, check FGT compatibility
					#if compatible:
					#	keep lastNode and examine next rec	
					print(rec.samples)
					
			# if not this_chrom:
			# 	this_chrom = rec.CHROM
			# 	start = 1
			# 	stop = 1
			# 	count = 0
			# #If we entered new chromosome, submit old break
			# elif this_chrom != rec.CHROM:
			# 	t = tuple([this_chrom, start, contigs[this_chrom]])
			# 	regions.append(t)
			# 	this_chrom = rec.CHROM
			# 	start = 1
			# 	stop = 1
			# 	tree = IntervalTree()
			# 	k_lookup = dict()
			# 	nodes = list()
			# 	#resolveTree()
			# 	sys.exit(0)
		


'''
Processing algorithm:
	k = order of overlap (number of SNP nodes encompassed)
	for each SNP i, explore right until finding minimum-k conflict
	Add interval to data structure, increment i and continue
	
Data structure:
	Interval tree:
	Interval tree, but intervals indexed by k
	Solving goes from k=1 -> kmax

check out kerneltree for a prepackaged option: https://github.com/biocore-ntnu/kerneltree
or: https://github.com/bxlab/bx-python

Or nested containment list might be better:https://academic.oup.com/bioinformatics/article/23/11/1386/199545

	
	
Solving algorithm:
	For each layer from k=1 to kmax:
		for each interval in layer k
			query center point of interval to get interval depth
			if overlaps exist:
				query each SNP-node centerpoint between interval start and end
				place break at centerpoint which maximizes depth, or center of maximum depth region
				delete all intervals intersecting with breakpoint

'''


#Function to write list of regions tuples, in GATK format
def write_regions(f, r):

	with open(f, 'w') as fh:
		try:
			for reg in r:
				ol = str(reg[0]) + ":" + str(reg[1]) + "-" + str(reg[2]) + "\n"
				fh.write(ol)
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()



#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'v:r:c:h', \
			[])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.rule=1
		self.chrom=None

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 'v':
				self.vcf = arg
			elif opt == 'r':
				self.rule=int(arg)
			elif opt == 'c':
				self.chrom = arg
			elif opt in ('h'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Must provide VCF file <-v,--vcf>")

		if self.rule not in [1, 2, 3]:
			self.display_help("Value for <-r> must be one of: 1, 2 or 3")

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nFGTpartitioner.py\n")
		print ("Contact:Tyler K. Chafin, tylerkchafin@gmail.com")
		print ("\nUsage: ", sys.argv[0], "-v <input.vcf> -r <1|2|3> [-c chr1]\n")
		print ("Description: Computes minimal breakpoints to partition chromosomes into recombination-free blocks")

		print("""
	Arguments:
		-v	: VCF file for parsing
		-r	: Strategy to treat heterozygotes. Options:
			   1: Randomly haploidize heterozygous sites
			   2: Fail FGT if heterozygote might be incompatible
			   3: Pass FGT if heterozygote might be compatible
		-c	: Chromosome or contig in VCF to partition
		-h	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
