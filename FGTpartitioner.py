#!/usr/bin/python

import re
import sys
import os
import getopt
import vcf
import pysam
import random
import math
import collections
import operator
import intervaltree
from SNPcall import SNPcall
from intervaltree import Interval, IntervalTree
from collections import OrderedDict
from functools import partial
import multiprocessing

def main():
	params = parseArgs()

	print("\nOpening VCF file:",params.vcf)
	vfh = vcf.Reader(filename=params.vcf)

	#grab contig sizes
	contigs = dict()
	if params.chrom:
		print("Only reading chromosome:",params.chrom)
	else:
		print("Reading all chromosomes from VCF...")

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

	print("Diploid resolution strategy: ", end="")
	if params.rule==1:
		print("Random", end="")
	elif params.rule==2:
		print("Pessimistic", end="")
	elif params.rule==3:
		print("Optimistic",end="")
	else:
		print("Invalid!",end="")
		sys.exit(1)
	print(" (change with -r)\n")

	breakpoints = collections.OrderedDict()

	print("Starting FGT sweeps.....\n")
	#for each chromosome
	for this_chrom in contigs:
		print("Performing FGT check on:", this_chrom)
		#initialize data structures
		tree = IntervalTree()
		k_lookup = dict()
		nodes = list()
		start = 0
		stop =1
		index=1 #keys for intervals
		breaks = list()

		#Gather relevant SNP calls from the VCF
		records = vfh.fetch(this_chrom)
		if not records:
			print("Not enough records found for chromosome:",this_chrom)
		else:
			nodes = fetchNodes(records, this_chrom)

		#Traverse node list to find FGT conflicts
		if len(nodes) > 2:

			#call parallel worker runs here if params.threads>1
			print("\tSeeking intervals across",this_chrom,"using",params.threads,"threads.")
			if params.threads > 1:
				findFGTs_parallel(tree, nodes, params, k_lookup)
			else:
				findFGTs(tree,nodes,params,k_lookup)

			print("\tFound ",len(tree),"intervals.")
			#print(tree)

			#order k_lookup by k
			#NOTE: Over-rode the __lt__ function for Intervals: see __main__ below
			#otherwise this would sort on start position!
			sorted_k = sorted(k_lookup.items(), key=operator.itemgetter(1)) #gets ordered tuples
			#print(sorted_k)

			#start resolving from lowest k
			print("\tResolving FGT incompatibilities...")
			breaks = resolveFGTs(tree, sorted_k, nodes)

		else:
			print("\tNo passing variants found for chromosome",this_chrom,"")

		#submit all selection breakpoints
		breakpoints[this_chrom] = breaks
		print("\tFound",len(breaks),"most parsimonious breakpoints.")

	print("\nWriting regions to file (format: chromosome:start-end)")
	regions = getRegions(breakpoints, contigs)

	if len(regions) > 0:
		write_regions(params.out, regions)
	else:
		print("No regions found.")

	print("Done\n")



'''
Skip sites with <4 genotyped individuals
Skip sites with >2 alleles

Processing algorithm:
	k = order of overlap (number of SNP nodes encompassed)
	for each SNP i, explore right until finding minimum-k conflict
	Add interval to data structure, increment i and continue

Data structure:
	Interval tree:
	Interval tree, but intervals indexed by k
	Solving goes from k=1 -> kmax

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


def resolveFGTs(tree, sorted_k, nodes):
	breaks = list()
	for i in sorted_k:
		this_interval = i[1]
		if this_interval in tree:
			#get all intervals overlapping range(start:end)
			#overlaps = tree[this_interval.begin: this_interval.end]
			#query each centerpoint
			global_center = (nodes[this_interval.data.start].position + nodes[this_interval.data.end].position) / 2
			max_centerpoint = global_center
			#max_intersects = tree[global_center]
			max_depth = len(tree[global_center])
			#print("K=",this_interval.data.k)
			#find centerpoint of greatest depth in range of target interval
			for check_start in range(this_interval.data.start,this_interval.data.end-1):
				#print(check_start)
				centerpoint = (nodes[check_start].position + nodes[check_start+1].position) / 2
				intersects = tree[centerpoint] #all intervals intersecting with current centerpoint
				local_depth = len(intersects) #depth at centerpoint
				#if current depth highest seen, keep it
				if local_depth > max_depth:
					max_depth = local_depth
					max_centerpoint = centerpoint
					#max_intersects = intersects
			#remove all intervals overlapping with centerpoint at greatest depth
			del tree[max_centerpoint]

			#add to breakpoints
			breaks.append(max_centerpoint)
			#print("Selected breakpoint:",max_centerpoint)
	return(breaks)


#TODO:try to speed this up. 13% of total runtime
def findFGTs(tree, nodes, params, k_lookup):
	start = 0
	end = 1
	count=0
	index = 0
	while start <= len(nodes):
		#print("Start=",start)
		#print("End=",end)
		if end >= len(nodes):
			start = start + 1
			end= start + 1
			continue
		#print("Comparing:",nodes[start].position,"and",nodes[end].position)
		#Check if start and end are compatible
		compat = nodes[start].FGT(nodes[end], params.rule)
		if compat == True: #if compatible, increment end and continue
			#print("Compatible! Checking next SNP")
			end+=1
			continue
		else: #if FGT fails, submit interval to IntervalTree, and increment start
			#print("Not compatible!")
			interval = Interval(nodes[start].position, nodes[end].position, IntervalData(start, end, index))
			k_lookup[index] = interval #k-layer for this interval
			tree.add(interval) #add interval from start.position to end.position
			#print(interval)
			index +=1 #increment key, so all will be unique
			start = start+1 #move start to next SNP
			end = start+1 #reset end to neighbor of new start

#findFGTs function for the parallel call
def findFGTs_worker(local_nodes, params, k_lookup, proc_number):
	print("proc")
	try:
		start = 0 + proc_number #everyprocess starts at an offset
		skip = params.threads #every process checks for FGTs at an interval
		end = 1
		count=0
		index = 0 + proc_number
		local_tree = IntervalTree()
		maximum = len(local_nodes)

		#initialize local random number seed
		random.seed(random.randrange(sys.maxsize)+proc_number)

		while start <= maximum:
			#print("Start=",start)
			#print("End=",end)
			if end >= maximum:
				start = start + skip
				end= start + 1
				continue
			#print("Comparing:",nodes[start].position,"and",nodes[end].position)
			#Check if start and end are compatible
			compat = local_nodes[start].FGT(local_nodes[end], params.rule)
			if compat == True: #if compatible, increment end and continue
				#print("Compatible! Checking next SNP")
				end+=1
				continue
			else: #if FGT fails, submit interval to IntervalTree, and increment start
				#print("Not compatible!")
				interval = Interval(local_nodes[start].position, local_nodes[end].position, IntervalData(start, end, index))
				k_lookup[index] = interval #k-layer for this interval
				local_tree.add(interval) #add interval from start.position to end.position
				print(interval)
				index +=skip #increment key, so all will be unique
				start = start+skip #move start to next SNP
				end = start+1 #reset end to neighbor of new start

		#return local_tree
		return(local_tree)

	except Exception as e:
		raise Exception(e)


def findFGTs_parallel(tree, nodes, params, k_lookup):
	print("parallel")

	#calculate skip sizes for each process
	#e.g. process 1 of 4 runs FGT chackes for SNP 1, 5, 9, etc

	#for each process call, generate: deep copies of nodes
	#local_nodes = list()
	proc_numbers = list()
	for i in range(0,params.threads):
		proc_numbers.append(i)
		#local_nodes.append()

	#multiprocess call
	print("multiprocess call for processes:",proc_numbers)
	try:
		print("try")
		with multiprocessing.Pool(processes=params.threads) as pool:
			print("with")
			func = partial(findFGTs_worker, nodes, params, k_lookup)
			results = pool.map(func, proc_numbers)
			print(results)
	except Exception as e:
		pool.close()
		raise Exception(e)
		sys.exit(e)

	pool.close()
	pool.join()

	sys.stdout.flush()


	sys.exit(1)



def fetchNodes(records, this_chrom):
	nodes = list()
	miss_skips = 0
	allel_skips = 0
	c=0
	for rec in records:
		#if this SNP
		if rec.is_snp:
			if rec.num_called < 2:
				miss_skips +=1
			elif len(rec.alleles) > 2:
				allel_skips +=1
			else:
				#print(rec.samples)
				samps = [s.gt_type for s in rec.samples]
				nodes.append(SNPcall(rec.POS, samps))
				c+=1
		if c>5000:
			break
	if miss_skips > 0:
		print("\tChromosome",this_chrom,"skipped",str(miss_skips),"sites for too much missing data.")
	if allel_skips > 0:
		print("\tChromosome",this_chrom,"skipped",str(allel_skips),"sites for >2 alleles.")
	return(nodes)


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

#breaks is a dict where key = chromosome and value = list of breakpoints
#lengths is a dict where key = chromosome and value = total length
#function to calculate region bounds from breakpoints
#returns list of region tuples
def getRegions(breaks, lengths):
	ret = list()
	for chrom in breaks.keys():
		sorted_breaks = sorted(breaks[chrom])
		if len(breaks[chrom]) == 1:
			ret.append(tuple([chrom, 1, int(round(breaks[chrom][0]))]))
			ret.append(tuple([chrom, int(math.ceil(sorted_breaks[-1])), int(lengths[chrom])]))
		elif len(breaks[chrom]) > 1:
			first = tuple([chrom, 1, int(round(sorted_breaks[0]))])
			ret.append(first)

			for idx, br in enumerate(sorted_breaks[1:-1]):
				ret.append(tuple([chrom, int(round(sorted_breaks[idx])), int(math.ceil(sorted_breaks[idx+1]))]))

			last = tuple([chrom, int(math.ceil(sorted_breaks[-1])), int(lengths[chrom])])
			ret.append(last)
		else:
			ret.append([chrom, 1, int(lengths[chrom])])
	return(ret)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'v:r:c:o:t:h', \
			[])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.rule=1
		self.chrom=None
		self.out="regions.out"
		self.threads=1

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
			elif opt == "o":
				self.out = arg
			elif opt == "t":
				self.threads=int(arg)
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
		print ("Description: Computes parsimonious breakpoints to partition chromosomes into recombination-free blocks")

		print("""
	Arguments:
		-v	: VCF file for parsing
		-r	: Strategy to treat heterozygotes. Options:
			   1: Randomly haploidize heterozygous sites
			   2: Fail FGT if heterozygote might be incompatible
			   3: Pass FGT if heterozygote might be compatible
		-c	: Chromosome or contig in VCF to partition
		-o	: Output file name [default: regions.out]
		-t	: Number of threads for parallel execution
		-h	: Displays help menu
""")
		print()
		sys.exit()

###overriding __lt__ method for Interval
def IntervalSort(self,other):
	return(self.data < other.data)

class IntervalData():
	def __init__(self, start, end, index):
		self.start = start
		self.end = end
		self.index=index
		self.k = end-start

	def __lt__(self, other):
		return(self.k < other.k)

	def getK(self):
		return(self.k)

	def getStart(self):
		return(self.start)

	def getEnd(self):
		return(self.end)

	def getIndex(self):
		return(self.index)

	def __repr__(self):
		return("IntervalData()")

	def __str__(self):
		s="IntervalData.k="+str(self.k)
		return(s)


#Call main function
if __name__ == '__main__':
	try:
		Interval.__lt__ = IntervalSort #override sort for Interval
		main()
	except KeyboardInterrupt:
		sys.exit(1)
