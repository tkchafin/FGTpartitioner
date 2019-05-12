#!/usr/bin/python

import random

def rebuild(pos, calls):
	sc = SNPcall(pos, calls)
	return(sc)

cdef class SNPcall(object):
	cdef public:
		int position
		list calls

	def __init__(self, pos, samps):
		self.position = int(pos)
		self.calls = list(samps)
	
	# def __lt__(self, other):
	# 	return(self.position < other.position)
	
	def __richcmp__(self, other, op):
		if op == 2:#Py_EQ
			return(self.position == other.position)
		elif op == 3:#Py_NE
			return(not self.position == other.position)
		elif op == 0:
			return(self.position < other.position)
		elif op == 1:
			return(self.position <= other.position)
		elif op == 4:
			return(self.position > other.position)
		elif op == 5:
			return(self.position >= other.position)
		else:
			assert False

	#trying to make it pickle-able
	@staticmethod
	def rebuild(pos, calls):
		sc = SNPcall(pos, calls)
		return(sc)
		
	def __reduce__(self):
		return(rebuild, (self.position, self.calls))


##############################################################################

	#TODO:try to speed this up. 32% of runtime currently
	def FGT(self, other, rule):
		#print("Four gamete test for",self.position,"and",other.position)
		cdef list gametes = [0,0,0,0] #00, 01, 10, 11
		cdef list hets = list()
		cdef int gamete
		
		#TODO: This line takes a long time. 21% of total runtime 
		valid =set([0, 1, 2])
		#cdef list genotypes = [[gt, other.calls[i]] for i, gt in enumerate(self.calls) if gt and other.calls[i] in valid]
		cdef list genotypes = [[x,other.calls[i]] for i,x in enumerate(self.calls) if set([x, other.calls[i]]).issubset(valid)]

		#print("geno is:",genotypes)

		for geno in genotypes:
			#print(geno)
			gamete = self.hapCheck(geno)
			if gamete <= 3: #9 is designated 'missing data' placeholder
				gametes[gamete] = 1
			else:
				if 1 in geno:
					if rule == 1:
						if geno[0] == 1:
							geno[0] = random.choice([0,2])
						if geno[1] == 1:
							geno[1] = random.choice([0,2])
						#print(geno)
						gametes[self.hapCheck(geno)] = 1
					elif rule == 2:
						possible1 = list()
						possible2 = list()
						if geno[0] == 1:
							possible1 = [0,2]
						else:
							possible1 = [geno[0]]
						if geno[1] == 1:
							possible2 = [0,2]
						else:
							possible2 = [geno[1]]
						for i in possible1:
							for j in possible2:
								gametes[self.hapCheck([i,j])] = 1
					elif rule == 3:
						hets.append(geno)
		#print("found gametes:",gametes)
		if sum(gametes) == 4:
			return(False) #return False if not compatible
		elif hets:
			if not self.optimisticFGT(gametes, hets):
				return(False) #return false if not compatible 
			else:
				return(True)
		else:
			return(True)
############################################################################
	
	@staticmethod
	#TODO: Optimize; currently 11% of runtime after 2X speedup
	def hapCheck(geno):
		#print(geno)
		if geno[0] == 0:
			if geno[1] == 0:
				return(0)
			elif geno[1] ==2:
				return(1)
			else:
				return(9)
		elif geno[0] == 2:
			if geno[1] == 0:
				return(2)
			elif geno[1] ==2:
				return(3)
			else:
				return(9)
		else:
			return(9)
	
	def optimisticFGT(self, seen, hets):
		cdef list possibilities = list()
		cdef list locals = list()
		cdef possible1 = list()
		cdef possible2 = list()
		for het in hets:
			locals = list()
			possible1 = list()
			possible2 = list()
			if het[0] == 1:
				possible1 = [0,2]
			else:
				possible1 = [het[0]]
			if het[1] == 1:
				possible2 = [2,0]
			else:
				possible2 = [het[1]]
			for i in possible1:
				for j in possible2:
					copy = seen[:] #deep copy
					copy[self.hapCheck([i,j])] = 1
					locals.append(copy)
			possibilities = locals[:]
		#print(possibilities)
		for opt in possibilities:
			#print(opt)
			if sum(opt) != 4: #if ANY possibilities 
				return True #return False if not compatible
		return(False)