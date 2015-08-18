"""
T7 model
========

A model of T7 life cycle in E.coli host.  The model uses the Stochastic Simulation Algorithm from Gibson & Bruck (2000) [1].

[1] Gibson M.A and Bruck J. (2000), "Exact stochastic simulation of Chemical Systems with Many Species and Many Channels", J.Phys. Chem. 104:1876-1889

"""
import os
import sys
import pandas as pd
import argparse
import numpy as np
import copy as copy
from itertools import *
import pickle


class Species(object):

	def __init__(self, name, count=False):
		"""Initiate the species class.

		Parameters
		----------
		name : str
			   The name of the species, i.e. DNA, mRNA, etc.

		count : int, optional
				The quantity of the species, defualts to zero.

		Returns
		-------
		Species object
		"""
		assert type(name) == str
		self.name = name
		if count:
			self.count = count
		else:
			self.count = 0

	
	@property
	def count(self):
		return self.c
	
	@count.setter
	def count(self, count):
		self.c = count

	def molar(self):
		self.molar_conc = (self.count/float(6.02*1e23) * (1/float(1e-15)))


class DNA(Species):
	def __init__(self, name, tag_dna, dna_length=False, count=False, nut_sequence=False, termination_sequence=False):
		 


		"""Initiate the species class.

		Parameters
		----------
		name : str
			   The name of the species, i.e. DNA, mRNA, etc.

		count : int, optional
				The quantity of the species, defualts to zero.

		tag_dna : int
				Integer representing promoter id.

		dna_length : int
				Total nucleotide length (default=False)

		nut_sequence : list
				First and Second elements nucleotide position. Third element tuple of reaction constants. Ex. [1000, 1060, (0..1,0.2,0.3,0.4)] 

		termination_sequence : list
				First and Second elements nucletide position. Third element tuple of reaction constants. Ex. [1000, 1060, (0.1,0.2,0.7)]

				

		Returns
		-------
		DNA object
		"""


		assert type(tag_dna) == int
		if nut_sequence:
			assert type(nut_sequence) == list
		if termination_sequence:
			assert type(termination_sequence) == list
		Species.__init__(self, name, count=False)
		self.tag_dna = tag_dna
		self.dna_length = dna_length
		self.count = count
		self.nut_sequence = nut_sequence
		self.termination_sequence = termination_sequence


class RNA(Species):
	def __init__(self, name, tag_rna, rna_length, promoter_affiliation, pre_termination=False, count=False):

		"""Initiate the RNA class.

		Parameters
		----------
		name : str
			   The name of the species, i.e. DNA, mRNA, etc.

		count : int, optional
				The quantity of the species, defualts to zero.

		tag_rna : int
				Int corresponding to Protein Object.

		rna_length : int
				Total ribonucletide length of segment.

		promoter_affiliation : list
				List which includes tag_dna as element.

		pre_termination : bool
				If true, increase count of RNA species in the event of termination
		Returns
		-------
		Species object
		"""
		
		assert type(tag_rna) == int
		assert type(rna_length) == int
		assert type(promoter_affiliation) == list
		assert type(pre_termination) == bool
			
		Species.__init__(self, name, count=False)
		self.tag_rna = tag_rna
		self.rna_length = rna_length
		self.promoter_affiliation = promoter_affiliation
		self.pre_termination = pre_termination
		self.count = count
		

class Protein(Species):
	
	def __init__(self, name, tag_rna, count=False):
		"""Initiate the RNApol class.

		Parameters
		----------
		name : str
			   The name of the species, i.e. DNA, mRNA, etc.

		tag_rna : int
				Corresponds to the int value of tag_rna for RNA object 

		count : int, optional
				The quantity of the species, defualts to zero.

		Returns
		-------
		Protein object

		See Also
		--------
		Species : This class inherits from the species class

		"""
		assert type(tag_rna) == int
		Species.__init__(self, name, count=False)
		self.tag_rna = tag_rna

class Reaction(object):

	def __init__(self, reactants, products, ks):

		"""Initiate the RNApol class.

		Parameters
		----------
		reactants : list
			   list elements include Species Objects corresponding to reactants. Stoichiometry represented by multiple entries

		products : list
				list elements include Species Objects corresponding to products 

		ks : tuple
				forward and reverse reaction constants.

		Returns
		-------
		Reaction object

		See Also
		--------
		Species : This class inherits from the species class

		"""
		#assert type(reactants[0]) == Species
		#assert type(products[0]) == Species
		assert type(ks) == tuple
		self.reactants = reactants
		self.products = products
		self.ks = ks
		self.tau = 0 #time for reaction to occur
		self.tau_old = 0 #time that is recorded if self.tau not equals to 0 
		self.prop_old = 0
	

	
	@property
	def time(self):
		return self.t

	@time.setter
	def time(self,time):
		self.t = time #time is the tau for the reaction that has just been executed"""
   
	def get_propensity(self):

		""" 
		calculates likelihood of reaction occurring.
		
		"""

		multiplicity_dict = {}
		unique_elements = set(self.reactants)
		for i in unique_elements:
			elements = []
			for j in self.reactants:
				if j == i:
					elements.append(j)
			multiplicity_dict.update({i:elements})

		"""
		The dictionary has been created. The value of each key represents a list with each element identical to the key.
		If the number of elements inside the list is greater than one (corresponding to a stoichiometry greater than one),
		then the Combinatorics formulae (nCr) is applied.  
		
		"""   

		penultimate = []

		for i in multiplicity_dict.keys():
			if len(multiplicity_dict[i]) > 1:
				value = 1
				for j in range(len(multiplicity_dict[i])):
					value *= multiplicity_dict[i][j].c - j 

				extend = [j for j in range(len(multiplicity_dict[i])) if j]
				extend.append(len(multiplicity_dict[i]))
				factorial = 1
				for j in extend:
					factorial *= j
				intermediate = value/float(factorial) #application of the combinatoric formulae

				penultimate.append(intermediate)

			else:
				penultimate.append(multiplicity_dict[i][len(multiplicity_dict[i]) - 1].c)

		""" 
		At the end of the for loop, the penultimate list will contain values pertaining to the adjusted stoichiometic count for  
		each species. The next step is to iterate through this list and multiply all values in the list with the rate constant. 

		"""
	   
		a = 1
		for i in penultimate:
			a *= i
		a *= self.ks[0]
		self.prop = a

		""" 
		At the end of the procedure, the propensity value is stored in self.prop. If the value is non-zero, it is also 
		stored in self.prop_old (this is necessary to retrive a propensity value for a reaction that was previously disabled
		according to the Next Reaction Method)

		"""


	def get_tau(self, time): #this method is to recalculate the tau for the reaction that has just been executed

		"""
		This method is necessary to recalculate tau for the reaction that has just been executed according to the
		Next Reaction Method. If the reaction that has just been executed has an updated propensity value of zero, then
		the tau value is stored in self.tau_old to await recall when the propensity becomes non-zero again. The moment the 
		value becomes non-zero, a random number is drawn from a uniform distribution and the tau value is recalculated. This
		delayed random draw is not explicitly stated in the Gibson and Brucks Next Reaction Method, but is inferred based
		upon the principles of the Next Reaction Method.  
		"""
		if self.prop == 0:
			self.tau_old = self.tau
			self.t = time
			self.tau = np.inf
		else:
			 self.tau = ((-1)*(np.log(np.random.random())))/float(self.prop) + time
			 self.tau_old = self.tau
			 self.prop_old = self.prop

		#self.tau = tau

	def get_det_tau(self, time): #this method is to recalculate the tau for reactions that are dependent on the just executed reaction
		"""
		This method calculates the taus of reactions that are affected by the execution of the previous reaction according to the NRM.
		The "deterministic" calculation can be found in Gibson and Brucks Next Reaction Method. For reactions that were previously
		set to zero, an alternate formulae is used to get the tau value.
		"""
		if self.prop == 0:

			if self.tau != np.inf: 
				self.tau_old = self.tau
				self.t = time
				self.tau = np.inf
		else:
			if self.tau == np.inf:
				if self.prop_old == 0:
					self.tau = ((-1)*(np.log(np.random.random())))/float(self.prop) + time
				
				elif self.tau_old == self.t:
					self.tau = ((-1)*(np.log(np.random.random())))/float(self.prop) + time
				
				else:
					self.tau = (self.prop_old/float(self.prop))*(self.tau_old - self.t) + time
			else:
				self.tau = ((self.prop_old/float(self.prop))*((self.tau) - time)) + time
			self.prop_old = self.prop
			self.tau_old = self.tau

	

class NextReactionMethod(object):
	
	def __init__(self, k_elong=False, k_translation_elong=False):
		"""
		k_elong : float

		k_translation_elong : float

		"""
		self.k_elong = k_elong
		self.k_translation_elong = k_translation_elong


	def create_rec_list(self, rec_list): #creates the reaction list
		assert type(rec_list[0]) == Reaction
		self.temp_reactions = rec_list


		
	def create_species_list(self): #creates a list of species objects - this does not include transcription and translation reactions
	   
		#Creacting a list of species
		
		arch_list = []

		for i in self.temp_reactions:
			temp = i.reactants + i.products
			arch_list.append(temp)

		reference = arch_list[0]

		for i in arch_list:
			if arch_list.index(i) > 0:
				catch = i + reference
				reference = catch

		species_set = set(reference)

		self.species = species_set

		store_react = []
		for i in self.temp_reactions:
			hit = 0 
			for j in i.reactants:
				if type(j) == DNA:
					hit += 1
			if hit == 0:
				store_react.append(i)

		self.reactions = store_react


	   

	def create_transcription_elongation_species(self):
		"""
		Creating the transcription elongation consists of the following steps. The sequential termination and NUT 
		sites have to be modeled along the DNA. Two properties of the object include the termination sequence and the NUT sequence.
		Therefore, the species objects must include RNApolymerase on NUT sites, RNApolymerase + N protein on Nut sites,
		RNApolymerase on termination sites, RNApolymerase + N protein on termination sites. This method creates a dictionary
		of species objects (DNA species objects in this simulation are actually promoters) as keys and values as a list of the
		individual nucleotides associated with the DNA. These nucleotides include all possible configuration of proteins on that 
		nucleic acid (RNApolymerase + N protein for example). If there is no need to model NUT and termination sites, then the
		number of species objects in the value list is simply the length of the transcript. If species and nut sites are included,
		then the value for that DNA object is a list of list with each list containing species objects of different nucleotide
		length and protein make-up (one list will contain all the species that convey RNApolymerase + N protein translocating on the 
		NUT sites). The other list will contain species objects pertaining to the termination sites including the possible combinations
		of RNApolymerase and N protein on these sites.
		"""
		elong_species_dict = {} 
		for i in self.species:
			if type(i) == DNA:
				if i.nut_sequence and i.termination_sequence: #not all DNA objects will have NUT sites and termination sites
					elong_list = []
					for j in range(i.dna_length):
						if j < i.nut_sequence[0]:
							name = str(i.name) + '-' + str(j)
							elong_list.append(Species(name,count=0))
						elif i.nut_sequence[0] <= j <= i.nut_sequence[1]:
							name = str(i.name) + '-' + 'Nut' + '-' + str(j)
							elong_list.append(Species(name,count=0))
						elif i.termination_sequence[0] <= j <= i.termination_sequence[1]:
							name = str(i.name) + '-' + 'TR' + '-' + str(j)
							elong_list.append(Species(name,count=0))
						else:
							name = str(i.name) + '-' + str(j)
							elong_list.append(Species(name,count=0))
					elong_species_dict.update({i:[elong_list]})
				else: 
					elong_list = []
					for j in range(i.dna_length):
						name = str(i.name) + '-' + str(j)
						elong_list.append(Species(name,count=0))
					elong_species_dict.update({i:elong_list})
		self.transcription_elong_species = elong_species_dict


		for i in self.transcription_elong_species.keys():
			if i.nut_sequence and i.termination_sequence:
				nut_list = []
				term_list = []
				for j in range(i.nut_sequence[0],i.nut_sequence[1]+1):
					if j == i.nut_sequence[0]:
						name = str(i.name) + '-' + 'Nut' + '-' + str(j)
						nut_list.append(Species(name,count=0))
					else:
						name = str(i.name) + '-' + 'N' + '-' + 'Nut' + '-' + str(j)
						nut_list.append(Species(name,count=0))
				self.transcription_elong_species[i].append(nut_list)

				term_list.append(nut_list[len(nut_list) - 1])
				for j in range(i.termination_sequence[0],i.dna_length):
					if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
						name = str(i.name) + '-' + 'N' '-' + 'TR' + '-' + str(j)
						term_list.append(Species(name,count=0))
					else:
						name = str(i.name) + '-' + 'N' + '-' + str(j)
						term_list.append(Species(name,count=0))
				self.transcription_elong_species[i].append(term_list)

	def create_trancription_elongation_reactions(self):
		"""
		The transcription elongation reactions are represented in a dictionary. The keys in the dictionary represent the DNA species objects
		and the values mirror the list of list for the species odbjects representing the different possible configurations of RNApolymerase
		on the DNA. At the end of the code, the dictionary is exhaustively compiled into a single list of reactions representing all
		possible reactions associated with a given DNA object (including NUT sites and termination). This list (along with the
		translation elongation reactions) are then appended to the self.reactions list and ready for the Next Reaction Method Execution. 
		"""
		elong_reactions_dict = {}
		for i in self.transcription_elong_species.keys():
			if i.nut_sequence and i.termination_sequence:
				for j in range(len(self.transcription_elong_species[i])):
					Y = self.transcription_elong_species[i][j]
					if j == 0:
						elong_react = []
						for k in self.species:
							if k.name.split('-')[0] == 'Open':
								if k.name.split('-')[1] == i.name.split('-')[0]:
									for l in self.species:
										if type(l) == DNA:
											if l.tag_dna == i.tag_dna:
												elong_react.append(Reaction([k],[l,Y[0]],(self.k_elong,self.k_elong)))
						for k in range(len(Y)):
							if k < i.nut_sequence[0]:
								elong_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
							elif i.nut_sequence[0] <= k <= i.nut_sequence[1]:
								Z = self.transcription_elong_species[i][1]
								for l in range(len(Z)):
									if 0 <= l < (len(Z) - 1):
										if Z[l].name.split('-')[len(Z[l].name.split('-')) - 1] == Y[k].name.split('-')[len(Y[k].name.split('-')) - 1]:
											for m in self.species:
												if m.name == 'N':
													elong_react.append(Reaction([Y[k],m],[Z[l+1]],(i.nut_sequence[2][1],i.nut_sequence[2][1])))
													elong_react.append(Reaction([Z[l+1]],[Y[k],m],(i.nut_sequence[2][2],i.nut_sequence[2][2])))
								elong_react.append(Reaction([Y[k]],[Y[k+1]],(i.nut_sequence[2][0],i.nut_sequence[2][0])))
							elif i.termination_sequence[0] <= k <= i.termination_sequence[1]:
								for l in self.species:
									if l.name == 'RNAP':
										for m in self.species:
											if type(m) == RNA:
												for n in m.promoter_affiliation:
													if n == i.tag_dna:
														if m.pre_termination == True:
															elong_react.append(Reaction([Y[k]],[Y[k+1]],(i.termination_sequence[2][0],i.termination_sequence[2][0])))
															elong_react.append(Reaction([Y[k]],[l,m],(i.termination_sequence[2][2],i.termination_sequence[2][2])))							
							elif i.termination_sequence[1] <= k < (i.dna_length - 1):
								elong_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
							else:
								products = []
								for l in self.species:
									if l.name == 'RNAP':
										products.append(l)
									elif type(l) == RNA:
										for m in l.promoter_affiliation:
											if m == i.tag_dna:
												products.append(l)
								elong_react.append(Reaction([Y[k]],products,(self.k_elong,self.k_elong)))
						elong_reactions_dict.update({i:[elong_react]})
					elif j == 1:
						elong_nut_react = []
						for k in range(len(Y)):
							if 0 < k < (len(Y) - 1):
								elong_nut_react.append(Reaction([Y[k]],[Y[k+1]],(i.nut_sequence[2][3],i.nut_sequence[2][3])))
						elong_reactions_dict[i].append(elong_nut_react)
					elif j == 2:
						elong_term_react = []
						for k in range(len(Y)):
							if k < ((i.termination_sequence[1] - i.termination_sequence[0]) + 1):
								elong_term_react.append(Reaction([Y[k]],[Y[k+1]],(i.termination_sequence[2][1],i.termination_sequence[2][1])))
							elif ((i.termination_sequence[1] - i.termination_sequence[0]) + 1) <= k < (len(Y) - 1):
								elong_term_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
							elif k == (len(Y) - 1):
								products = []
								for l in self.species:
									if l.name == 'RNAP':
										products.append(l)
									elif l.name == 'N':
										products.append(l)
									elif type(l) == RNA:
										for m in l.promoter_affiliation:
											if m == i.tag_dna:
												products.append(l)
								elong_term_react.append(Reaction([Y[k]],products,(self.k_elong,self.k_elong)))
						elong_reactions_dict[i].append(elong_term_react)
			else:
				Y = self.transcription_elong_species[i]
				elong_reactions = []
				for j in range(len(Y)):
					if j == 0:
						for k in self.species:
							if k.name.split('-')[0] == 'Open':
								if k.name.split('-')[1] == i.name.split('-')[0]:
									for l in self.species:
										if type(l) == DNA:
											if l.tag_dna == i.tag_dna:
												elong_reactions.append(Reaction([k],[l,Y[j]],(self.k_elong,self.k_elong)))
					elif 0 < j < (len(Y) - 1):
						elong_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_elong,self.k_elong)))
					elif j == (len(Y) - 1):
						products = []
						for k in self.species:
							if k.name == 'RNAP':
								products.append(k)
							elif type(k) == RNA:
								for l in k.promoter_affiliation:
									if l == i.tag_dna: 
										products.append(k)
						elong_reactions.append(Reaction([Y[j]],products,(self.k_elong,self.k_elong)))
				elong_reactions_dict.update({i:[elong_reactions]})
			self.elong_reactions_dict = elong_reactions_dict

		temp_list = []
		for i in elong_reactions_dict.keys():
			if i.nut_sequence and i.termination_sequence:
				aggregate_list = np.sum(elong_reactions_dict[i])
				temp_list.append(aggregate_list)
			else:
				temp_list.append(elong_reactions_dict[i][0])

		for i in range(len(temp_list)):
			hold = temp_list[0]
			if i > 0:
				hold += temp_list[i]
		self.transcription_elong_reactions = hold



	def create_translation_elongation_species(self):
		translation_species_dict = {}
		for i in self.species:
			if type(i) == RNA:
				translation_species_list = []
				for j in range(i.rna_length):
					name = str(i.name) + '-' + str(j)
					translation_species_list.append(Species(name,count=0))
				translation_species_dict.update({i:translation_species_list})
		self.translation_elong_species = translation_species_dict

	def create_translation_elongation_reactions(self):
		translation_react_dict = {}
		for i in self.translation_elong_species.keys():
			translation_reactions = []
			Y = self.translation_elong_species[i]
			for j in range(len(self.translation_elong_species[i])):
				if j == 0:
					temp = []
					valor = []
					for k in self.species:
						if type(k) == RNA:
							if k.tag_rna == i.tag_rna:
								valor.append(k)
						elif k.name == 'Ribosome':
							temp.append(k)
						elif k.name == 'ElongationRibosome':
							if k.name.split()[2] == i.tag_rna:
								temp.append(k)
					valor.append(Y[j])
					translation_reactions.append(Reaction(temp,valor,(self.k_translation_elong,self.k_translation_elong)))
					translation_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_translation_elong,self.k_translation_elong)))
				elif 0 < j < (len(Y) - 1):
					translation_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_translation_elong,self.k_translation_elong)))
				elif j == (len(Y) - 1):
					for k in self.species:
						if type(k) == Protein:
							if k.tag_rna == i.tag_rna:
								for l in self.species:
									if l.name == 'Ribosome':
										translation_reactions.append(Reaction([Y[j]],[k,l],(self.k_translation_elong,self.k_translation_elong)))
			translation_react_dict.update({i:translation_reactions})

		Y = translation_react_dict.values()
		for i in range(len(Y)):
			hold = Y[0]
			if i > 0:
				hold += Y[i]
		self.translation_elong_reactions = hold


	def update_reaction_list(self):
		self.total_reactions = self.reactions + self.transcription_elong_reactions + self.translation_elong_reactions


	def generate_dep_graph(self):
		self.dep_graph = {}
		for i in self.reactions:
			temp = i.reactants + i.products
			dep_rec = []
			for j in self.reactions:
				set_value = [val for val in temp if val in j.reactants]
				if len(set_value) != 0:
					dep_rec.append(j)
			self.dep_graph[i] = dep_rec
		#return self.dep_graph


	def PR_PRM_model_dependencies1(self):
		#Use the positions of RNAP on the list(either position 1 or 3 - zero based indexing - to determine which reactions are affected
		config = np.array([[1, 0, 0, 0, 0.0], [2, 0, 0, 'R', -11.7], [3, 0, 'R', 0, -10.1], [4, 'R', 0, 0, -10.1], [5, 0, 0, 'C', -10.8],
		[6, 0, 'C', 0, -10.8], [7, 'C', 0, 0, -12.1], [8,'RNAP', 0, 0, -11.5], [9, 0, 0, 'RNAP', -12.5], 
		[10, 0, 'R', 'R', -23,7], [11, 'R',  0, 'R', -21.8], [12, 'R', 'R', 0, -22.2], [13, 0, 'C', 'C', -21.6],
		[14, 'C', 0, 'C', -22.9], [15, 'C', 'C', 0, -22.9], [16, 'RNAP', 0, 'RNAP', -24], [17, 0, 'C', 'R', -22.5],
		[18, 0, 'R', 'C', -20.9], [19, 'R', 0, 'C', -20.9], [20, 'C', 0, 'R', -23.8], [21, 'R', 'C', 0, -20.9], 
		[22, 'C', 'R', 0, -22.2], [23, 'R', 0, 'RNAP', -22.6], [24, 'RNAP', 'R', 0, -21.6], [25, 'RNAP', 0, 'R', -23.2],
		[26, 'C', 0, 'RNAP', -24.6], [27, 'RNAP', 'C', 0, -22.3], [28, 'RNAP', 0, 'C', -22.3], [29, 'R', 'R', 'R', -33.8], 
		[30, 'C', 'C', 'C', -33.7], [ 31, 'C', 'R', 'R', -35.8], [32, 'R', 'C', 'R', -32.6], [33, 'R', 'R', 'C', -33.0], 
		[34, 'R', 'C', 'C', -31.7], [35, 'C', 'R', 'C', -33.0], [36, 'C', 'C', 'R', -34.6], [37, 'RNAP', 'R', 'R', -35.2], 
		[38, 'RNAP', 'C', 'C', -33.1], [39, 'RNAP', 'C', 'R', -34.0], [40, 'RNAP', 'R', 'C', -32.4]])
		
		occupancy_species_dict = {}
		"""
		Creates a dictionary which relates a particular occupancy state with a list of relevant species.  
		"""
		for i in range(len(config)):
			occupancy_decrease = []
			occupancy_increase = []
			for j in range(len(config[i])):
				if config[i][j] == 'RNAP':
					if j == 1:
						for k in self.species:
							if k.name == 'PRM':
								occupancy_decrease.append(k)
							elif k.name == 'RNAP':
								occupancy_decrease.append(k)
							elif k.name.split('-')[0] == 'Closed':
								if k.name.split('-')[1] == 'PRM':
									occupancy_increase.append(k)

					elif j == 3:
						for k in self.species:
							if k.name == 'PR':
								occupancy_decrease.append(k)
							elif k.name == 'RNAP':
								occupancy_decrease.append(k)
							elif k.name.split('-')[0] == 'Closed':
								if k.name.split('-')[1] == 'PR':
									occupancy_increase.append(k)
				else:
					if type(config[i][j]) != float:
						for k in self.species:
							if k.name == config[i][j]:
								occupancy_decrease.append(k)

			if len(occupancy_increase) > 0:
				occupancy_species_dict.update({i:(occupancy_increase,occupancy_decrease)})
			else:
				occupancy_species_dict.update({i:occupancy_decrease})
		self.occupancy_species1 = occupancy_species_dict

		
		
		occupancy_reaction_dict = {}
		"""
		Creates a dictionary which relates a particular occupancy state with the Gillespie reactions (self.reactions) that are affected. 
		There are no Gillespie reactions explicitly associated with RNApolymerase (RNAP) or DNA objects(promoters).
		These are 'implicit' reactions as they are alluded to in the stat_thermodynamic model for operator/promoter binding.
		"""
		
		X = copy.deepcopy(occupancy_species_dict)

		for i in range(len(X)):
			if type(X[i]) == tuple:
				final_list = X[i][0]
				for j in range(len(X[i])):
					if j > 0:
						final_list += X[i][j]
				reaction_occupancy_list = list(set(final_list))
			else:
				reaction_occupancy_list = list(set(X[i]))

			reaction_list = []
			for j in reaction_occupancy_list:
				for k in self.reactions:
					for l in k.reactants:
						if j.name == l.name:
							reaction_list.append(k)
			occupancy_reaction_dict.update({i:reaction_list})
		self.occupancy_reaction1 = occupancy_reaction_dict
		

	def PR_PRM_stat_energy_model_selection1(self):
		config = np.array([[1, 0, 0, 0, 0.0], [2, 0, 0, 'R', -11.7], [3, 0, 'R', 0, -10.1], [4, 'R', 0, 0, -10.1], [5, 0, 0, 'C', -10.8],
		[6, 0, 'C', 0, -10.8], [7, 'C', 0, 0, -12.1], [8,'RNAP', 0, 0, -11.5], [9, 0, 0, 'RNAP', -12.5], 
		[10, 0, 'R', 'R', -23,7], [11, 'R',  0, 'R', -21.8], [12, 'R', 'R', 0, -22.2], [13, 0, 'C', 'C', -21.6],
		[14, 'C', 0, 'C', -22.9], [15, 'C', 'C', 0, -22.9], [16, 'RNAP', 0, 'RNAP', -24], [17, 0, 'C', 'R', -22.5],
		[18, 0, 'R', 'C', -20.9], [19, 'R', 0, 'C', -20.9], [20, 'C', 0, 'R', -23.8], [21, 'R', 'C', 0, -20.9], 
		[22, 'C', 'R', 0, -22.2], [23, 'R', 0, 'RNAP', -22.6], [24, 'RNAP', 'R', 0, -21.6], [25, 'RNAP', 0, 'R', -23.2],
		[26, 'C', 0, 'RNAP', -24.6], [27, 'RNAP', 'C', 0, -22.3], [28, 'RNAP', 0, 'C', -22.3], [29, 'R', 'R', 'R', -33.8], 
		[30, 'C', 'C', 'C', -33.7], [ 31, 'C', 'R', 'R', -35.8], [32, 'R', 'C', 'R', -32.6], [33, 'R', 'R', 'C', -33.0], 
		[34, 'R', 'C', 'C', -31.7], [35, 'C', 'R', 'C', -33.0], [36, 'C', 'C', 'R', -34.6], [37, 'RNAP', 'R', 'R', -35.2], 
		[38, 'RNAP', 'C', 'C', -33.1], [39, 'RNAP', 'C', 'R', -34.0], [40, 'RNAP', 'R', 'C', -32.4]])
		
		boltzman = (1.9872041 * (10**(-3)))
		state_energy_list = []
		probability_configuration_list = []

		for i in config:
			R_power = 0
			C_power = 0
			RNAP_power = 0
			for j in range(len(i)):
				if i[j] == 'C':
					C_power += 1
				elif i[j] == 'R':
					R_power += 1
				elif i[j] == 'RNAP':
					RNAP_power += 1
			for j in self.species:
				if j.name == 'C':
					for k in self.species:
						if k.name == 'R':
							for l in self.species:
								if l.name == 'RNAP': #when the volume function is finished, change count to concentration
									state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(k.molar_conc**R_power)*(j.molar_conc**C_power)*(l.molar_conc**RNAP_power))
			state_energy_list.append(state_energy)



		for i in range(len(config)):
			for j in range(len(config[i])):
				if config[i][j] == 'RNAP':
					if j == 1:
						for k in self.species:
							if k.name == 'PRM':
								if k.count == 0:
									state_energy_list[i] == 0
					elif j == 3:
						for k in self.species:
							if k.name == 'PR':
								if k.count == 0:
									state_energy_list[i] == 0
		for i in state_energy_list:
			prob = i/float(np.sum(state_energy_list))
			probability_configuration_list.append(prob)

		Y = np.random.multinomial(1,probability_configuration_list)
		X = np.argmax(Y)

		self.current_PR_PRM_config1 = X


	def PR_PRM_model_config_update1(self):
		
		Z = self.current_PR_PRM_config1
		
		if type(self.occupancy_species1[Z]) == tuple:
			for i in self.occupancy_species1[Z][0]:
				i.count += 1
			for i in self.occupancy_species1[Z][1]:
				i.count -= 1
		else:
			for i in self.occupancy_species1[Z]:
				i.count -= 1
		for i in self.occupancy_reaction1[Z]:
			i.get_propensity()

	def PRE_model_dependencies2(self):
		"This is the statistical binding model for the PRE promoter"
		config2 = np.array([[0,0,0.0],[0,'RNAP',-9.9],['CII',0,-9.7],['CII','RNAP',-21.5]])

		species_config2_dict = {}
		for i in range(len(config2)):
			species_increase = []
			species_decrease = []
			for j in config2[i]:
				if j == 'RNAP':
					for k in self.species:
						if k.name.split('-')[0] == 'Closed':
							if k.name.split('-')[1] == 'PRE':
								species_increase.append(k)
						elif k.name == 'PRE':
							species_decrease.append(k)
						elif k.name == j:
							species_decrease.append(k)
				else:
					if type(j) != float:
						for k in self.species:
							if type(k) == Protein:
								if k.name == j:
									species_decrease.append(k)
			if len(species_increase) > 0:
				species_config2_dict.update({i:(species_increase,species_decrease)})
			else:
				species_config2_dict.update({i:species_decrease})
		self.occupancy_species2 = species_config2_dict

		X = copy.deepcopy(species_config2_dict)

		reaction_config2_dict = {}
		for i in X.keys():
			if type(X[i]) == tuple:
				final_list = X[i][0]
				for j in range(len(X[i])):
					if j > 0:
						final_list += X[i][j]
				reaction_occupancy_list = list(set(final_list))
			else:
				reaction_occupancy_list = list(set(X[i]))

			reaction_list = []
			for j in reaction_occupancy_list:
				for k in self.reactions:
					for l in k.reactants:
						if l.name == j.name:
							reaction_list.append(k)
			reaction_config2_dict.update({i:reaction_list})
		self.occupancy_reaction2 = reaction_config2_dict


	def PRE_stat_energy_model_selection2(self):
		config2 = np.array([[0,0,0.0],[0,'RNAP',-9.9],['CII',0,-9.7],['CII','RNAP',-21.5]])

		boltzman = (1.9872041 * (10**(-3)))
		state_energy_list = []
		probability_configuration_list = []

		for i in config2:
			RNAP_power = 0
			CII_power = 0
			for j in range(len(i)):
				if i[j] == 'RNAP':
					RNAP_power += 1
				elif i[j] == 'CII':
					CII_power += 1
			for j in self.species:
				if j.name == 'CII':
					for k in self.species:
						if k.name == 'RNAP':
							state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(j.molar_conc**CII_power)*(k.molar_conc**RNAP_power))
			state_energy_list.append(state_energy)


		for i in range(len(config2)):
			for j in range(len(config2[i])):
				if config2[i][j] == 'RNAP':
					for k in self.species:
						if k.name == 'PRE':
							if k.count == 0:
								state_energy_list[i] = 0

		for i in state_energy_list:
			prob = i/float(np.sum(state_energy_list))
			probability_configuration_list.append(prob)

		Y = np.random.multinomial(1,probability_configuration_list)
		X = np.argmax(Y)

		self.current_PRE_config2 = X

	def PRE_model_config_update2(self):
		Z = self.current_PRE_config2

		if type(self.occupancy_species2[Z]) == tuple:
			for j in self.occupancy_species2[Z][0]:
				j.count += 1
			for j in self.occupancy_species2[Z][0]:
				j.count -= 1
		else:
			for j in self.occupancy_species2[Z]:
				j.count -= 1

		for i in self.occupancy_reaction2[Z]:
			i.get_propensity()


	def PL_model_dependencies3(self):
		config3 = np.array([[0,0,0.0],['C',0,-10.9],[0,'C',-12.1],['R',0,-11.7],[0,'R',-10.1],[0,'RNAP',-12.5],
		['C','C',-22.9],['C','R',-20.9],['R','C',-22.8],['R','R',-23.7]])

		species_config3_dict = {}
		for i in range(len(config3)):
			species_increase = []
			species_decrease = []
			for j in config3[i]:
				if j == 'RNAP':
					for k in self.species:
						if k.name.split('-')[0] == 'Closed':
							if k.name.split('-')[1] == 'PL':
								species_increase.append(k)
						elif k.name  == 'RNAP':
							species_decrease.append(k)
						elif k.name == 'PL':
							species_decrease.append(k)

				else:
					if type(j) != float:
						for k in self.species:
							if k.name == j:
								species_decrease.append(k)

			if len(species_increase) > 0:
				species_config3_dict.update({i:(species_increase,species_decrease)})
			else:
				species_config3_dict.update({i:species_decrease})
		self.occupancy_species3 = species_config3_dict

		X = copy.deepcopy(species_config3_dict)

		reaction_config3_dict = {}
		for i in X.keys():
			if type(X[i]) == tuple:
				for j in range(len(X[i])):
					final_list = X[i][0]
					if j > 0:
						final_list += X[i][j]
					reaction_occupancy_list = list(set(final_list))
			else:
				reaction_occupancy_list = list(set(X[i]))

			reaction_list = []
			for j in reaction_occupancy_list:
				for l in self.reactions:
					for k in l.reactants:
						if j.name == k.name:
							reaction_list.append(l)
			reaction_config3_dict.update({i:reaction_list})
		self.occupancy_reaction3 = reaction_config3_dict


	def PL_stat_energy_model_selection3(self):
		config3 = np.array([[0,0,0.0],['C',0,-10.9],[0,'C',-12.1],['R',0,-11.7],[0,'R',-10.1],[0,'RNAP',-12.5],
		['C','C',-22.9],['C','R',-20.9],['R','C',-22.8],['R','R',-23.7]])

		boltzman = (1.9872041 * (10**(-3)))
		state_energy_list = []
		probability_configuration_list = []

		for i in config3:
			R_power = 0
			C_power = 0
			RNAP_power = 0
			for j in range(len(i)):
				if j == 'C':
					C_power += 1
				elif j == 'R':
					R_power += 1
				elif j == 'RNAP':
					RNAP_power += 1
			for j in self.species:
				if j.name == 'R':
					for k in self.species:
						if k.name == 'C':
							for l in self.species:
								if l.name == 'RNAP':
									state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(k.molar_conc**C_power)*(j.molar_conc**R_power)*(l.molar_conc**RNAP_power))
			state_energy_list.append(state_energy)
		
		for i in range(len(config3)):
			for j in config3[i]:
				if j == 'RNAP':
					for k in self.species:
						if k.name == 'PL':
							if k.count == 0:
								state_energy_list[i] = 0

		for i in state_energy_list:
			prob = i/float(np.sum(state_energy_list))
			probability_configuration_list.append(prob)

		Y = np.random.multinomial(1,probability_configuration_list)

		X = np.argmax(Y)

		self.current_PL_config3 = X

	def PL_model_config_update3(self):

		Z = self.current_PL_config3

		if type(self.occupancy_species3[Z]) == tuple:
			for j in self.occupancy_species3[Z][0]:
				j.count += 1
			for j in self.occupancy_species3[Z][1]:
				j.count -= 1
		else:
			for j in self.occupancy_species3[Z]:
				j.count -= 1

		for i in self.occupancy_reaction3[Z]:
			i.get_propensity()

	def stat_occup_change(self, prior_PR_PRM_configuration, prior_PRE_configuration, prior_PL_configuration, tau_list, system_time):

		
		L = [self.current_PR_PRM_config1, self.current_PRE_config2, self.current_PL_config3]

		J = [prior_PR_PRM_configuration, prior_PRE_configuration, prior_PL_configuration]

		M = [self.occupancy_species1, self.occupancy_species2, self.occupancy_species3]

		N = [self.occupation_reaction1, self.occupation_reaction2, self.occupation_reaction3]

		Get = [val for val in zip(L,J,M,N) if val[0] != val[1]]

		for i in Get:
			for j in i[2][i[1]]:
				if type(i[2][i[1]]) == tuple:
					for k in i[2][i[1]][0]:
						if k.count != 0:
							k.count -= 1
							for l in i[2][i[1]][1]:
								l.count += 1
				else:
					for k in i[2][i[1]]:
						k.count += 1
			for j in i[3][i[1]]:
				j.get_propensity()
				for k in range(len(self.reactions)):
					if j == self.reactions[k]:
						j.get_det_tau(system_time)
						tau_list[k] = j.tau


	def tau_current_occup_update(self, tau_list, system_time):


		X = self.current_PR_PRM_config1 

		Y = self.current_PRE_config2 

		Z = self.current_PL_config3 

		alpha = [X, Y, Z]

		beta = [self.occupancy_reaction1, self.occupancy_reaction2, self.occupancy_reaction3]

		
		Get = [val for val in zip(alpha, beta)]

		for i in Get:
			for j in i[1][i[0]]:
				for k in range(len(self.reactions)):
					if j == self.reactions[k]:
						j.get_det_tau(system_time)
						print j.prop, j.tau
						tau_list[k] = j.tau



	def NRM_execution(self, step):

		self.create_transcription_elongation_species()

		self.create_trancription_elongation_reactions()

		self.create_translation_elongation_species()

		self.create_translation_elongation_reactions()

		self.update_reaction_list()

		system_time = 0

		K0 = (4.76)*(1e-18)

		#First Step, create a species dictionary from the species list

		#keys will be species, values will be a list of counts at each iteration of the stochastic algorithm

		species_dict = {}
		for i in self.species:
			i.molar()
			species_dict.update({i.name:[i.molar_conc]})

		species_dict.update({'time':[system_time]})

		tau_list = []
		for m in self.total_reactions:
			m.get_propensity()
			m.get_det_tau(system_time)
			tau_list.append(m.tau)

		for i in self.total_reactions:
			print i.prop

		print tau_list

		self.generate_dep_graph()

		self.PR_PRM_model_dependencies1()

		self.PRE_model_dependencies2()

		self.PL_model_dependencies3()


		for i in range(step):

			self.PR_PRM_stat_energy_model_selection1()

			self.PRE_stat_energy_model_selection2()

			self.PL_stat_energy_model_selection3()

			if i > 0:
				stat_occup_change(prior_PR_PRM_configuration, prior_PRE_configuration, prior_PL_configuration, tau_list, system_time)

			self.PR_PRM_model_config_update1()

			self.PRE_model_config_update2()

			self.PL_model_config_update3()

			prior_PR_PRM_configuration = self.current_PR_PRM_config1

			prior_PRE_configuration = self.current_PRE_config2

			prior_PL_configuration = self.current_PL_config3

			self.tau_current_occup_update(tau_list, system_time)

			for i in self.total_reactions:
				print i.prop

			print tau_list

			reaction_index = np.argmin(tau_list)

			dependency_list = self.dep_graph[self.reactions[reaction_index]]

			system_time = tau_list[reaction_index] #this will give us the time variable input necessary for get_det_tau calculation


			if self.reactions[reaction_index].reactants == self.reactions[reaction_index].products:
				for j in self.total_reactions[reaction_index].reactants:
					new = j.count - 1
					j.count = new

			else:
				for j in self.total_reactions[reaction_index].reactants:
					new = j.count - 1
					j.count = new
				for j in self.total_reactions[reaction_index].products:
					new = j.count + 1
					j.count = new



			for j in dependency_list:
				if j == self.total_reactions[reaction_index]:       
					j.get_propensity()
					
					j.get_tau(system_time)

					tau_list[reaction_index] = j.tau

				else:
					j.get_propensity()

					j.get_det_tau(system_time)

					tau_list[self.reactions.index(j)] = j.tau

			cell_volume = (1 + (K0*system_time))*(1e-15)


			for j in self.species:
				j.molar_conc = (j.count)*(1/float(6.02*1e-23))*(1/float(cell_volume))
				species_dict[j.name].append(j.molar_conc)

			species_dict['time'].append(system_time)

		return species_dict    

		"""
		file = open("home\\geurenam\\simoutput.py",'w')

		pickle.dump(species_dict,file)

		file.close()"""




	def multiple_simulation(self,trajectories,step):

		simulation_dict = {}

		self.create_species_list()

		for i in range(trajectories):

			shadow_reaction = copy.deepcopy(self)

			Y = shadow_reaction.NRM_execution(step)

			simulation_dict.update({i:Y})

		self.stats = simulation_dict
