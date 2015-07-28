"""
T7 model
========

A model of T7 life cycle in E.coli host.  The model uses the Stochastic Simulation Algorithm from Gibson & Bruck (2000) [1].

[1] Gibson M.A and Bruck J. (2000), "Exact stochastic simulation of Chemical Systems with Many Species and Many Channels", J.Phys. Chem. 104:1876-1889

"""
import pandas as pd
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

class DNA(Species):
    def __init__(self, name, tag_dna, dna_length=False, count=False, nut_sequence=False, termination_sequence=False):
        """ All genes will be associated with an id tag comprising of an integer 
        to facillitate generation of reactions associated with transcription elongation

        DNA species will have an associated transcript length and a tag for identification

        nut_sequence represents the N - utilization sites present in the gene. Is represented with a list containing the 
        start length and end length followed by a tuple containing the rate constants. The tuples are important because there are 
        multiple reactions associated with the binding of N-utilization sites to the respective DNA elements (therefore multiple
        rate constants can be accessed via the indice)

		termination_sequence represents the DNA elements in which the POSSIBILITY for dissociation of RNApolymerase from the DNA is
		made available in the set of reactions. If the RNApolymerase + DNA complex is joined by an N-protein, then the possibility
		of dissociation is abrogated.

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
        """
        if count:
            self.count = count
        else:
            self.count = 0
            """

class RNA(Species):
    def __init__(self, name, tag_rna, rna_length, promoter_affiliation, polycistronic=False, count=False):
        """ All mRNA transcript will be associated with an id tag of integer 
        corresponding to expressed gene. 

        tag_rna does NOT correspond to the DNA segment it comes from (tag_dna). The tag_rna is a unique integer assigned
        to a particular RNA that is present in the viral polycistronic RNA transcript. It directly corresponds to the tag_protein
        which equates the RNA segment (Cro for example) to the protein that results from translation of that RNA transcript.


        rna_length is an integer that denotes the total number of ribonucleotides for the particular segment of the polycistronic
        transcript.

        promoter_affiliation is a list that contains all the DNA/promoter integer tags the RNA is transcribed from.
        Certain genes on the lambda phage genome (such as CI) are read/ transcribed from multiple promoters. This
        list allows the RNA count to keep track of transcription of the gene product from multiple promoters.

        polycistronic is an integer number (starting from one) corresponding to the the location of the RNA segment on the
        polycistronic RNA transcript. This allows us to keep track of which portion of the polycistronic RNA has been transcribed 
        prior to a termination site. For example, a polycistronic number of 1 stipulates that the RNA segment is before the
        termination_sequence and the count of that associated RNA is increased if there is a termination event.

        """

        assert type(tag_rna) == int
        assert type(rna_length) == int
        assert type(promoter_affiliation) == list
        if polycistronic:
        	assert type(polycistronic) == int
        Species.__init__(self, name, count=False)
        self.tag_rna = tag_rna
        self.rna_length = rna_length
        self.promoter_affiliation = promoter_affiliation
        self.polycistronic = polycistronic
        self.count = count
        """
        if count:
            self.count = count
        else:
            self.count = 0
            """

class RNApol(Species):
    
    def __init__(self, name, rate_of_decay, count=False):
        """Initiate the RNApol class.

        Parameters
        ----------
        name : str
               The name of the species, i.e. DNA, mRNA, etc.

        rate_of_decay : float
                        The rate of decay for the species of interest.

        count : int, optional
                The quantity of the species, defualts to zero.

        Returns
        -------
        RNApol object

        See Also
        --------
        Species : This class inherits from the species class

        """
        assert type(rate_of_decay) == float
        Species.__init__(self, name, count=False)
        self.r_decay = rate_of_decay
        
class DNApol(Species):
    
    def __init__(self, name, rate_of_decay, count=False):
        """Initiate the RNApol class.

        Parameters
        ----------
        name : str
               The name of the species, i.e. DNA, mRNA, etc.

        rate_of_decay : float
                        The rate of decay for the species of interest.

        count : int, optional
                The quantity of the species, defualts to zero.

        Returns
        -------
        DNApol object

        See Also
        --------
        Species : This class inherits from the species class

        """
        assert type(rate_of_decay) == float
        Species.__init__(self, name, count=False)
        self.r_decay = rate_of_decay

class Protein(Species):
    
    def __init__(self, name, tag_protein, count=False):
        """Initiate the RNApol class.

        Parameters
        ----------
        name : str
               The name of the species, i.e. DNA, mRNA, etc.

        rate_of_decay : float
                        The rate of decay for the species of interest.

        count : int, optional
                The quantity of the species, defualts to zero.

        Returns
        -------
        Protein object

        See Also
        --------
        Species : This class inherits from the species class

        """
        assert type(tag_protein) == int
        Species.__init__(self, name, count=False)
        self.tag_protein = tag_protein

class Reaction(object):

    """ Reaction object will be a list containing 2 lists and 1 tuple. The first list contains the reactant species. If a reactant species appears
        more than once, it is included in the list. (the propensity function can take into account stoichiometry- more of that in 
        the propensity function). The second list contains the product species. If a product species occurs more than once it is included.
        The tuple contains the rate constants for the forward and backward reactions (reactants to products and products to reactants)


    The self.tau method calculates the time for that reaction to occur (absolute time) in accordance with Gillespie's Next Reaction Method.

    The self.get_propensity method calculates the propensity."""

    def __init__(self, reactants, products, ks):
        #assert type(reactants[0]) == Species
        #assert type(products[0]) == Species
        assert type(ks) == tuple
        self.reactants = reactants
        self.products = products
        self.ks = ks
        self.tau = 0 #time for reaction to occur
        self.tau_old = 0 #time that is recorded if self.tau not equals to 0 
        self.prop_old = 0
        self.get_propensity()
    

    
    @property
    def time(self):
        return self.t

    @time.setter
    def time(self,time):
        self.t = time #time is the tau for the reaction that has just been executed
   
    def get_propensity(self):

    	""" 
    	
    	The get_propensity function must be able to handle reactions of varying stoichiometry. This code takes that into
    	account. The central part of the code consists of identifying degeneracies (multiple inputs of one species) in the
    	reactant list of the reaction object (For more information of how the reaction object is structured, please consult the
    	Reaction class). A dictionary is created for each reactant species with keys as the species, value as a list
		with each element as the same species. The number of entries in the list corresponds to the stoichiometry of 
		that reactant species.
    	
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
                    value *= multiplicity_dict[i][j].count - j 

                extend = [j for j in range(len(multiplicity_dict[i])) if j]
                extend.append(len(multiplicity_dict[i]))
                factorial = 1
                for j in extend:
                    factorial *= j
                intermediate = value/float(factorial) #application of the combinatoric formulae

                penultimate.append(intermediate)

            else:
                penultimate.append(multiplicity_dict[i][len(multiplicity_dict[i]) - 1].count)

        """ 
        At the end of the for loop, the penultimate list will contain values pertaining to the adjusted stoichiometic count for  
        each species. The next step is to iterate through this list and multiply all values in the list with the rate constant. 

        """
       
        a = 1
        for i in penultimate:
            a *= i
        a *= self.ks[0]
        if a:
            self.prop_old = a
        self.prop = a

        """ 
        At the end of the procedure, the propensity value is stored in self.prop. If the value is non-zero, it is also 
        stored in self.prop_old (this is necessary to retrive a propensity value for a reaction that was previously disabled
        according to the Next Reaction Method)

        """


    def get_tau(self, time): #this method is to recalculate the tau for the reaction that has just been executed
        if self.prop == 0:
            self.tau_old = self.tau
            self.t = time
            self.tau = np.inf
        else:
             self.tau = ((-1)*(np.log(np.random.random())))/float(self.prop) + time
             self.tau_old = self.tau
        #self.tau = tau

    def get_det_tau(self, time): #this method is to recalculate the tau for reactions that are dependent on the just executed reaction
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
            self.tau_old = self.tau

    

class NextReactionMethod(object):
    
    def __init__(self, directory=False, num_elong=False, k_elong=False, k_translation_elong=False):
        self.k_elong = k_elong
        self.k_translation_elong = k_translation_elong
        """
        if num_elong:
            self.num_elong = num_elong
            self.k_elong = k_elong
            self.create_elong_reactions()
        if directory:
            self.read_from_file(directory)
            """

    def create_rec_list(self, rec_list):
        assert type(rec_list[0]) == Reaction
        self.reactions = rec_list


        
    def create_species_list(self):
       
        #Creacting a list of species
        
        arch_list = []

        for i in self.reactions:
            temp = i.reactants + i.products
            arch_list.append(temp)

        reference = arch_list[0]

        for i in arch_list:
            if arch_list.index(i) > 0:
                catch = i + reference
                reference = catch

        species_set = set(catch)

        self.species = species_set

        #return self.species
       

    def create_transcription_elongation_species(self):
        elong_species_dict = {}
        for i in self.species:
	        if type(i) == DNA:
	            if i.nut_sequence and i.termination_sequence:
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
														if m.polycistronic == 1:
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
    						if k.tag_protein == i.tag_rna:
    							for l in self.species:
    								if l.name == 'Ribosome':
    									translation_reactions.append([Y[j]],[k,l],(self.k_translation_elong,self.k_translation_elong))
    		translation_react_dict.update({i:translation_reactions})

    	Y = translation_react_dict.values()
    	for i in range(len(Y)):
    		hold = Y[0]
    		if i > 0:
    			hold += Y[i]
    	self.translation_elong_reactions = hold


    def update_reaction_list(self):
    	self.reactions += self.transcription_elong_reactions + self.translation_elong_reactions


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

    def stat_reaction_list(self, stat_list):
    	self.stat_reaction_list = stat_list

    def operator_stat_model_dependencies(self):
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
		
		stat1_config_dependency = {}
		for i in self.transcription_elong_species.keys():
			for j in config:
				for k in range(len(j)):
					dep_list = []
					for l in self.stat_reaction_list:
						if j[k] == 1:
							for m in l.reactants:
								if m.name.split('-')[0] == 'Open':
									if m.name.split('-')[1] == 'PRM':
										dep_list.append(l)
					stat1_config_dependency.update({i:dep_list})





    def operator_stat_model(self):
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
			state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(R**r_power)*(C**c_power)*(RNAP**rnap_power))
			state_energy_list.append(state_energy)

		for i in state_energy_list:
			prob = i/float(np.sum(state_energy_list))
			probability_configuration_list.append(prob)

		config_dependency_dict = {}
		for i in range(len(config)):
			dependency_config = []
			for j in range(len(config[i])):
				if type(config[i][j]) == str:
					for k in self.reactions:
						for l in k.reactants:
							if config[i][j] == l.name:
								dependency_config.append(k)
								break
			config_dependency_dict.update({i,dependency_config})


		Y =np.argmax(np.random.multinomial(1,probability_configuration_list))

		for i in config[Y]:
			if type(i) == str:
				for j in self.species:
					if j.name == i:
						j.count -= 1


    def NRM_execution(self, step):
             
        system_time = 0

        #First Step, create a species dictionary from the species list

        #keys will be species, values will be a list of counts at each iteration of the stochastic algorithm

        species_dict = {}
        for i in self.species:
            species_dict.update({i.name:[i.count]})

        species_dict.update({'time':[system_time]})

        tau_list = []
        for m in self.reactions:
            m.get_tau(system_time) 
            tau_list.append(m.tau)

        self.generate_dep_graph()
        

        for i in range(step):
            reaction_index = np.argmin(tau_list)
            
            dependency_list = self.dep_graph[self.reactions[reaction_index]]

            system_time = tau_list[reaction_index] #this will give us the time variable input necessary for get_det_tau calculation
            
            """
            prop_list = []
            for h in self.reactions:
                prop_list.append(h.prop)
            print species_dict,prop_list,tau_list,reaction_index
            """
            
            if self.reactions[reaction_index].reactants == self.reactions[reaction_index].products:
                for j in self.reactions[reaction_index].reactants:
                    new = j.count - 1
                    j.count = new
            
            else:
                for j in self.reactions[reaction_index].reactants:
                    new = j.count - 1
                    j.count = new
                for j in self.reactions[reaction_index].products:
                    new = j.count + 1
                    j.count = new



            for j in dependency_list:
                if j == self.reactions[reaction_index]:       
                    j.get_propensity()
                    
                    j.get_tau(system_time)
                        
                    tau_list[reaction_index] = j.tau
                
                else:
                    j.get_propensity()

                    j.get_det_tau(system_time)

                    tau_list[self.reactions.index(j)] = j.tau


            for j in self.species:
                species_dict[j.name].append(j.count)

            species_dict['time'].append(system_time)
                
        return species_dict


    def multiple_simulation(self,trajectories,step):

        simulation_dict = {}

        self.create_species_list()

        for i in range(trajectories):

            shadow_reaction = copy.deepcopy(self)

            Y = shadow_reaction.NRM_execution(step)

            simulation_dict.update({i:Y})

        self.stats = simulation_dict