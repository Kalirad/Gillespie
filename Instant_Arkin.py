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
        self.count = count
        self.tag_dna = tag_dna
        self.dna_length = dna_length
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
        self.count = count
        self.tag_rna = tag_rna
        self.rna_length = rna_length
        self.promoter_affiliation = promoter_affiliation
        self.pre_termination = pre_termination
        

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
        self.count = count
        self.tag_rna = tag_rna

class IsomerComplex(Species):

    def __init__(self, name, promoter_state, tag_complex, dna, count=False):

        """Specify object with parameters special for initital transcription complex.

        Parameters
        ----------

        name : str

        promoter_state : int (0 or 1)
                Determines whether the transcribing complex is a 'Closed' config (0) or an 'Open' config (1)

        tag_complex : int
                Relates a Closed-Complex corresponding to a particular promoter with its Open complex and DNA copy number. An int number

        DNA : str
                The name of the DNA object (Promoter) the initial transcribed complex is associated with

        count : int
            defaults to zero
         """
        Species.__init__(self, name, count=False)
        self.count = count
        self.promoter_state = promoter_state
        self.tag_complex = tag_complex
        self.dna = dna

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
   
    def get_propensity(self,const=False):

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
        if const:
            a *= self.ks[1]
        else:
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
            hit = 0
            for j in i.reactants:
                if type(j) == DNA:
                    hit += 1
            if hit == 0:
                temp = i.reactants + i.products
            else:
                temp = i.reactants
            arch_list.append(temp)

        reference = arch_list[0]

        for i in arch_list:
            if arch_list.index(i) > 0:
                catch = i + reference
                reference = catch

        species_set = list(set(reference))

        self.species = species_set

        # the next bit of code removes explicit DNA-binding reactions from the list of Gillespie's reactions

        store_react = []
        for i in self.temp_reactions:
            hit = 0 
            for j in i.reactants:
                if type(j) == DNA:
                    hit += 1
            if hit == 0:
                store_react.append(i)

        self.intermediate_reactions = store_react

    def create_isomerization_species(self):

        for i in self.species:
            if type(i) == DNA:
                copy_number = i.count
                for j in range(copy_number):
                    name_closed = 'Closed' + '-' + i.name + '-' + str(j+1)
                    name_open = 'Open' + '-' + i.name + '-' + str(j+1)
                    self.species.append(IsomerComplex(name_closed, 0, j+1, i.name, count=False))
                    self.species.append(IsomerComplex(name_open, 1, j+1, i.name, count=False))

    def create_isomerization_reaction(self,reaction_constants):
        """ 
        reaction_constants : dict
                A dictionary of keys as the First entry before the hyphen of the IsomerComplex object and values as a tuple of reaction constants
                must be provided """

        isomer_react = []
        for i in self.species:
            products = []
            if type(i) == IsomerComplex:
                if i.promoter_state == 0:
                    for j in self.species:
                        if type(j) == IsomerComplex:
                            if j.promoter_state == 1:
                                if j.dna == i.dna:
                                    if j.tag_complex == i.tag_complex:
                                        products.append(j)
                    values = reaction_constants[i.name.split('-')[1]] #values will be a tuple of reactions constants
                    isomer_react.append(Reaction([i],products,values))

        self.reactions = self.intermediate_reactions + isomer_react


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
                inter_list = []
                term_list = []
                for j in range(i.nut_sequence[0],i.termination_sequence[1]):
                    if j == i.nut_sequence[0]:
                        name = str(i.name) + '-' + 'Nut' + '-' + str(j)
                        nut_list.append(Species(name,count=0))
                    elif i.nut_sequence[0] < j < i.nut_sequence[1]+1:
                        name = str(i.name) + '-' + 'N' + '-' + 'Nut' + '-' + str(j)
                        nut_list.append(Species(name,count=0))
                        if j == i.nut_sequence[1]:
                            inter_list.append(nut_list[-1])
                    elif i.nut_sequence[1]+1 <= j < i.termination_sequence[0]:
                        name = str(i.name) + '-' + 'N' + '-' + str(j)
                        inter_list.append(Species(name,count=0))
                self.transcription_elong_species[i].append(nut_list)
                self.transcription_elong_species[i].append(inter_list)

                if len(inter_list) != 0:
                    term_list.append(inter_list[-1])
                    for j in range(i.termination_sequence[0],i.dna_length):
                        if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                            name = str(i.name) + '-' + 'N' '-' + 'TR' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                        else:
                            name = str(i.name) + '-' + 'N' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                    self.transcription_elong_species[i].append(term_list)

                else:
                    term_list.append(nut_list[-1])
                    for j in range(i.termination_sequence[0],i.dna_length):
                        if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                            name = str(i.name) + '-' + 'N' '-' + 'TR' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                        else:
                            name = str(i.name) + '-' + 'N' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                    self.transcription_elong_species[i].append(term_list)

    def create_transcription_elongation_reactions(self):
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
                            if type(k) == IsomerComplex:
                                if k.promoter_state == 1:
                                    if k.dna == i.name:
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
                                                if type(m) == Protein:
                                                    if m.name == 'N':
                                                        elong_react.append(Reaction([Y[k],m],[Z[l+1]],(i.nut_sequence[2][1],i.nut_sequence[2][1])))
                                                        elong_react.append(Reaction([Z[l+1]],[Y[k],m],(i.nut_sequence[2][2],i.nut_sequence[2][2])))
                                elong_react.append(Reaction([Y[k]],[Y[k+1]],(i.nut_sequence[2][0],i.nut_sequence[2][0])))

                            elif i.nut_sequence[1] < k < i.termination_sequence[0]:
                                if i.nut_sequence[3] == k+1:
                                    for l in self.species:
                                        if type(l) == RNA:
                                            for n in l.promoter_affiliation:
                                                if n == i.tag_dna:
                                                    if l.pre_termination == True:
                                                        elong_react.append(Reaction([Y[k]],[Y[k+1],l],(self.k_elong,self.k_elong)))
                                else:
                                    elong_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
                            elif i.termination_sequence[0] <= k <= i.termination_sequence[1]:
                                for l in self.species:
                                    if l.name == 'RNAP':
                                        elong_react.append(Reaction([Y[k]],[Y[k+1]],(i.termination_sequence[2][0],i.termination_sequence[2][0])))
                                        elong_react.append(Reaction([Y[k]],[l],(i.termination_sequence[2][2],i.termination_sequence[2][2])))                          
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
                                                if not l.pre_termination:
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
                        if len(Y) != 0:
                            elong_inter_react = []
                            for k in range(len(Y)):
                                if 0 <= k < (len(Y) - 1):
                                    if int(Y[k+1].name.split('-')[2]) == i.nut_sequence[3]:
                                        for l in self.species:
                                            if type(l) == RNA:
                                                for m in l.promoter_affiliation:
                                                    if m == i.tag_dna:
                                                        if l.pre_termination:
                                                            elong_inter_react.append(Reaction([Y[k]],[Y[k+1],l],(self.k_elong,self.k_elong)))
                                    else:
                                        elong_inter_react.append(Reaction([Y[k]],[Y[k+1],],(self.k_elong,self.k_elong)))
                            elong_reactions_dict[i].append(elong_inter_react)
                    elif j == 3:
                        elong_term_react = []
                        for k in range(len(Y)):
                            if k < ((i.termination_sequence[1] - i.termination_sequence[0]) + 1):
                                if int(Y[k+1].name.split('-')[3]) == i.nut_sequence[3]:
                                        for l in self.species:
                                            if type(l) == RNA:
                                                for m in l.promoter_affiliation:
                                                    if m == i.tag_dna:
                                                        if l.pre_termination:
                                                            elong_term_react.append(Reaction([Y[k]],[Y[k+1],l],(self.k_elong,self.k_elong)))
                                else:
                                    elong_term_react.append(Reaction([Y[k]],[Y[k+1]],(i.termination_sequence[2][1],i.termination_sequence[2][1])))
                            elif ((i.termination_sequence[1] - i.termination_sequence[0]) + 1) <= k < (len(Y) - 1):
                                elong_term_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
                            elif k == (len(Y) - 1):
                                products = []
                                for l in self.species:
                                    if l.name == 'RNAP':
                                        products.append(l)
                                    elif type(l) == Protein:
                                        if l.name == 'N':
                                            products.append(l)
                                    elif type(l) == RNA:
                                        for m in l.promoter_affiliation:
                                            if m == i.tag_dna:
                                                if not l.pre_termination:
                                                    products.append(l)
                                elong_term_react.append(Reaction([Y[k]],products,(self.k_elong,self.k_elong)))
                        elong_reactions_dict[i].append(elong_term_react)
            else:
                Y = self.transcription_elong_species[i]
                elong_reactions = []
                for j in range(len(Y)):
                    if j == 0:
                        for k in self.species:
                            if type(k) == IsomerComplex:
                                if k.promoter_state == 1:
                                    if k.dna == i.name:
                                        for l in self.species:
                                            if type(l) == DNA:
                                                if l.tag_dna == i.tag_dna:
                                                    elong_reactions.append(Reaction([k],[l,Y[j]],(self.k_elong,self.k_elong)))
                        elong_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_elong,self.k_elong)))
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
            for j in range(len(Y)):
                if j == 0:
                    temp = []
                    valor = []
                    for k in self.species:
                        if type(k) == Species:
                            if k.name.split('-')[0] == i.name:
                                if k.name.split('-')[1] == 'Ribosome':
                                    temp.append(k)
                        elif type(k) == RNA:
                            if k.tag_rna == i.tag_rna:
                                valor.append(k)
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
        for i in self.total_reactions:
            temp = i.reactants + i.products
            dep_rec = []
            for j in self.total_reactions:
                set_value = [val for val in temp if val in j.reactants]
                if len(set_value) != 0:
                    dep_rec.append(j)
            self.dep_graph[i] = dep_rec
        #return self.dep_graph

    
    def RNAP_data_store(self):
        RNAP_elong_exit = []
        for i in range(len(self.total_reactions)):
            R = self.total_reactions[i]
            for j in R.reactants:
                if type(j) != IsomerComplex:
                    for k in R.products:
                        if type(k) == Species:
                            if k.name == 'RNAP':
                                RNAP_elong_exit.append(i)
                else:
                    if j.promoter_state == 1:
                        RNAP_elong_entry.append(i)
        self.RNAP_elong_store = {'RNAP_entry':RNAP_elong_entry, 'RNAP_exit':RNAP_elong_exit}

    def config_prom(self):

        config_dict = {}

        config1 = np.array([[1, 0, 0, 0, 0.0], [2, 0, 0, 'R', -11.7], [3, 0, 'R', 0, -10.1], [4, 'R', 0, 0, -10.1], [5, 0, 0, 'C', -10.8],
        [6, 0, 'C', 0, -10.8], [7, 'C', 0, 0, -12.1], [8,'RNAP', 0, 0, -11.5], [9, 0, 0, 'RNAP', -12.5], 
        [10, 0, 'R', 'R', -23.7], [11, 'R',  0, 'R', -21.8], [12, 'R', 'R', 0, -22.2], [13, 0, 'C', 'C', -21.6],
        [14, 'C', 0, 'C', -22.9], [15, 'C', 'C', 0, -22.9], [16, 'RNAP', 0, 'RNAP', -24.0], [17, 0, 'C', 'R', -22.5],
        [18, 0, 'R', 'C', -20.9], [19, 'R', 0, 'C', -20.9], [20, 'C', 0, 'R', -23.8], [21, 'R', 'C', 0, -20.9], 
        [22, 'C', 'R', 0, -22.2], [23, 'R', 0, 'RNAP', -22.6], [24, 'RNAP', 'R', 0, -21.6], [25, 'RNAP', 0, 'R', -23.2],
        [26, 'C', 0, 'RNAP', -24.6], [27, 'RNAP', 'C', 0, -22.3], [28, 'RNAP', 0, 'C', -22.3], [29, 'R', 'R', 'R', -33.8], 
        [30, 'C', 'C', 'C', -33.7], [ 31, 'C', 'R', 'R', -35.8], [32, 'R', 'C', 'R', -32.6], [33, 'R', 'R', 'C', -33.0], 
        [34, 'R', 'C', 'C', -31.7], [35, 'C', 'R', 'C', -33.0], [36, 'C', 'C', 'R', -34.6], [37, 'RNAP', 'R', 'R', -35.2], 
        [38, 'RNAP', 'C', 'C', -33.1], [39, 'RNAP', 'C', 'R', -34.0], [40, 'RNAP', 'R', 'C', -32.4]])

        config2 =  np.array([[0,0,0.0],[0,'RNAP',-9.9],['CII',0,-9.7],['CII','RNAP',-21.5]])

        config3 = np.array([[0,0,0.0],['C',0,-10.9],[0,'C',-12.1],['R',0,-11.7],[0,'R',-10.1],[0,'RNAP',-12.5],['C','C',-22.9],['C','R',-20.9],['R','C',-22.8],['R','R',-23.7]])

        config_dict.update({'PR_PRM':config1, 'PRE':config2, 'PL':config3})

        self.config_dict = config_dict




    def PR_PRM_model_dependencies1(self):
        #Use the positions of RNAP on the list(either position 1 or 3 - zero based indexing - to determine which reactions are affected
        
        for i in self.species:
            if type(i) == DNA:
                if i.name == 'PRM' or i.name == 'PR':
                    copy_number = []
                    for j in range(i.count):
                        copy_number.append(j+1)

                    break


        config = np.array([[1, 0, 0, 0, 0.0], [2, 0, 0, 'R', -11.7], [3, 0, 'R', 0, -10.1], [4, 'R', 0, 0, -10.1], [5, 0, 0, 'C', -10.8],
        [6, 0, 'C', 0, -10.8], [7, 'C', 0, 0, -12.1], [8,'RNAP', 0, 0, -11.5], [9, 0, 0, 'RNAP', -12.5], 
        [10, 0, 'R', 'R', -23.7], [11, 'R',  0, 'R', -21.8], [12, 'R', 'R', 0, -22.2], [13, 0, 'C', 'C', -21.6],
        [14, 'C', 0, 'C', -22.9], [15, 'C', 'C', 0, -22.9], [16, 'RNAP', 0, 'RNAP', -24.0], [17, 0, 'C', 'R', -22.5],
        [18, 0, 'R', 'C', -20.9], [19, 'R', 0, 'C', -20.9], [20, 'C', 0, 'R', -23.8], [21, 'R', 'C', 0, -20.9], 
        [22, 'C', 'R', 0, -22.2], [23, 'R', 0, 'RNAP', -22.6], [24, 'RNAP', 'R', 0, -21.6], [25, 'RNAP', 0, 'R', -23.2],
        [26, 'C', 0, 'RNAP', -24.6], [27, 'RNAP', 'C', 0, -22.3], [28, 'RNAP', 0, 'C', -22.3], [29, 'R', 'R', 'R', -33.8], 
        [30, 'C', 'C', 'C', -33.7], [ 31, 'C', 'R', 'R', -35.8], [32, 'R', 'C', 'R', -32.6], [33, 'R', 'R', 'C', -33.0], 
        [34, 'R', 'C', 'C', -31.7], [35, 'C', 'R', 'C', -33.0], [36, 'C', 'C', 'R', -34.6], [37, 'RNAP', 'R', 'R', -35.2], 
        [38, 'RNAP', 'C', 'C', -33.1], [39, 'RNAP', 'C', 'R', -34.0], [40, 'RNAP', 'R', 'C', -32.4]])
        
        occupancy_species_dict = {}

        PR_PRM_isomer_const = []

        PR_PRM_transcription_states = []

        for i in copy_number:
            single_dict = {}
            for j in range(len(config)):
                counter = 0
                occupancy_increase = []
                occupancy_decrease = []
                for k in range(len(config[j])):
                    if config[j][k] == 'RNAP':
                        counter += 1
                        if k == 1:
                            for l in self.species:
                                if type(l) == DNA:
                                    if l.name == 'PRM':
                                        occupancy_decrease.append(l)
                                elif type(l) == Species:
                                    if l.name == 'RNAP':
                                        occupancy_decrease.append(l)
                                elif type(l) == IsomerComplex:
                                    if l.promoter_state == 0:
                                        if l.dna == 'PRM':
                                            if l.tag_complex == i:
                                                occupancy_increase.append(l)
                            if i == 1:
                                for l in range(len(config[j])): #The stipulation of stimulated isomerization rate constant should not be included in this section. Consider Removing !
                                    if l == 2:
                                        if config[j][l] == 'R':
                                            PR_PRM_isomer_const.append(j)
                        elif k == 3:
                            for l in self.species:
                                if type(l) == DNA:
                                    if l.name == 'PR':
                                        occupancy_decrease.append(l)
                                elif type(l) == Species:
                                    if l.name == 'RNAP':
                                        occupancy_decrease.append(l)
                                elif type(l) == IsomerComplex:
                                    if l.promoter_state == 0:
                                        if l.dna == 'PR':
                                            if l.tag_complex == i:
                                                occupancy_increase.append(l)

                if len(occupancy_increase) > 0: #This allows us to distinguish between multiple promoter activation, promoter activation, and no promoter activation
                    hit = 0
                    for k in config[j]:
                        if k == 'RNAP':
                            hit += 1
                    if hit > 1:
                        single_dict.update({j:(occupancy_increase,list(set(occupancy_decrease)))})
                    else:
                        single_dict.update({j:(occupancy_increase,occupancy_decrease)})
                else:
                    single_dict.update({j:occupancy_decrease})

                if counter != 0:
                    V = [val for val in PR_PRM_transcription_states if val == j]
                    if len(V) == 0:
                        PR_PRM_transcription_states.append(j)
            occupancy_species_dict.update({i:single_dict})
        self.occupancy_species1 = occupancy_species_dict
        self.PR_PRM_isomer_const = PR_PRM_isomer_const 


        """
        Creates a dictionary which relates a particular occupancy state with the Gillespie reactions (self.reactions) that are affected. 
        There are no Gillespie reactions explicitly associated with RNApolymerase (RNAP) or DNA objects(promoters).
        These are 'implicit' reactions as they are alluded to in the stat_thermodynamic model for operator/promoter binding.
        """
        X = copy.deepcopy(self.occupancy_species1)

        occupancy_reaction_dict = {}

        for i in X.keys():
            temp_dict = {}
            for j in X[i].keys():
                reaction_list = []
                if type(X[i][j]) == tuple:
                    reaction_occupancy_list = list(set(X[i][j][0]))
                    for k in reaction_occupancy_list:
                        for l in self.reactions:
                            for m in l.reactants:
                                if k.name == m.name:
                                    reaction_list.append(l)
                    temp_dict.update({j:reaction_list})
            occupancy_reaction_dict.update({i:temp_dict})
        self.occupancy_reaction1 = occupancy_reaction_dict
        self.PR_PRM_transcription_states = PR_PRM_transcription_states

    def PR_PRM_stat_energy_model_selection1(self,copy_number):
        config = np.array([[1, 0, 0, 0, 0.0], [2, 0, 0, 'R', -11.7], [3, 0, 'R', 0, -10.1], [4, 'R', 0, 0, -10.1], [5, 0, 0, 'C', -10.8],
        [6, 0, 'C', 0, -10.8], [7, 'C', 0, 0, -12.1], [8,'RNAP', 0, 0, -11.5], [9, 0, 0, 'RNAP', -12.5], 
        [10, 0, 'R', 'R', -23.7], [11, 'R',  0, 'R', -21.8], [12, 'R', 'R', 0, -22.2], [13, 0, 'C', 'C', -21.6],
        [14, 'C', 0, 'C', -22.9], [15, 'C', 'C', 0, -22.9], [16, 'RNAP', 0, 'RNAP', -24.0], [17, 0, 'C', 'R', -22.5],
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
                                    if j.count < C_power:
                                        state_energy = 0
                                    elif k.count < R_power:
                                        state_energy = 0
                                    elif l.count < RNAP_power:
                                        state_energy = 0
                                    else:
                                        state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(k.molar_conc**R_power)*(j.molar_conc**C_power)*(l.molar_conc**RNAP_power))
            state_energy_list.append(state_energy)

        for i in state_energy_list:
            prob = i/float(np.sum(state_energy_list))
            probability_configuration_list.append(prob)


        Y = np.random.multinomial(1,probability_configuration_list)

        current_PR_PRM_config1 = np.argmax(Y)

        return current_PR_PRM_config1


    def PR_PRM_model_config_update1(self, copy_number, val, system_time, tau_list):

        #Remeber to change value of Z to val input
            
        Z = val

        X = self.occupancy_species1[copy_number]

        Y = self.occupancy_reaction1[copy_number]

        if type(X[Z]) == tuple:
            for i in X[Z][0]:
                i.count += 1
            for i in X[Z][1]:
                if type(i) == Species:
                    if i.name == 'RNAP':
                        i.count -= len(X[Z][0])
                    else:
                        i.count -= 1
                elif type(i) == DNA:
                    if i.name == 'PR' or i.name == 'PRM':
                        i.count -= 1

            for i in Y[Z]:
                hit = 0
                for j in self.PR_PRM_isomer_const:
                    if j == Z:
                        for k in i.reactants:
                            if type(k) == IsomerComplex:
                                if k.promoter_state == 0:
                                    hit += 1

                if hit == 0:
                    i.get_propensity()
                else:
                    i.get_propensity(const=True)

                for j in range(len(self.reactions)):
                    if self.reactions[j] == i:
                        i.get_det_tau(system_time)
                        tau_list[j] = i.tau

    def PRE_model_dependencies2(self):
        "This is the statistical binding model for the PRE promoter"
        
        for i in self.species:
            if type(i) == DNA:
                if i.name == 'PRE':
                    copy_number = []
                    for j in range(i.count):
                        copy_number.append(j+1)
                    break

        config2 = np.array([[0,0,0.0],[0,'RNAP',-9.9],['CII',0,-9.7],['CII','RNAP',-21.5]])

        occupancy_species_dict2 = {}

        PRE_isomer_const = []

        PRE_transcription_states = []

        for i in copy_number:
            temp_dict = {}
            for j in range(len(config2)):
                counter = 0
                occupancy_decrease = []
                occupancy_increase = []
                for k in range(len(config2[j])):
                    if config2[j][k] == 'RNAP':
                        counter += 1
                        for l in self.species:
                            if type(l) == DNA:
                                if l.name == 'PRE':
                                    occupancy_decrease.append(l)
                            elif type(l) == Species:
                                if l.name == 'RNAP':
                                    occupancy_decrease.append(l)
                            elif type(l) == IsomerComplex:
                                if l.promoter_state == 0:
                                    if l.dna == 'PRE':
                                        if l.tag_complex == i:
                                            occupancy_increase.append(l)
                        if i == 1:
                            for l in config2[j]:
                                if l == 'CII':
                                    PRE_isomer_const.append(j)
                if len(occupancy_increase) > 0:
                    temp_dict.update({j:(occupancy_increase,occupancy_decrease)})
                else:
                    temp_dict.update({j:occupancy_decrease})
                if counter != 0:
                    V = [val for val in PRE_transcription_states if val == j]
                    if len(V) == 0:
                        PRE_transcription_states.append(j)
            occupancy_species_dict2.update({i:temp_dict})
        self.occupancy_species2 = occupancy_species_dict2
        self.PRE_isomer_const = PRE_isomer_const

        X = copy.deepcopy(self.occupancy_species2)

        reaction_dict = {}
        for i in X.keys():
            temp_react = {}
            for j in X[i].keys():
                reaction_occupancy = []
                if type(X[i][j]) == tuple:
                    species = list(set(X[i][j][0]))
                    for k in species:
                        for l in self.reactions:
                            for m in l.reactants:
                                if k.name == m.name:
                                    reaction_occupancy.append(l)
                    temp_react.update({j:reaction_occupancy})
            reaction_dict.update({i:temp_react})
        self.occupancy_reaction2 = reaction_dict

        self.PRE_transcription_states = PRE_transcription_states

    def PRE_stat_energy_model_selection2(self, copy_number):
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
                if type(j) == Protein:
                    if j.name == 'CII':
                        for k in self.species:
                            if k.name == 'RNAP':
                                if j.count < CII_power:
                                    state_energy = 0
                                elif k.count < RNAP_power:
                                    state_energy = 0
                                else:
                                    state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(j.molar_conc**CII_power)*(k.molar_conc**RNAP_power))
            state_energy_list.append(state_energy)

        for i in state_energy_list:
            prob = i/float(np.sum(state_energy_list))
            probability_configuration_list.append(prob)


        Y = np.random.multinomial(1,probability_configuration_list)
        
        self.current_PRE_config2 = np.argmax(Y)

        return self.current_PRE_config2

    def PRE_model_config_update2(self ,copy_number, val, system_time, tau_list):
        Z = val

        X = self.occupancy_species2[copy_number]

        Y = self.occupancy_reaction2[copy_number]

        if type(X[Z]) == tuple:
            for j in X[Z][0]:
                j.count += 1
            for j in X[Z][1]:
                j.count -= 1

            for i in Y[Z]:
                hit = 0
                for j in self.PRE_isomer_const:
                    if j == Z:
                        for k in i.reactants:
                            if type(k) == IsomerComplex:
                                if k.promoter_state == 0:
                                    hit += 1
                if hit == 0:
                    i.get_propensity()
                else:
                    i.get_propensity(const=True)

                for j in range(len(self.reactions)):
                    if self.reactions[j] == i:
                        i.get_det_tau(system_time)
                        tau_list[j] = i.tau

    def PL_model_dependencies3(self):

        for i in self.species:
            if type(i) == DNA:
                if i.name == 'PL':
                    copy_number = []
                    for j in range(i.count):
                        copy_number.append(j+1)
                    break

        config3 = np.array([[0,0,0.0],['C',0,-10.9],[0,'C',-12.1],['R',0,-11.7],[0,'R',-10.1],[0,'RNAP',-12.5],['C','C',-22.9],['C','R',-20.9],['R','C',-22.8],['R','R',-23.7]])

        PL_transcription_states = []
        occupancy_species3_dict = {}
        for i in copy_number:
            temp_dict = {}
            for j in range(len(config3)):
                counter = 0
                occupancy_increase = []
                occupancy_decrease = []
                for k in range(len(config3[j])):
                    if config3[j][k] == 'RNAP':
                        counter += 1
                        for l in self.species:
                            if type(l) == DNA:
                                if l.name == 'PL':
                                    occupancy_decrease.append(l)
                            elif type(l) == Species:
                                if l.name == 'RNAP':
                                    occupancy_decrease.append(l)
                            elif type(l) == IsomerComplex:
                                if l.promoter_state == 0:
                                    if l.dna == 'PL':
                                        if l.tag_complex == i:
                                            occupancy_increase.append(l)
                if len(occupancy_increase) > 0:
                    temp_dict.update({j:(occupancy_increase,occupancy_decrease)})
                else:
                    temp_dict.update({j:occupancy_decrease})
                if counter != 0:
                    V = [val for val in PL_transcription_states if val == j]
                    if len(V) == 0:
                        PL_transcription_states.append(j)
            occupancy_species3_dict.update({i:temp_dict})
        self.occupancy_species3 = occupancy_species3_dict

        X = copy.deepcopy(self.occupancy_species3)

        reaction_occupancy_dict = {}

        for i in X.keys():
            reaction_dict = {}
            for j in X[i].keys():
                reaction_list = []
                if type(X[i][j]) == tuple:
                    species = list(set(X[i][j][0]))
                    for k in species:
                        for l in self.reactions:
                            for m in l.reactants:
                                if k.name == m.name:
                                    reaction_list.append(l)
                    reaction_dict.update({j:reaction_list})
            reaction_occupancy_dict.update({i:reaction_dict})
        self.occupancy_reaction3 = reaction_occupancy_dict
        self.PL_transcription_states = PL_transcription_states

    def PL_stat_energy_model_selection3(self, copy_number):
        config3 = np.array([[0,0,0.0],['C',0,-10.9],[0,'C',-12.1],['R',0,-11.7],[0,'R',-10.1],[0,'RNAP',-12.5],['C','C',-22.9],['C','R',-20.9],['R','C',-22.8],['R','R',-23.7]])

        boltzman = (1.9872041 * (10**(-3)))
        state_energy_list = []
        probability_configuration_list = []

        for i in config3:
            R_power = 0
            C_power = 0
            RNAP_power = 0
            for j in i:
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
                                    if j.count < R_power:
                                        state_energy = 0
                                    elif k.count < C_power:
                                        state_energy = 0
                                    elif l.count < RNAP_power:
                                        state_energy = 0
                                    else:
                                        state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(k.molar_conc**C_power)*(j.molar_conc**R_power)*(l.molar_conc**RNAP_power))
            state_energy_list.append(state_energy)

        for i in state_energy_list:
            prob = i/float(np.sum(state_energy_list))
            probability_configuration_list.append(prob)


        Y = np.random.multinomial(1,probability_configuration_list)

        current_PL_config3 =  np.argmax(Y)

        return current_PL_config3

    def PL_model_config_update3(self, copy_number, val, system_time, tau_list):
 
        Z = val

        X = self.occupancy_species3[copy_number]

        Y = self.occupancy_reaction3[copy_number]

        if type(X[Z]) == tuple:
            for j in X[Z][0]:
                j.count += 1
            for j in X[Z][1]:
                j.count -= 1

            for i in Y[Z]:
                i.get_propensity()
                for j in range(len(self.reactions)):
                    if self.reactions[j] == i:
                        i.get_det_tau(system_time)
                        tau_list[j] = i.tau

    def promoter_maintenance(self):
        prom_maint = {}
        for i in self.species:
            if type(i) == DNA:
                temp = {}
                hit = 0
                V = 'unbound'
                if i.name == 'PRM' or i.name == 'PR':
                    if hit == 0:
                        count = 1
                        for j in range(i.count):
                            temp.update({count:V})
                            count += 1
                        hit += 1
                        prom_maint.update({'PR_PRM':temp})
                else:
                    count = 1
                    for j in range(i.count):
                        temp.update({count:V})
                        count += 1
                    prom_maint.update({i.name:temp})
        self.prom_maint = prom_maint


    def promoter_snapshot(self):

        prom_snapshot = {}
        for l in range(4):
            temp2 = {}

            for i in self.prom_maint.keys():
                temp = {}
                for j in self.config_dict.keys():
                    if i == j:
                        for k in self.prom_maint[i].keys():
                            prom_list = np.zeros(len(self.config_dict[j]))
                            temp.update({k:prom_list})
                temp2.update({i:temp})

            prom_snapshot.update({l:temp2})

        self.prom_snapshot = prom_snapshot


    def promoter_occupancy_stat(self, cell_volume,system_time,tau_list):

        bound_elem = {}
        for i in self.prom_maint.keys():
            L = []
            for j in self.prom_maint[i].keys():
                if self.prom_maint[i][j] == 'unbound':
                    L.append(j)
            np.random.shuffle(L)
            bound_elem.update({i:L})

        M = [self.PR_PRM_stat_energy_model_selection1,self.PRE_stat_energy_model_selection2,self.PL_stat_energy_model_selection3]
        N = [self.PR_PRM_model_config_update1,self.PRE_model_config_update2,self.PL_model_config_update3]
        W = ['PR_PRM','PRE','PL']

        V = [val for val in zip(W,M,N)]

        PR_PRM = bound_elem['PR_PRM']
        PRE = bound_elem['PRE']
        PL = bound_elem['PL']
        
        temp_PR_PRM = {}
        for i in PR_PRM:
            temp_PR_PRM.update({i:0})
        temp_PRE = {}
        for i in PRE:
            temp_PRE.update({i:0})
        temp_PL = {}
        for i in PL:
            temp_PL.update({i:0})
        prom_status = {'PR_PRM':temp_PR_PRM,'PRE':temp_PRE,'PL':temp_PL}


        while np.sum(bound_elem.values()) > 0:

            F = range(3)
            np.random.shuffle(F)
            for i in F:
                Y = bound_elem[V[i][0]]
                if len(Y) > 0:
                    temp = []
                    val = Y[0]  
                    X = V[i][1](val)
                    V[i][2](val,X,system_time,tau_list)
                    for j in self.species:
                        j.molar_conc = (j.count)*(1/float(6.02*1e23))*(1/float(cell_volume))
                    prom_status[V[i][0]][val] = X
                    for j in range(len(Y)):
                        if j != 0:
                            temp.append(Y[j])
                    bound_elem[V[i][0]] = temp


                    self.prom_snapshot[self.time_window][V[i][0]][val][X] += 1

        self.prom_status = prom_status


    def stat_occup_change(self, system_time, tau_list):

        F = ['PR_PRM','PRE','PL']

        G = [self.occupancy_species1, self.occupancy_species2, self.occupancy_species3]

        I = [self.occupancy_reaction1,self.occupancy_reaction2,self.occupancy_reaction3]

        M = [val for val in zip(F,G,I)]

        for j in M:
            for i in self.prom_maint.keys():
                if j[0] == i:
                    for k in self.prom_maint[i].keys():
                        if type(self.prom_maint[i][k]) == str:
                            Y = self.prom_status[i][k]
                            if type(j[1][k][Y]) == tuple:
                                for l in j[1][k][Y][0]:
                                    l.count -= 1
                                for l in j[1][k][Y][1]:
                                    if type(l) == Species:
                                        if l.name == 'RNAP':
                                            l.count += len(j[1][k][Y][0])
                                    elif type(l) == DNA:
                                        l.count += 1

                    for k in self.prom_maint[i].keys():
                        if type(self.prom_maint[i][k]) == str:
                            Y = self.prom_status[i][k]
                            if type(j[1][k][Y]) == tuple:
                                for l in j[2][k][Y]:
                                    l.get_propensity()
                                    for m in range(len(self.reactions)):
                                        if self.reactions[m] == l:
                                            l.get_det_tau(system_time)
                                            tau_list[m] = l.tau     


    def check_isomerization(self, reaction_index, system_time, tau_list):

        # For multiple RNAP binding situations (as in PR and PRM) it is important to ensure only one of the transcribing complex is isomerized
        # Therefore, the 

        ref = [('PR_PRM',self.occupancy_species1,self.occupancy_reaction1),('PRE',self.occupancy_species2,self.occupancy_reaction2),('PL',self.occupancy_species3,self.occupancy_reaction3)]
        
        X = self.total_reactions[reaction_index] 
        if type(X.reactants[0]) == IsomerComplex:
            if type(X.products[0]) == IsomerComplex:
                Y = (X.reactants[0].dna, X.reactants[0].tag_complex)
                V = [val for val in self.prom_maint.keys() if Y[0] in val.split('_')]
                self.prom_maint[V[0]][Y[1]] = self.prom_status[V[0]][Y[1]]
                state = self.prom_maint[V[0]][Y[1]]

                for i in ref:
                    if i[0] == V[0]:
                        W = i[1][Y[1]][state]
                        assert type(W) == tuple
                        if len(W[0]) > 1:
                            for j in W[0]:
                                if j.dna != Y[0]:
                                    j.count -= 1
                            for j in W[1]:
                                if type(j) == Species:
                                    if j.name == 'RNAP':
                                        j.count += (len(W[0]) - 1)
                                elif type(j) == DNA:
                                    if j.name != V[0]:
                                        j.count += 1

                        M = i[2][Y[1]][state]
                        for j in M:
                            j.get_propensity()
                            for k in range(len(self.reactions)):
                                if j == self.reactions[k]:
                                    j.get_det_tau(system_time)
                                    tau_list[k] = j.tau

    def check_elong_react(self, reaction_index): #This function detects is an inititation elongation reaction has taken place. Updates the status of promoter occupancy
        """This function should be done prior to stat_occup change methods.
        If elongation reaction has occured, it will update the promoter status to unbound, permitting future initiation of transcription"""
        R = self.total_reactions[reaction_index]
        V = [val for val in R.products if type(val) == DNA]

        if len(V) == 1:

            M = [val for val in self.prom_maint.keys() if V[0].name in val.split('_')]
            X = R.reactants[0].tag_complex
            self.prom_maint[M[0]][X] = 'unbound'
            self.RNAP_elong_check += 1


    def check_transcription_elongation_complete(self, reaction_index):
        R = self.total_reactions[reaction_index]
        V = [val for val in R.products if val.name == 'RNAP']
        if len(V) == 1:
            self.RNAP_elong_check -= 1

    def check_translation_elongation_initiation(self, reaction_index):
        R = self.total_reactions[reaction_index]
        V = [val for val in R.reactants if len(val.name.split('-')) == 2 and val.name.split('-')[1] == 'Ribosome']
        if len(V) == 1:
            self.Ribosome_elong_check += 1

    def check_translation_elongation_complete(self, reaction_index):
        R = self.total_reactions[reaction_index]
        V = [val for val in R.products if val.name == 'Ribosome']
        if len(V) == 1:
            self.Ribosome_elong_check -= 1


    def increment_initiation(self, cell_volume):

        """ This function updates the counts of the E coli host proteins 
            to remain at similar concentrations throughout cell growth """

        for i in self.species:
            if type(i) == Species:
                if i.name == 'RNAP':
                    i.count = int(np.round((30*1e-9)*(6.02*1e23)*cell_volume))
                elif i.name == 'P1':
                    i.count = int(np.round((35*1e-9)*(6.02*1e23)*cell_volume)) 
                elif i.name == 'P2':
                    i.count = int(np.round((140*1e-9)*(6.02*1e23)*cell_volume))
                elif i.name == 'Ribosome':
                    i.count = int(np.round((500*1e-9)*(6.02*1e23)*cell_volume))


    def increment_calc(self, cell_volume):
    
        for i in self.species:
            if type(i) == Species:
                if i.name == 'RNAP':
                    num = int(np.round((30*1e-9)*(6.02*1e23)*cell_volume))
                    count = 0
                    hit = 0
                    V_open = [val for val in self.species if val.name.split('-')[0] == 'Open' and val.count != 0]
                    V_closed = [val for val in self.species if val.name.split('-')[0] == 'Closed' and val.count != 0]
                    elong = self.RNAP_elong_check
                    for j in V_open:
                        count += j.count
                    for j in V_closed:
                        hit += j.count
                    assert num >= (i.count + count + hit + elong)
                    if num > (i.count + count + hit + elong):
                        value = num - (i.count + count + hit + elong)
                        i.count += value

                elif i.name == 'Ribosome':
                    num = int(np.round((500*1e-9)*(6.02*1e23)*cell_volume))
                    count = 0
                    V = [val for val in self.species if len(val.name.split('-')) == 2 and val.name.split('-')[1] == 'Ribosome']
                    elong = self.Ribosome_elong_check
                    for j in V:
                        count += j.count
                    assert num >= (i.count + count + elong)
                    if num > (i.count + count + elong):
                        value = num - (i.count + count + elong)
                        i.count += value
                
                elif i.name == 'P1':
                    num = int(np.round((35*1e-9)*(6.02*1e23)*cell_volume))
                    count = 0
                    V = [val for val in self.species if len(val.name.split('-')) > 1 and val.name.split('-')[0] == 'P1']
                    for j in V:
                        count += j.count
                    if num > (i.count + count):
                        assert num > (i.count + count)
                        value = num - (i.count + count)
                        i.count += value

                elif i.name == 'P2':
                    num = int(np.round((140*1e-9)*(6.02*1e23)*cell_volume))
                    count = 0
                    V = [val for val in self.species if len(val.name.split('-')) > 1 and val.name.split('-')[0] == 'P2']
                    for j in V:
                        count += j.count 
                    if num > (i.count + count):
                        assert num > (i.count + count)
                        value = num - (i.count + count)
                        i.count += value

    
    def plot_dict(self,cell_volume,system_time):
        self.species_dict = {}
        for i in self.species:
            i.molar()
            self.species_dict.update({i:[i.count]})
        self.species_dict.update({'time':[system_time]})
        self.species_dict.update({'cell_vol':[cell_volume]})

        self.plot_dict = {}
        for i in self.species:
            if type(i) == Protein:
                if i.name == 'Cro':
                    self.plot_dict.update({i.name:[i.count]})
                elif i.name == 'CI':
                    self.plot_dict.update({i.name:[i.count]})
                elif i.name == 'CII':
                    self.plot_dict.update({i.name:[i.count]})
                elif i.name == 'CIII':
                    self.plot_dict.update({i.name:[i.count]})
            elif type(i) == Species:
                if i.name == 'C':
                    self.plot_dict.update({i.name:[i.count]})
                elif i.name == 'R':
                    self.plot_dict.update({i.name:[i.count]})
        self.plot_dict.update({'time':[system_time]})
        self.plot_dict.update({'cell_vol':[cell_volume]})
        
        P1 = 0
        for i in self.P1:
            P1 += i.count
        P1_m = (P1/float(6.02e23)) * (1/float(cell_volume))
        P1_tot = 0
        for i in self.P1_tot:
            P1_tot += i.count
        P1_tot_m = (P1_tot/float(6.02e23)) * (1/float(cell_volume))
        P2 = 0
        for i in self.P2:
            P2 += i.count
        P2_m = (P2/float(6.02e23)) * (1/float(cell_volume))
        P2_tot = 0
        for i in self.P2_tot:
            P2_tot += i.count
        P2_tot_m = (P2_tot/float(6.02e23)) * (1/float(cell_volume))


        prot_rate = (((0.002)*(P1_m)/float(P1_tot_m)) + ((0.6)*(P2_m)/float(P2_tot_m)))/float(0.002+0.6)
        self.plot_dict.update({'prot_lyt':[prot_rate]})



    def proteolytic_species(self):

        M = [val for val in self.species if type(val) == Species and val.name.split('-')[0] == 'P1' or val.name.split('-')[0] == 'P2']

        P1_tot = []
        P2_tot = []
        for i in M:
            if i.name.split('-')[0] == 'P1':
                P1_tot.append(i)
            elif i.name.split('-')[0] == 'P2':
                P2_tot.append(i)

        P1 = []
        P2 = []
        for i in P1_tot:
            V = [val for val in i.name.split('-') if val == 'CIII']
            if len(V) == 0:
                P1.append(i)
        for i in P2_tot:
            V = [val for val in i.name.split('-') if val == 'CIII']
            if len(V) == 0:
                P2.append(i)
        self.P1 = P1
        self.P2 = P2
        self.P1_tot = P1_tot
        self.P2_tot = P2_tot



    def prom_status_time(self, system_time):
        if system_time >= self.check_time:
            if self.check_time < 1800:
                self.check_time += 600
            elif self.check_time == 1800:
                self.check_time += 300
            
            self.time_window += 1

    def NRM_initialization_protocol(self,reaction_constants):

        self.create_isomerization_species()

        self.create_isomerization_reaction(reaction_constants)
        
        self.create_transcription_elongation_species()

        self.create_transcription_elongation_reactions()

        self.create_translation_elongation_species()

        self.create_translation_elongation_reactions()

        self.update_reaction_list()

        self.config_prom()

        self.promoter_snapshot()

    def promoter_decoupling(self, reaction_index, tau_list, system_time):

        self.check_transcription_elongation_complete(reaction_index)

        self.check_translation_elongation_initiation(reaction_index)

        self.check_translation_elongation_complete(reaction_index)

        self.check_isomerization(reaction_index, system_time, tau_list)

        self.stat_occup_change(system_time,tau_list)
        
        self.check_elong_react(reaction_index)
           
    def NRM_execution(self, step, leap, end, reaction_constants):

        self.time_window = 0

        self.check_time = 600

        system_time = 0

        self.RNAP_elong_check = 0

        self.Ribosome_elong_check = 0

        self.promoter_maintenance()

        self.NRM_initialization_protocol(reaction_constants)

        k = leap

        cell_volume = 1e-15

        K0 = 1/float(end)

        #First Step, create a species dictionary from the species list

        #keys will be species, values will be a list of counts at each iteration of the stochastic algorithm

        self.increment_initiation(cell_volume)

        self.proteolytic_species()

        self.plot_dict(cell_volume,system_time)

        tau_list = []
        for m in self.total_reactions:
            m.get_propensity()
            m.get_tau(system_time)
            tau_list.append(m.tau)


        self.generate_dep_graph()

        self.PR_PRM_model_dependencies1()

        self.PRE_model_dependencies2()

        self.PL_model_dependencies3()

        for i in range(step):

            while system_time >= k:
                Z = self.species_dict['time']
                temp = {}
                for j in range(len(Z)):
                    if Z[j] <= k:
                        diff = k - Z[j]
                        temp.update({j:diff})
                diff_values = temp.values()
                min_diff = np.min(diff_values)
                for j in temp.keys():
                    if temp[j] == min_diff:
                        update_index = j

                for j in self.species_dict.keys():
                    X = self.species_dict[j]
                    if type(j) == Protein:
                        if j.name == 'Cro':
                            self.plot_dict[j.name].append(X[update_index])
                        elif j.name == 'CI':
                            self.plot_dict[j.name].append(X[update_index])
                        elif j.name == 'CII':
                            self.plot_dict[j.name].append(X[update_index])
                        elif j.name == 'CIII':
                            self.plot_dict[j.name].append(X[update_index])
                    elif type(j) == Species:
                        if j.name == 'C':
                            self.plot_dict[j.name].append(X[update_index])
                        elif j.name == 'R':
                            self.plot_dict[j.name].append(X[update_index])
                
                self.plot_dict['time'].append(k)
                self.plot_dict['cell_vol'].append(self.species_dict['cell_vol'][update_index])

                P1 = 0
                for i in self.P1:
                    P1 += i.count
                P1_m = (P1/float(6.02e23)) * (1/float(cell_volume))
                P1_tot = 0
                for i in self.P1_tot:
                    P1_tot += i.count
                P1_tot_m = (P1_tot/float(6.02e23)) * (1/float(cell_volume))
                P2 = 0
                for i in self.P2:
                    P2 += i.count
                P2_m = (P2/float(6.02e23)) * (1/float(cell_volume))
                P2_tot = 0
                for i in self.P2_tot:
                    P2_tot += i.count
                P2_tot_m = (P2_tot/float(6.02e23)) * (1/float(cell_volume))


                prot_rate = (((0.002)*(P1_m)/float(P1_tot_m)) + ((0.6)*(P2_m)/float(P2_tot_m)))/float(0.002+0.6)
                self.plot_dict['prot_lyt'].append(prot_rate)


                k += leap

                if k >= float(end):
                    break
                            
            if system_time >= float(end): #ensures that lists are at uniform lenghts based on final time point (end) and interval values (leap) 
                key = [val for val in self.plot_dict.keys() if len(self.plot_dict[val]) != (end/float(leap) + 1)]
                if len(key) != 0:
                    for j in key:
                        self.plot_dict[j].append(self.plot_dict[j][-1])

                self.plot_dict.update({'prom_activation':self.prom_snapshot})
                break

            self.prom_status_time(system_time)

            if i > 0:

                self.promoter_decoupling(reaction_index, tau_list, system_time)

                self.increment_calc(cell_volume)
                
                for j in self.species:
                    j.molar_conc = (j.count)*(1/float(6.02e23))*(1/float(cell_volume))
            

            self.promoter_occupancy_stat(cell_volume,system_time,tau_list)

            reaction_index = np.argmin(tau_list)

            dependency_list = self.dep_graph[self.total_reactions[reaction_index]]

            system_time = tau_list[reaction_index] #this will give us the time variable input necessary for get_det_tau calculation


            if self.total_reactions[reaction_index].reactants == self.total_reactions[reaction_index].products:
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

                    tau_list[self.total_reactions.index(j)] = j.tau

            cell_volume = (1 + (K0*system_time))*(1e-15)

            for j in self.species:
                j.molar_conc = (j.count)*(1/float(6.02e23))*(1/float(cell_volume))
                self.species_dict[j].append(j.count)
            self.species_dict['time'].append(system_time)
            self.species_dict['cell_vol'].append(cell_volume)





