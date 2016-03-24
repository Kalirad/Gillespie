"""
Lambda Phage Model
========

mRNA tracking during simulation.  The model uses the Stochastic Simulation Algorithm from Gibson & Bruck (2000) [1].

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
    def __init__(self, name, tag_dna, mol_number, dna_length=False, count=False, nut_sequence=False, termination_sequence=False):
         


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
        self.mol_number = mol_number
        self.dna_length = dna_length
        self.nut_sequence = nut_sequence
        self.termination_sequence = termination_sequence

class RNA(Species):
    def __init__(self, name, tag_rna, rna_length, promoter_affiliation, mol_numb, pre_termination=False, count=False):

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
        assert type(mol_numb) == int
            
        Species.__init__(self, name, count=False)
        self.count = count
        self.tag_rna = tag_rna
        self.rna_length = rna_length
        self.promoter_affiliation = promoter_affiliation
        self.mol_numb = mol_numb
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

    def __init__(self, name, promoter_state, dna, count=False):

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
        self.dna = dna

class Translate_Elong(Species):
    
    def __init__(self, name, rna, nucleotide, count=False):
        """Initiate Translation elongation species. This object
        will be a great boon when dealing with mRNA methods.

        Parameters

        ----------

        name : str

        rna : The name of the RNA object it is associated with 

        nucleotide : int, Corresponds to position on mRNA transcript

        count : int
        """
        Species.__init__(self, name, count=False)
        self.rna = rna
        self.nucleotide = nucleotide

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
   
    def get_propensity(self, const=False, transcript=False):

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
        if transcript:
            a *= (self.ks[0] * transcript)
        else:
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
    
    def __init__(self, copy_number=False, k_elong=False, k_translation_elong=False):
        """
        k_elong : float

        k_translation_elong : float

        """
        self.k_elong = k_elong
        self.k_translation_elong = k_translation_elong
        self.copy_number = copy_number

    def create_rec_list(self, rec_list): #creates the reaction list
        assert type(rec_list[0]) == Reaction
        self.temp_reactions = rec_list

    def create_species_list(self): #creates a list of species objects - this does not include transcription and translation reactions
       
        #Creacting a list of species
        react_spec = []
        prod_spec = []
        
        for i in self.temp_reactions:
            for j in i.reactants:
                react_spec.append(j)
            for j in i.products:
                prod_spec.append(j)
        X = list(set(react_spec))
        Y = list(set(prod_spec))

        spec_list = X + Y
        self.species = list(set(spec_list))

    def create_DNA_objects(self, DNA_init):
        for i in range(self.copy_number):
            for j in DNA_init.keys():
                name = j + '-' + 'CN' + '-' +str(i+1)
                self.species.append(DNA(name, DNA_init[j][0], i+1, dna_length=DNA_init[j][1], count=1, nut_sequence=DNA_init[j][2], termination_sequence=DNA_init[j][3]))
    
    def append_RNA_objects(self, RNA_list):
        for i in RNA_list:
            self.species.append(i)

    def create_Isomer_species(self):
        self.species.append(Species('RNAP',count=0))
        V = [val for val in self.species if type(val) == Species and val.name == 'RNAP']
        assert len(V) == 1
        for i in self.species:
            if type(i) == DNA:
                name1 =  'Closed' + '-' + i.name + '-' + str(i.mol_number)
                self.species.append(IsomerComplex(name1, 0, i.name, count=False))
                name2 = 'Open' + '-' + i.name + '-' + str(i.mol_number)
                self.species.append(IsomerComplex(name2, 1, i.name, count=False))

    def create_transcript_activation_species(self):
        activation = []
        V = [val for val in self.species if type(val) == DNA]
        for i in V:
            name = i.name.split('-')[0] + '-' + 'transcript_activate' + '-' + str(i.mol_number)
            activation.append(Species(name,count=False))
        self.transcript_activation = activation


    def create_transcribe_reactions(self):
        transcribe_reactions = []
        X = [val for val in self.species if type(val) == Species and val.name == 'RNAP']
        for i in self.species:
            if type(i) == DNA:
                V = [val for val in self.species if type(val) == IsomerComplex and val.dna == i.name]
                Q = [val for val in V if val.promoter_state == 0]
                assert len(Q) == 1
                for j in self.transcript_activation:
                    Y = j.name.split('-')
                    if Y[0] == i.name.split('-')[0]:
                        if int(Y[-1]) == i.mol_number:
                            transcribe_reactions.append(Reaction([X[0],i,j], [Q[0]], (1,1)))

        self.reactions = self.temp_reactions + transcribe_reactions

    def create_isomerization_reactions(self,reaction_constants):
        X = [val for val in self.species if type(val) == IsomerComplex and val.promoter_state == 0]
        Y = [val for val in self.species if type(val) == IsomerComplex and val.promoter_state == 1]
        isomer_reaction = []
        for i in X:
            for j in Y:
                if j.dna == i.dna:
                    isomer_reaction.append(Reaction([i],[j],reaction_constants[j.dna.split('-')[0]]))
        self.reactions += isomer_reaction

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
        DNA_specs = []
        for i in self.species:
            if type(i) == DNA:
                if i.mol_number == 1:
                    DNA_specs.append(i)

        assert len(DNA_specs) == 4 #There should be four promoters
        
        elong_species_dict = {} 
        for i in DNA_specs:
            ref = i.name.split('-')[0]
            if i.nut_sequence and i.termination_sequence: #not all DNA objects will have NUT sites and termination sites
                elong_list = []
                for j in range(i.dna_length):
                    if j < i.nut_sequence[0]:
                        name = str(ref) + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                    elif i.nut_sequence[0] <= j <= i.nut_sequence[1]:
                        name = str(ref) + '-' + 'Nut' + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                    elif i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                        name = str(ref) + '-' + 'TR' + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                    else:
                        name = str(ref) + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                elong_species_dict.update({i:[elong_list]})
            else: 
                elong_list = []
                for j in range(i.dna_length):
                    name = str(ref) + '-' + str(j)
                    elong_list.append(Species(name,count=0))
                elong_species_dict.update({i:elong_list})

        self.transcription_elong_species = elong_species_dict


        for i in self.transcription_elong_species.keys():
            ref = i.name.split('-')[0]
            if i.nut_sequence and i.termination_sequence:
                nut_list = []
                inter_list = []
                term_list = []
                for j in range(i.nut_sequence[0],i.termination_sequence[1]):
                    if j == i.nut_sequence[0]:
                        name = str(ref) + '-' + 'Nut' + '-' + str(j)
                        nut_list.append(Species(name,count=0))
                    elif i.nut_sequence[0] < j < i.nut_sequence[1]+1:
                        name = str(ref) + '-' + 'N' + '-' + 'Nut' + '-' + str(j)
                        nut_list.append(Species(name,count=0))
                        if j == i.nut_sequence[1]:
                            inter_list.append(nut_list[-1])
                    elif i.nut_sequence[1]+1 <= j < i.termination_sequence[0]:
                        name = str(ref) + '-' + 'N' + '-' + str(j)
                        inter_list.append(Species(name,count=0))
                self.transcription_elong_species[i].append(nut_list)
                self.transcription_elong_species[i].append(inter_list)

                if len(inter_list) != 0:
                    term_list.append(inter_list[-1])
                    for j in range(i.termination_sequence[0],i.dna_length):
                        if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                            name = str(ref) + '-' + 'N' '-' + 'TR' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                        else:
                            name = str(ref) + '-' + 'N' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                    self.transcription_elong_species[i].append(term_list)

                else:
                    term_list.append(nut_list[-1])
                    for j in range(i.termination_sequence[0],i.dna_length):
                        if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                            name = str(ref) + '-' + 'N' '-' + 'TR' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                        else:
                            name = str(ref) + '-' + 'N' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                    self.transcription_elong_species[i].append(term_list)

    def create_transcription_elongation_reactions(self):
        """
        The transcription elongation reactions are represented in a dictionary. The keys in the dictionary represent the DNA species objects
        and the values mirror the list of list for the species objects representing the different possible configurations of RNApolymerase
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
                            if type(k) == DNA:
                                if i.name.split('-')[0] == k.name.split('-')[0]:
                                    for l in self.species:
                                        if type(l) == IsomerComplex:
                                            if l.dna == k.name:
                                                if l.promoter_state == 1:
                                                    elong_react.append(Reaction([l],[k,Y[0]],(self.k_elong,self.k_elong)))
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
                            if type(k) == DNA:
                                if i.name.split('-')[0] == k.name.split('-')[0]: 
                                    for l in self.species:
                                        if type(l) == IsomerComplex:
                                            if l.dna == k.name:
                                                if l.promoter_state == 1:
                                                    elong_reactions.append(Reaction([l],[k,Y[j]],(self.k_elong,self.k_elong)))
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
                    translation_species_list.append(Translate_Elong(name, i.name, j, count=0))
                translation_species_dict.update({i:translation_species_list})
        self.translation_elong_species = translation_species_dict

    def create_translation_elongation_reactions(self):
        translation_react_dict = {}
        for i in self.translation_elong_species.keys():
            translation_reactions = []
            Y = self.translation_elong_species[i]
            for j in range(len(Y)):
                if j == 0:
                    translation_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_translation_elong,0)))
                elif 0 < j < (len(Y) - 1):
                    translation_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_translation_elong,0)))
                elif j == (len(Y) - 1):
                    for k in self.species:
                        if type(k) == Protein:
                            if k.tag_rna == i.tag_rna:
                                for l in self.species:
                                    if l.name == 'Ribosome':
                                        translation_reactions.append(Reaction([Y[j]],[k,l],(self.k_translation_elong,0)))
            translation_react_dict.update({i:translation_reactions})

        translation_elong_dict = {}
        for i in translation_react_dict.keys():
            translation_elong_dict.update({i.name:translation_react_dict[i]})

        M = []
        for i in translation_react_dict.keys():
            for j in range(len(translation_react_dict[i])):
                M.append(translation_react_dict[i][j])

        self.translation_react_dict = translation_react_dict
        self.translation_elong_reactions = M
        self.translation_elong_dict = translation_elong_dict
               
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

    def unique_reactions_promoter(self):
        unique_prom = {}
        for i in self.transcript_activation:
            for j in self.reactions:
                if len(j.products) == 1:
                    if type(j.products[0]) == IsomerComplex:
                        if j.products[0].promoter_state == 0:
                            Y = j.products[0]
                            if i.name.split('-')[-1] == Y.name.split('-')[-1]:
                                if Y.dna.split('-')[0] == i.name.split('-')[0]:
                                    unique_prom.update({i:j})
        self.unique_prom = unique_prom

    def isomerization_stim_activity(self):
        isomer_react_dict = {}
        V = [val for val in self.reactions if len(val.reactants) == len(val.products) and type(val.reactants[0]) == IsomerComplex]
        for i in V:
            isomer_react_dict.update({i:False})
        self.isomer_react_dict = isomer_react_dict
    
    def PR_PRM_model_dependencies1(self):
        #Use the positions of RNAP on the list(either position 1 or 3 - zero based indexing - to determine which reactions are affected
        copy_number = [(val+1) for val in range(self.copy_number)]
        
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
        
        PR_PRM_isomer_const = []

        PR_PRM_transcription_states = {}
        PRM_states = []
        for l in copy_number:
            temp = {}
            for i in range(len(config)):
                reaction_list = []
                for j in range(len(config[i])):
                    if config[i][j] == 'RNAP':
                        if j == 1:
                            for k in self.species:
                                if type(k) == DNA:
                                    if k.name.split('-')[0] == 'PRM':
                                        if k.mol_number == l:
                                            V = [val for val in self.reactions if type(val.products[0]) == IsomerComplex and len(val.reactants) != 1 and val.products[0].dna == k.name]
                                            assert len(V) == 1
                                            reaction_list.append(V[0])
                            if l == 1:
                                PRM_states.append(i)
                                if config[i][j+1] == 'R':
                                    PR_PRM_isomer_const.append(i)
                        elif j == 3:
                            for k in self.species:
                                if type(k) == DNA:
                                    if k.name.split('-')[0] == 'PR':
                                        if k.mol_number == l:
                                            V = [val for val in self.reactions if type(val.products[0]) == IsomerComplex and len(val.reactants) != 1 and val.products[0].dna == k.name]
                                            assert len(V) == 1
                                            reaction_list.append(V[0])
                if len(reaction_list) != 0:
                    temp.update({i:reaction_list})
            PR_PRM_transcription_states.update({l:temp})
        self.occupancy_reaction1 = PR_PRM_transcription_states
        self.PR_PRM_isomer_const = PR_PRM_isomer_const
        self.PRM_states = PRM_states

    def PR_PRM_stat_energy_model_selection1(self, copy):
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
                                    state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(k.molar_conc**R_power)*(j.molar_conc**C_power)*(l.molar_conc**RNAP_power))
            state_energy_list.append(state_energy)

        for i in state_energy_list:
            prob = i/float(np.sum(state_energy_list))
            probability_configuration_list.append(prob)


        Y = np.random.multinomial(1,probability_configuration_list)

        current_PR_PRM_config1 = np.argmax(Y)

        prob_value = probability_configuration_list[current_PR_PRM_config1]

        X = (current_PR_PRM_config1, prob_value)

        Q = [val for val in self.PR_PRM_isomer_const if val == current_PR_PRM_config1]

        if len(Q) == 1:
            V = [val for val in self.PRM_states if val == current_PR_PRM_config1]
            if len(V) == 1:
                for j in self.isomer_react_dict.keys():
                    if j.reactants[0].dna.split('-')[0] == 'PRM':
                        if int(j.reactants[0].name.split('-')[-1]) == copy:
                            self.isomer_react_dict.update({j:True})
        return X

    def PR_PRM_model_config_update1(self, copy, val, system_time, tau_list):

        #Remeber to change value of Z to val input

        V = [i for i in self.occupancy_reaction1[copy].keys() if i == val[0]]
    
        if len(V) == 1:
            for i in range(len(self.occupancy_reaction1[copy][val[0]])):
                R = self.occupancy_reaction1[copy][val[0]][i]
                for j in self.unique_prom.keys():
                    if self.unique_prom[j] == R:
                        j.count += 1
                R.get_propensity(transcript=val[1])
                for j in range(len(self.reactions)):
                    if self.reactions[j] == R:
                        R.get_det_tau(system_time)
                        tau_list[j] = R.tau
    
    def PRE_model_dependencies2(self):
        "This is the statistical binding model for the PRE promoter"
        copy_number = [(val+1) for val in range(self.copy_number)]

        config2 = np.array([[0,0,0.0],[0,'RNAP',-9.9],['CII',0,-9.7],['CII','RNAP',-21.5]])

        PRE_isomer_const = []
        PRE_transcription_states = {}
        for l in copy_number:
            temp = {}
            for i in range(len(config2)):
                reaction_list = []
                for j in range(len(config2[i])):
                    if config2[i][j] == 'RNAP':
                        for k in self.species:
                            if type(k) == DNA:
                                if k.name.split('-')[0] == 'PRE':
                                    if k.mol_number == l:
                                        V = [val for val in self.reactions if type(val.products[0]) == IsomerComplex and len(val.reactants) != 1 and val.products[0].dna == k.name]
                                        assert len(V) == 1
                                        reaction_list.append(V[0])  
                        if l == 1:     
                            if config2[i][j-1] == 'CII':
                                PRE_isomer_const.append(i)
                if len(reaction_list) != 0:
                    temp.update({i:reaction_list})
            PRE_transcription_states.update({l:temp})
        self.PRE_isomer_const = PRE_isomer_const
        self.occupancy_reaction2 = PRE_transcription_states
  
    def PRE_stat_energy_model_selection2(self, copy):
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
                                state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(j.molar_conc**CII_power)*(k.molar_conc**RNAP_power))
            state_energy_list.append(state_energy)

        for i in state_energy_list:
            prob = i/float(np.sum(state_energy_list))
            probability_configuration_list.append(prob)


        Y = np.random.multinomial(1,probability_configuration_list)
        
        current_PRE_config2 = np.argmax(Y)

        prob_value = probability_configuration_list[current_PRE_config2]

        X = (current_PRE_config2, prob_value)

        Q = [val for val in self.PRE_isomer_const if val == current_PRE_config2]
        if len(Q) == 1:
            for j in self.isomer_react_dict.keys():
                if j.reactants[0].dna.split('-')[0] == 'PRE':
                    if int(j.reactants[0].name.split('-')[-1]) == copy:
                        self.isomer_react_dict.update({j:True})
        return X 

    def PRE_model_config_update2(self, copy, val, system_time, tau_list):

        V = [i for i in self.occupancy_reaction2[copy].keys() if i == val[0]]

        if len(V) == 1:
           for i in range(len(self.occupancy_reaction2[copy][val[0]])):
                R = self.occupancy_reaction2[copy][val[0]][i]
                for j in self.unique_prom.keys():
                    if self.unique_prom[j] == R:
                        j.count += 1
                R.get_propensity(transcript=val[1])
                for j in range(len(self.reactions)):
                    if self.reactions[j] == R:
                        R.get_det_tau(system_time)
                        tau_list[j] = R.tau

    def PL_model_dependencies3(self):

        copy_number = [(val + 1) for val in range(self.copy_number)]

        config3 = np.array([[0,0,0.0],['C',0,-10.9],[0,'C',-12.1],['R',0,-11.7],[0,'R',-10.1],[0,'RNAP',-12.5],['C','C',-22.9],['C','R',-20.9],['R','C',-22.8],['R','R',-23.7]])

        PL_transcript_states = {}
        for l in copy_number:
            temp = {}
            for i in range(len(config3)):
                reaction_list = []
                for j in range(len(config3[i])):
                    if config3[i][j] == 'RNAP':
                        for k in self.species:
                            if type(k) == DNA:
                                if k.name.split('-')[0] == 'PL':
                                    if k.mol_number == l:
                                        V = [val for val in self.reactions if type(val.products[0]) == IsomerComplex and len(val.reactants) != 1 and val.products[0].dna == k.name]
                                        assert len(V) == 1
                                        reaction_list.append(V[0])
                if len(reaction_list) != 0:
                    temp.update({i:reaction_list})
            PL_transcript_states.update({l:temp})
        self.occupancy_reaction3 = PL_transcript_states

    def PL_stat_energy_model_selection3(self,copy):
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
                                    state_energy = (np.exp(((-1)*float(i[-1]))/float(boltzman*310.15))*(k.molar_conc**C_power)*(j.molar_conc**R_power)*(l.molar_conc**RNAP_power))
            state_energy_list.append(state_energy)

        for i in state_energy_list:
            prob = i/float(np.sum(state_energy_list))
            probability_configuration_list.append(prob)


        Y = np.random.multinomial(1,probability_configuration_list)

        current_PL_config3 =  np.argmax(Y)

        prob_value = probability_configuration_list[current_PL_config3] 

        X = (current_PL_config3, prob_value)

        return X

    def PL_model_config_update3(self, copy, val, system_time, tau_list):

        V = [i for i in self.occupancy_reaction3[copy].keys() if i == val[0]]
        if len(V) == 1:
           for i in range(len(self.occupancy_reaction3[copy][val[0]])):
                R = self.occupancy_reaction3[copy][val[0]][i]
                for j in self.unique_prom.keys():
                    if self.unique_prom[j] == R:
                        j.count += 1
                R.get_propensity(transcript=val[1])
                for j in range(len(self.reactions)):
                    if self.reactions[j] == R:
                        R.get_det_tau(system_time)
                        tau_list[j] = R.tau

    def promoter_occupancy_states(self, system_time, tau_list):
        switch_list = ['PR_PRM','PRE','PL']
        switch = {}
        for i in switch_list:
            L = []
            for j in range(self.copy_number):
                L.append(j+1)
            np.random.shuffle(L)
            switch.update({i:L})

        M = [self.PR_PRM_stat_energy_model_selection1,self.PRE_stat_energy_model_selection2,self.PL_stat_energy_model_selection3]
        N = [self.PR_PRM_model_config_update1,self.PRE_model_config_update2,self.PL_model_config_update3]
        V = [val for val in zip(switch_list,M,N)]

        while np.sum(switch.values()) > 0:
            F = range(3)
            np.random.shuffle(F)
            for i in F:
                Y = switch[V[i][0]]
                if len(Y) > 0:
                    temp = []
                    val = Y[0]
                    X = V[i][1](val)
                    V[i][2](val,X,system_time,tau_list)
                    for j in range(len(Y)):
                        if j != 0:
                            temp.append(Y[j])
                    switch[V[i][0]] = temp

    def reverse_initiation_constant(self, system_time, tau_list):
        V = [val for val in self.transcript_activation if val.count != 0]
        for i in V:
            i.count -= 1
            for j in self.unique_prom.keys():
                if j == i:
                    R = self.unique_prom[j]
                    R.get_propensity()
                    for k in range(len(self.reactions)):
                        if self.reactions[k] == R:
                            R.get_det_tau(system_time)
                            tau_list[k] = R.tau

    def reverse_stim_constant(self):
        V = [val for val in self.isomer_react_dict.keys() if val.prop == 0]
        for i in V:
            self.isomer_react_dict.update({i:False})

    def engage_isomer_reactions(self, system_time, tau_list, reaction_object):
        R = reaction_object
        rate_const = self.isomer_react_dict[R]
        if rate_const:
            R.get_propensity(const=True)
        else:
            R.get_propensity()
        for i in range(len(self.reactions)):
            if self.reactions[i] == R:
                R.get_det_tau(system_time)
                tau_list[i] = R.tau

    def check_transcription_elong_initiation(self, reaction_index): #This function detects if an inititation elongation reaction has taken place. Updates the status of promoter occupancy
        """This function should be done prior to stat_occup change methods.
        If elongation reaction has occured, it will update the promoter status to unbound, permitting future initiation of transcription"""
        R = self.total_reactions[reaction_index]
        V = [val for val in R.products if type(val) == DNA]

        if len(V) == 1:
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

    def find_first_translate_elongation(self):
    #Makes a dictionary with the RNA name as keys and the first elongation reaction object 0 -> 1 as value
        find_elongate = {}
        for i in self.translation_elong_dict.keys():
            V = [val for val in self.translation_elong_dict[i] if len([j for j in val.reactants if type(j) == Translate_Elong and j.nucleotide == 0]) == 1]
            assert len(V) == 1
            find_elongate.update({i:V[0]})
        self.find_first_elongate = find_elongate

    def find_last_translate_elongation(self):
        find_last_elongate = {}
        for i in self.translation_react_dict.keys():
            V = [val for val in self.translation_react_dict[i] if len([j for j in val.reactants if type(j) == Translate_Elong and j.nucleotide == (i.rna_length - 1)]) == 1]
            assert len(V) == 1
            find_last_elongate.update({i.name:V[0]})#remember, the translation_react_dict has an rna object as its key, not a string correspondong to i.name
        self.find_last_elongate = find_last_elongate

    def RNA_synthesis(self, V):
        for i in self.species:
            if type(i) == RNA:
                if i.name == V[0].name:
                    R = RNA(name = i.name, tag_rna = i.tag_rna, rna_length = i.rna_length, promoter_affiliation = i.promoter_affiliation, mol_numb = i.count, \
                        pre_termination = i.pre_termination, count = 1)
                    self.RNA_output[R.name].update({R.mol_numb:0})
        return R

    def RNA_reaction_init(self, rna_obj, RNA_reaction_init, system_time, tau_list):
        name1 = rna_obj.name + '-' +  str(rna_obj.mol_numb) + '-' + 'Ribosome'
        for i in self.species:
            if type(i) == Species and i.name == 'Ribosome':
                spec_obj = Species(name = name1, count=False)
                R1 = Reaction([rna_obj,i],[spec_obj],RNA_reaction_init[0])
        R2 = Reaction([rna_obj],[rna_obj],RNA_reaction_init[1])

        R3 = self.find_first_elongate[rna_obj.name]
        
        R4 = Reaction([spec_obj],[rna_obj,R3.reactants[0]],(self.k_translation_elong,self.k_translation_elong))
        X = [R1,R2,R4]
        for i in X:
            self.total_reactions.append(i)
            i.get_propensity()
            i.get_tau(system_time)
            tau_list.append(i.tau)
        Y = [R1,R2,R3,R4]

        return Y

    def update_dep_graph(self, reaction_list):
        #First way of updating dependency graph - Cunning way - reaction list is the Y output from previous method
        temp_graph = {}
        for i in reaction_list:
            hold = i.reactants + i.products
            dep_rec = []
            for j in reaction_list:
                set_value = [val for val in hold if val in j.reactants]
                if len(set_value) != 0:
                    dep_rec.append(j)
            temp_graph[i] = dep_rec

        for i in self.RNA_synth_react:
            temp_graph[reaction_list[0]].append(i)
            self.dep_graph[i].append(reaction_list[0])

        for i in temp_graph.keys():
            if i != reaction_list[2]: # do not reform the self.find_first_elongate reaction dependency entry
                self.dep_graph.update({i:temp_graph[i]})

        for i in self.find_last_elongate.keys():
            rxn_obj = self.find_last_elongate[i]
            self.dep_graph[rxn_obj].append(reaction_list[0])
        
        self.RNA_synth_react.append(reaction_list[0])
    
    def initialize_mRNA_times(self, rna_object, system_time):
        self.RNA_times[rna_object.name].update({rna_object.mol_numb:[system_time]})
        self.initiation_tracker[rna_object.name].update({rna_object.mol_numb:[]})

    def decay_mRNA_times(self, reaction_obj, system_time):
        if type(reaction_obj.reactants[0]) == RNA:
            X = reaction_obj.reactants[0]
            self.RNA_times[X.name][X.mol_numb].append(system_time)

    def RNA_translate_init_time(self, reaction_obj, system_time):
        V = [val for val in reaction_obj.reactants if (len(val.name.split('-')) > 1 and val.name.split('-')[-1] == 'Ribosome')]
        if len(V) == 1:
            X = [val for val in reaction_obj.products if type(val) == RNA]
            assert len(X) == 1
            self.initiation_tracker[X[0].name][X[0].mol_numb].append(system_time)
            self.RNA_output[X[0].name][X[0].mol_numb] += 1
     
    def check_initialize_RNA_objects(self, V, RNA_reaction_init, system_time, tau_list):
        X = self.RNA_synthesis(V)
        Y = self.RNA_reaction_init(X, RNA_reaction_init, system_time, tau_list)
        self.update_dep_graph(Y)
        self.initialize_mRNA_times(X, system_time)

    def create_RNA_data_structures(self):   
        RNA_times = {}
        RNA_output = {}
        initiation_tracker = {}
        for i in self.mRNA_list:
            RNA_times.update({i:{}})
            RNA_output.update({i:{}})
            initiation_tracker.update({i:{}})
        self.RNA_times = RNA_times
        self.RNA_output = RNA_output
        self.initiation_tracker = initiation_tracker

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
                    V = [val for val in self.species if type(val) == IsomerComplex]
                    elong = self.RNAP_elong_check
                    for j in V:
                        count += j.count
                    assert num >= (i.count + count + elong)
                    if num > (i.count + count + elong):
                        value = num - (i.count + count + elong)
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

    def NRM_initialization_protocol(self,DNA_init,reaction_constants):

        self.create_DNA_objects(DNA_init)

        self.create_Isomer_species()

        self.create_transcript_activation_species()

        self.create_transcribe_reactions()

        self.create_isomerization_reactions(reaction_constants)

        self.isomerization_stim_activity()
        
        self.create_transcription_elongation_species()

        self.create_transcription_elongation_reactions()

        self.create_translation_elongation_species()

        self.create_translation_elongation_reactions()

        self.update_reaction_list()

        self.unique_reactions_promoter()

    def promoter_decoupling(self, reaction_index, tau_list, system_time):

        self.check_transcription_elong_initiation(reaction_index)

        self.check_transcription_elongation_complete(reaction_index)

        self.check_translation_elongation_initiation(reaction_index)

        self.check_translation_elongation_complete(reaction_index)

        self.reverse_initiation_constant(system_time, tau_list)

        self.reverse_stim_constant()
       
    def NRM_execution(self, step, leap, end, DNA_init, reaction_constants, RNA_reaction_init):

        system_time = 0

        self.RNAP_elong_check = 0

        self.Ribosome_elong_check = 0

        self.NRM_initialization_protocol(DNA_init,reaction_constants)

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

        self.mRNA_list = ['Cro','CI','CII','CIII','N']

        self.create_RNA_data_structures()

        self.find_first_translate_elongation()

        self.find_last_translate_elongation()
        
        self.RNA_synth_react = []

        self.reaction_check = 0

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
                for j in self.P1:
                    P1 += j.count
                P1_m = (P1/float(6.02e23)) * (1/float(cell_volume))
                P1_tot = 0
                for j in self.P1_tot:
                    P1_tot += j.count
                P1_tot_m = (P1_tot/float(6.02e23)) * (1/float(cell_volume))
                P2 = 0
                for j in self.P2:
                    P2 += j.count
                P2_m = (P2/float(6.02e23)) * (1/float(cell_volume))
                P2_tot = 0
                for j in self.P2_tot:
                    P2_tot += j.count
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
                break

            if i > 0:

                self.promoter_decoupling(reaction_index, tau_list, system_time)

            self.promoter_occupancy_states(system_time,tau_list)

            reaction_index = np.argmin(tau_list)

            dependency_list = self.dep_graph[self.total_reactions[reaction_index]]

            reaction_obj = self.total_reactions[reaction_index]

            system_time = tau_list[reaction_index] #this will give us the time variable input necessary for get_det_tau calculation

            self.RNA_translate_init_time(self.total_reactions[reaction_index], system_time)

            if self.total_reactions[reaction_index].reactants == self.total_reactions[reaction_index].products:
                for j in self.total_reactions[reaction_index].reactants:
                    new = j.count - 1
                    j.count = new
                
                self.decay_mRNA_times(self.total_reactions[reaction_index], system_time) # method screens to see if reaction involves protein or RNA decay

            else:
                for j in self.total_reactions[reaction_index].reactants:
                    new = j.count - 1
                    j.count = new
                for j in self.total_reactions[reaction_index].products:
                    new = j.count + 1
                    j.count = new

            for j in dependency_list:
                if j == self.total_reactions[reaction_index]:
                    V = [val for val in j.products if type(val) == RNA and val.mol_numb == 0]

                    j.get_propensity()
                    
                    j.get_tau(system_time)

                    tau_list[reaction_index] = j.tau
                    
                    if len(V) == 1:
                        self.check_initialize_RNA_objects(V, RNA_reaction_init, system_time, tau_list)

                else:
                    Q = [val for val in self.isomer_react_dict.keys() if val == j]
                
                    if len(Q) == 1:
                        self.engage_isomer_reactions(system_time, tau_list, j)

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

            self.reaction_check += 1

    def save_stat_method(self):
        tup = (self.plot_dict,self.RNA_times,self.RNA_output)
        file_obj = open("save_file","wb")
        pickle.dump(tup, file_obj)

        
