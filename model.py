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
        Protein object

        See Also
        --------
        Species : This class inherits from the species class

        """
        assert type(rate_of_decay) == float
        Species.__init__(self, name, count=False)
        self.r_decay = rate_of_decay

class Reaction(object):

    def __init__(self, reactants, products, ks):
        assert type(reactants[0]) == Species
        assert type(products[0]) == Species
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
        if len(self.reactants) > 1:
            for i in self.reactants:
                storage = []
                if i != self.reactants[0]:
                    storage.append(i)

            if len(storage) == 0:
                val = self.reactants[0].count
                for j in range(len(self.reactants)):
                    if j > 0:
                        val *= (self.reactants[j].count - j)
                val *= self.ks[0]

                factorial = [k for k in range(len(self.reactants)) if k]
                factorial.append(len(self.reactants))
                a = val/float(np.product(factorial))
                if a:
                    self.prop_old = a
                self.prop = a
            else:
                a = 1
                for m in self.reactants:
                    a *= m.count

                a *= self.ks[0]
                if a:
                    self.prop_old = a
                self.prop = a

        else:
            a = self.reactants[0].count * self.ks[0]
            if a:
                self.prop_old = a
            self.prop = a

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
    
    def __init__(self, directory=False, num_elong=False, k_elong=False):
        if num_elong:
            self.num_elong = num_elong
            self.k_elong = k_elong
            self.create_elong_reactions()
        if directory:
            self.read_from_file(directory)

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
        


        """
        reference_list = []

        if len(self.reactions) == 1:
            for i in self.reactions:
                intersect = [val for val in i.products if val not in i.reactants]
                if len(intersect) != 0:
                    temp = i.reactants + intersect
                else:
                    temp = i.reactants
                reference_list = temp
            species_list = reference_list

        else:
            for i in self.reactions:
                intersect = [val for val in i.products if val not in i.reactants]
                if len(intersect) != 0:
                    temp = i.reactants + intersect
                else:
                    temp = i.reactants
                reference_list.append(temp)


            reference = reference_list[0]
            for i in reference_list:
                if reference_list.index(i) > 0:
                    catch = [val for val in i if val not in reference]
                    if catch != 0:
                        species_list = reference + catch
                        reference = species_list  

        self.species = species_list """

    def create_elongation_species(self):
        for i in self.species:
            if i.name.split()[0] == 'DNA':
                nucleotide_length = int(i.name.split()[2])

        temp = []
        for i in range(nucleotide_length):
            temp.append(Species(self,'NL-i',count=0))



    def create_elong_reactions(self):
        temp = []
        for i in range(self.num_elong - 1):
            rec = Reaction([i], [i + 1], (self.k_elong, self.k_elong))
            temp.append(rec)
        temp.append(Reaction([i + 1], ['RNA'] , (self.k_elong, self.k_elong)))
        self.reactions = self.reactions + temp
        

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
        
    def read_from_file(self, directory):
        '''Read a file and create a list of reactions object
        '''
        self.reactions = reactions_list
        #self.species = species


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
            

            for j in dependency_list:
                
                if j == self.reactions[reaction_index]:
                    for k in j.reactants:
                        new = k.count - 1
                        k.count = new
               
                    for k in j.products:
                        new = k.count + 1
                        k.count = new
                       
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









"""def Immigration_Death(self,stop):
    def Immigration_Death_Model(self,step):
        
        system_time = np.ndarray(step)

        track_species =np.ndarray(step)


        propensity_list []
        
        for i in self.reactions:
            propensity_list.append(i.get_propensity())


        for i in range(step):
            tau = (((-1)*(math.log(np.random.random())))/float(np.sum(propensity_list)))

            system_time[i] = tau

            weighted_sum = np.random.random() * np.sum(propensity_list)

            reaction_index = [j for j in range(len(propensity_list)) if np.sum(propensity_list[0:(j + 1)]) > weighted_sum][0]

            for j in propensity_list:

                if j == propensity_list[reaction_index]:

                    for k in j.reactants:
                        k.count -= 1

                    for k in j.products:
                        k.count +"""    