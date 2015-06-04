"""
T7 model
========

A model of T7 life cycle in E.coli host.  The model uses the Stochastic Simulation Algorithm from Gibson & Bruck (2000) [1].

[1] Gibson M.A and Bruck J. (2000), "Exact stochastic simulation of Chemical Systems with Many Species and Many Channels", J.Phys. Chem. 104:1876-1889

"""
import pandas as pd
import numpy as np
from itertools import *
from copy import *
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
        a = 1
        for i in self.reactants:
            a *= i.count
        a *= self.ks[0]
        if a:
            self.prop_old = a
        self.prop = a

    def get_tau(self, time): #this method is to recalculate the tau for the reaction that has just been executed
        if self.prop == 0:
            self.tau_old = self.tau
            self.t = time
            self.tau = np.inf
        else:
             self.tau = ((-1)*(np.log(np.random.random())))/float(self.prop)
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
                    self.tau = ((-1)*(np.log(np.random.random())))/float(self.prop)  
                else:
                    self.tau = (self.prop_old/float(self.prop))*(self.tau_old - self.t) + time
            else:
                self.tau = ((self.prop_old/float(self.prop))*((self.tau) - time)) + time

    

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
        

        
    def generate_dep_graph(self):
        self.dep_graph = {}
        for i in self.reactions:
            temp = i.reactants + i.products
            dep_rec = []
            for j in self.reactions:
                V = [val for val in temp if val in j.reactants]
                if len(V) != 0:
                    dep_rec.append(j)
            self.dep_graph[i] = dep_rec
        #return self.dep_graph
        
    def create_elong_reactions(self):
        temp = []
        for i in range(self.num_elong - 1):
            rec = Reaction([i], [i + 1], (self.k_elong, self.k_elong))
            temp.append(rec)
        temp.append(Reaction([i + 1], ['RNA'] , (self.k_elong, self.k_elong)))
        self.reactions = self.reactions + temp
            
    
    def read_from_file(self, directory):
        '''Read a file and create a list of reactions object
        '''
        self.reactions = reactions_list
        #self.species = species


    def NRM_execution(self, step):
        system_time = 0
        history= []

        tau_list = []
        for j in self.reactions:
            j.get_tau(system_time) 
            tau_list.append(j.tau)

        self.generate_dep_graph()

        for i in range(step):
            reaction_index = np.argmin(tau_list)

            dependency_list = self.dep_graph[self.reactions[reaction_index]]

            system_time = tau_list[reaction_index] #this will give us the time variable input necessary for get_det_tau calculation
            for j in dependency_list:
                for k in j.reactants:
                    new = k.count - 1
                    k.count = new
                    print k.count, k.name
                for k in j.products:
                    new = k.count + 1
                    k.count = new
                j.get_propensity()
                if j == self.reactions[reaction_index]:
                    j.get_tau(system_time)
                    tau_list[reaction_index] = j.tau
                else:
                    j.get_det_tau(system_time) 
                    tau_list[self.reactions.index(j)] = j.tau



























         
    