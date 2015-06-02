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
        self.tau_old = 0
        self.prop_old = 0
    
    @property
    def time(self):
        return self.time

    @time.setter
    def time(self,time):
        self.time = time

    @property    
    def propensity(self):
        a = 1
        for i in self.reactants:
            a *= i.count
        a *= self.ks[0]
        if a:
            self.prop_old = a
        return a

    @property
    def tau(self):
        if self.propensity == 0:
            self.tau = inf
        else:
             self.tau = ((-1)*(math.log(np.random.random())))/float(self.propensity)
             self.tau_old = tau
        return self.tau
    

class NextReactionMethod(object):
    
    def __init__(self, num_elong, k_elong, directory):
        self.num_elong = num_elong
        self.k_elong = k_elong
        self.read_from_file(directory)
        self.create_elong_reactions()
        
    def generate_dep_graph(self):
        self.dep_graph = {}
        for i in self.reactions:
            temp = i.reactants + i.products
            dep_rec = []
            for j in rec_list:
                if set(temp) in set(j.reactants):
                    dep_rec.append(j)
            self.dep_graph[i] = dep_rec
        
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
        self.species = species


    def NRM_singlegene_model(self, step):
         
        total_react_prop = []
        for i in reactions:
            total_react_prop.append(i.propensity)
        
        """ Propensity and time storage vectors """
        
        
        last_propensity_value = []
        for i in range(len(total_react_prop)): #gaskjsdkasud
            last_propensity_value.append(total_react_prop[i])
        
        new_propensity_value = np.ndarray(len(total_react_prop))
        
        tau_store = np.ndarray(len(total_react_prop))
        propensity_store = np.ndarray(len(total_react_prop))
        time_store = np.ndarray(len(total_react_prop))
        
        """ Generate putative times for each reaction to occur """
        
        for u,p in enumerate(total_react_prop): #Calculating taus associated with the non-elongation propensities
            if p == 0:
                tau[u] = inf
                tau_store[u] = 0
                propensity_store[u] = 0
                time_store[u] = 0
                
        """ Initiate the SSA - n refers to the number of reactios """
        for i in range(n):
            B = tau.argmin()
            Y = reactant_dict[B]
            Z = product_dict[B]
            
            """ the following conditional statement is important due to the presence of 
            mRNA and protein decay reactions that have Y and Z equal to each other """
            if Y != Z:
                for r in Y:
                    species_list[r] -= 1
                for rr in Z:
                    species_list[rr] += 1
            else:
                for r in Y:
                    species_list[r] -= 1
                
            q = tau[B]
            time[i] = q
            
            dependency = self.generate_dep_graph()
            
            """ to recalculate the propensities, all one needs to do is iterate through
                the reactant list stored in the reactant dictionary - Reactant_dict - and multiply all items in the 
                list with the reaction constant- C[9] """
            
            for g in dependency:
                if g == B:
                    """ the following for loop goes through the reactant list and multiplies all the reactants
                    using variable M initially set tp 1. The propensity function can then be calculated
                    generically - without specifying the exact reactant species involved """
                    M = 1
                    for k in range(len(Y)):
                        M *= species_list[Y[k]]
                    total_react_prop[g] = C[g]*M
                    new_propensity_value[g] = total_react_prop[g] 
                    if new_propensity_value[g] == 0:
                        tau_store[g] = tau[g]
                        propensity_store[g] = last_propensity_value[g]
                        time_store[g] = time[i]
                        tau[g] = inf
                    else:
                        tau[g] = (((-1)*(math.log(np.random.random())))/float(new_propensity_value[g])) + time[i]
                    last_propensity_value[g] = new_propensity_value[g]
                        
                else:
                    YY = reactant_dict[g]
                    MM = 1
                    for k in range(len(YY)):
                        MM *= species_list[YY[k]]
                    total_react_prop[g] = C[g]*MM
                    new_propensity_value[g] = total_react_prop[g]
                    if new_propensity_value[g] == 0: 
                        if tau[g] != inf:
                            tau_store[g] = tau[g]
                            propensity_store[g] = last_propensity_value[g]
                            time_store[g] = time[i]
                            tau[g] = inf
                    else:
                        if tau[g] == inf:
                            if propensity_store[g] == 0:
                                tau[g] = (((-1)*(math.log(np.random.random())))/float(new_propensity_value[g])) + time[i]
                            elif tau_store[g] == time_store[g]:
                                tau[g] = (((-1)*(math.log(np.random.random())))/float(new_propensity_value[g])) + time[i]
                            else:
                                tau[g] = (((propensity_store[g]) / float(new_propensity_value[g]))*(tau_store[g] - time_store[g])) + time[i]  
                        else:
                            tau[g] = ((last_propensity_value[g]/float(new_propensity_value[g]))*(tau[g] - time[i])) + time[i]
                    last_propensity_value[g] = new_propensity_value[g]
                    
        return species_list