"""
T7 model
========

A model of T7 life cycle in E.coli host.  The model uses the Stochastic Simulation Algorithm from Gillespie (1977) [1].

[1] Gillespie D.T (1977), "Exact stochastic simulation of coupled chemical reactions", J.Phys. Chem. 81:2340-2361

"""
import pandas as pd
import numpy as np
from itertools import *
from copy import *
import pickle


class Species(object):
    
    def __init__(self, name, reaction, count=False):
        
        assert type(name) == str
        assert type(reaction) == list
        self.name = name
        self.reaction = reaction
        if count:
            self.count = count
        else:
            self.count = 0

class RNApol(Species):
    
    def __init__(self, name, reaction, rate_of_decay, count=False):
        Species.__init__(self, name, reaction, count=False)
        self.r_decay = rate_of_decay
        
class DNApol(Species):
    
    def __init__(self, name, reaction, rate_of_decay, count=False):
        Species.__init__(self, name, reaction, count=False)
        self.r_decay = rate_of_decay

class Protein(Species):
    
    def __init__(self, name, reaction, rate_of_decay, count=False):
        Species.__init__(self, name, reaction, count=False)
        self.r_decay = rate_of_decay

class Reaction(object):
    
    def __init__(self, reactants, products, ks):
        self.reactants = reactants
        self.products = products
        self.ks = ks

class NextReactionMethod(object):
    
    def __init__(self, num_elong, k_elong, directory):
        self.num_elong = num_elong
        self.k_elong = k_elong
        self.read_from_file(directory)
        self.create_elong_reactions()
        
    
    #@property    
    #def propensity(self):
        
        
    def generate_dep_graph(self):
        self.dep_graph = {}
        for i in self.reactions:
            temp = i.reactants + i.products
            dep_rec = []
            for j in rec_list:
                if set(temp) and set(j.reactants):
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