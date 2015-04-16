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
        self.r_decay = r_decay
        
class DNApol(Species):
    
    def __init__(self, name, reaction, rate_of_decay, count=False):
        Species.__init__(self, name, reaction, count=False)
        self.r_decay = r_decay

class Protein(Species):
    
    def __init__(self, name, reaction, rate_of_decay, count=False):
        Species.__init__(self, name, reaction, count=False)
        self.r_decay = r_decay

class reaction(object):
	def __init__(self, reactants, products, ks):
	    self.reactants = reactants
	    self.products = products
	    self.ks = ks

class SSA(object):
	"""
	I dont know whta I am doing. Have different list of reactions for each set of related biological reactions, e.g. Transcription, Translation, and etc. 
	"""
