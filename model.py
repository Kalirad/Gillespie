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
	return obj

class RNApol(Species):
	return obj

class DNApol(Species):
	return obj

class Protein(Species):
	return obj

class reaction(object):
	return obj

class SSA(object):
	"""
	Have different list of reactions for each set of related biological reactions, e.g. Transcription, Translation, and etc. 
	"""