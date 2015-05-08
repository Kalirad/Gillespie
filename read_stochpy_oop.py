#!/usr/bin/python
import re
import itertools



class ReadFile(object):
    def __init__(self):
	self.filename = open('input_sample.txt', 'r')

    def list_species(self):
        """ Create a list of the species from the input file
        """
        file =self.filename.readlines()
        rx_species = []
        s_no_duplicate = []
        for lines in file:
            if lines.find(">") > -1:
                x= re.split("\W+", lines)
                rx_species.append(x)
        all_species = list(itertools.chain(*rx_species))
        filter_obj  = filter(None, all_species)    #filter the empty objects in the list
        [s_no_duplicate.append(item) for item in filter_obj if item not in s_no_duplicate]

        print "The list of species is:"
        print s_no_duplicate #list of species without duplicates 
        print "\n"

    def species_values(self):
        """Create a dictionary of the Species initial values from the input file
        """
        species = {}
        file = open ('input_sample.txt', 'r')
        for line in file:
            if (line.find("k")== -1) and (line.find("=")> -1):
                (key, value) = line.split("=")
                species[str(key)] = value

        species_values = {key.strip():value.strip() for key, value in species.iteritems()}
        
        print "The species values are:"
        print species_values
        print "\n"

    def constant_rates(self):
        """Create a dictionary of the constant rates and values
        """
        k_rates = {}
        file = open ('input_sample.txt', 'r')
        for line in file:
            if (line.find("k") > -1) and (line.find("=")> -1):
                (key, value) = line.split("=")
                k_rates[str(key)] = value

        rates_clean_spaces = {key.strip():value.strip() for key, value in k_rates.iteritems()}
        print "The rates constants and values are:"
        print rates_clean_spaces ##rates without the spaces    


####################################################################################
# Use the function

test= read()
test.list_species()
test.species_values()
test.constant_rates()
