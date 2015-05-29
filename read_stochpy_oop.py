import re
import itertools

class read(object):
    
    def __init__(self):
        self.species = {} 
        self.k_rates = {}
        
    def get_file(self):
        file = open('input_sample.txt', 'r')
        return file

    def get_list_species(self):
        """ Create a list of the species from the input file                                                 
        """
        #file = self.filename.readlines()
        rx_species = []
        s_no_duplicate = []
        species = []
        for lines in self.get_file():
            if lines.find(">") > -1:
                x= re.split("\W+", lines)
                
                rx_species.append(x)
                #print "the rx", rx_species                                                                  
                all_species = list(itertools.chain(*rx_species))
                #print "after itertools", all_species                                                        
        filter_obj  = filter(None, all_species)    #filter the empty objects in the list
        
        [s_no_duplicate.append(item) for item in filter_obj if item not in s_no_duplicate]
        
        for i in s_no_duplicate:
            y = re.sub('^[0-9]+', '', i)
            species.append(y)
        filter_obj2 = filter(None,species)
        
        return  s_no_duplicate, filter_obj2
        
    def get_species_values(self):
        """Create a dictionary of the Species initial values from the input file                             
        """
       
        for lines in self.get_file():
            if (lines.find("k")== -1) and (lines.find("=")> -1):
                (key, value) = lines.split("=")
                self.species[str(key)] = value

        species_values = {key.strip():value.strip() for key, value in self.species.iteritems()}
        return species_values

    def get_constant_rates(self):
        """Create a dictionary of the constant rates and values                                              
        """
        
        for lines in self.get_file():
            if (lines.find("k") > -1) and (lines.find("=")> -1):
                (key, value) = lines.split("=")
                self.k_rates[str(key)] = value

        rates_clean_spaces = {key.strip():value.strip() for key, value in self.k_rates.iteritems()}
        return rates_clean_spaces                                                