import re
import itertools

class read_second(object):
    def __init__(self):
        self.species = {} 
        self.k_rates = {}
        
    def get_file(self):
        file = open('input_sample_2.txt', 'r')
        return file

    def get_reactants(self):
        react_all = []
        k_rates = []
        for line in self.get_file():
            if line.find(">") > -1:
                rx_species = line.split(">")
                reactants0 = rx_species[0]
                reactants = re.split("\W+", reactants0)                      
                reactant_filter = filter(None, reactants)
                react_all.append(reactant_filter)
        return react_all
        
    def get_products(self):
        prod_all = []
        for line in self.get_file():
            if line.find(">") > -1:
                rx_species = line.split(">")
                products0 = rx_species[1]
                products = re.split("\W+", products0)
                product_filter= filter(None, products)
                prod_all.append(product_filter)
        return prod_all
                
    def rx_list_of_list(self):
        react_all = []
        prod_all = []
        k_rates = []
        for line in self.get_file():
            if line.find(">") > -1:
                rx_species = line.split(">")

                reactants0 = rx_species[0]
                reactants = re.split("\W+", reactants0)                    
                reactant_filter = filter(None, reactants)

                products0 = rx_species[1]
                products = re.split("\W+", products0)
                product_filter= filter(None, products)

                react_all.append(reactant_filter)
                prod_all.append(product_filter)

            if (line.find("k") > -1) and (line.find("=") == -1):
                constant = line
                prueba = re.split("\W+", constant)
                prueba_filter = filter(None, prueba)
                k_rates.append(prueba_filter)
        solve = zip(react_all, prod_all,k_rates)
        return solve
    
    def get_k_rates(self):
        """NOT working"""
        
        k_rates = []
        for line in self.get_file():
            if (line.find("k") > -1) and (line.find("=") == -1):
                constant = line
                prueba = re.split("\W+", constant)
                prueba_filter = filter(None, prueba)
                k_rates.append(prueba_filter)
            return k_rates
            
    def get_list_species(self):
        """ Create a list of the species from the input file                                                 
        """
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
        return  filter_obj2
        
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
    
    def dict_react(self):
        number_of_rx = 3 
        rx_index = [i+1 for i in range(number_of_rx)]
        react_index = test2.get_reactants()
        dict_list = zip(rx_index, react_index)
        react_dict = dict(dict_list)
        print react_dict


    def dict_prod(self):
        number_of_rx = 3
        rx_index = [i+1 for i in range(number_of_rx)]
        prod_index = test2.get_products()
        dict_list2 = zip(rx_index, prod_index)
        prod_dict = dict(dict_list2)
        print prod_dict