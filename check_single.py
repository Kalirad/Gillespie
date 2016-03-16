import numpy as np

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

class Next_Reaction_Method(object):

    def __init__(self,k1_count,k2_count,k3_count,k4_count,k5_count,k6_count,const_list):

        assert type(const_list) == list
        self.k1_count = k1_count
        self.k2_count = k2_count
        self.k3_count = k3_count
        self.k4_count = k4_count
        self.k5_count = k5_count
        self.k6_count = k6_count
        self.const_list = const_list
    
    def species_list(self):
        species_list = []
        self.species = species_list

    def create_k1_species_objects(self):
        name = 'k1'
        spec_obj = Species(name,count=self.k1_count)
        self.species.append(spec_obj)

    def create_k3_species_objects(self):
        name = 'k3' 
        spec_obj = Species(name,count=self.k3_count)
        self.species.append(spec_obj)

    def create_other_species(self):
        name1 = 'k2'
        name2 = 'k4'
        name3 = 'k5'
        name4 = 'k6'
        name5 = 'pool'
        spec_obj1 = Species(name1,count=self.k2_count)
        spec_obj2 = Species(name2,count=self.k4_count)
        spec_obj3 = Species(name3,count=self.k5_count)
        spec_obj4 = Species(name4,count=self.k6_count)
        spec_obj5 = Species(name5,count=False)
        self.species.append(spec_obj1)
        self.species.append(spec_obj2)
        self.species.append(spec_obj3)
        self.species.append(spec_obj4)
        self.species.append(spec_obj5)

    def reaction_list(self):
        reactions = []
        self.reactions = reactions
        self.k6_reactions = []

    def create_k2_syth_reactions(self):
        products = [val for val in self.species if val.name == 'k2']
        reactants = []
        for i in self.species:
            if i.name == 'k1':
                reactants.append(i)
                reactants.append(i)
        self.reactions.append(Reaction(reactants,products,self.const_list[0]))

    def create_k4_synth_reactions(self):
        products = [val for val in self.species if val.name == 'k4']
        reactants = []
        for i in self.species:
            if i.name == 'k3':
                reactants.append(i)
                reactants.append(i)
        self.reactions.append(Reaction(reactants,products,self.const_list[1]))

    def create_k2_deg_reactions(self):
        reactants = [val for val in self.species if val.name == 'k2']
        products = []
        for i in self.species:
            if i.name == 'k1':
                products.append(i)
                products.append(i)
        self.reactions.append(Reaction(reactants,products,self.const_list[2]))

    def create_k4_deg_reactions(self):
        reactants = [val for val in self.species if val.name == 'k4']
        products = []
        for i in self.species:
            if i.name == 'k3':
                products.append(i)
                products.append(i)
        self.reactions.append(Reaction(reactants,products,self.const_list[3]))

    def create_k5_synth_reactions(self):
        reactants = []
        for i in self.species:
            if i.name == 'k2':
                reactants.append(i)
            elif i.name == 'k3':
                reactants.append(i)
        products = [val for val in self.species if val.name == 'k5']

        self.reactions.append(Reaction(reactants,products,self.const_list[4]))

    def create_k6_synth_reactions(self):
        reactants = []
        for i in self.species:
            if i.name == 'k1':
                reactants.append(i)
            elif i.name == 'k4':
                reactants.append(i)
        products = [val for val in self.species if val.name == 'k6']
        self.reactions.append(Reaction(reactants,products,self.const_list[5]))

    def create_k5_deg(self):
        reactants = [val for val in self.species if val.name == 'k5']
        products = [val for val in self.species if val.name == 'pool']

        self.reactions.append(Reaction(reactants,products,self.const_list[6]))

    def create_k6_deg(self):
        reactants = [val for val in self.species if val.name == 'k6']
        products = [val for val in self.species if val.name == 'pool']

        self.k6_reactions.append(Reaction(reactants,products,self.const_list[7]))

    def create_k5_k6_dim(self):
        reactants = []
        for i in self.species:
            if i.name == 'k5':
                reactants.append(i)
            elif i.name == 'k6':
                reactants.append(i)
        products = []
        for i in self.species:
            if i.name == 'k1':
                products.append(i)
            elif i.name == 'k3':
                products.append(i)
        self.k6_reactions.append(Reaction(reactants,products,self.const_list[8]))

    def create_dep_graph(self):
        temp = {}
        for i in self.reactions:
            sub = []
            hold = i.reactants + i.products
            for j in self.reactions:
                V = [val for val in hold if val in j.reactants]
                if len(V):
                    sub.append(j)
            temp.update({i:sub})
        self.dep_graph = temp

    def create_k6_obj(self, spec_list):
        assert len(spec_list) == 1
        ref = spec_list[0]
        name = ref.name + '-' + str(ref.count)
        spec_obj = Species(name,count=1)
        return spec_obj

    def create_k6_reactions(self, spec_obj):
        V = [val for val in self.k6_reactions if len([i for i in val.products if i.name == 'k6'])]
        
        assert not len(V)

        rxn_list = []
        for i in self.k6_reactions:
            Q = [val for val in i.reactants if val.name != 'k6']
            Q.append(spec_obj)
            reaction_obj = Reaction(Q,i.products,i.ks)
            rxn_list.append(reaction_obj)

        return rxn_list

    def create_dep_entries(self,rxn_list):
        temp = {}
        for i in rxn_list:
            hold = i.reactants + i.products
            sub  = []
            for j in rxn_list:
                V = [val for val in j.reactants if val in hold]
                if len(V):
                    sub.append(j)
            temp.update({i:sub})

        ulti = {}
        for i in rxn_list:
            react = [val for val in i.reactants if len(val.name.split('-')) == 1]
            prod = [val for val in i.products if len(val.name.split('-')) == 1]
            tot = react + prod
            tot = set(tot)
            ulti.update({i:list(tot)})

        omni = {}
        for i in ulti.keys():
            keep = []
            V = [val for val in self.reactions if len([k for k in val.reactants if k in ulti[i]])]
            if len(V):
                for j in V:
                    keep.append(j)
            omni.update({i:keep})

        #update newly created reaction entries in temp graph
        for i in omni.keys():
            for j in omni[i]:
                temp[i].append(j)

        #self.focal_reactions keeps all reactions involving the k6 + k5 > k1 + k3 reactions
        for i in temp.keys():
            if len(i.reactants) > 1:
                for j in self.focal_reactions:
                    temp[i].append(j)
                    self.dep_graph[j].append(i)
                self.focal_reactions.append(i)

        #update the dependency graph
        for i in temp.keys():
            self.dep_graph.update({i:temp[i]})

    def update_reaction_list(self, rxn_list, system_time, tau_list):
        assert len(self.reactions) == len(tau_list)
        for i in rxn_list:
            i.get_propensity()
            i.get_tau(system_time)
            self.reactions.append(i)
            tau_list.append(i.tau)

    def initialize_method(self):

        self.species_list()

        self.create_k1_species_objects()

        self.create_k3_species_objects()

        self.create_other_species()

        self.reaction_list()

        self.create_k2_syth_reactions()

        self.create_k4_synth_reactions()

        self.create_k2_deg_reactions()

        self.create_k4_deg_reactions()

        self.create_k5_synth_reactions()

        self.create_k6_synth_reactions()

        self.create_k5_deg()

        self.create_k6_deg()

        self.create_k5_k6_dim()

        self.create_dep_graph()

    def create_new_species_reactions(self, spec_list, system_time, tau_list):

        Y = self.create_k6_obj(spec_list)

        X = self.create_k6_reactions(Y)

        self.create_dep_entries(X)

        self.update_reaction_list(X, system_time, tau_list)

    #for trouble-shooting    
    def create_focal_reaction_list(self):
        self.focal_reactions = []

    def NRM(self,steps,leap,end):

        k = leap

        self.initialize_method()
        
        system_time = 0
       
        tau_list = []
        for i in self.reactions:
            i.get_propensity()
            i.get_tau(system_time)
            tau_list.append(i.tau)

        plot_dict = {}
        for i in self.species:
            plot_dict.update({i.name:[i.count]})

        plot_dict.update({'time':[0]})

        self.plot_dict = plot_dict

        self.focal_reactions = []

        time = []

        for i in range(steps):

            while system_time >= k:
                Z = time
                temp = {}
                for j in range(len(Z)):
                    if Z[j] <= k:
                        diff = k - Z[j]
                        temp.update({j:diff})
                diff_values = temp.values()
                if not len(diff_values):
                    for j in self.species:
                        self.plot_dict[j.name].append(j.count)
                    self.plot_dict['time'].append(k)

                    k += leap
                else:
                    min_diff = np.min(diff_values)
                    for j in temp.keys():
                        if temp[j] == min_diff:
                            update_index = j

                    for j in self.species:
                        self.plot_dict[j.name].append(j.count)
                    self.plot_dict['time'].append(k)

                    k += leap

                    if k >= float(end):
                        break

            if system_time >= float(end): #ensures that lists are at uniform lenghts based on final time point (end) and interval values (leap) 
                key = [val for val in self.plot_dict.keys() if len(self.plot_dict[val]) != (end/float(leap) + 1)]
                if len(key) != 0:
                    for j in key:
                        self.plot_dict[j].append(self.plot_dict[j][-1])
                break

            reaction_index = np.argmin(tau_list)

            system_time = tau_list[reaction_index]

            time.append(system_time)

            react_obj = self.reactions[reaction_index]

            dep_slot = self.dep_graph[react_obj]

            for j in react_obj.reactants:
                j.count -= 1
            for j in react_obj.products:
                j.count += 1

            for j in dep_slot:
                if j == react_obj:
                    j.get_propensity()
                    j.get_tau(system_time)
                    tau_list[reaction_index] = j.tau

                else:
                    index = self.reactions.index(j)
                    j.get_propensity()
                    j.get_tau(system_time)
                    tau_list[index] = j.tau

            V = [val for val in react_obj.products if val.name == 'k6']
            if len(V):
                self.create_new_species_reactions(V, system_time, tau_list)

        self.iterations = i








