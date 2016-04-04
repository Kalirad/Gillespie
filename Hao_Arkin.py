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
            
class DNA(Species):
    def __init__(self, name, tag_dna, mol_number, dna_length=False, count=False, nut_sequence=False, termination_sequence=False):
         


        """Initiate the species class.

        Parameters
        ----------
        name : str
               The name of the species, i.e. DNA, mRNA, etc.

        count : int, optional
                The quantity of the species, defualts to zero.

        tag_dna : int
                Integer representing promoter id.

        dna_length : int
                Total nucleotide length (default=False)

        nut_sequence : list
                First and Second elements nucleotide position. Third element tuple of reaction constants. Ex. [1000, 1060, (0..1,0.2,0.3,0.4)] 

        termination_sequence : list
                First and Second elements nucletide position. Third element tuple of reaction constants. Ex. [1000, 1060, (0.1,0.2,0.7)]

                

        Returns
        -------
        DNA object
        """


        assert type(tag_dna) == int
        if nut_sequence:
            assert type(nut_sequence) == list
        if termination_sequence:
            assert type(termination_sequence) == list
        Species.__init__(self, name, count=False)
        self.count = count
        self.tag_dna = tag_dna
        self.mol_number = mol_number
        self.dna_length = dna_length
        self.nut_sequence = nut_sequence
        self.termination_sequence = termination_sequence

class RNA(Species):

    def __init__(self, name, tag_rna, rna_length,promoter_affiliation, time, pool_num, curr_pool, absolute, pre_termination=False, count=False): 
        Species.__init__(self, name, count)
        self.tag_rna = tag_rna
        self.rna_length = rna_length
        self.promoter_affiliation = promoter_affiliation
        self.time = time
        self.pool_num = pool_num
        self.curr_pool = curr_pool
        self.absolute = absolute
        self.pre_termination = pre_termination

    def create_rna_output(self):
        temp = {}
        for i in range(self.pool_num+1):
            if i:
                temp.update({i:0})
        self.RNA_init = temp

    def create_rna_times(self):
        self.RNA_times = {}
        self.RNA_t_entry = []

    def initialize_time(self, system_time):
        self.RNA_t_entry.append(system_time)
        
    def end_elong(self, system_time):
        val = self.RNA_t_entry[0]
        del_t = system_time - val
        self.absolute += 1
        self.RNA_times.update({self.absolute:del_t})

        self.RNA_t_entry = self.RNA_t_entry[1:len(self.RNA_t_entry)]

    def calculate_rate_init(self):
        self.RNA_init[self.curr_pool] += 1

    def change_pool_attributes(self, system_time, track_leap):
        if system_time > self.time:
            self.time += track_leap
            self.curr_pool = int(self.time//track_leap)


class IsomerComplex(Species):

    def __init__(self, name, promoter_state, dna, count=False):

        """Specify object with parameters special for initital transcription complex.

        Parameters
        ----------

        name : str

        promoter_state : int (0 or 1)
                Determines whether the transcribing complex is a 'Closed' config (0) or an 'Open' config (1)

        tag_complex : int
                Relates a Closed-Complex corresponding to a particular promoter with its Open complex and DNA copy number. An int number

        DNA : str
                The name of the DNA object (Promoter) the initial transcribed complex is associated with

        count : int
            defaults to zero
         """
        Species.__init__(self, name, count=False)
        self.count = count
        self.promoter_state = promoter_state
        self.dna = dna

class Translate_Elong(Species):

    def __init__(self, name, nucleotide, rna, count=False):
        Species.__init__(self, name, count)
        assert type(rna) == str
        assert type(nucleotide) == int
        self.rna = rna
        self.nucleotide = nucleotide

class Protein(Species):

    def __init__(self, name, rna, count=False):
        Species.__init__(self, name, count)
        self.rna = rna

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

class StochSim(object):

    def __init__(self, init_const, elong_const, copy_number, k_elong, rna_deg_const, species_list):
        self.init_const = init_const
        self.elong_const = elong_const #translation_elong
        self.k_elong = k_elong #transcription elong
        self.copy_number = copy_number 
        self.species = species_list
        self.copy_number = copy_number
        self.rna_deg_const = rna_deg_const

    def create_DNA_objects(self, DNA_init):
        for i in range(self.copy_number):
            for j in DNA_init.keys():
                name = j + '-' + 'CN' + '-' +str(i+1)
                self.species.append(DNA(name, DNA_init[j][0], i+1, dna_length=DNA_init[j][1], count=1, nut_sequence=DNA_init[j][2], termination_sequence=DNA_init[j][3]))

    def create_Isomer_species(self):
        V = [val for val in self.species if type(val) == Species and val.name == 'RNAP']
        assert len(V) == 1
        for i in self.species:
            if type(i) == DNA:
                name1 =  'Closed' + '-' + i.name + '-' + str(i.mol_number)
                self.species.append(IsomerComplex(name1, 0, i.name, count=False))
                name2 = 'Open' + '-' + i.name + '-' + str(i.mol_number)
                self.species.append(IsomerComplex(name2, 1, i.name, count=False))

    def create_initiation_reaction(self):
        rxn_list = []
        rna_list = [val for val in self.species if type(val) == RNA]
        for i in rna_list:
            reactants = [val for val in self.species if type(val) == Species and val.name == 'Ribosome']
            reactants.append(i)
            name = i.name + '-' + 'Rib'
            spec_obj = Species(name,count=False)
            self.species.append(spec_obj)
            rxn_obj = Reaction(reactants,[spec_obj],self.init_const)
            rxn_list.append(rxn_obj)
        self.reactions = rxn_list

    def create_DNA_init_reactions(self):
        self.DNA_reactions = []
        Q = [val for val in self.species if type(val) == DNA]
        ISO = [val for val in self.species if type(val) == IsomerComplex and not val.promoter_state]
        for i in Q:
            for j in ISO:
                if i.name == j.dna:
                    reactants = [val for val in self.species if val.name == 'RNAP']
                    reactants.append(i)
                    rxn_obj = Reaction(reactants,[j],(0,0))
                    self.reactions.append(rxn_obj)
                    self.DNA_reactions.append(rxn_obj)

    def create_isomer_reactions(self, reaction_constants):
        Open = [val for val in self.species if type(val) == IsomerComplex and val.promoter_state]
        Closed = [val for val in self.species if type(val) == IsomerComplex and not val.promoter_state]

        for i in Closed:
            for j in Open:
                if j.dna == i.dna:
                    rxn_obj = Reaction([i],[j],reaction_constants[j.dna.split('-')[0]])
                    self.reactions.append(rxn_obj)

    def create_transcription_elongation_species(self):
        """
        Creating the transcription elongation consists of the following steps. The sequential termination and NUT 
        sites have to be modeled along the DNA. Two properties of the object include the termination sequence and the NUT sequence.
        Therefore, the species objects must include RNApolymerase on NUT sites, RNApolymerase + N protein on Nut sites,
        RNApolymerase on termination sites, RNApolymerase + N protein on termination sites. This method creates a dictionary
        of species objects (DNA species objects in this simulation are actually promoters) as keys and values as a list of the
        individual nucleotides associated with the DNA. These nucleotides include all possible configuration of proteins on that 
        nucleic acid (RNApolymerase + N protein for example). If there is no need to model NUT and termination sites, then the
        number of species objects in the value list is simply the length of the transcript. If species and nut sites are included,
        then the value for that DNA object is a list of list with each list containing species objects of different nucleotide
        length and protein make-up (one list will contain all the species that convey RNApolymerase + N protein translocating on the 
        NUT sites). The other list will contain species objects pertaining to the termination sites including the possible combinations
        of RNApolymerase and N protein on these sites.
        """
        DNA_specs = []
        for i in self.species:
            if type(i) == DNA:
                if i.mol_number == 1:
                    DNA_specs.append(i)

        assert len(DNA_specs) == 4 #There should be four promoters
        
        elong_species_dict = {} 
        for i in DNA_specs:
            ref = i.name.split('-')[0]
            if i.nut_sequence and i.termination_sequence: #not all DNA objects will have NUT sites and termination sites
                elong_list = []
                for j in range(i.dna_length):
                    if j < i.nut_sequence[0]:
                        name = str(ref) + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                    elif i.nut_sequence[0] <= j <= i.nut_sequence[1]:
                        name = str(ref) + '-' + 'Nut' + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                    elif i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                        name = str(ref) + '-' + 'TR' + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                    else:
                        name = str(ref) + '-' + str(j)
                        elong_list.append(Species(name,count=0))
                elong_species_dict.update({i:[elong_list]})
            else: 
                elong_list = []
                for j in range(i.dna_length):
                    name = str(ref) + '-' + str(j)
                    elong_list.append(Species(name,count=0))
                elong_species_dict.update({i:elong_list})

        self.transcription_elong_species = elong_species_dict


        for i in self.transcription_elong_species.keys():
            ref = i.name.split('-')[0]
            if i.nut_sequence and i.termination_sequence:
                nut_list = []
                inter_list = []
                term_list = []
                for j in range(i.nut_sequence[0],i.termination_sequence[1]):
                    if j == i.nut_sequence[0]:
                        name = str(ref) + '-' + 'Nut' + '-' + str(j)
                        nut_list.append(Species(name,count=0))
                    elif i.nut_sequence[0] < j < i.nut_sequence[1]+1:
                        name = str(ref) + '-' + 'N' + '-' + 'Nut' + '-' + str(j)
                        nut_list.append(Species(name,count=0))
                        if j == i.nut_sequence[1]:
                            inter_list.append(nut_list[-1])
                    elif i.nut_sequence[1]+1 <= j < i.termination_sequence[0]:
                        name = str(ref) + '-' + 'N' + '-' + str(j)
                        inter_list.append(Species(name,count=0))
                self.transcription_elong_species[i].append(nut_list)
                self.transcription_elong_species[i].append(inter_list)

                if len(inter_list) != 0:
                    term_list.append(inter_list[-1])
                    for j in range(i.termination_sequence[0],i.dna_length):
                        if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                            name = str(ref) + '-' + 'N' '-' + 'TR' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                        else:
                            name = str(ref) + '-' + 'N' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                    self.transcription_elong_species[i].append(term_list)

                else:
                    term_list.append(nut_list[-1])
                    for j in range(i.termination_sequence[0],i.dna_length):
                        if i.termination_sequence[0] <= j <= i.termination_sequence[1]:
                            name = str(ref) + '-' + 'N' '-' + 'TR' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                        else:
                            name = str(ref) + '-' + 'N' + '-' + str(j)
                            term_list.append(Species(name,count=0))
                    self.transcription_elong_species[i].append(term_list)

    def create_transcription_elongation_reactions(self):
        """
        The transcription elongation reactions are represented in a dictionary. The keys in the dictionary represent the DNA species objects
        and the values mirror the list of list for the species objects representing the different possible configurations of RNApolymerase
        on the DNA. At the end of the code, the dictionary is exhaustively compiled into a single list of reactions representing all
        possible reactions associated with a given DNA object (including NUT sites and termination). This list (along with the
        translation elongation reactions) are then appended to the self.reactions list and ready for the Next Reaction Method Execution. 
        """
        elong_reactions_dict = {}
        for i in self.transcription_elong_species.keys():
            if i.nut_sequence and i.termination_sequence:
                for j in range(len(self.transcription_elong_species[i])):
                    Y = self.transcription_elong_species[i][j]
                    if j == 0:
                        elong_react = []
                        for k in self.species:
                            if type(k) == DNA:
                                if i.name.split('-')[0] == k.name.split('-')[0]:
                                    for l in self.species:
                                        if type(l) == IsomerComplex:
                                            if l.dna == k.name:
                                                if l.promoter_state == 1:
                                                    elong_react.append(Reaction([l],[k,Y[0]],(self.k_elong,self.k_elong)))
                        for k in range(len(Y)):
                            if k < i.nut_sequence[0]:
                                elong_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
                            elif i.nut_sequence[0] <= k <= i.nut_sequence[1]:
                                Z = self.transcription_elong_species[i][1]
                                for l in range(len(Z)):
                                    if 0 <= l < (len(Z) - 1):
                                        if Z[l].name.split('-')[len(Z[l].name.split('-')) - 1] == Y[k].name.split('-')[len(Y[k].name.split('-')) - 1]:
                                            for m in self.species:
                                                if type(m) == Protein:
                                                    if m.name == 'N':
                                                        elong_react.append(Reaction([Y[k],m],[Z[l+1]],(i.nut_sequence[2][1],i.nut_sequence[2][1])))
                                                        elong_react.append(Reaction([Z[l+1]],[Y[k],m],(i.nut_sequence[2][2],i.nut_sequence[2][2])))
                                elong_react.append(Reaction([Y[k]],[Y[k+1]],(i.nut_sequence[2][0],i.nut_sequence[2][0])))

                            elif i.nut_sequence[1] < k < i.termination_sequence[0]:
                                if i.nut_sequence[3] == k+1:
                                    for l in self.species:
                                        if type(l) == RNA:
                                            for n in l.promoter_affiliation:
                                                if n == i.tag_dna:
                                                    if l.pre_termination == True:
                                                        elong_react.append(Reaction([Y[k]],[Y[k+1],l],(self.k_elong,self.k_elong)))
                                else:
                                    elong_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
                            elif i.termination_sequence[0] <= k <= i.termination_sequence[1]:
                                for l in self.species:
                                    if l.name == 'RNAP':
                                        elong_react.append(Reaction([Y[k]],[Y[k+1]],(i.termination_sequence[2][0],i.termination_sequence[2][0])))
                                        elong_react.append(Reaction([Y[k]],[l],(i.termination_sequence[2][2],i.termination_sequence[2][2])))                          
                            elif i.termination_sequence[1] <= k < (i.dna_length - 1):
                                elong_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
                            else:
                                products = []
                                for l in self.species:
                                    if l.name == 'RNAP':
                                        products.append(l)
                                    elif type(l) == RNA:
                                        for m in l.promoter_affiliation:
                                            if m == i.tag_dna:
                                                if not l.pre_termination:
                                                    products.append(l)
                                elong_react.append(Reaction([Y[k]],products,(self.k_elong,self.k_elong)))
                        elong_reactions_dict.update({i:[elong_react]})
                    elif j == 1:
                        elong_nut_react = []
                        for k in range(len(Y)):
                            if 0 < k < (len(Y) - 1):
                                elong_nut_react.append(Reaction([Y[k]],[Y[k+1]],(i.nut_sequence[2][3],i.nut_sequence[2][3])))
                        elong_reactions_dict[i].append(elong_nut_react)
                    elif j == 2:
                        if len(Y) != 0:
                            elong_inter_react = []
                            for k in range(len(Y)):
                                if 0 <= k < (len(Y) - 1):
                                    if int(Y[k+1].name.split('-')[2]) == i.nut_sequence[3]:
                                        for l in self.species:
                                            if type(l) == RNA:
                                                for m in l.promoter_affiliation:
                                                    if m == i.tag_dna:
                                                        if l.pre_termination:
                                                            elong_inter_react.append(Reaction([Y[k]],[Y[k+1],l],(self.k_elong,self.k_elong)))
                                    else:
                                        elong_inter_react.append(Reaction([Y[k]],[Y[k+1],],(self.k_elong,self.k_elong)))
                            elong_reactions_dict[i].append(elong_inter_react)
                    elif j == 3:
                        elong_term_react = []
                        for k in range(len(Y)):
                            if k < ((i.termination_sequence[1] - i.termination_sequence[0]) + 1):
                                if int(Y[k+1].name.split('-')[3]) == i.nut_sequence[3]:
                                        for l in self.species:
                                            if type(l) == RNA:
                                                for m in l.promoter_affiliation:
                                                    if m == i.tag_dna:
                                                        if l.pre_termination:
                                                            elong_term_react.append(Reaction([Y[k]],[Y[k+1],l],(self.k_elong,self.k_elong)))
                                else:
                                    elong_term_react.append(Reaction([Y[k]],[Y[k+1]],(i.termination_sequence[2][1],i.termination_sequence[2][1])))
                            elif ((i.termination_sequence[1] - i.termination_sequence[0]) + 1) <= k < (len(Y) - 1):
                                elong_term_react.append(Reaction([Y[k]],[Y[k+1]],(self.k_elong,self.k_elong)))
                            elif k == (len(Y) - 1):
                                products = []
                                for l in self.species:
                                    if l.name == 'RNAP':
                                        products.append(l)
                                    elif type(l) == Protein:
                                        if l.name == 'N':
                                            products.append(l)
                                    elif type(l) == RNA:
                                        for m in l.promoter_affiliation:
                                            if m == i.tag_dna:
                                                if not l.pre_termination:
                                                    products.append(l)
                                elong_term_react.append(Reaction([Y[k]],products,(self.k_elong,self.k_elong)))
                        elong_reactions_dict[i].append(elong_term_react)
            else:
                Y = self.transcription_elong_species[i]
                elong_reactions = []
                for j in range(len(Y)):
                    if j == 0:
                        for k in self.species:
                            if type(k) == DNA:
                                if i.name.split('-')[0] == k.name.split('-')[0]: 
                                    for l in self.species:
                                        if type(l) == IsomerComplex:
                                            if l.dna == k.name:
                                                if l.promoter_state == 1:
                                                    elong_reactions.append(Reaction([l],[k,Y[j]],(self.k_elong,self.k_elong)))
                        elong_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_elong,self.k_elong)))
                    elif 0 < j < (len(Y) - 1):
                        elong_reactions.append(Reaction([Y[j]],[Y[j+1]],(self.k_elong,self.k_elong)))
                    elif j == (len(Y) - 1):
                        products = []
                        for k in self.species:
                            if k.name == 'RNAP':
                                products.append(k)
                            elif type(k) == RNA:
                                for l in k.promoter_affiliation:
                                    if l == i.tag_dna: 
                                        products.append(k)
                        elong_reactions.append(Reaction([Y[j]],products,(self.k_elong,self.k_elong)))
                elong_reactions_dict.update({i:[elong_reactions]})
            self.elong_reactions_dict = elong_reactions_dict

        temp_list = []
        for i in elong_reactions_dict.keys():
            if i.nut_sequence and i.termination_sequence:
                aggregate_list = np.sum(elong_reactions_dict[i])
                temp_list.append(aggregate_list)
            else:
                temp_list.append(elong_reactions_dict[i][0])

        for i in range(len(temp_list)):
            hold = temp_list[0]
            if i > 0:
                hold += temp_list[i]
        self.transcription_elong_reactions = hold

    def create_initiation_reaction(self):
        rxn_list = []
        rna_list = [val for val in self.species if type(val) == RNA]
        for i in rna_list:
            reactants = [val for val in self.species if type(val) == Species and val.name == 'Ribosome']
            reactants.append(i)
            name = i.name + '-' + 'Rib'
            spec_obj = Species(name,count=False)
            self.species.append(spec_obj)
            rxn_obj = Reaction(reactants,[spec_obj],self.init_const)
            rxn_list.append(rxn_obj)
        self.reactions = rxn_list

    def create_rna_deg_reactions(self):
        rna_list = [val for val in self.species if type(val) == RNA]

        for i in rna_list:
            rxn_obj = Reaction([i],[i],self.rna_deg_const)
            self.reactions.append(rxn_obj)

    def create_translate_elong_species(self):
        RNA_spec_list = [val for val in self.species if type(val) == RNA]
        rna_dict = {}
        for i in RNA_spec_list:
            elong_list = []
            for j in range(i.rna_length):
                spec_name = i.name + '-' + str(j)
                spec_obj = Translate_Elong(name = spec_name, nucleotide = j, rna = i.name, count=False)
                elong_list.append(spec_obj)
            rna_dict.update({i:elong_list})
        self.elong_dict = rna_dict

    def create_translate_elong_reactions(self):
        for i in self.elong_dict.keys():
            Y = self.elong_dict[i]
            for j in range(len(Y)):
                if not j:
                    for k in self.species:
                        X = k.name.split('-')
                        if len(X) > 1:
                            if X[0] == i.name:
                                for l in self.species:
                                    if type(l) == RNA:
                                        if i == l:
                                            rxn_obj1 = Reaction([k],[Y[j],l],self.elong_const)
                                            rxn_obj2 = Reaction([Y[j]],[Y[j+1]],self.elong_const)
                                            self.reactions.append(rxn_obj1)
                                            self.reactions.append(rxn_obj2)
                elif 0 < j < len(Y) - 1:
                    rxn_obj = Reaction([Y[j]],[Y[j+1]],self.elong_const)
                    self.reactions.append(rxn_obj)

                elif j == len(Y) - 1:
                    for k in self.species:
                        if type(k) == Protein:
                            if k.rna == i.tag_rna:
                                Q = [val for val in self.species if val.name == 'Ribosome']
                                rxn_obj = Reaction([Y[j]],[k,Q[0]],self.elong_const)
                                self.reactions.append(rxn_obj)

    def update_reaction_list(self):
        self.total_reactions = self.reactions + self.transcription_elong_reactions

    def create_dependency_graph(self):
        dep_dict = {}
        for i in self.total_reactions:
            hold = i.reactants + i.products
            trap = []
            for j in self.total_reactions:
                V = [val for val in j.reactants if val in hold]
                if len(V):
                    trap.append(j)
            dep_dict.update({i:trap})
        self.dep_graph = dep_dict

    def create_plot_dict(self, system_time):
        self.plot_dict = {}
        for i in self.species:
            self.plot_dict.update({i.name:[i.count]})
        self.plot_dict.update({'time':[system_time]})


    def promoter_selection(self, system_time, tau_list):
        for i in self.DNA_reactions:
            const = np.random.random()
            i.ks = (const,const)

            i.get_propensity()
            i.get_det_tau(system_time)
            index = self.total_reactions.index(i)
            tau_list[index] = i.tau

    def initialize_stochsim(self, DNA_init, reaction_constants):

        self.create_DNA_objects(DNA_init)

        self.create_Isomer_species()

        self.create_initiation_reaction()

        self.create_rna_deg_reactions()

        self.create_DNA_init_reactions()

        self.create_isomer_reactions(reaction_constants)

        self.create_transcription_elongation_species()

        self.create_transcription_elongation_reactions()

        self.create_translate_elong_species()

        self.create_translate_elong_reactions()

        self.update_reaction_list()

        self.create_dependency_graph()

    def NRM(self, steps, leap, track_leap, end, DNA_init, reaction_constants):
        
        self.initialize_stochsim(DNA_init, reaction_constants)

        k = leap

        tau_list = []

        system_time = 0

        self.create_plot_dict(system_time)

        for i in self.total_reactions:
            i.get_propensity()
            i.get_tau(system_time)
            tau_list.append(i.tau)

        for i in self.species:
            if type(i) == RNA:
                i.create_rna_output()
                i.create_rna_times()

        for i in range(steps):
            
            while system_time >= k:
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


            self.promoter_selection(system_time, tau_list)

            reaction_index = np.argmin(tau_list)

            dependency_list = self.dep_graph[self.total_reactions[reaction_index]]

            react_obj = self.total_reactions[reaction_index]

            system_time = tau_list[reaction_index]

            if self.total_reactions[reaction_index].reactants == self.total_reactions[reaction_index].products:
                for j in self.total_reactions[reaction_index].reactants:
                    new = j.count - 1
                    j.count = new

            else:
                for j in self.total_reactions[reaction_index].reactants:
                    new = j.count - 1
                    j.count = new
                for j in self.total_reactions[reaction_index].products:
                    new = j.count + 1
                    j.count = new

            for j in dependency_list:
                if j == self.total_reactions[reaction_index]:   
                    hit = 0
                    if len(j.reactants) == 1:
                        if type(j.reactants[0]) == IsomerComplex:
                            hit += 1
                    if hit:
                        deal = np.random.random()
                        if deal > 0.5:
                            j.get_propensity(const=True)
                        else:
                            j.get_propensity()
                    else:
                        j.get_propensity()

                    j.get_tau(system_time)

                    tau_list[reaction_index] = j.tau

                else:
                    hit = 0
                    if len(j.reactants) == 1:
                        if type(j.reactants[0]) == IsomerComplex:
                            hit += 1
                    if hit:
                        deal = np.random.random()
                        if deal > 0.5:
                            j.get_propensity(const=True)
                        else:
                            j.get_propensity()
                    else:
                        j.get_propensity()

                    j.get_det_tau(system_time)

                    tau_list[self.total_reactions.index(j)] = j.tau

            V = [val for val in react_obj.products if val.name == 'Ribosome']

            if len(V):
                for j in self.species:
                    if type(j) == RNA:
                        prot = [val for val in react_obj.products if type(val) == Protein]
                        if prot[0].rna == j.tag_rna: #protein class has attribute rna
                            j.end_elong(system_time)

            Q = [val for val in react_obj.products if type(val) == Translate_Elong and not val.nucleotide]    

            if len(Q):
                for j in self.species:
                    if type(j) == RNA:
                        if Q[0].rna == j.name:
                            j.initialize_time(system_time)

            Z = [val for val in react_obj.products if len(val.name.split('-')) > 1 and val.name.split('-')[-1] == 'Rib']

            if len(Z):
                for j in self.species:
                    if type(j) == RNA:
                        if Z[0].name.split('-')[0] == j.name:
                            j.calculate_rate_init()

            for j in self.species:
                if type(j) == RNA:
                    j.change_pool_attributes(system_time, track_leap)






















    



    
    



    