
"""
Data Structures for the RNA methods:

    1) self.RNA_synth_react :  A list of all created reactions involving translation initiation - mRNA + Rib - mRNA*Rib.
    
    2) self.RNA_times : A dictionary with keys as numerals 1 through N and values as a list of lenght 2. First entry is time of synthesis. Second entry is time of decay.

    3) self.translate_times : A dictionary of times which signify which RNA molecule will have produced the  newly synthesized protein. 

    4) self.RNA_output : A dictionary which has copy number of specific RNAs as keys and values as integers corresponding to #of proteins synthesized
    
    5) self.elong_dep_graph : A dictionary that serves as a dep graph within the larger overall dep graph. This dep graph specifically caters to translation elongation reactions.

"""
def RNA_synthesis(self, reaction_obj):
    V = [val for val in reaction_obj.products if len(reaction_obj.products) != 1 and type(val) == RNA]
    if len(V) == 1:
        for i in self.species:
            if i.name == V[0].name:
                name1 = i.name + '-' + str(i.count)
                R = RNA(name = i.name, tag_rna = i.tag_rna, rna_length = rna_length, promoter_affiliation = i.promoter_affiliation, mol_numb = i.count, \
                    pre_termination = i.pre_termination, count = 1)
                self.RNA_output[R.name].update({R.mol_numb:0})
    return R

def RNA_reaction_init(self, rna_object, RNA_reaction_init, system_time, tau_list):
    # RNA_reaction_init a list of tuples of rate constants 
    name1 = rna_object.name + '-' + 'Ribosome'
    for i in self.species:
        if type(i) == Species and i.name == 'Ribosome':
            spec_obj = Species(name = name1, count=False)
            R1 = Reaction([rna_object,i],[spec_obj],RNA_reaction_init[0])
    R2 = Reaction([rna_object],[rna_object],RNA_reaction_init[1])
    
    V = [val for val in self.translation_elong_dict[rna_object.name] if \
    len([i for i in val.reactants if i.name.split('-')[0] == rna_object.name.split('-')[0] and i.name.split('-')[-1] == str(0)]) != 0]

    R3 = V[0]
        
    R4 = Reaction([spec_obj],[rna_object,R3[0].reactants[0]],RNA_reaction_init[2])
    X = [R1,R2,R4]
    for i in X:
        self.total_reactions.append(i)
        i.get_propensity()
        i.get_tau(system_time)
        tau_list.append(i.tau)
    Y = [R1,R2,R3,R4]
    return Y

def create_elong_dep_graph(self):
    for i in self.transaltion_elong_dict.keys():
        for j,p in enumerate(self.translation_elong[i]):
            if j == 0:
                self.elong_dep_graph[i].update({p:[]})
            else:
                self.elong_dep_graph[i].update({p:[self.translation_elong[i][j-1]]})

def update_elong_dep(self,reaction_list):
    V = [val for val in reaction_list if len(val.products) > 1 and len([i for i in val.products if type(i) == RNA]) != 0]
    assert len(V) == 1
    R = [val for val in V[0].products if type(val) == RNA]
    self.elong_dep_graph[R[0].name][self.translation_elong_dict[R[0].name][0]].append(V[0])

def update_dep_graph(self, reaction_list):
    #First way of updating dependency graph - Cunning way - reaction list is the Y output from previous method
    temp_graph = {}
    for i in reaction_list:
        hold = i.reactants + i.products
        dep_rec = []
        for j in reaction_list:
            set_value = [val for val in hold if val in j.reactants]
            if len(set_value) != 0:
                dep_rec.append(j)
        temp[i] = dep_rec

    for i in self.RNA_synth_react:
        temp_graph[reaction_list[0]].append(i)
        self.dep_graph[i].append(reaction_list[0])

    for i in temp_graph.keys():
        self.dep_graph.update({i:temp_graph[i]})

    #Second way of updating dependency graph - Brute Force
    """
    self.dep_graph.update({reaction_list[0]:[reaction_list[1],reaction_list[3]]})
    self.dep_graph.update({reaction_lsit[1]:[reaction_list[0],reaction_list[2],reaction_list[3]]})
    self.dep_graph.update({reaction_list[2]:[reaction_list[0]})

    for i in self.RNA_synth_react:
        self.dep_graph[i].append(reacction_list[0])
        self.dep_graph[reaction_list[0]].append(i)
    """
    self.RNA_synth_react.append(reaction_list[0])

#Initialize mRNA function gets rna_object from output of RNA_synthesis method
def initialize_mRNA_times(self, rna_object, system_time):
    self.RNA_times[rna_object.name].update({ran_object.mol_numb:[system_time]})

def decay_mRNA_times(self, reaction_obj, system_time):
    if reaction_obj.reactants == reaction_obj.products:
        if len(reaction_obj.reactants) == 1:
            if type(reaction_obj.reactants[0]) == mRNA:
                X = reaction_obj.reactants[0]
                self.RNA_times[X.name][X.mol_numb].append(system_time)

def check_elong_translate_(self, reaction_obj, system_time, tau_list):
    V = [val for val in reaction_obj.products if len(val.name.split('-')) > 1 and len([i for i in val.name.split('-') if i in mRNA_list]) == 1 and val.name.split('-')[-1] != 'Ribosome']
    if len(V) == 1:
        if V[0].count:
            reaction_obj.get_propensity()
            reaction_obj.get_det_tau(system_time)
            tau_list[self.total_reactions.index(reaction_obj)] = reaction_obj.tau

def engage_elong_dep(self, reaction_obj, system_time, tau_list):
    QQ = [val for val in reaction_obj.products if len(val.name.split('-')) > 1 and len([i for i in val.name.split('-') if i in self.mRNA_list]) == 1 and val.name.split('-')[-1] != 'Ribosome']
    if len(QQ) == 1:
        if len(reaction_obj.reactants) == len(reaction_obj.products):
            Y = reaction_obj.products[0].name.split)('-')
            for i in self.elong_dep_graph[Y[0]][reaction_obj]:
                if i.products[0].count:
                    i.get_propensity()
                    i.get_det_tau(system_time)
                    tau_list[self.total_reactions.index(i)] = i.tau

def RNA_translate_init_time(self, reaction_obj, system_time):
    V = [val for val in reaction_obj.reactants if len(val.name.split('-')[0]) > 1 and val.name.split('-')[-1] == 'Ribosome']
    if len(V) == 1:
        X = [val for val in reaction_obj.products if type(val) == RNA]
        if len(X) == 1:
            self.translate_times[X.name].update({X.mol_numb:system_time})

def RNA_translate_finish_time(self,reaction_obj):
    V = [val for val in reaction_obj.products if type(val) == Protein]
    if len(V) == 1:
        Q = [val for val in self.translate_times[V[0].name].values() if val == np.max(self.translate_times[V[0].name].values())]
        assert len(Q) == 1
        temp = {}
        for i in self.translate_times[V[0].name].keys():
            if self.translate_times[V[0].name][i] == Q[0]:
                self.RNA_output[V[0].name][i] += 1
            else:
                temp.update({i:self.transalte_times[V[0].name][i]})
        self.translate_times[V[0].name] = temp

def check_initialize_RNA_objects(self, reaction_obj, RNA_reaction_init, system_time, tau_list):
    X = self.RNA_synthesis(reaction_obj)
    Y = self.RNA_reaction_init(X, RNA_reaction_init, system_time, tau_list)
    self.update_dep_graph(Y)
    self.initialize_mRNA_times(rna_object, system_time)