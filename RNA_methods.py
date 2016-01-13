"""
Data Structures for the RNA methods:

    1) self.RNA_synth_react :  A list of all created reactions involving translation initiation - mRNA + Rib - mRNA*Rib.
    
    2) self.RNA_times : A dictionary with keys as numerals 1 through N and values as a list of lenght 2. First entry is time of synthesis. Second entry is time of decay.

    3) self.translate_times : A dictionary of times which signify which RNA molecule will have produced the  newly synthesized protein. 

    4) self.RNA_output : A dictionary which has ke  

"""
def RNA_synthesis(self, reaction_object, RNA_species):
    V = [val for val in self.total_reactions[reaction_object].products if len(self.total_reaction[reaction_object].products) != 1 and type(val) == RNA]
    if len(V) == 1:
        for i in RNA_species:
            if i.name == V[0].name:
                name1 = i.name + '-' + str(i.count)
                R = RNA(name = i.name, tag_rna = i.tag_ran, rna_length = rna_length, promoter_affiliation = i.promoter_affiliation, mol_numb = i.count, \
                    pre_termination = i.pre_termination, count = 1)
                self.species.append(R)
    return R

def RNA_reaction_init(self, rna_object, RNA_reaction_init):
    name1 = rna_object.name + '-' + 'Ribosome'
    for i in self.species:
        if type(i) == Species and i.name == 'Ribosome':
            spec_obj = Species(name = name1, count=False)
            self.species.append(spec_obj)
            R1 = Reaction([reaction_object,i],[spec_obj],(RNA_reaction_init[0],RNA_reaction_init[0]))
    R2 = Reaction([rna_object],[rna_object],(RNA_reaction_init[1], RNA_reaction_init[1]))
    
    R3 = [val for val in self.translation_elong_reactions if len([i for i in val.reactants if i.name.split('-')[0] == rna_object.name.split('-')[0] and  i.name.split('-')[-1] == 0]) != 0]
    
    R4 = Reaction([spec_obj],[rna_object,R3[0].reactants[0]],(RNA_reaction_init[2],RNA_reaction_init[2]))
    X = [R1,R2,R3]
    for i in X:
        self.total_reactions.append(i)
    Y = [R1,R2,R3,R4]
    return Y

def update_dep_graph(self, reaction_list):
    #First way of updating dependency graph - Cunning way
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


def check_elong_translate(self, reaction_obj):
    mRNA_list = ['Cro','CI','CII','CIII','N']
    V = [val for val in reaction_obj.products if len(val.name.split('-')) > 1 and len([i for i in val.name.split('-') if i in mRNA_list]) == 1 and val.name.split('-')[-1] != 'Ribosome']
    if len(V) == 1:
        Y = V[0].name.split('-')
        for i in self.translation_elong_dict.keys():
            if i.name == Y[0]:
                if self.translation_elong_dict[i][int(Y[-1]) + 1].prop != 0:
                    reaction_obj.get_propensity(const=True)
                else:
                    reaction_obj.get_propensity()

def initialize_mRNA_times(self, rna_object,system_time):
    self.RNA_times.update({rna_object.name:{ran_object.mol_numb:[system_time]}})

def decay_mRNA_times(self, reaction_obj, system_time):
    if reaction_obj.reactants == reaction_obj.products:
        if len(reaction_obj.reactants) == 1:
            if type(reaction_obj.reactants[0]) == mRNA:
                X = reaction_obj.reactants[0]
                self.RNA_times[X.name][X.mol_numb].append(system_time)

def RNA_translate_init_time(self, reaction_obj, system_time):
    V = [val for val in reaction_obj.reactants if len(val.name.split('-')[0]) > 1 and val.name.split('-')[-1] == 'Ribosome']
    if len(V) == 1:
        X = reaction_obj.products[0]
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



