def NRM_singlegene_model(S,nucleotide_count,C,n):
    
    """ S - Species involved in non-elongation/transcription reactions.
        V - Elongation species. Initalized as a vector of zeros
        C - Constant list. Must be equal in length to the total number of reactions and 
            must match the indices of the concatenated list X
        n - number of reactions/iterations """
    elong_step = np.zeros(nucleotide_count)
    
    time = np.ndarray(n)
    
    species_list = S + list(elong_step) 
    
    reactant_dict = {0:[0,1],1:[2],2:[2],3:[3],4:[len(species_list) - 1],5:[4],6:[4,5],7:[6],8:[7],9:[8]} 
    species_dict = {0:[0,1,2],1:[2,3],2:[0,1,2],3:[0,3,9],4:[1,4,(len(species_list) - 1)],5:[4],6:[5,6,4],7:[6,7],8:[4,5,7,8],9:[8]}
    
    """ Create Elongation Species Reaction """
    G = range(len(species_list))
    for i in range(len(G)):
        if i > 9:
            species_dict.update({i:[i-1,i]})
            reactant_dict.update({i:[i - 1]})
            
    species_dict_val = species_dict.values()
    reactant_dict_val = reactant_dict.values()
    
    """ Create Dependency Graph """
    
    dep_graph = {}
    
    for i in species_dict.keys():
        temp = []
        for r,p in enumerate(reactant_dict_val):
            M = [val for val in species_dict_val[i] if val in p]
            if len(M) != 0:
                temp.append(r)
                dep_graph.update({i:temp})
                
    """ Create Dictionary of Products and Reactants """
                
    
    product_dict = {}
    
    for i in species_dict.keys():
        temp = [val for val in species_dict[i] if val not in reactant_dict[i]]
        product_dict.update({i:temp})
        
    product_dict.update({5:[4],9:[8]})# necassry t 
        
        
    """ Define propensity functions """
    
    a0 = C[0]*S[0]*S[1]
    a1 = C[1]*S[2]
    a2 = C[2]*S[2]
    a3 = C[3]*S[3]
    a4 = C[4]*elong_step[len(elong_step) - 1]
    a5 = C[5]*S[4]
    a6 = C[6]*S[4]*S[5]
    a7 = C[7]*S[6]
    a8 = C[8]*S[7]
    a9 = C[9]*S[8]
    non_elong_prop = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9]
    
    elong_prop_func = []
    for i,p in enumerate(elong_step): 
        if i < (len(elong_step) - 1):
            elong_propensity = C[3] * p
            elong_prop_func.append(elong_propensity)
            
    total_react_prop = non_elong_prop + elong_prop_func
    
    """ Propensity and time storage vectors """
    
    
    last_propensity_value = []
    for i in range(len(total_react_prop)):
        last_propensity_value.append(total_react_prop[i])
    
    new_propensity_value = np.ndarray(len(total_react_prop))
    
    tau_store = np.ndarray(len(total_react_prop))
    propensity_store = np.ndarray(len(total_react_prop))
    time_store = np.ndarray(len(total_react_prop))
    
    """ Generate putative times for each reaction to occur """
    
    tau = np.ndarray(len(total_react_prop)) # time axis should be equal to length of propensity array
    for u,p in enumerate(total_react_prop): #Calculating taus associated with the non-elongation propensities
        if p == 0:
            tau[u] = inf
            tau_store[u] = 0
            propensity_store[u] = 0
            time_store[u] = 0
        else:
            tau[u] = ((-1)*(math.log(np.random.random())))/float(p)
            
    """ Initiate the SSA - n refers to the number of reactios """
    for i in range(n):
        B = tau.argmin()
        Y = reactant_dict[B]
        Z = product_dict[B]
        
        """ the following conditional statement is important due to the presence of 
        mRNA and protein decay reactions that have Y and Z equal to each other """
        if Y != Z:
            for r in Y:
                species_list[r] -= 1
            for rr in Z:
                species_list[rr] += 1
        else:
            for r in Y:
                species_list[r] -= 1
            
        q = tau[B]
        time[i] = q
        
        dependency = dep_graph[B]
        
        """ to recalculate the propensities, all one needs to do is iterate through
            the reactant list stored in the reactant dictionary - Reactant_dict - and multiply all items in the 
            list with the reaction constant- C[9] """
        
        for g in dependency:
            if g == B:
                """ the following for loop goes through the reactant list and multiplies all the reactants
                using variable M initially set tp 1. The propensity function can then be calculated
                generically - without specifying the exact reactant species involved """
                M = 1
                for k in range(len(Y)):
                    M *= species_list[Y[k]]
                total_react_prop[g] = C[g]*M
                new_propensity_value[g] = total_react_prop[g] 
                if new_propensity_value[g] == 0:
                    tau_store[g] = tau[g]
                    propensity_store[g] = last_propensity_value[g]
                    time_store[g] = time[i]
                    tau[g] = inf
                else:
                    tau[g] = (((-1)*(math.log(np.random.random())))/float(new_propensity_value[g])) + time[i]
                last_propensity_value[g] = new_propensity_value[g]
                    
            else:
                YY = reactant_dict[g]
                MM = 1
                for k in range(len(YY)):
                    MM *= species_list[YY[k]]
                total_react_prop[g] = C[g]*MM
                new_propensity_value[g] = total_react_prop[g]
                if new_propensity_value[g] == 0: 
                    if tau[g] != inf:
                        tau_store[g] = tau[g]
                        propensity_store[g] = last_propensity_value[g]
                        time_store[g] = time[i]
                        tau[g] = inf
                else:
                    if tau[g] == inf:
                        if propensity_store[g] == 0:
                            tau[g] = (((-1)*(math.log(np.random.random())))/float(new_propensity_value[g])) + time[i]
                        elif tau_store[g] == time_store[g]:
                            tau[g] = (((-1)*(math.log(np.random.random())))/float(new_propensity_value[g])) + time[i]
                        else:
                            tau[g] = (((propensity_store[g]) / float(new_propensity_value[g]))*(tau_store[g] - time_store[g])) + time[i]  
                    else:
                        tau[g] = ((last_propensity_value[g]/float(new_propensity_value[g]))*(tau[g] - time[i])) + time[i]
                last_propensity_value[g] = new_propensity_value[g]
                
    return species_list