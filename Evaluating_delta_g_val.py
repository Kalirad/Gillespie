""" Partition function for R_conc = 200nm and P_conc = 30nm """

#PArtition coefficient
L = [0,11.7,10.1,23.7,22.2,21.8,10.1,11.7,10.1,33.8,23.7]

M = []
for i in L:
    val = np.exp(i/float(310.5))
    M.append(val)

part_coef = []

R_conc = [1,200,200,(200**2),(200**2),(200**2),200,200,200,(200**3),(200**2)]
for i in range(len(R_conc)):
    R_conc[i] *= 1e-9

for i,p in enumerate(M):
    if i == 0:
        val = 1
    else:
        val = p*R_conc[i]
    part_coef.append(val)

coef_num = np.sum(part_coef)
print coef_num

#partition variable coeeficient
Q = [0,0,0,10.1,11.7,10.1,23.7]
R = [1,1,1,200,200,200,(200**2)]
P = [30,30,(30**2),30,30,30,30]

for i in range(len(R)):
    R[i] *= 1e-9
for i in range(len(P)):
    P[i] *= 1e-9

part_varia = []
for i,p in enumerate(Q):
    val = np.exp(p/float(310.5))*R[i]*P[i]
    part_varia.append(val)
    
varia_list = ['x1','x2','x3','x2','x1','x1','x1']
print part_varia

m_varia = []
for i,p in enumerate(part_varia):
    val = str(p) + ',' + varia_list[i]
    m_varia.append(val)


part_list = []
while len(m_varia) > 0:
    ref = m_varia[0].split(',')
    val0 = ref[1]
    val1 = float(ref[0])
    L = []
    for i in range(len(m_varia)):
        if i != 0:
            comp = m_varia[i].split(',')
            if comp[1] == val0:
                L.append(i)
                val1 += float(comp[0])
    tot = str(val1) + val0
    part_list.append(tot)
    if len(L) == 0:
        if len(m_varia) > 1:
            sub = []
            for i in range(len(m_varia)):
                if i != 0:
                    sub.append(m_varia[i])
            m_varia = sub
        else:
            break
    else:
        V = [val for val in range(len(m_varia)) if val != 0 and val not in L]
        sub = []
        for i in V:
            sub.append(m_varia[i])
        m_varia = sub        


part_list.append(str(coef_num))

'R_conc = 200nm'
part_funct_200nm = part_list


'System of Linear Equations'

'Equation A:'

'left-side of equation'

kprm1 = 0.011
left_side = []
for i in PRM_200:
    T = []
    ref = i.split(',')
    val0 = ref[0]
    val1 = ref[1]
    for j,p in enumerate(part_funct_0nm):
        if j != (len(part_funct_0nm) - 1):
            check = p.split(',')
            check0 = check[0]
            check1 = check[1]
            item0 = float(val0) * float(check0)*kprm1
            item1 = val1 + ',' + check1
            tot = str(item0) + ',' + item1
        else:
            check = float(p)
            item = float(val0) * check
            tot = str(item)
        T.append(tot)
    left_side.append(T)
    
v = left_side[0]
for j,p in enumerate(left_side):
    if j != 0:
        v += p
        
tot_left = v

'right-side of equation - multiplied by ratio'

master_right = []
for x in np.arange(7,11,1):
    kprm2 = 0.001
    right_side = []
    for i in PRM_0:
        T = []
        ref = i.split(',')
        val0 = ref[0]
        val1 = ref[1]
        for j,p in enumerate(part_funct_200nm):
            if j != (len(part_funct_200nm) - 1):
                check = p.split(',')
                check0 = check[0]
                check1 = check[1]
                item0 = float(val0) * float(check0)*kprm2*x
                item1 = val1 + ',' + check1
                tot = str(item0) + ',' + item1
            else:
                check = float(p)
                item = float(val0) * check
                tot = str(item)
            T.append(tot)
        right_side.append(T)

    v = right_side[0]
    for j,p in enumerate(right_side):
        if j != 0:
            v += p

    tot_right = v
    master_right.append(tot_right)
