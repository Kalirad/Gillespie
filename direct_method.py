import random as rng
import math
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

def Direct_Method(t, runs):
    total = np.zeros(runs)
    DNA= 1
    Pol = 700
    c1 = 0
    c2 = 0
    c3 = 0
    RBS = 0
    Rib = 10000
    c4 = 0
    c5 = 0
    Prot = 0

    kon = 4e7
    koff = 4
    kinit = 1.2
    kelon = 0.23
    krecy = 0
    krbs = 20
    krib = 1.15e4
    kclea = 0.14
    kprot = 0.645
    krnadeg = 2e-3
    kpdeg = 7e-4
    max_t = t   # simulation time
    for i in range(1,runs+1):
        #Reactions
        def react1(molec):

            molec[0] -=1  
            molec[1] -=1
            molec[2] +=1
            # k on
            return molec


        def react2(molec):

            molec[2] -=1  
            molec[0] +=1
            molec[1] +=1
            # k off
            return molec


        def react3(molec):
    
            molec[2] -=1  
            molec[3]  +=1
            #k initiation
            return molec

        def react4(molec):
    
            molec[3] -=1  
            molec[2] +=1
            # k recycle
            return molec


        def react5(molec):
    
            molec[3] -=1
            molec[4] +=1
            molec[0] +=1
            # k elongation
            return molec


        def react6(molec):
    
            molec[4] -=1
            molec[1] +=1
            molec[5] +=1
            # k rbs ??
            return molec

        def react7(molec):
            molec[5] -=1
            # k degradation
            return molec

        def react8(molec):
    
            molec[5] -=1
            molec[6] -=1
            molec[7] +=1
            # k ribon
            return molec

        def react9(molec):
    
            molec[7] -=1
            molec[5] +=1
            molec[8] +=1
            # k clear
            return molec

        def react10(molec):
    
            molec[8] -=1
            molec[9] +=1
            molec[6] +=1
            # k 
            return molec

        def react11(molec):
    
            molec[9] -=1
            # k protein-degradation
            return molec

        def compute_prob(molec,rates):
    
            return  (rates[0]*molec[0]*molec[1], 
                 rates[1]*molec[2],
                 rates[2]*molec[2],
                 rates[4]*molec[3], 
                 rates[3]*molec[3],
                 rates[5]*molec[4],
                 rates[9]*molec[5],
                 rates[6]*molec[5]*molec[6],
                 rates[7]*molec[7],
                 rates[8]*molec[8],
                 rates[10]*molec[9]) 

    ### Gillespie
        
        def main():
    
            time       = 0.0
            iteration  = 0                         
            molec  = [DNA, Pol, c1, c2, c3, RBS, Rib, c4, c5, Prot]
            rates      = [kon, koff, kinit, kelon, krecy, krbs, krib, kclea, kprot, krnadeg, kpdeg]            
            updaters   = [react1,react2,react3,react4, react5, react6, react7, react8, react9, react10, react11]
            
            # simulate to max_t
            while time < max_t:  #and molec[9] < 1: ##this code in # is part of the TFP old code
        
                a_i = compute_prob(molec, rates)
                a_0 = sum(a_i)

                #print "this is ai", a_i
                #print  "this is a0", a_0

                # first random number
                rand_1 = rng.random()

                tau    = (1.0/a_0) * math.log(1.0/rand_1)
                time  += tau
                
                rand_2    = rng.random()
                threshold = a_0 * rand_2
                #print "number 1 is", rand_1, "number 2 is", rand_2
                suma = 0
                count     = 0
               
                while threshold > suma:
                    suma += a_i[count]
                    count += 1
            
                # actualization of the values
                molec = updaters[count-1](molec)
              
                
                if (iteration % 1) == 0:
                    results = [time, molec[9]] #if we wanna use this to some array
                    print("iteration %d   time %4.3g protein %4.3g" % (iteration, time, molec[9])) 
                
                iteration += 1
                
            #plt.plot_date(x=time, y=molec[9])
              
            #print "Total time of the run", float(time) #Show the time for each run 

            total[i-1] = time #time of every run
        
        if __name__ == "__main__":
            main()
        