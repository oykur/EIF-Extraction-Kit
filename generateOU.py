#!/usr/bin/env python
# coding: utf-8

# In[ ]:


######### General script for generating a OU process#########

#imported numpy and math module from python
import numpy as np


#definition of the function that generate a OU process with zero mean and unitary variance

def generateOU (dt, tau, T):
    
    x = np.zeros(int(T/dt)) # vector used to store the OU process
    
# x is the OU process, the number of iterations id defined by the ratio of the simulation total time over the time step.  
    
    x[0] = 0    
    #loop on all the lenght of the vector to generate the OU process
    for i in range (1, len(x)):
        epsy = np.random.normal(0, 1)# indipendent Gaussian  pseudo-random number with unitary variance and 0 mean
        x[i] = (1- dt / tau) * x[i-1] + np.sqrt( 2 * dt / tau) * epsy #generate an OU process with 0 mean e unitary variance
    #print(x)
    return x

