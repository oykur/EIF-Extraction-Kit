######### General script for generating a OU process#########
#for this function the library numpy is needed
import numpy as np

# definition of the function that generate a OU process with zero mean and unitary variance
# epsy generate an OU process with 0 mean e unitary variance
# tau is autocorrelation time of the process
def generateOU (dt, tau, T):
    np.random.seed(10)              #seed has been fixed to do all the process 
    x = np.zeros(int(T/dt)) # vector used to store the OU process, its length depends on the time parameters
    
# x is the OU process, the number of iterations id defined by the ratio of the simulation total time over the time step.  
    
    x[0] = 0    #fixed starting point
    #loop on all the lenght of the vector to generate the OU process
    for i in range (1, len(x)):
        epsy = np.random.normal(0, 1)
        x[i] = (1- dt / tau) * x[i-1] + np.sqrt( 2 * dt / tau) * epsy       #generate an OU process with 0 mean e unitary variance
    return x

