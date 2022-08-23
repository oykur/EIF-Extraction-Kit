
#Sometimes smaller or a bigger interval of membrane potential for I_d gives a fit
#with better performance.
#Purpose of this code is to find the best fit.

path = "/home/stack/L5_TTPC1_cADpyr232_1"
import numpy as np  

from scipy.optimize import curve_fit
from scipy.stats import chisquare
from matplotlib import pyplot as plt  
import ipympl
from ipywidgets import interactive
import os
#import sys

os.chdir(path)
trial_no = 1
name_files = "exp_array_TTPC1_Trial_1_exp_array_2_cr.npy"  ## Insert the output file from extraction here
experiment_array_2 = np.load(str(name_files))

Cap = np.load(path+"/capacitance_cr.npz") # Insert the capacitance found by the extraction (written in .txt)
C = Cap["Capacitance"]
print("\n ", C, "\n")
experiment_array_2 = experiment_array_2[(experiment_array_2[:,4] > -17000)] ## Filter outliers


## Rest is the same steps with EIF_Extraction. Only difference is made by changing
## min_max_range and range_V intervals according to the performance of fit. 
I_d = np.array([])
error_I = np.array([])
min_value = -90
max_value = -52
min_max_range = np.arange(min_value, max_value, 0.5)  ## Change this interval
# I5d = [] ## store I_d values, found by formula I_d = <I_m> (mean of I_m for little intervals such as V-1.5 mV to V+1.5 mV)
i = 0
while i < len(min_max_range)-1:
    V_interval = experiment_array_2[(experiment_array_2[:,0] > min_max_range[i]) & (experiment_array_2[:,0] < min_max_range[i+1])]
    I_d = np.append(I_d, np.mean(V_interval[:,4])) 
    error_I = np.append(error_I,3.0*np.std(V_interval[:,4])/np.sqrt(len(V_interval[:,4]))) #
    i += 1
range_V = np.arange((min_value+0.25), (max_value-0.25), 0.5) ## Change this interval according to min_max_rangeec
range_V = range_V[~np.isnan(I_d)]
I_d = I_d[~np.isnan(I_d)]
error_I = error_I[~np.isnan(I_d)]

plt.figure() 
plt.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 2)
plt.yticks(np.arange(-1500, 2000, 250))
plt.ylim(-1500,2000)
plt.ylim(min_value,max_value)
plt.axhline(y=0, color='k', linestyle='dotted')
plt.plot(range_V, I_d, markersize = 5)
plt.legend(loc='best')
plt.title("I_ion As a Function of V(t)")
plt.xlabel('V (mV)')
plt.ylabel("I_ion (pA)")
plt.show()
plt.savefig("Iion.png")

F_actual = -I_d/C        ## F(V) = -I_d/C from the article, find each F value for V_ms
err_F_actual = error_I/C
#print("\n the F(V) error is: ",err_F_actual,"\n")
## Find the EIF Fit for F(V)
popt, pcov = curve_fit(lambda V, tau, E_m, Vt, delta_t: (1/tau)*(E_m - V + delta_t * np.exp((V-Vt)/delta_t)),   np.array(range_V),  np.array(F_actual).flatten(), method='trf', p0 = [10,-72,-55,1.0],  bounds=((0,-80, -65,0), (25,-68,-45, 3)),ftol=1e-15,xtol=1e-10,gtol=1e-15,maxfev=10000, loss='cauchy') #,
tau, E_m, Vt, delta_t = popt        ## Parameters are saved
F_fit = (1/tau)*(E_m-range_V+delta_t*np.exp((range_V - Vt)/delta_t)) ## F(V) = -I_d/C from the article, find each F_fit value for V_m's
#print("the function:\n",F_fit, "\n" )
## Calculating r_square of the fit
F_array = np.array(F_actual)
#F_array = np.reshape(F_array, (len(F_array),))
F_fit_array = np.array(F_fit)

print( "the square difference point by point ",(F_array - F_fit_array)**2, "\n the denominator: ",err_F_actual**2,"\n" )
chi_square = np.sum(((F_array - F_fit_array)/err_F_actual)**2);
print("the residual:\n", ((F_array - F_fit_array)/err_F_actual)**2)
#sum2 = sum((F_array - np.mean(F_fit_array))**2)
#r_sq = 1- (sum1/sum2) #chisquare(F_fit_array, F_array, ddof=len(F_fit_array)-4)

## Plot I_ion-V, I_d-V, F(V)-V, and EIF Fit together
fig = plt.figure(figsize=(10,10))
plt.legend()
#, (ax1, ax2) = plt.subplots(1,2)
#plt.subplots_adjust(hspace=0.5)
# fig.suptitle("I_m Curve and F(V) Curve with EIF Fit For \n NEURON Cell = " + (name_files.replace("_Trial_1.npy","")).replace("exp_array_","") +"\n")
# ax1.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 1)
# ax1.plot(range_V, I_d,linestyle='-', color='k', markersize = 5)
# ax1.set_yticks(np.arange(-1500, 2000, 250))
# ax1.set_ylim(-1500,2000)
# ax1.axhline(y=0, color='k', linestyle='dotted')
# ax1.legend(loc='best')
# ax1.set_title("I_ion As a Function of V(t)")
# ax1.set_xlabel('V (mV)',fontsize=15)
# ax1.set_ylabel("I_ion (pA)",fontsize=15)
# ax2.plot(range_V,F_actual, marker='o', color='b', label="Dynamic I-V Data from NEURON Model")
# ax2.plot(range_V,F_fit, linestyle='-', color='k', label="Fit to EIF Model")
# ax2.set_title("Quantification of Dynamic I-V Curve")
# ax2.axhline(y=0, color='k', linestyle='dotted')
# ax2.set_xlabel("V(mV)",fontsize=15)
# ax2.set_ylabfl("F(V)(mV/ms)",fontsize=15)
# ax2.legend()
#plt.subplot(1,2,1)
plt.plot(range_V, I_d,linestyle='-', color='k', markersize =4, label='best')
plt.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 1)
plt.xlabel('V (mV)',fontsize=15)
plt.ylim(-1500,2000)
plt.ylabel("I_ion (pA)",fontsize=15)
plt.ylabel("I_ion (pA)",fontsize=15)
plt.axhline(y=0.5, xmin=0.0, xmax=1.0, color='r')

plt.title("I_ion As a Function of V(t)")
fig = plt.figure(figsize=(10,10))
plt.legend()
#plt.subplot(1,2,2)
plt.errorbar(range_V,F_actual, yerr= err_F_actual,  fmt='o',color="blue", markersize=3,ecolor = 'blue', label="Dynamic I-V Data from NEURON Model")
plt.plot(range_V,F_actual,"bo", markersize=3, label="Dynamic I-V Data from NEURON Model")
plt.plot(range_V,F_fit, linestyle='-', color='k', label="Fit to EIF Model")
plt.title("Quantification of Dynamic I-V Curve")
plt.xlabel("V(mV)",fontsize=18)
plt.ylabel("F(V)(mV/ms)",fontsize=18)
plt.axhline(y=0, color='k', linestyle='dotted')
plt.show()
plt.savefig("F_I_accurate.png")

cell = path.split("/")
print(cell[3])

###################### save the dT#####

# np.savez(cell[3]+".npz", dT=popt[3], parameters = popt, errors = pcov)
# np.save(cell[3], popt)
# ########à#######################################
print("Parameters of the fit: tau= "+str(popt[0])+" ms, V_rest= " + str(popt[1])+
       " mV, V_thres= "+str(popt[2])+" mV, delta_t= "+ str(popt[3])+" mV"+
 " error for each parameter is" + str(np.sqrt(np.diag(pcov)))) # std over my parameter, is the error over them
print("\n the DF of freedoom are: ",len(F_fit_array)-4, " chi square: ",str(chi_square), "\n" )  #" and error of the fit is "+str(chi_square)+ 