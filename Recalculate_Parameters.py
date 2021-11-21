"""
Sometimes smaller or a bigger interval of membrane potential for I_d gives a fit
with better performance.
Purpose of this code is to find the best fit.
"""
import numpy as np                     
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt  

trial_no = 1
name_files = "exp_array_h.cADpyr232_L5_TTPC2_a28017c6c7(0)_Trial_1.npy"  ## Insert the output file from extraction here
experiment_array_2 = np.load(str(name_files))

C = 435  ## Insert the capacitance found by the extraction (written in .txt)
experiment_array_2 = experiment_array_2[(experiment_array_2[:,4] > -17000)] ## Filter outliers

## Rest is the same steps with EIF_Extraction. Only difference is made by changing
## min_max_range and range_V intervals according to the performance of fit. 
I_d = np.array([])
min_max_range = np.arange(-90, -53.5, 0.5)  ## Change this interval
# I_d = [] ## store I_d values, found by formula I_d = <I_m> (mean of I_m for little intervals such as V-1.5 mV to V+1.5 mV)
i = 0
while i < len(min_max_range)-1:
    V_interval = experiment_array_2[(experiment_array_2[:,0] > min_max_range[i]) & (experiment_array_2[:,0] < min_max_range[i+1])]
    I_d = np.append(I_d, np.mean(V_interval[:,4])) 
    i += 1
range_V = np.arange((-89.75), (-53.75), 0.5) ## Change this interval according to min_max_range
range_V = range_V[~np.isnan(I_d)]
I_d = I_d[~np.isnan(I_d)]

plt.figure() 
plt.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 2)
plt.yticks(np.arange(-1500, 2000, 250))
plt.ylim(-1500,2000)
plt.axhline(y=0, color='k', linestyle='dotted')
plt.plot(range_V, I_d, markersize = 5)
plt.legend(loc='best')
plt.title("I_ion As a Function of V(t)")
plt.xlabel('V (mV)')
plt.ylabel("I_ion (pA)")
plt.show()
plt.savefig("Iion.png")

F_actual = -I_d/C        ## F(V) = -I_d/C from the article, find each F value for V_ms
## Find the EIF Fit for F(V)
popt, pcov = curve_fit(lambda V, tau, E_m, Vt, delta_t: (1/tau)*(E_m - V + delta_t * np.exp((V-Vt)/delta_t)),  np.array(range_V),  np.array(F_actual).flatten(),p0 = [10,-72,-55,1.5])
tau, E_m, Vt, delta_t = popt        ## Parameters are saved
F_fit = (1/tau)*(E_m-range_V+delta_t*np.exp((range_V - Vt)/delta_t)) ## F(V) = -I_d/C from the article, find each F_fit value for V_m's

## Calculating r_square of the fit
F_array = np.array(F_actual)
F_array = np.reshape(F_array, (len(F_array),))
F_fit_array = np.array(F_fit)
sum1 = sum((F_array - F_fit_array)**2)

sum2 = sum((F_array - np.mean(F_fit_array))**2)

r_sq = (sum1/sum2)


## Plot I_ion-V, I_d-V, F(V)-V, and EIF Fit together
fig, (ax1, ax2) = plt.subplots(1,2)
fig.subplots_adjust(hspace=0.5)
fig.suptitle("I_m Curve and F(V) Curve with EIF Fit For \n NEURON Cell = " + (name_files.replace("_Trial_1.npy","")).replace("exp_df_","") +"\n")
ax1.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 1)
ax1.plot(range_V, I_d,linestyle='-', color='k', markersize = 5)
ax1.set_yticks(np.arange(-1500, 2000, 250))
ax1.set_ylim(-1500,2000)
ax1.axhline(y=0, color='k', linestyle='dotted')
ax1.legend(loc='best')
ax1.set_title("I_ion As a Function of V(t)")
ax1.set_xlabel('V (mV)',fontsize=15)
ax1.set_ylabel("I_ion (pA)",fontsize=15)
ax2.plot(range_V,F_actual, marker='o', color='b', label="Dynamic I-V Data from NEURON Model")
ax2.plot(range_V,F_fit, linestyle='-', color='k', label="Fit to EIF Model")
ax2.set_title("Quantification of Dynamic I-V Curve")
ax2.axhline(y=0, color='k', linestyle='dotted')
ax2.set_xlabel("V(mV)",fontsize=15)
ax2.set_ylabel("F(V)(mV/ms)",fontsize=15)
ax2.legend()

print("Parameters of the fit: tau= "+str(popt[0])+" ms, V_rest= " + str(popt[1])+
      " mV, V_thres= "+str(popt[2])+" mV, delta_t= "+ str(popt[3])+" mV"+
      " and error of the fit is "+str(r_sq)+ 
      " error for each parameter is" + str(np.diag(pcov)))