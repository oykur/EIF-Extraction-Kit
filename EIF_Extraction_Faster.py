def extraction_EIF_faster(cell_info, alfa_a, trial_no):
    """
    A function to do simulations on BBP NEURON cells, fit the output data
    to EIF model, and obtain parameters (tau, V_rest, V_thres, delta_T) for 
    different cell types. 
    
    INPUTS:
    cell_info: Gives the information about the NEURON cell. Its name and copy
    number are used to obtain code of the cell to run the simulation
    
    alfa_a: Used to multiply input values to generate an injected current that 
    generates firing rate within the range 1-15 Hz
    
    trial_no: Number of experiment done with the cell. To save the outputs in 
    separate files than the files of previous experiments
    
    RETURNS:
    popt: Parameters of fit obtained from EIF fit of F(V). time constant/tau (ms),
    resting potential/V_rest (mV), threshold potential/V_thres (mV), and steepness
    of the curve/delta_T (ms) are stored in popt in this order. 
    
    pcov: Returns parameter covariance matrix of popt
    
    r_sq: Error of the fit 
    
    np.diag(pcov): Error of each parameter 
    
    freqs: Firing frequencies for each simulation. If out of 1-15 Hz, experiment
    can be run again.
    
    """
    import time as tm
    start_time = tm.time()
    from neuron import h
    import numpy as np
    from OU_Process import ou_proc
    from pandas import read_pickle
    from scipy.optimize import curve_fit
    import warnings
    warnings.filterwarnings("ignore")
    
    cells_df = read_pickle("cells_df.pkl")

    h.load_file("stdrun.hoc")         ## Import the run function
    h.load_file("import3d.hoc")       ## 3D morphology
    
    ## Specific files, related to the BBP models, are loaded
    h.load_file("constants.hoc")
    h.load_file("morphology.hoc")
    h.load_file("biophysics.hoc")
    h.load_file("template.hoc")
    #####################################################################################
    cell_name = "h." + cells_df[cell_info[0]][int(cell_info[1])-1] + "(0)"
    cell= eval(cell_name) ## specify the cell here,
    
    T, transient, h.dt = 40000, 300, 0.025      ## ms - Simulation lifetime, same as in article
                                                ## ms - Transitory time, same as article, used to be sure to be far away from the first 300 ms that are "unstable"
    dt = h.dt                                   ## ms - Integration time step, the highest as possible. It is the highest as possible using NEURON

    ## 2 OU processes one slow, one fast are produced by using the code from
    ## OU_Process.py with the sigmas and taus indicated below
    total_time_exp = T + transient*2
    datapoints = int(total_time_exp/0.025)   
    tau_f, tau_s = 3, 10  #fast
    exp_num = 0                 ## To count the simulations in while loop 
    ## 8 types of inputs same as in article. DC and sigma values for fast and slow OU processes
    DC_sigma = [[0,.18,.18],[.02,.18,.18],[.03,.18,.18],[.06,.18,.18],
                [0,.25,.36],[.02,.25,.36],[.03,.25,.36],[.06,.25,.36]]
    ## Define zero arrays to store outputs of 8 simulation alltogether in one vector
    Iinj = np.zeros(int((T+transient*2)/dt))      
    time_total = np.zeros(datapoints*len(DC_sigma)) 
    membrane_potential = np.zeros(datapoints*len(DC_sigma))
    currents = np.zeros(datapoints*len(DC_sigma))
    spiketimes = np.array([])   ## Array stores spike counts of each simulation
    freqs = np.array([])        ## Array to store firing frequency of each simulation (should be in 1-15 Hz range)
    while exp_num < len(DC_sigma): ##Continue loop until simulation for 8 input is done
        apc        = h.APCount(cell.soma[0](0.5)) ## introduce the Action potentials (APs) counter, and it is attached to the soma compartment
        apc.thresh = 0                            ## the voltage threshold that if it is ovecome by the voltage membrane, then an AP is counted
        spkTimes   = h.Vector()                   ## hoc vetor for saving the APs that will be stored in the Python variable spkTimes
        apc.record(spkTimes,dt)                   ## the APs have been recorded in this vector
        exp_num += 1                              ## Increase count for each new experiment
        time_exp = np.linspace(0,int(total_time_exp),int((total_time_exp)/dt)) + total_time_exp * (exp_num-1)      ## Time vector, it contains all the simulation lifetime dt by dt
        ## Setting up a IClamp mechanism, for somatic current injection as the techinique of current clamp
        stim = h.IClamp(cell.soma[0](0.5))
        
        ## Parameters of input defined: 
        ## DC_inj adds a magnitude to input, sigma defines deviation of OU process from 0, and alfa is to obtain a firing frequency in the wanted range 1-15 Hz
        DC_inj, sigma_s, sigma_f, alfa = DC_sigma[exp_num-1][0], DC_sigma[exp_num-1][1], DC_sigma[exp_num-1][2], alfa_a  
        stim.dur = T+transient*2                    ## Simulation lifetime
        x = ou_proc(dt, tau_s, T, sigma_s, DC_inj)  ## Slow OU process generated
        y = ou_proc(dt, tau_f, T, sigma_f, DC_inj)  ## Fast OU process generated
        ## Total OU input is multiplied with the alfa and put inside current vector
        Iinj[int((transient)/dt):int((T+transient)/dt)] = alfa * (x+y)
        
        inp_vec  = h.Vector(Iinj)          ## Saving the injected current vector in a NEURON vector in order to use later
        inp_vec.play(stim._ref_amp, dt)    ## Injection into the model cell the stimulus previously generated as a current
        somavec = h.Vector()                          ## It charges the soma membrane voltage values in a vector
        somavec.record(cell.soma[0](0.5)._ref_v,dt)   ## Record  the membrane values of voltage

        h.tstop  = T+transient*2           ## Lifetime of the simulation
        h.run()                            ## The simulation is finally launched
    
        ## Translating the hoc vector in Python vector and storing output inside arrays
        currents[(exp_num-1)*datapoints:((exp_num)*datapoints)] = np.asarray(inp_vec, dtype=np.float64)
        membrane_potential[(exp_num-1)*datapoints:((exp_num)*datapoints)] = np.asarray(somavec, dtype= np.float64)
        apccount = np.asarray(spkTimes, dtype= np.float64) ## AP times are saved in array to be counted
        time_total[(exp_num-1)*datapoints:((exp_num)*datapoints)] = time_exp
        if np.size(membrane_potential) != np.size(currents):  ## Sometimes an error occurs and potential vector has 1 more data, this part is to prevent this problem
            membrane_potential = np.delete(membrane_potential,-1)
        spikes = apccount + (exp_num-1)*total_time_exp
        spiketimes = np.hstack((spiketimes, spikes))
        freqs = np.hstack((freqs, np.size(spiketimes)/(T/1000))) ## Firing frequencies are stored
    
    spikes_store = spiketimes.copy()       ## Spikes stored in another array since it will be modified
    dV_dt = np.diff(membrane_potential)/dt ## Voltage change over time is calculate for each dt 
    dV_dt = np.append(dV_dt,0)             ## Vector has 1 value less than output vectors, change is converging to 0 in the end, thus 0 can be added to the end of array
    experiment_array = np.vstack((membrane_potential, currents, time_total, dV_dt)) ## Output array is generated
    V_maxes = np.array([])                 ## Array to store the time of max points of APs

    ## Filtering of data 200 ms after each spike 
    ## Delete transient part from array
    experiment_array = np.delete(experiment_array,np.s_[:int(transient/0.025)],axis=1)
    spiketimes -= transient
    ## Delete data after the last maximum point, final transient part
    experiment_array = np.delete(experiment_array,np.s_[-int(transient/0.025):],axis=1)
    
    spike = 0                               ## Count the spike inside the while loop
    while spike < len(spiketimes)-1:
        V_s = int(spiketimes[spike]/dt)-100 ## 100 ms before the AP
        V_max = np.argmax(experiment_array[0][V_s:(V_s+300)]) + V_s ## Use V_s to find where AP has peaked
        V_s_next = int(spiketimes[spike+1]/dt)-100 ## 100 ms before the next peak
        V_max_next = np.argmax(experiment_array[0][V_s_next:(V_s_next+300)]) + V_s_next ## Use V_s_next to find where next AP has peaked
        V_maxes = np.append(V_maxes, V_max) ## Add current peak to V_max vector
    
        if V_max*dt+200 >= V_max_next*dt:   ## If peak of next spike is in the interval 200 ms
                                            ## after the peak of current AP, next spike cannot be used,
                                            ## so data from after the current spike to 200 ms after
                                            ## the next spike is deleted from the experiment array
            experiment_array = np.delete(experiment_array,np.s_[V_max:int(V_max_next+(200/dt))],axis=1)
            spike += 1                      ## Helps loop to skip the next spike
            shift = int(V_max_next+(200/dt) - V_max)*0.025 ## With deletion, place of next spikes are shifted, find the value of shiftment 
            V_maxes = np.append(V_maxes, V_max_next) ## Add the next spike into V_max vector
            print("This is the "+ str(spike) + "th spike. Not keeping the next spike.")
        else:                               ## If there is no spike 200 ms after the current
            experiment_array = np.delete(experiment_array,np.s_[V_max:int(V_max+(200/dt))],axis=1) ## Delete data in array that is 200 ms after the current spike
            shift = 200                     ## Index for each spike in potential vector is shifted 200 ms
            print("This is the "+ str(spike) + "th spike.")                                   
        spike += 1                          ## Increase spike count 
        spiketimes -= shift                 ## Apply the shiftment 
    
    experiment_array = experiment_array.T[np.argsort(experiment_array.T[:, 0])] ## Array sorted according to V values
    experiment_array_2 = experiment_array[(experiment_array[:,0] < -60)]        ## Subthreshold values are taken for calculation of capacitance
    ## Calculation of capacitance by variance minimization method
    C_e = np.arange(0,1,0.001)              ## nF, candidates to be the C, use each and calculate varience
                                            ## As the increment gets smaller (0.001), accuracy increases
    var_list = np.array([])                 ## Store the variances to find the minimum
    ce = 0                                  ## For applying method inside a loop for each candidate C_e value
    while ce < np.size(C_e):
        inj_ce = experiment_array_2[:,1]/C_e[ce]      ## to find var[I_inj/C-dVdt], I_inj/C_e values put in a vector
        difference = inj_ce - experiment_array_2[:,3] ## [I_inj/C-dVdt] 
        var_list = np.hstack((var_list,np.var(difference))) ## Variance of [I_inj/C-dVdt] for C_e[ce]
        ce += 1                    ## Continue the loop with next C_e
    C_e = C_e[~np.isnan(var_list)] ## Delete values where variance vector has nan values
    var_list = var_list[~np.isnan(var_list)] ## Delete nan values
    minimum = np.argmin(var_list)  ## Find where is the minimum of variance
    C = C_e[minimum] *1000         ## pF, capacitance of cell
    print("Capacitance of cell is " + str(C)+ " pF")
    
    ## Visualize the variance minimization process and C value in the graph
    # plt.figure()
    # plt.plot(C_e, var_list)
    # plt.ylim(min(var_list)-5,min(var_list)+5)
    # plt.xlim(C_e[minimum]-0.3,C_e[minimum]+0.2)
    # plt.scatter(C_e[minimum], var_list[minimum], color="red", )
    # plt.annotate("minimum", (C_e[minimum], var_list[minimum]))
    # plt.ylabel("Var [I_in/C_e - dV/dt]",fontsize=15)
    # plt.xlabel('C_e (nF)',fontsize=15)
    # plt.savefig("Find_c.png")
    
    ## Finding the Transmembrane Current (I_m) and I_d (mean of I_m within small intervals)
    experiment_array_2 = experiment_array[(experiment_array[:,0] > -90) & (experiment_array[:,0] < -46)]
    inj_pA = experiment_array_2[:,1]*1000   ## I_inj converted from nA to pA
    dvdt_C = experiment_array_2[:,3]*C      ## for calculating I_ion = I_inj - dvdt*C
    I_ion = inj_pA - dvdt_C                 ## Transmembrane current is calculated
    experiment_array_2 = np.column_stack((experiment_array_2, I_ion))
                                            ## add I_ion as a column to the matrix
    I_d = np.array([]) ## Array to store mean of I_m
    min_max_range = np.arange(-90, -49, 0.5)  ## From min to max, every V_m value with increment 0.5 mV
    i = 0              ## For the while loop to be applied to each value in the range
    while i < len(min_max_range)-1:
        V_interval = experiment_array_2[(experiment_array_2[:,0] > min_max_range[i]) & (experiment_array_2[:,0] < min_max_range[i+1])]
        ## V_m values in the interval V-0.5 mV to V+0.5mV
        I_d = np.append(I_d, np.mean(V_interval[:,4])) ## Mean of the interval is found and added as a point of I_d
        i += 1         ## Next interval is called 
    
    range_V = np.arange((-89.75), (-49.25), 0.5) ## Points of means for the graph
    range_V = range_V[~np.isnan(I_d)]            ## Delete values where I_d has nans 
    I_d = I_d[~np.isnan(I_d)]                    ## Delete nan values 
    
    ## Graph I-m-V and I_d-V Together 
    # plt.figure() 
    # plt.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 2)
    # # plt.yticks(np.arange(-1500, 2000, 250))
    # # plt.ylim(-1500,2000)
    # plt.axhline(y=0, color='k', linestyle='dotted')
    # # plt.axvline(x=xreq, color='k', linestyle='-')
    # plt.plot(range_V, I_d, markersize = 5)
    # plt.legend(loc='best')
    # plt.title("I_ion As a Function of V(t)")
    # plt.xlabel('V (mV)')
    # plt.ylabel("I_ion (pA)")
    # plt.show()
    # plt.savefig("Iion.png")
    
    ## Find function of voltage (F(V))
    F_actual = -I_d/C   ## F(V) = -I_d/C from the article, find each F value for V_ms
    ## Find fit of F(V) by EIF model    
    popt, pcov = curve_fit(lambda V, tau, E_m, Vt, delta_t: (1/tau)*(E_m - V + delta_t * np.exp((V-Vt)/delta_t)),  np.array(range_V),  np.array(F_actual).flatten(),p0 = [10,-72,-55,1.5])
    tau, E_m, Vt, delta_t = popt      ## Parameters are saved
    F_fit = (1/tau)*(E_m-range_V+delta_t*np.exp((range_V - Vt)/delta_t)) ## Find each F_fit value for voltages in the range
        
    ## Calculating r_square of the fit
    F_array = np.array(F_actual)
    F_array = np.reshape(F_array, (len(F_array),))
    F_fit_array = np.array(F_fit)
    sum1 = sum((F_array - F_fit_array)**2)
    sum2 = sum((F_array - np.mean(F_fit_array))**2)
    r_sq = (sum1/sum2)
    
    ## Graph I_m-V with I_d-V and F(V)-V with EIF Fit-V Together
    # fig, (ax1, ax2) = plt.subplots(1,2)
    # # make a little extra space between the subplots
    # fig.subplots_adjust(hspace=0.5)
    # fig.suptitle("I_m Curve and F(V) Curve with EIF Fit For \n NEURON Cell = " + str(cells_df[cell_info[0]][int(cell_info[1])-1]))
    # ax1.plot(experiment_array_2[:,0], experiment_array_2[:,4], 'o', markersize = 2)
    # ax1.plot(range_V, I_d,linestyle='-', color='k')
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
    # ax2.set_ylabel("F(V)(mV/ms)",fontsize=15)
    # ax2.legend()
    # fig.savefig("Iion_Id_F_EIF.png")

    end_time =tm.time()  ## Save the ending time of simulation to see the time performance of the code
    # Save the results in a .txt file. 
    name_files = "exp_array_" + str(cell_name) + "_Trial_" +str(trial_no)
    name_files_txt = name_files +"_output.txt"
    with open(name_files_txt, "a") as o:
        o.write("\n Faster EIF Extraction Kit. "+"Experiment with " + str(exp_num)+ " types of input, DC and sigma = " + str(DC_sigma)
                 + " stimulation time= " + str(T) + " alfa= " + str(alfa_a))
        o.write("Results are: firing frequency in each stimulation= "+str(freqs)+" Capacitance = "+str(C)+
                ", Parameters of the fit: tau= "+str(popt[0])+" ms, V_rest= " + str(popt[1])+" mV, V_thres= "+str(popt[2])+" mV, delta_t= "+ str(popt[3])+" mV"+
                " and error of the fit is "+str(r_sq)+ " error for each parameter is" + str(np.diag(pcov)) + " Total time passed: " 
                + str(end_time - start_time) + "\n")
    ## Save the outputs
    ## Saving just the experiment_array_2 is enough since one can reconstruct
    ## the same or a better result with a different voltage interval
    ## by using Recalculate_Parameters.py code within seconds (see ReadMe.txt).
    name_files = name_files +".npy"
    np.save(str(name_files), experiment_array_2)
    return popt, pcov, r_sq, np.diag(pcov), freqs

