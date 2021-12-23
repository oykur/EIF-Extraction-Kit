def extraction_EIF(cell_info, cell_path, alfa, trial_no):
    """
    A function to do simulations on BBP NEURON cells, fit the output data
    with the EIF model, and obtain parameters (tau, V_rest, V_thres, delta_T) 
    of EIF for different cell types. 
    
    Steps are from: Badel L, Lefort S, Brette R, Petersen CC, Gerstner W, Richardson MJ. 
    Dynamic I-V curves are reliable predictors ofnaturalistic pyramidal-neuron voltage traces. 
    J Neurophysiol. 2008 Feb;99(2):656-66. doi: 10.1152/jn.01107.2007. Epub2007 Dec 5. PMID: 18057107.
    where extraction is done on real L5 Pyramidal neurons. 
    
    INPUTS:
    cell_info: Gives the information about the NEURON cell. Its name and copy
    number are used to obtain code of the cell to run the simulation
    
    alfa_a: Used to multiply input values to generate an injected current that 
    generates firing rate within the range 1-15 Hz
    
    trial_no: Number of simulations done on the cell. To save the outputs of  
    each new simulation in new files. 
    
    RETURNS:
    popt: Parameters of fit obtained from EIF fit on F(V). Parameters are:
    time constant/tau (ms), resting potential/V_rest (mV), 
    threshold potential/V_thres (mV), and steepness of the curve/delta_T (ms) 
    They are stored in popt with this order. 
    
    pcov: Returns parameter covariance matrix of popt
    
    r_sq: Error of the fit 
    
    np.diag(pcov): Error of each parameter 
    
    freqs: Firing frequencies for each simulation. If out of 1-15 Hz, experiment
    should be run again.
    
    """
    import time as tm
    start_time = tm.time()
    import os
    from neuron import h
    import numpy as np
    from OU_process import ou_proc
    from pandas import read_pickle
    from scipy.optimize import curve_fit
    import warnings

    warnings.filterwarnings("ignore")
    cells_df = read_pickle("cells_df.pkl")       ## Read the cell code from the data frame
    os.getcwd()                                  ## Check current directory's path
    os.chdir(str(cell_path))                     ## Navigate to directory of cell
    
    h.load_file("stdrun.hoc")                    ## Import the run function
    h.load_file("import3d.hoc")                  ## 3D morphology
    
    ## Specific files, related to the BBP models, are loaded
    h.load_file("constants.hoc")
    h.load_file("morphology.hoc")
    h.load_file("biophysics.hoc")
    h.load_file("template.hoc")
    #####################################################################################
    cell_name = "h." + cells_df[cell_info[0]][int(cell_info[1])-1] + "(0)"
    cell= eval(cell_name)                        ## Specify the cell here by the code obtained from the data frame
                                                 ## Cell code can also be found from run.py file of each cell
    
    T, transient, h.dt = 40000, 300, 0.025       ## T: ms - Simulation lifetime, same as in the article
                                                 ## transient: ms - Transitory time, same as article, used to be sure to be far away from the first 300 ms that are "unstable"
                                                 ## h.dt: ms - Integration time step. 
    dt = h.dt
    

    total_time_exp = T + transient*2            
    datapoints = int(total_time_exp/0.025)   
    ## Generation of injected current. 
    ## 2 OU processes one slow, one fast are produced by using the code from
    ## OU_process.py with the sigmas and taus indicated below
    Iinj = np.zeros(int((T+transient*2)/dt))      ## Vector where the generated current is saved
    DC_sigma = [[0,.18,.18],[.02,.18,.18], [.03,.18,.18],[.06,.18,.18], ## First element indicates, DC content, 
                [0,.25,.36],[.02,.25,.36],[.03,.25,.36],[.06,.25,.36]]  ## second and third indicates standart deviation of OU processes
                
    tau_f, tau_s = 3, 10                          ## Time constant for fast and slow OU processes respectively
    
    ## Stimulating BBP neuron with 8 types of inputs. 
    
    exp_num = 0                                   ## Will count the experiment numbers
    ## Define zero arrays to store outputs of 8 simulation alltogether in one vector
    time_total = np.zeros(datapoints*len(DC_sigma)) 
    membrane_potential = np.zeros(datapoints*len(DC_sigma))
    currents = np.zeros(datapoints*len(DC_sigma))
    spiketimes = np.array([])                     ## Array stores spike counts of each simulation
    freqs = np.array([])                          ## Array to store firing frequency of each simulation (should be in 1-15 Hz range)
    while exp_num < len(DC_sigma):                ## Continue loop until simulation for 8 input is done
        apc        = h.APCount(cell.soma[0](0.5)) ## introduce the Action potentials (APs) counter, and it is attached to the soma compartment
        apc.thresh = 0                            ## the voltage threshold that if it is ovecome by the voltage membrane, then an AP is counted
        spkTimes   = h.Vector()                   ## hoc vetor for saving the APs that will be stored in the Python variable spkTimes
        apc.record(spkTimes,dt)                   ## the APs have been recorded in this vector
        exp_num += 1                              ## Increase count for each new experiment
        time_exp = np.linspace(0,int(total_time_exp),int((total_time_exp)/dt)) + total_time_exp * (exp_num-1) ## time vector, it contains all the simulation lifetime dt by dt
        ## Setting up a IClamp mechanism, for somatic current injection as the techinique of current clamp
        stim = h.IClamp(cell.soma[0](0.5))
        ## Parameters of input defined: 
        ## DC_inj adds a magnitude to input, sigma defines deviation of OU process from 0, and alfa is to obtain a firing frequency in the wanted range 1-15 Hz
        DC_inj, sigma_s, sigma_f, alfa = DC_sigma[exp_num-1][0], DC_sigma[exp_num-1][1], DC_sigma[exp_num-1][2], 1  
        stim.dur = T+transient*2                    ## Simulation lifetime
        x = ou_proc(dt, tau_s, T, sigma_s, DC_inj)  ## Slow OU process generated
        y = ou_proc(dt, tau_f, T, sigma_f, DC_inj)  ## Fast OU process generated
        ## Total OU input is multiplied with the alfa and put inside current vector
        Iinj[int((transient)/dt):int((T+transient)/dt)] = alfa * (x+y)
    
        inp_vec  = h.Vector(Iinj)          ## Saving the injected current vector in a NEURON vector in order to use later
        inp_vec.play(stim._ref_amp, dt)    ## Injection into the model cell the stimulus previously generated as a current
        somavec = h.Vector()                        ## It charges the soma membrane voltage values in a vector
        somavec.record(cell.soma[0](0.5)._ref_v,dt) ## Record  the membrane values of voltage
    
        h.tstop  = T+transient*2           ## Lifetime of the simulation
        h.run()                            ## The simulation is finally launched
    
        ## Translating the hoc vector in python vector in order to manipulate them
        currents[(exp_num-1)*datapoints:((exp_num)*datapoints)] = np.asarray(inp_vec, dtype=np.float64)
        membrane_potential[(exp_num-1)*datapoints:((exp_num)*datapoints)] = np.asarray(somavec, dtype= np.float64)
        apccount = np.asarray(spkTimes, dtype= np.float64)    ## AP times are saved in array to be counted and plotted
        time_total[(exp_num-1)*datapoints:((exp_num)*datapoints)] = time_exp
        if np.size(membrane_potential) != np.size(currents):  ## Sometimes an error occurs and potential vector has 1 more data, this part is to prevent this
            membrane_potential = np.delete(membrane_potential,-1)
        spikes = apccount + (exp_num-1)*total_time_exp
        spiketimes = np.hstack((spiketimes, spikes))
        freqs = np.hstack((freqs, np.size(spikes)/(T/1000)))  ## Firing frequencies are stored
    ## A function to create specific names while saving each file. 
    cell_type = str(cell_info[0])
    def filename(keyword):
        name = "exp_array_" + str(cell_type) + "_Trial_" +str(trial_no) + "_" + keyword + ".npy"
        return name

    dV_dt = np.diff(membrane_potential)/dt             ## Voltage change over time is calculate for each dt
    dV_dt = np.append(dV_dt,0)                         ## Vector has 1 value less than output vectors, change is converging to 0 in the end, thus 0 can be added to the end of array
    np.save(filename("v_mem"), membrane_potential)     ## Outputs are stored in .npy files
    np.save(filename("spiketimes"),spiketimes)
    ## Output array is generated
    experiment_array = np.vstack((membrane_potential, currents, time_total, dV_dt))
    np.save(filename("exp_array"),experiment_array)
    
    V_maxes = np.array([])                              ## store the max point of APs in a vector
    V_maxes_ind = np.array([])                          ## store the indexes of max point in the experiment array 
    
    V_s = int(spiketimes[0])                            ## time of first spike
    V_s_index = np.argmax(experiment_array[2,:] > V_s)  ## index of time point of V_s in the experiment array
    V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)])) ## V_max should be V_s_plus elements after the V_s_index
    V_max = V_s_plus*0.025+ V_s                         ## find the time point of peak of spike
    V_maxes = np.append(V_maxes, V_max)                 ## store the time of peak
    V_max_array_ind = int(V_s_index+V_s_plus)           ## find the index of peak in the experiment array
    V_max_last_ind = V_max_array_ind                    ## to be used in filtering, store this spike as the last spike
    V_last_spike = V_s
    
    not_keeping_next = 0 ## if there one or more spikes in the 200 ms after a spike, in the interval that will be filtered, 
                         ## we should fing the last spike in this interval and filter from peak of first spike to 200 ms after 
                         ## after the last spike. Thus, if the current spike is still in the 200 ms interval, not_keeping_next = 1
                         ## and changed back to 0 when the filtering is completed with the last spike. 
    spike = 1   ## start from the second spike
    
    while spike < len(spiketimes):                      ## till the end of the spikes, continue filtering
        V_s = int(spiketimes[spike]) 
        V_s_index = np.argmax(experiment_array[2,:] > V_s) 
        V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
        V_max = V_s_plus*0.025+ V_s                     ## use V_s to find where AP has peaked
        V_maxes = np.append(V_maxes, V_max)
        V_max_array_ind = int(V_s_index+V_s_plus)
    
        if V_last_spike + 200 >= V_s and not_keeping_next ==0:
            print("This is the "+str(spike)+"th spike. We will skip that.")
            V_old = V_last_spike                        ## store the first spike in the 200 ms period as the old spike
            V_old_ind = V_max_last_ind  
            not_keeping_next = 1                        ## this spike will not be filtered, yet
            spike += 1                                  ## continue with the next spike
        elif V_last_spike + 200 >= V_s and not_keeping_next ==1:  ## current spike is still in the 200 ms period, no filtering, skip
            print("This is the "+str(spike)+"th spike. We will skip that.")
            spike += 1
        elif V_last_spike + 200 < V_s and not_keeping_next ==1:      ## current spike is out of the 200 ms period, filtering step
            V_max_last_200 = (V_max_last_ind+int(200/dt))            ## 200 ms after the last spike
            delete_range = np.arange(V_old_ind, V_max_last_200+1,1)  ## delete from max of old spike to the 200 ms after the last spike
            experiment_array = np.delete(experiment_array, delete_range,1)   
            V_maxes_ind = np.append(V_maxes_ind, V_old_ind)          ## since the index of current spike is changed, store the new index
            V_s_index = np.argmax(experiment_array[2,:] > V_s)
            V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
            V_max_array_ind = int(V_s_index+V_s_plus)        
            not_keeping_next = 0                                     ## go back to 0 since this 200 ms is completed
            print("This is the "+str(spike)+"th spike. Finally filtered.")
            spike += 1
        elif V_last_spike + 200 < V_s and not_keeping_next ==0:      ## there aren't any extra spikes in the 200 ms period, filter directly
            V_max_last_200 = (V_max_last_ind+int(200/dt))    
            delete_range = np.arange(V_max_last_ind, V_max_last_200+1,1) 
            experiment_array = np.delete(experiment_array, delete_range,1) 
            V_maxes_ind = np.append(V_maxes_ind, V_max_last_ind)
            V_s_index = np.argmax(experiment_array[2,:] > V_s)
            V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
            V_max_array_ind = int(V_s_index+V_s_plus) 
            print("This is the "+str(spike)+"th spike. Filtered.")
            spike += 1
        V_max_last_ind = V_max_array_ind                             ## always store the last spike since it can be used in the next filtering
        V_last_spike = V_s
    if not_keeping_next==1:                                          ## sometimes there are a lot of spike in the very last 200 ms 
        V_max_last_200 = (V_max_last_ind+int(200/dt))                ## and loop is completed without them filtered, filter them if this is the case 
        delete_range = np.arange(V_old_ind, V_max_last_200+1,1)
        experiment_array = np.delete(experiment_array, delete_range,1) 
        print("This is the last spike. Filtering is completed.")
    
    ## Calculation of Capacitance With Variance Minimization Method
    
    experiment_array = experiment_array.T[np.argsort(experiment_array.T[:, 0])]  ## sort the array by increasing values of membrane potential
    experiment_array_2 = experiment_array[(experiment_array[:,0] < -60)]         ## store subthreshold data into another array
    C_e = np.arange(0,1,0.005)                                ## nF, candidates to be the C, each on is used to calculate corresponding varience
                                                              ## as the increment (0.005 nF) gets smaller, accuracy increases
    var_list = np.array([])                                   ## store the variances in an array to find the minimum
    ce = 0
    while ce < np.size(C_e):                         
        inj_ce = experiment_array_2[:,1]/C_e[ce]              ## to find [I_inj/C-dVdt], I_inj/C_e values put in a vector
        difference = inj_ce - experiment_array_2[:,3]         ## [I_inj/C-dVdt] values for the current C_e is found
        var_list = np.hstack((var_list,np.var(difference)))   ## calculated var[I_inj/C-dVdt] is stored 
        ce += 1                                               ## continue with the next C_e
    C_e = C_e[~np.isnan(var_list)]                            
    var_list = var_list[~np.isnan(var_list)]
    minimum = np.argmin(var_list)                             ## find where is the minimum of variance
    C = C_e[minimum] *1000
    print("Capacitance of cell is " + str(C)+ " pF")   
    
    ## Plot the variance calculated for each C_e. Graph is magnified to show minimum point. 
    # plt.figure()
    # plt.plot(C_e, var_list)
    # plt.ylim(min(var_list)-5,min(var_list)+5)
    # plt.xlim(C_e[minimum]-0.3,C_e[minimum]+0.2)
    # plt.scatter(C_e[minimum], var_list[minimum], color="red", )
    # plt.annotate("minimum", (C_e[minimum], var_list[minimum]))
    # plt.ylabel("Var [I_in/C_e - dV/dt]",fontsize=15)
    # plt.xlabel('C_e (nF)',fontsize=15)
    # plt.savefig("find_c.png")
    
    experiment_array_2 = experiment_array[(experiment_array[:,0] > -105) & (experiment_array[:,0] < -45)]
    
    inj_pA = experiment_array_2[:,1]*1000              ## I_inj converted from nA to pA
    dvdt_C = experiment_array_2[:,3]*C                 ## for calculating I_ion = I_inj - dvdt*C
    I_ion = inj_pA - dvdt_C                                  ## Transmembrane current is calculated
    experiment_array_2 = np.column_stack((experiment_array_2, I_ion))
                                                       ## add I_ion to last dataframe
    interval_dictionary = {"STPC": [-105,-50], "UTPC": [-105, -47], "TTPC2": [-105,-48]} ## best EIF fits for these cells are found in these specific volage intervals
    interval_V = interval_dictionary[cell_type]
    I_d = np.array([])
    min_max_range = np.arange(interval_V[0], interval_V[1], 0.5)  ## from min to max, every V_m value with increment 0.5 mV
    # I_d = [] ## store I_d values, found by formula I_d = <I_m> (mean of I_m for little intervals such as V-1.5 mV to V+1.5 mV)
    i = 0
    while i < len(min_max_range)-1:
        V_interval = experiment_array_2[(experiment_array_2[:,0] > min_max_range[i]) & (experiment_array_2[:,0] < min_max_range[i+1])]
                                                       ## V_m values in the interval V-0.5 mV to V+0.5mV
        I_d = np.append(I_d, np.mean(V_interval[:,4])) ## mean of the interval is found and added as a point of I_d
        i += 1
    range_V = np.arange((interval_V[0]+0.25), (interval_V[1]-0.25), 0.5) ## points of means
    range_V = range_V[~np.isnan(I_d)]
    I_d = I_d[~np.isnan(I_d)]
    
    ## plot I_m-V and I_d-V graphs together
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
    
    ## Calculation of Function of Voltage (F(V)) and Finding the Parameters for the EIF Model
    F_actual = -I_d/C                   ## F(V) = -I_d/C from the article, find each F value for V_ms    
    popt, pcov = curve_fit(lambda V, tau, E_m, Vt, delta_t: (1/tau)*(E_m - V + delta_t * np.exp((V-Vt)/delta_t)),  np.array(range_V),  np.array(F_actual).flatten(),p0 = [10,-72,-55,1.5])
    tau, E_m, Vt, delta_t = popt        ## parameters are saved
        
    F_fit = (1/tau)*(E_m-range_V+delta_t*np.exp((range_V - Vt)/delta_t)) ## F(V) = -I_d/C from the article, find each F_fit value for V_m's
        
    ## Calculating r_square of the fit
    F_array = np.array(F_actual)
    F_array = np.reshape(F_array, (len(F_array),))
    F_fit_array = np.array(F_fit)
    
    sum1 = sum((F_array - F_fit_array)**2)
    sum2 = sum((F_array - np.mean(F_fit_array))**2)
    r_sq = 1 - (sum1/sum2)  
    
    ## Plot the I_m-V and I_d-V, together with F(V)-V and EIF Fit-V
    # fig, (ax1, ax2) = plt.subplots(1,2)
    # # make a little extra space between the subplots
    # fig.subplots_adjust(hspace=0.5)
    
    # fig.suptitle("I_m Curve and F(V) Curve with EIF Fit For \n NEURON Cell = cADpyr232_L5_UTPC_5e3840b51e")
    
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



    ## Save the results in a .txt file.
    with open("output_params.txt", "a") as o:
        o.write("Experiment with " + str(exp_num)+ " types of input, DC and sigma = " + str(DC_sigma)
                 + " stimulation time= " + str(T))
        o.write("Results are: firing frequency in each stimulation= "+str(freqs)+" Capacitance = "+str(C)+
                ", Parameters of the fit: tau= "+str(popt[0])+" ms, V_rest= " + str(popt[1])+" mV, V_thres= "+str(popt[2])+" mV, delta_t= "+ str(popt[3])+" mV"+
                " and error of the fit is "+str(r_sq)+ " error for each parameter is" + str(np.diag(pcov)) + " Total time passed: " 
                + str(tm.time() - start_time))
    # Save the outputs
    # Saving just the experiment_df_2 is important to reconstruct each result using the lines
    # below the construction of the dataframe within seconds. Graphs can be easily and quickly 
    ## generated with these lines. For this job, Recalculate_Parameters_With_Outputs.py can be used.
    np.save(filename("exp_array_2"), experiment_array_2)
    return popt, pcov, r_sq, np.diag(pcov), freqs

