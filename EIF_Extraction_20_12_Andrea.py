def extraction_EIF(cell_info, cell_path, alfa, trial_no):
    import time as tm
    start_time = tm.time()
    # def extraction_EIF_faster(cell_info, cell_path, alfa, trial_no):
    import os
    
    import numpy as np
    from OUprocess_simple import ou_proc
    #from pandas import read_pickle
    from scipy.optimize import curve_fit
    import warnings
    import BBP_type
    import os
    warnings.filterwarnings("ignore")
 #########################################################

    os.chdir(cell_path)
    print(os.getcwd)
    from neuron import h
    cell = BBP_type.BBP_type(cell_path)
    #cell= eval(cell_name) ## specify the cell here,
                                              ## can be found from run.py file of each cell
    
    T, transient, h.dt = 40000, 300, 0.025      #ms - Simulation lifetime, same as in article
                                                #ms - Transitory time, same as article, used to be sure to be far away from the first 300 ms that are "unstable" it should have been 3000s, with no noise in the expertiment I need to measure the background noise; in a simulation I don't need since I injected the noise, and the cell doesn't need to relax after a stimulation. Anyway there are 300s of silence at the beginning and the end of each injection since the cell can reach the V_rest before the new injection
                                                #ms - Integration time step, the highest as possible. It is the highest as possible using NEURON
    dt = h.dt
    
    ## 2 OU processes one slow, one fast are produced by using the code from
    ## OUprocess_simple.py with the sigmas and taus indicated below
    total_time_exp = T + transient*2
    datapoints = int(total_time_exp/0.025)   
    tau_f, tau_s = 3, 10  #fast
    
    Iinj = np.zeros(int((T+transient*2)/dt)) #vector where the generated current is saved
    
    ## lists for saving the values
    
    
    exp_num = 0
    DC_sigma = [[0,.18,.18],[.02,.18,.18], [.03,.18,.18],[.06,.18,.18],
                [0,.25,.36],[.02,.25,.36],[.03,.25,.36],[.06,.25,.36]]
    time_total = np.zeros(datapoints*len(DC_sigma))
    membrane_potential = np.zeros(datapoints*len(DC_sigma))
    currents = np.zeros(datapoints*len(DC_sigma))
    spiketimes = np.array([])
    freqs = np.array([])
    while exp_num < len(DC_sigma):
        apc        = h.APCount(cell.soma[0](0.5)) # introduce the Action potentials (APs) counter, and it is attached to the soma compartment
        apc.thresh = 0                            # the voltage threshold that if it is ovecome by the voltage membrane, then an AP is counted
        spkTimes   = h.Vector()                   # hoc vetor for saving the APs that will be stored in the Python variable spkTimes
        apc.record(spkTimes,dt)                   # the APs have been recorded in this vector
        exp_num += 1
    
        time_exp = np.linspace(0,int(total_time_exp),int((total_time_exp)/dt)) + total_time_exp * (exp_num-1)      #time vector, it contains all the simulation lifetime dt by dt
    
        # Setting up a IClamp mechanism, for somatic current injection as the techinique of current clamp
    
        stim = h.IClamp(cell.soma[0](0.5))
        DC_inj, sigma_s, sigma_f, alfa = DC_sigma[exp_num-1][0], DC_sigma[exp_num-1][1], DC_sigma[exp_num-1][2], 1  # alfa: to obtain a firing frequency in the wanted range 1-15 Hz
                                                                            # total OU input is multiplied with the alfa
        stim.dur = T+transient*2                    #simulation lifetime
        x = ou_proc(dt, tau_s, T, sigma_s, DC_inj)  # slow OU process
        y = ou_proc(dt, tau_f, T, sigma_f, DC_inj)  # fast ou process
    
        Iinj[int((transient)/dt):int((T+transient)/dt)] = alfa * (x+y)
    
        inp_vec  = h.Vector(Iinj)          #saving the injected current vector in a NEURON vector in order to use later
        inp_vec.play(stim._ref_amp, dt)    #injection into the model cell the stimulus previously generated as a current
        somavec = h.Vector()                          #it charges the soma membrane voltage values in a vector
        somavec.record(cell.soma[0](0.5)._ref_v,dt)   #record  the membrane values of voltage
    
        h.tstop  = T+transient*2  #lifetime of the simulation
        h.run()      # The simulation is finally launched
    
        #translating the hoc vector in python vector in order to manipulate them
        currents[(exp_num-1)*datapoints:((exp_num)*datapoints)] = np.asarray(inp_vec, dtype=np.float64)
        membrane_potential[(exp_num-1)*datapoints:((exp_num)*datapoints)] = np.asarray(somavec, dtype= np.float64)
        apccount = np.asarray(spkTimes, dtype= np.float64) ## AP times are saved in array to be counted and plotted
        time_total[(exp_num-1)*datapoints:((exp_num)*datapoints)] = time_exp
        if np.size(membrane_potential) != np.size(currents):  ## sometimes an error occurs and potential vector has 1 more data, this part is to prevent this
            membrane_potential = np.delete(membrane_potential,-1)
        spikes = apccount + (exp_num-1)*total_time_exp
        spiketimes = np.hstack((spiketimes, spikes))
        freqs = np.hstack((freqs, np.size(spikes)/(T/1000)))
    cell_type = str(cell_info[0])
    def filename(keyword):
        name = "exp_array_" + str(cell_type) + "_Trial_" +str(trial_no) + "_" + keyword + ".npy"
        return name

    # membrane_potential, currents, time_total, freqs, spiketimes, mem_pot_exp = extraction_EIF_faster(['TTPC1', '1'], "C:/L5_TTPC1_cADpyr232_1", 1, 5)
    dV_dt = np.diff(membrane_potential)/dt ## voltage change over time
    dV_dt = np.append(dV_dt,0)
    np.save(filename("v_mem"), membrane_potential)
    np.save(filename("spiketimes"),spiketimes)
    experiment_array = np.vstack((membrane_potential, currents, time_total, dV_dt))
    np.save(filename("exp_array"),experiment_array)
    # membrane_potential, currents, time_total, freqs, spiketimes, mem_pot_exp = extraction_EIF_faster(['TTPC1', '1'], "C:/L5_TTPC1_cADpyr232_1", 1, 5)
    dV_dt = np.diff(membrane_potential)/dt ## voltage change over time
    dV_dt = np.append(dV_dt,0)
    experiment_array = np.vstack((membrane_potential, currents, time_total, dV_dt))
    
    V_maxes = np.array([]) ## store the max point of APs in a vector
    V_maxes_ind = np.array([]) ## store the max point of APs in a vector
    
    
    V_s = int(spiketimes[0]) ## 100 ms before the AP
    V_s_index = np.argmax(experiment_array[2,:] > V_s)
    V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
    V_max = V_s_plus*0.025+ V_s## use V_s to find where AP has peaked
    V_maxes = np.append(V_maxes, V_max)
    V_max_array_ind = int(V_s_index+V_s_plus)
    V_max_last_ind = V_max_array_ind
    V_last_spike = V_s
    
    
    not_keeping_next = 0
    spike = 1
    
    while spike < len(spiketimes):
    # while spike < 6:
        V_s = int(spiketimes[spike]) ## 100 ms before the AP
        V_s_index = np.argmax(experiment_array[2,:] > V_s)
        V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
        V_max = V_s_plus*0.025+ V_s## use V_s to find where AP has peaked
        V_maxes = np.append(V_maxes, V_max)
        V_max_array_ind = int(V_s_index+V_s_plus)
    
        if V_last_spike + 200 >= V_s and not_keeping_next ==0:
            print("This is the "+str(spike)+"th spike. We will skip that.")
            V_old = V_last_spike
            V_old_ind = V_max_last_ind
            not_keeping_next = 1
            spike += 1
        elif V_last_spike + 200 >= V_s and not_keeping_next ==1:
            spike += 1
            print("This is the "+str(spike)+"th spike. We will skip that.")
        elif V_last_spike + 200 < V_s and not_keeping_next ==1:
            V_max_last_200 = (V_max_last_ind+int(200/dt))    
            delete_range = np.arange(V_old_ind, V_max_last_200+1,1)
            experiment_array = np.delete(experiment_array, delete_range,1)   
            V_maxes_ind = np.append(V_maxes_ind, V_old_ind)
            V_s_index = np.argmax(experiment_array[2,:] > V_s)
            V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
            V_max_array_ind = int(V_s_index+V_s_plus)     ## helps loop to skip the next spike     
            not_keeping_next = 0
            print("This is the "+str(spike)+"th spike. Finally filtered.")
            spike += 1
        elif V_last_spike + 200 < V_s and not_keeping_next ==0:
            V_max_last_200 = (V_max_last_ind+int(200/dt))    
            delete_range = np.arange(V_max_last_ind, V_max_last_200+1,1)
            experiment_array = np.delete(experiment_array, delete_range,1) 
            V_maxes_ind = np.append(V_maxes_ind, V_max_last_ind)
            V_s_index = np.argmax(experiment_array[2,:] > V_s)
            V_s_plus = (np.argmax(experiment_array[0][V_s_index:(V_s_index+100)]))
            V_max_array_ind = int(V_s_index+V_s_plus) 
            spike += 1
            print("This is the "+str(spike)+"th spike. Filtered.")
        V_max_last_ind = V_max_array_ind
        V_last_spike = V_s
    if not_keeping_next==1:
        V_max_last_200 = (V_max_last_ind+int(200/dt))    
        delete_range = np.arange(V_old_ind, V_max_last_200+1,1)
        experiment_array = np.delete(experiment_array, delete_range,1) 
        print("This is the last spike. Finally filtered.")
    
    
    # transient = 300
    # experiment_array = np.delete(experiment_array,np.s_[:int(transient/0.025)],axis=1)
    # ## Delete data after the last maximum point 
    # experiment_array = np.delete(experiment_array,np.s_[-int(transient/0.025):],axis=1)

    experiment_array = experiment_array.T[np.argsort(experiment_array.T[:, 0])]
    experiment_array_2 = experiment_array[(experiment_array[:,0] < -60)]
    
    C_e = np.arange(0,1,0.005) ## nF, candidates to be the C, use each and calculate varience
                                ## as the increment gets smaller, accuracy increases
    var_list = np.array([])            ## store the variances to find the minimum
    ce = 0
    while ce < np.size(C_e):
        inj_ce = experiment_array_2[:,1]/C_e[ce] ## to find var[I_inj/C-dVdt], I_inj/C_e values put in a vector
        difference = inj_ce - experiment_array_2[:,3]
        var_list = np.hstack((var_list,np.var(difference)))
        ce += 1
    C_e = C_e[~np.isnan(var_list)]
    var_list = var_list[~np.isnan(var_list)]
    minimum = np.argmin(var_list)     ## where the minimum of variance
    C = C_e[minimum] *1000
    print("Capacitance of cell is " + str(C)+ " pF")
    np.savez(cell_path+"/capacitance.npz", Capacitance =C)
    # plt.figure()
    # plt.plot(C_e, var_list)
    # plt.ylim(min(var_list)-5,min(var_list)+5)
    # plt.xlim(C_e[minimum]-0.3,C_e[minimum]+0.2)
    # plt.scatter(C_e[minimum], var_list[minimum], color="red", )
    # plt.annotate("minimum", (C_e[minimum], var_list[minimum]))
    # plt.ylabel("Var [I_in/C_e - dV/dt]",fontsize=15)
    # plt.xlabel('C_e (nF)',fontsize=15)
    # plt.savefig("find_c.png")
    
    experiment_array_2 = experiment_array[(experiment_array[:,0] > -95) & (experiment_array[:,0] < -40)]
    
    inj_pA = experiment_array_2[:,1]*1000           ## I_inj converted from nA to pA
    dvdt_C = experiment_array_2[:,3]*C              ## for calculating I_ion = I_inj - dvdt*C
    I_ion = inj_pA - dvdt_C                                  ## Transmembrane current is calculated
    experiment_array_2 = np.column_stack((experiment_array_2, I_ion))
                                                  ## add I_ion to last dataframe
        
    I_d = np.array([])
    min_max_range = np.arange(-90, -39, 0.5)  ## from min to max, every V_m value with increment 0.5 mV
        # I_d = [] ## store I_d values, found by formula I_d = <I_m> (mean of I_m for little intervals such as V-1.5 mV to V+1.5 mV)
    i = 0
    while i < len(min_max_range)-1:
        V_interval = experiment_array_2[(experiment_array_2[:,0] > min_max_range[i]) & (experiment_array_2[:,0] < min_max_range[i+1])]
                                           ## V_m values in the interval V-0.5 mV to V+0.5mV
        I_d = np.append(I_d, np.mean(V_interval[:,4])) ## mean of the interval is found and added as a point of I_d
        i += 1
    
    
    
    range_V = np.arange((-89.75), (-39.25), 0.5) ## points of means
    range_V = range_V[~np.isnan(I_d)]
    I_d = I_d[~np.isnan(I_d)]
    
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
    
    F_actual = -I_d/C        ## F(V) = -I_d/C from the article, find each F value for V_ms
    # F_actual = F_actual[~np.isnan(F_actual)]
    
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



    # Save the results in a .txt file.
    
    with open("output_params.txt", "a") as o:
        o.write("Experiment with " + str(exp_num)+ " types of input, DC and sigma = " + str(DC_sigma)
                 + " stimulation time= " + str(T))
        o.write("Results are: firing frequency in each stimulation= "+str(freqs)+" Capacitance = "+str(C)+
                ", Parameters of the fit: tau= "+str(popt[0])+" ms, V_rest= " + str(popt[1])+" mV, V_thres= "+str(popt[2])+" mV, delta_t= "+ str(popt[3])+" mV"+
                " and error of the fit is "+str(r_sq)+ " error for each parameter is" + str(np.diag(pcov)) + " Total time passed: " 
                + str(tm.time() - start_time))
    # Save the outputs
    # Saving just the experiment_df is enough since one can construct each result using the lines
    # below the construction of the dataframe. Graphs can be easily and quickly generated with these lines.
    
    # In addition, V_maxes and mem_pot are saved to be compare the spikes with the EIF spikes.
    np.save(filename("exp_array_2"), experiment_array_2)

    #data = np.load(str(name_files))

    return popt, pcov, r_sq, np.diag(pcov), freqs

