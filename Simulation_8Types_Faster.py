from EIF_Extraction_Faster import extraction_EIF_faster

print("Before starting, please check if you downloaded, extracted and compiled the NEURON cell."+ 
      "Then check if the EIF_Extraction_faster, Simulation_8Types_Faster, and OU_Process .py codes,"+
      "also the cells_df.pkl are all in the same file which is the NEURON cell file.")

alfa_a = 1   ## To multiply input values to generate an injected current that 
             ## generates firing rate within the range 1-15 Hz

trial_no = 1 ## Number of experiment done with the cell. 
             ## To save the outputs in separate files than the files of previous experiments
cell_path = "/home/ookur/newcells/L5_TTPC2_cADpyr232_2/"  ## (e.g. /home/ookur/L5_TTPC2_cADpyr232_2/) : ")

## Obtains the cell information (cell name, copy number) from the path
split_path = cell_path.split("/")  
while("" in split_path) :
    split_path.remove("")
split_name = split_path[-1].split("_")
cell_type = split_name[1]
cell_copy = split_name[3]
cell_info = [cell_type, cell_copy]

print("Simulation is starting... It will be completed approximately in an hour.")
print(" Go have some time with your friends and family, computer will be working hard for you... Beep Beep!")

## Start the simulation and return the necessary outputs with the function 
popt, pcov, r_sq, paramerr, freqs = extraction_EIF_faster(cell_info, alfa_a, trial_no)

print("Frequencies for 8 Trials: " + str(freqs) + " EIF parameters are " +str(popt)+ " so, delta_T you were looking for is "+ str(popt[3]))
print("Error of fit is " + str((r_sq)) + ". Error of the parameters are " + str(paramerr))
print("To find the best fit, sometimes manually changing the range for the F(V) function is necessary.")
print("If you don't like the result you can use Recalculate_Parameters.py to manually change the range and calculate the parameters within seconds.")
print("Thanks for using Oyku Okur's EIF Extraction Kit. See you in the next simulation!")



