from EIF_Extraction_20_12 import extraction_EIF
import sys
print("Before starting, please check if you downloaded, extracted and compiled the NEURON cell. Then check if the EIF_Extraction, Simulation_8Types, and OUProcess_Simple codes, also the cells_df.pkl are all in the same file.")
alfa = 1
trial_no = 1
cell_path = sys.argv[1] ## (e.g. /home/ookur/L5_TTPC2_cADpyr232_2/) : ")

split_path = cell_path.split("/")
while("" in split_path) :
    split_path.remove("")
split_name = split_path[-1].split("_")
cell_type = split_name[1]
cell_copy = split_name[3]
cell_info = [cell_type, cell_copy]
print("To arrange the firing frequency, amplitude (alfa) to magnify the noise is set to 1, you can change it from EIF_Extraction.py")

print("Simulation is starting... If you didn't write anything wrong, simulation will be completed approximately in 4 hours.")
print(" Go have some time with your friends and family, computer will be working hard for you... Beep Beep!")

popt, pcov, r_sq, paramerr, freqs = extraction_EIF(cell_info, cell_path, alfa, trial_no)

print("Frequencies for 8 Trials: " + str(freqs) + " EIF parameters are " +str(popt)+ " so, delta_T you were looking for is "+ str(popt[3]))
print("Error of fit is " + str((1-r_sq)) + ". Error of the parameters are " + str(paramerr))
print("To find the best fit, sometimes manually changing the range for the F(V) function is necessary.")
print("If you don't like the result you can use Recalculate_Parameters_With_Outputs.py to manually change the range and calculate the parameters within seconds.")
print("Thanks for using Oykur's EIF Extraction Kit. See you in the next simulation!")



