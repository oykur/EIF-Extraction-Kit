# Oyku-Okur-s-EIF-Extraction-Kit-Faster-Version
Faster version of this kit has the same principles and structure as before. With the improvements and new comments to make process faster, more simple, and easy to understand.
Speed has been significantly improved, thus kit named "faster".

This EIF Extraction Kit has 5 files that enables extracting parameters of EIF model (tau, V_resting, V_threshold, delta_t) from a BBP NEURON model.

1) Before starting the simulation, cell should be downloaded, extracted and compiled. NEURON should be downloaded. 
Follow the steps in https://neuron.yale.edu/neuron/download

2) Download the NEURON cell of your choice from https://bbp.epfl.ch/nmc-portal/downloads.html 
Remember this kit contains L5 - TTPC1, TTPC2, UTPC, and STPC cells in its dictionary. If you want to investigate another cell,
you should add it inside the dataframe "cells_df.pkl"

3) Extract the cell from the .zip file. You can see the possible compilation methods at https://www.neuron.yale.edu/neuron/download/compile_linux
on the left of the webpage Windows and Mac compilation methods can be seen. 
Compilation should be done inside the "mechanisms" file of the cell file.
Compiled .mod files should be taken outside the mechanisms file. Placed where the .py files are.

4) Files found in this kit should be placed inside the cell file, together with "run.py" or .mod files that are transferred in the previous step.
5) Open the Simulation_8Types_Faster.py. Insert the path of the cell, experiment number and alfa. Run the simulation. 
6) Wait for an hour for the simulation to finish. Have a cup of tea with your family and/or friends. 
7) Check the .txt output file 
   - Is the firing rate is between 1-15 Hz? If not change the alfa and experiment number, run the simulation again
   - Are the error for the fit and parameters smaller then 0.05? If not, use Recalculate_Parameters.py and calculate again by 
     changing the membrane potential interval. Find the best fit from with smallest error rate and save the graphs and outputs.

Files In This Kit:

*"cells_df.pkl" contains a dataframe where the cell names for the different types and copy numbers are stored. It is necessary for code to find specific
   cell code given to cell to start working with NEURON. e.g. Layer 5 (L5), Small Tufted Pyramidal Cell (STPC) Copy number 1 has the cell code > cADpyr232_L5_STPC_d16b0be14e

   It now has all 5 copies of L5 TTPC1, TTPC2, UTPC and STPC. You can add new cells and copies inside the dataframe by opening in Python. 

*OU_Process.py is needed to generate Ornstein Uhlenbeck processes for the input. 

   OU process is a noisy process where noise deviates from mean with standart deviation, but comes back to mean 
   with a time constant tau using the Brownian motion. 2 OU processes with mean 0 and standard deviations [[.18,.18],[.25,.36]] and tau = [3, 10 ms] 
   are used in this simulation following the work of Badel et al., 2008. 
   These processes are combinetd with 4 different DC currents [0, 0.02, 0.03, 0.06].
   So, 4 DC currents combined with 2 couples of standart deviations, 8 different OU process are given to cell. 

No changes will be made to this file.

*"EIF_Extraction_faster.py" has the function that extracts EIF parameters from the given neuron. 
   It follows the steps in the Badel et al., 2008. 
   Simulation injects generated OU process to NEURON cell for 40s after 300 ms transitory time. Injection is to the middle 
   of the soma and outputs obtained from the same spot. Output data is filtered (200 ms after each spike is discarded)
   and capacitance of cell is found with variance minimization method suggested in the article. 
   Transmembrane current (I_m), Dynamic I-V curve (I_d), and Function of Voltage F(V) is found. F(V) is fitted by 
   EIF model to obtain parameters. Function returns parameters for the EIF model, their error rates, error of the whole fit 
   and saves the output array of simulation for further ivestigation.

   One needs to specify the cell type, alfa (amplitude of input, for the firing frequency range of interest), and the trial number to use this function. 
   These parameters are changed inside the Simulation_8Types_Faster.py

*Cell type can be changed inside the "Simulation_8Types_Faster.py" file. Write the path of the cell by copying from cell folder directory.
   When directory is changed, code recognizes the name of the cell and copy number. Inputs the name of the cell into EIF_Extraction function.

   The directory should be changed manually inside the code in the format: home/user/other folders/layer_cell_type_copy/ e.g. home/ookur/newcells/L5_TTPC1_cADpyr232_2/
   It should direct the function inside the cell folder where other .py and .pkl files are.

  -Alfa is set to 1. 
   To increase or decrease the firing frequency it can be changed inside the "Simulation_8Types_spc.py" file.

  -Trial_no is set to 1. 
   To prevent overwriting or deletion of the output files, it is important to increase the trial number for each experiment inside "Simulation_8Types_spc.py"

*Outputs of the simulation can be observed in .txt file with the name of the cell, copy number and experiment number.

*I, V, dt, t, and I_m are stored as experiment array in .npy file with the name of the cell, copy number and experiment number.
   They can be used in Recalculate_Parameters.py 
   This will be useful since changing the range of the F(V) values can improve the accuracy of the fit. 
   If the performance of the fit is weak, then this code enables changing the range, seeing the graphs of interest (uncomment to visualize) 
   and obtaining new fits. Best fit can be used for further investigations.

Reference for steps: 
Badel L, Lefort S, Brette R, Petersen CC, Gerstner W, Richardson MJ. Dynamic I-V curves are reliable predictors of naturalistic pyramidal-neuron
voltage traces. J Neurophysiol. 2008 Feb;99(2):656-66. doi: 10.1152/jn.01107.2007. Epub2007 Dec 5. PMID: 18057107.

Enjoy your extraction!

Öykü Okur

For any questions contact > oyku.okur@ug.bilkent.edu.tr




