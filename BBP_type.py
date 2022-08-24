#!/usr/bin/env python
# coding: utf-8


# Neuron is imported, specifically its main module


import os
import subprocess

def BBP_type(path):
    os.chdir(path)

    command = 'grep \"cell = new\" ./$path/createsimulation.hoc | awk \'{print $4}\' | sed \'s/synapses_enabled/0/g\''
    cell_name = subprocess.run(command,shell= True, stdout=subprocess.PIPE).stdout
    cell_name = cell_name.decode().split(sep = "\r")
    cell_name = cell_name[0]
###############################################################################
    from neuron import h
    print("I'm in the folder " + path)
# Standard NEURON components are loaded

    h.load_file("stdrun.hoc") #funzione run
  
    h.load_file("import3d.hoc")#3D morphology

# Specific files, related to the BBP cell, are loaded
    h.load_file("constants.hoc")
    
    h.load_file("morphology.hoc")

    h.load_file("template.hoc")
    
    h.load_file("biophysics.hoc")

#####################################################################################
    cell_name = "h." +cell_name
#     print(cell_name)
    cell = eval(cell_name) ## specify the cell here,# can be found from run.py file of each cell
    #cell= h.cADpyr232_L5_TTPC2_a28017c6c7(0)
    return cell
# The cell is "created" into NEURON-PYTHON, it is a model cell from BBP

