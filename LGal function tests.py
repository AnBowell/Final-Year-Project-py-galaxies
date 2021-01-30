# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:20:49 2020

@author: Andrew

"""



import sys

#from PyGalaxies import Parameters

import PyGalaxies.L_Galaxies


# Init cooling functions. (Reads data from file)


# model_param_filepath='Input_Params/input_params.yml'
# verbosity = 1

# debug_flag = 0

# model_params = Parameters.ModelParams(model_param_filepath,verbosity,debug_flag)
# model_params.output_params()

# #L_Galaxies.C_read_cooling_functions() # Reset kernel everytime this is called += in code.

for i in range(-10,-9):

    temp = 100 * i
    log_z = i
    
    print('\n',L_Galaxies.C_get_metaldependent_cooling_rate(temp, log_z))