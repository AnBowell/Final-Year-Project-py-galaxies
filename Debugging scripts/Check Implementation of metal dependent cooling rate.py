# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:22:26 2021

@author: Andrew
"""


import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies')) # This is needed to import C routines.
import Parameters, Monitor, HDF_graph_properties, Halos
import numpy as np
import matplotlib.pyplot as plt
from LGalaxies.L_Galaxies import C_get_metaldependent_cooling_rate


model_param_filepath='../Input_Params/input_params.yml'
debug_flag = False
verbosity = 1 
time_code = True
model_params = Parameters.ModelParams(model_param_filepath,verbosity,debug_flag)
model_params.load_paramters_to_C()
model_params.read_in_data_tables_from_c()


H0 = model_params.H0

fig,ax1 = plt.subplots(figsize= (6,5))



for z in range(0,8, 2):

    # r0 in Mpc
    r0 = np.arange(0.01,10,0.001)
    
    circular_velocity = Halos.HaloProperties.calculate_virial_circular_velocity(z, H0, r0)
    
    virial_temps =  Halos.HaloProperties.calculate_virial_temperature(circular_velocity,
                                                                     model_params.mu,
                                                                     model_params.m_p,
                                                                     model_params.k_B)
    
    
    
    ax1.plot(r0,virial_temps,label='z = {}'.format(z))
ax1.legend()
# ax1.set_yscale('log')
# ax1.set_xscale('log')
ax1.set_xlabel('Halo radius / Mpc')
ax1.set_ylabel('Hot gas temperature / K')
plt.show()






temps = np.arange(1e4, 2e8, 1e3)



fig, ax1 = plt.subplots(figsize=(6,5))


metalicities = [0,0.012,0.001,0.00001,0.5]

for metalicity in metalicities:
    log_Z = np.log10(metalicity)
    metal_dependent_cooling_storage = np.empty(len(temps))
    for counter, log_temp in enumerate(np.log10(temps)):
        
        
    
        metal_dependent_cooling_storage[counter] = C_get_metaldependent_cooling_rate(
                                                                    log_temp,
                                                                    log_Z)
    
    ax1.plot(temps, metal_dependent_cooling_storage,label='Z= {}'.format(metalicity))
        
    
    
    
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('Cooling rate  / s^-1 erg cm^3')
ax1.set_xlabel('Gas temperature / K')
ax1.legend()
plt.savefig('output graphs/Cooling rates.png',dpi=600)
plt.show()