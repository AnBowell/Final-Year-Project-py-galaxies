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

plt.style.use('classic')


plt.rc('font', family='serif')

plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 

model_param_filepath='../Input_Params/input_params.yml'
debug_flag = False
verbosity = 1 
time_code = True
model_params = Parameters.ModelParams(model_param_filepath,verbosity,debug_flag)
model_params.load_paramters_to_C()
model_params.read_in_data_tables_from_c()


H0 = model_params.H0




temps = np.arange(1e4, 2e8, 1e3)

fig, ax1 = plt.subplots(figsize=(7,4))
ax1.grid(True,linestyle='--')
ax1.set_axisbelow(True)
metalicities = [0,0.001,0.012,0.1]

metalicities.reverse()

for metalicity in metalicities:
    log_Z = np.log10(metalicity)
    metal_dependent_cooling_storage = np.empty(len(temps))
    for counter, log_temp in enumerate(np.log10(temps)):
        
        
    
        metal_dependent_cooling_storage[counter] = C_get_metaldependent_cooling_rate(
                                                                    log_temp,
                                                                    log_Z)

    if metalicity == 0.012:
        label = 'Z = Z$_\odot$ = 0.012'
    else:
        label = 'Z = {}'.format(metalicity)
    
    ax1.plot(temps, metal_dependent_cooling_storage,label=label,
             lw = 3.5,zorder=2)
        

# ax1.grid(True,linestyle='--',zorder=0)
# ax1.set_axisbelow(True)
ax1.set_xlim(8.*10**3,1.5*10**8)
ax1.set_xscale('log')
ax1.tick_params(bottom=True, top=True, left=True, right=True)
ax1.set_yscale('log')
ax1.set_ylabel('Cooling Rate ($\Lambda$) / erg cm$^3$s$^{-1}$')#,fontsize=16)
ax1.set_xlabel('Gas Temperature (T) / K')#,fontsize=16)
ax1.legend()
plt.tight_layout()
plt.savefig('output graphs/Cooling_rates.eps')
plt.savefig('output graphs/Cooling_rates.png',dpi=400)
plt.show()