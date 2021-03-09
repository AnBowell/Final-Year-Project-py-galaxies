# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:28:25 2021

@author: Andrew
"""


import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies')) # This is needed to import C routines.
import Parameters, Monitor, HDF_graph_properties, Halos
import numpy as np
import matplotlib.pyplot as plt

model_param_filepath='../Input_Params/input_params.yml'
debug_flag = False
verbosity = 1
time_code = True
model_params = Parameters.ModelParams(model_param_filepath,verbosity,debug_flag)
model_params.reionize_model = 0
model_params.load_paramters_to_C()
model_params.read_in_data_tables_from_c()


if model_params.reionize_model == 0:
    title = 'Okamoto et al (2008)'
    
    masses = np.arange(10**8, 10**10, 10**6)
    redshifts = np.arange(0, 10, 0.01)
    
    extent = [10**8, 10**10,0,10]
    
else:
    title = 'Kravtsov et al (2004)'
    
    masses = np.arange(10**8, 10**12, 10**8)
    redshifts = np.arange(0, 10, 0.01)
    
    extent = [10 ** 8, 10 ** 12,0,10]


results_array = np.empty((len(masses),len(redshifts)))


for x_counter,mass in enumerate(masses):
    for y_counter, z in enumerate(redshifts):
        
        
        results_array[x_counter, y_counter] = Halos.HaloProperties.calculate_reionization_frac(
            mass,z)
        
        
    if x_counter % 100 == 0:
        print(x_counter, 'is done')
        
        


fig, ax1 = plt.subplots(figsize = (8,5))
img = ax1.imshow(results_array, extent= extent,aspect ='auto',
                 vmin=0,vmax=1 )
ax1.set_xscale('log')
ax1.set_ylabel('Redshift z')
ax1.set_xlabel('Halo Mass / Mo')
ax1.set_title('Reionization Fraction - {} method'.format(title))
fig.colorbar(img)
plt.savefig('output graphs/Reionization Fraction {}.png'.format(title), dpi = 600)
plt.show()