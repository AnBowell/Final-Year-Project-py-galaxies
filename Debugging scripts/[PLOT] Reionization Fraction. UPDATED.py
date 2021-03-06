# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:28:25 2021

@author: Andrew
"""


import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies')) # This is needed to import C routines.
import Parameters, Monitor, HDF_graph_properties, Halos
import numpy as np
import matplotlib.pyplot as plt


# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# plt.rc('axes', labelsize=18)  
# plt.rc('text', usetex=True) 
# plt.rc('xtick',labelsize=14)
# plt.rc('ytick',labelsize=14)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.style.use('classic')
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)
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
    
    # masses = np.arange(10**8, 10**11, 10**8)
    
    masses = np.logspace(np.log10(10**9), np.log10(10**11),1000)
    
    
    redshifts = np.arange(0, 10, 0.01)
    
    
    
    
    extent = [10**9, 10**11,0,10]
    
    
    
    
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
        
        


fig, ax1 = plt.subplots(figsize = (7,4))
img = ax1.imshow(results_array, extent= extent,aspect ='auto',
                 vmin=0,vmax=1,cmap='inferno')
ax1.set_xlim(10**9,10**11)

ax1.set_xscale('log')
ax1.set_ylabel('z')
ax1.set_xlabel(r'Halo Mass/M$_\odot$')
cbar = fig.colorbar(img)
cbar.ax.get_yaxis().labelpad = 17
cbar.set_label('Baryon Fraction Multiplier', rotation=270,fontsize=16)
#plt.savefig('output graphs/Reionization Fraction {}.png'.format(title), dpi = 600)
plt.tight_layout()
plt.savefig('output graphs/ReionizationFraction.eps'.format(title))
plt.savefig('output graphs/ReionizationFraction.png'.format(title),dpi=400)
plt.show()