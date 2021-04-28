# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 14:32:49 2021

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import os,sys 
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies'))
from PyGalaxies import Halos
plt.style.use('classic')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')
plt.style.use('classic')
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)

halo_masses = np.arange(10**10, 10**15, 10**11)

# redshifts = [0,1,2,3,4,5,6,7,8]
redshifts = [0,2,4,6,8]


stellar_masses = np.empty((len(halo_masses), len(redshifts)))


for z_counter, z in enumerate(redshifts):
    
    a = 1 / (1 + z)
    
    for mass_counter, halo_mass  in enumerate(halo_masses):
    
        stellar_masses[mass_counter,z_counter] = Halos.HaloProperties.Behroozi_formula(a, z, halo_mass)
        
        




fig, (ax1,ax2) = plt.subplots(2,1,figsize=(8,6.25),sharex = True,gridspec_kw={'hspace':0.025})

for counter in range(0,len(redshifts)):
    
    ax1.plot(halo_masses, stellar_masses[:,counter],lw=3,label = 'z = {:.1f}'.format(redshifts[counter]))
    
    
for counter in range(0,len(redshifts)):
    
    ax2.plot(halo_masses, stellar_masses[:,counter]/halo_masses,lw=3,label = 'z = {:.1f}'.format(redshifts[counter]))
    
    
    


ax1.set_xlim(10**10,10**15.5)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('M$_*$/M$_\odot$')
# ax1.legend(loc='lower right')
ax1.grid(True,linestyle='--')
ax2.set_ylim(0.00004,0.034)
ax1.set_ylim(10**6.5,10**11.5)
ax2.set_xlim(10**10,10**15)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel('M$_*$/M$_\mathrm{DM}$')
ax2.set_xlabel('M$_\mathrm{DM}$/M$_\odot$')
ax2.grid(True,linestyle='--')
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.9,0.624))#loc='center right')
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
# ax2.legend(loc='upper right')
plt.tight_layout()
plt.savefig('output graphs/Berhoozi.png', dpi = 400)
plt.savefig('output graphs/Berhoozi.eps')
plt.show()

