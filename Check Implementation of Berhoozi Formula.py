# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 14:32:49 2021

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
from PyGalaxies import Halos



halo_masses = np.arange(10**10, 10**15, 10**11)

redshifts = [0.1,1,2,3,4,5,6,7,8]

stellar_masses = np.empty((len(halo_masses), len(redshifts)))





for z_counter, z in enumerate(redshifts):
    
    a = 1 / (1 + z)
    
    for mass_counter, halo_mass  in enumerate(halo_masses):
    
        stellar_masses[mass_counter,z_counter] = Halos.HaloProperties.Behroozi_formula(a, z, halo_mass)/halo_mass
        
        




fig, ax1 = plt.subplots(figsize=(6,5))

for counter in range(0,len(redshifts)):
    
    ax1.plot(halo_masses, stellar_masses[:,counter],lw=3,label = 'z = {:.1f}'.format(redshifts[counter]))
    
    
    
ax1.set_ylim(0.0001,0.08)
ax1.set_xlim(10**10,10**15.5)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Stellar Mass [Mo]/Halo Mass [Mo]')
ax1.set_xlabel('Halo Mass [Mo]')
ax1.legend(loc='upper right')
plt.savefig('Berhoozi Formula check frac.png', dpi = 400)
plt.show()

