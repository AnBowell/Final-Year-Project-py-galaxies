# -*- coding: utf-8 -*-
""" Created on Wed Feb 10 09:28:25 2021

@author: Andrew
"""

import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies')) # This is needed to import C routines.
import Parameters, Monitor, HDF_graph_properties, Halos
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

save_plot = True

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# plt.rc('axes', labelsize=18)  
# plt.rc('text', usetex=True) 
# plt.rc('xtick',labelsize=14)
# # plt.rc('ytick',labelsize=14)

plt.style.use('classic')

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)



hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['subhalo_data']

SFR =  dataset['SFR'][:]

zero_mask = SFR > 1

dm_mass = dataset['DM_mass'][:]
redshifts = dataset['redshift'][:]
star_mass = dataset['C_stellar_mass'][:]



fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5),sharey=True,gridspec_kw={'wspace':0.025})

order = np.argsort(redshifts[zero_mask])
img = ax1.scatter(star_mass[zero_mask][order], np.log10(SFR[zero_mask][order]),
                  c = redshifts[zero_mask][order], s=4,edgecolor='none')
ax1.set_xscale('log')

ax1.set_ylabel('log$_{10}$(SFR / [M$_\odot$yr$^{-1}$])')
ax1.set_xlabel('M$_*$ / M$_\odot$') 

# ax1.set_yscale('log')
ax1.set_xlim(2*10**8,2*10**13)
ax2.set_xlim(10**10,4*10**14)
ax1.set_ylim(-0.05,2.5)
ax2.set_ylim(-0.05,2.5)


img = ax2.scatter(dm_mass[zero_mask][order], np.log10(SFR[zero_mask][order]),
                  c = redshifts[zero_mask][order], s=4,edgecolor='none')
ax2.set_xscale('log')
# ax1.set_yscale('log')
# ax1.set_xlim(-5,3.5*10**14)
# ax1.set_ylim(0,275)
ax2.set_xlabel('M$_\mathrm{DM}$ / M$_\odot$')




fig.subplots_adjust(top=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar_ax = fig.add_axes([0.1, 0.9, 0.8, 0.05])
cbar = fig.colorbar(img, cax=cbar_ax,orientation="horizontal")

cbar.set_label('z',fontsize=18,labelpad=-54)




plt.tight_layout()
plt.savefig('output graphs/SFR.png',dpi=400)
plt.show()
