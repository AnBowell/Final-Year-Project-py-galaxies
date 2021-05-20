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
import h5py as h5

save_plot = True

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# plt.rc('axes', labelsize=18)  
# plt.rc('text', usetex=True) 
# plt.rc('xtick',labelsize=14)
# plt.rc('ytick',labelsize=14)
plt.style.use('classic')

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)

PART_MASS = 1.3301059053364127E9

hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['halo_data']


DM_mass = dataset['mass_DM'][:]

Baryon_mass =  dataset['mass_baryon'][:]



bins = np.logspace(np.log10(1e10), np.log10(15e14), 300)

inds = np.digitize(DM_mass, bins)



bin_means = [(Baryon_mass[bin_idx==inds]/DM_mass[bin_idx == inds]).mean() for bin_idx in np.unique(inds)]
bin_medians = [np.median(Baryon_mass[bin_idx==inds]/DM_mass[bin_idx == inds]) for bin_idx in np.unique(inds)]







print(np.where((Baryon_mass/DM_mass) > 1))

fig, ax1 = plt.subplots(figsize = (7,4))
# fig, (ax1,ax2) = plt.subplots(2,1,figsize = (7,7))
ax1.scatter(DM_mass,Baryon_mass/DM_mass,s=1.5,color='navy')

ax1.plot(bins[np.unique(inds)], bin_means,color='red',zorder=1,lw=2, label='Mean')

ax1.plot(bins[np.unique(inds)], bin_medians,color='orange',zorder=1,lw=2,label='Median')

ax1.hlines(y=0.155,xmin=DM_mass.min(),xmax =DM_mass.max(),zorder=1,lw=2,color='lime',ls='--',label='Cosmic Mean')
ax1.set_yscale('log')
ax1.set_xlim(10**10,6*10**14)
ax1.set_xscale('log')
ax1.set_ylabel(r'M$_\mathrm{bary}$/M$_\mathrm{DM}$')
ax1.set_xlabel(r'M$_\mathrm{DM}$/M$_\odot$')
ax1.grid(True, linestyle ='--')
ax1.set_axisbelow(True)


ax1.legend()
plt.tight_layout()
plt.savefig('output graphs/Bary_frac.png',dpi=400)
# plt.savefig('output graphs/Bary_frac.eps',dpi=400, rasterized=True)
plt.show()

# ax2.scatter(DM_mass,Baryon_mass/DM_mass,zorder=-1,s=0.5)

# ax2.plot(bins[np.unique(inds)], bin_means,color='red',zorder=1,lw=2, label='Mean')

# ax2.plot(bins[np.unique(inds)], bin_medians,color='orange',zorder=1,lw=2,label='Median')




# ax2.hlines(y=0.155,xmin=DM_mass.min(),xmax =DM_mass.max(),zorder=1,color='black',ls='--')
# ax2.set_ylim(0.15,0.204)
# ax2.set_yscale('log')
# ax2.set_xscale('log')
# ax2.set_ylabel('Total Baryon Mass / Halo Mass')
# ax2.set_xlabel(r'Halo Mass/M$_\odot$')



# ax2.legend(title='Zoomed In')
#plt.savefig('output graphs/Reionization Fraction {}.png'.format(title), dpi = 600)
# plt.tight_layout()
# # plt.savefig('output graphs/BaryonFraction.eps')
# if save_plot:
#     plt.savefig('output graphs/BaryonFraction.png',dpi=400)
# plt.show()







