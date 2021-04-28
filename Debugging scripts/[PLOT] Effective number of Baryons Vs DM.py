# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 10:31:11 2021

@author: Andrew
"""
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

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('axes', labelsize=18)  

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

# fig, ax1 = plt.subplots(figsize = (5,3.5))
# fig, ax1 = plt.subplots(figsize = (7,4))
# ax1.scatter(DM_mass/PART_MASS,Baryon_mass/0.155/PART_MASS,zorder=-1,s=0.5)
# ax1.set_yscale('log')
# ax1.set_xscale('log')

# ax1.set_ylabel('Baryon_mass/BaryonFrac/PART_MASS')
# ax1.set_xlabel('DM_mass/PART_MASS')
# plt.savefig('output graphs/Effectiveparticlenumbers.png',dpi=400)
# plt.show()


bins = np.logspace(np.log10(1e10), np.log10(15e14), 300)

inds = np.digitize(DM_mass, bins)



bin_means = [(Baryon_mass[bin_idx==inds]/0.155/PART_MASS).mean() for bin_idx in np.unique(inds)]
bin_medians = [np.median(Baryon_mass[bin_idx==inds]/0.155/PART_MASS) for bin_idx in np.unique(inds)]




fig, (ax1,ax2) = plt.subplots(2,1,figsize = (7,7))
ax1.scatter(DM_mass/PART_MASS,Baryon_mass/0.155/PART_MASS,zorder=-1,s=0.4)
ax1.plot(bins[np.unique(inds)]/PART_MASS, bin_means,color='red',zorder=1,lw=2, label='Mean')

ax1.plot(bins[np.unique(inds)]/PART_MASS, bin_medians,color='orange',zorder=1,lw=2,label='Median')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Baryon_mass/BaryonFrac/PART_MASS')
ax1.set_xlabel('DM_mass/PART_MASS')
ax1.legend(loc='lower right',title='Normal')

ax2.scatter(DM_mass/PART_MASS,Baryon_mass/0.155/PART_MASS,zorder=-1,s=0.4)
ax2.set_yscale('log')
ax2.set_ylim(10e0,10e1)
ax2.set_xlim(10e0,10e1)
ax2.set_xscale('log')
ax2.set_ylabel('Baryon_mass/BaryonFrac/PART_MASS')
ax2.set_xlabel('DM_mass/PART_MASS')

ax2.plot(bins[np.unique(inds)]/PART_MASS, bin_means,color='red',zorder=1,lw=2, label='Mean')

ax2.plot(bins[np.unique(inds)]/PART_MASS, bin_medians,color='orange',zorder=1,lw=2,label='Median')
ax2.legend(loc='lower right',title='Zoomed')
plt.savefig('output graphs/Effectiveparticlenumbers_with_zoom.png',dpi=400)
plt.show()


# ax1.plot(bins[np.unique(inds)], bin_means,color='red',zorder=1,lw=2, label='Mean')

# ax1.plot(bins[np.unique(inds)], bin_medians,color='orange',zorder=1,lw=2,label='Median')

# ax1.hlines(y=0.155,xmin=DM_mass.min(),xmax =DM_mass.max(),zorder=1,color='black',ls='--')
# ax1.set_yscale('log')
# ax1.set_xscale('log')
# ax1.set_ylabel('Total Baryon Mass / Halo Mass')
# ax1.set_xlabel(r'Halo Mass/M$_\odot$')
# ax1.legend(title='Normal')



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
# #plt.savefig('output graphs/Reionization Fraction {}.png'.format(title), dpi = 600)
# plt.tight_layout()
# # plt.savefig('output graphs/BaryonFraction.eps')
# if save_plot:
#     plt.savefig('output graphs/BaryonFraction.png',dpi=400)
# plt.show()







