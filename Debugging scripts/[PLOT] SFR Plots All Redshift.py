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
from scipy import stats
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
plt.rc('legend', fontsize=10)



hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['subhalo_data']

SFR =  dataset['SFR'][:]

zero_mask = SFR > 1

dm_mass = dataset['DM_mass'][:]
redshifts = dataset['redshift'][:]
star_mass = dataset['C_stellar_mass'][:]



bins = np.logspace(np.log10(40e7), np.log10(10e14), 20)

redshift_bins = [2.5,5,7.5,10]


binned_stats = [stats.binned_statistic(star_mass[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
              (redshifts[zero_mask] < redshift_bin))], np.log10(SFR[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
              (redshifts[zero_mask] < redshift_bin))]), 'mean', bins) for redshift_bin in redshift_bins]



bins_DM = np.logspace(np.log10(15e9), np.log10(12e14), 20)

binned_stats_DM = [stats.binned_statistic(dm_mass[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
              (redshifts[zero_mask] < redshift_bin))], np.log10(SFR[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
              (redshifts[zero_mask] < redshift_bin))]), 'mean', bins_DM) for redshift_bin in redshift_bins]

# means = [[np.mean(star_mass[zero_mask][ind][((redshifts[zero_mask]>= redshift_bin - 2) & 
#              (redshifts[zero_mask] < redshift_bin))]) for ind in np.unique(inds)] for redshift_bin in redshift_bins]


norm = matplotlib.colors.Normalize(redshifts.min(), redshifts.max())



dm_masses_per_regime = [dm_mass[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
                (redshifts[zero_mask] < redshift_bin))] for redshift_bin in redshift_bins]


star_mass_per_regime = [star_mass[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
                (redshifts[zero_mask] < redshift_bin))] for redshift_bin in redshift_bins]


log10_sfr_per_regime = [np.log10(SFR[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
                             (redshifts[zero_mask] < redshift_bin))]) for redshift_bin in redshift_bins]


redshifts_per_regime = [redshifts[zero_mask][((redshifts[zero_mask]>= redshift_bin - 2) & 
                        (redshifts[zero_mask] < redshift_bin))] for redshift_bin in redshift_bins]

# plt.scatter(dm_masses_per_regime[0],log10_sfr_per_regime[0])

# plt.scatter(dm_masses_per_regime[1],log10_sfr_per_regime[1],color='red')
# plt.xscale('log')
# plt.show()


order = np.argsort(redshifts[zero_mask])
img = plt.scatter(star_mass[zero_mask][order], np.log10(SFR[zero_mask][order]),
                  c = redshifts[zero_mask][order], s=4,edgecolor='none')
plt.show()
fig, axes = plt.subplots(4,2,figsize=(10,20),sharey=True,sharex= True,gridspec_kw={'wspace':0.035,
                                                                                   'hspace':0.035})
colors = img.cmap(norm(redshift_bins))

# colors = ['blue','green','orange','red']

for counter, redshift in enumerate(reversed(redshift_bins)):
    
    ax1 = axes[counter][0]
    ax2 = axes[counter][1]
    
    ax1.set_ylim(-0.05,2.5)
    ax1.set_ylim(-0.05,2.5)
    
    order = np.argsort(redshifts_per_regime[counter])
    
    ax1.scatter(star_mass_per_regime[counter][order], log10_sfr_per_regime[counter][order],
                  c = redshifts_per_regime[counter][order], s=4,edgecolor='none')
    
    
    bin_centers = np.sqrt(binned_stats[counter][1][1:] * binned_stats[counter][1][:-1])                      
    ax1.plot(bin_centers,binned_stats[counter][0],linewidth=4,color='black', solid_capstyle='round')
    ax1.plot(bin_centers,binned_stats[counter][0],linewidth=3,
              label = '{} $\leq$ z $<$ {}'.format(redshift-2.5,redshift), solid_capstyle='round',
              color=colors[counter])
    ax1.set_xscale('log')


    ax2.scatter(dm_masses_per_regime[counter][order], log10_sfr_per_regime[counter][order],
                      c = redshifts_per_regime[counter][order], s=4,edgecolor='none')
        
    bin_centers_DM = np.sqrt(binned_stats_DM[counter][1][1:] * binned_stats_DM[counter][1][:-1])
    ax2.plot(bin_centers_DM,binned_stats_DM[counter][0],linewidth=4,color='black', solid_capstyle='round')
    ax2.plot(bin_centers_DM,binned_stats_DM[counter][0],linewidth=3,
              label = '{} $\leq$ z $<$ {}'.format(redshift-2.5,redshift), solid_capstyle='round', color=colors[counter])

    ax2.set_xscale('log')

# order = np.argsort(redshifts[zero_mask])

plt.show()









fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5),sharey=True,gridspec_kw={'wspace':0.035})

order = np.argsort(redshifts[zero_mask])
img = ax1.scatter(star_mass[zero_mask][order], np.log10(SFR[zero_mask][order]),
                  c = redshifts[zero_mask][order], s=4,edgecolor='none')

# colors = plt.cm.jet(redshift_bins)
colors = img.cmap(norm(redshift_bins))
for counter, redshift in enumerate(redshift_bins):#[::-1]:
    bin_centers = np.sqrt(binned_stats[counter][1][1:] * binned_stats[counter][1][:-1])                      
    ax1.plot(bin_centers,binned_stats[counter][0],linewidth=4,color='black', solid_capstyle='round')
    ax1.plot(bin_centers,binned_stats[counter][0],linewidth=3,
              label = '{} $\leq$ z $<$ {}'.format(redshift-2.5,redshift), solid_capstyle='round',
              color=colors[counter])
        

    
    bin_centers_DM = np.sqrt(binned_stats_DM[counter][1][1:] * binned_stats_DM[counter][1][:-1])
    ax2.plot(bin_centers_DM,binned_stats_DM[counter][0],linewidth=4,color='black', solid_capstyle='round')
    ax2.plot(bin_centers_DM,binned_stats_DM[counter][0],linewidth=3,
              label = '{} $\leq$ z $<$ {}'.format(redshift-2.5,redshift), solid_capstyle='round', color=colors[counter])

ax1.legend(loc='lower right',title='Mean')
ax2.legend(loc='lower right',title='Mean')


ax1.set_xscale('log')

ax1.set_ylabel('log$_{10}$(SFR / [M$_\odot$yr$^{-1}$])')
ax1.set_xlabel('M$_*$ / M$_\odot$') 

# ax1.set_yscale('log')
ax1.set_xlim(2*10**8,2*10**13)
ax2.set_xlim((10**10)-10**9,4*10**14)
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
# plt.savefig('output graphs/SFR.pdf',dpi=400,bbox_inches='tight')
# plt.savefig('output graphs/SFR.png',dpi=400,bbox_inches='tight')
plt.show()
