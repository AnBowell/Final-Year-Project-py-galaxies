# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 16:17:49 2021

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
import matplotlib
save_plot = True
plt.rcParams['agg.path.chunksize'] = 10000
plt.style.use('classic')

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=14)


hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['subhalo_data']

ber_stellar_mass = dataset['stellar_mass'][:]
c_stellar_mass = dataset['C_stellar_mass'][:]

redshifts = dataset['redshift'][:]


stellar_mask = ((ber_stellar_mass > 0) & (c_stellar_mass>0))




# ber_stellar_mass = ber_stellar_mass[stellar_mask]
# c_stellar_mass = c_stellar_mass[stellar_mask]



# bins = [10e6,10e7,10e8,10e9,10e10,10e11]

bins = np.logspace(np.log10(30e6), np.log10(10e11), 15)

redshift_bins = [2,4,6,8,10]

binned_stats = [stats.binned_statistic(ber_stellar_mass[stellar_mask][((redshifts[stellar_mask]>= redshift_bin - 2) & 
              (redshifts[stellar_mask] < redshift_bin))], c_stellar_mass[stellar_mask][((redshifts[stellar_mask]>= redshift_bin - 2) & 
              (redshifts[stellar_mask] < redshift_bin))], 'mean', bins) for redshift_bin in redshift_bins]



fig, ax1 = plt.subplots(figsize=(8,5))
ax1.grid(True, linestyle ='--')
ax1.set_axisbelow(True)
order = np.argsort(redshifts[stellar_mask])
scat = ax1.scatter(ber_stellar_mass[stellar_mask][order], c_stellar_mass[stellar_mask][order],s=4,zorder=-1,edgecolor='none', c=redshifts[stellar_mask][order])

# ax1.plot(ber_stellar_mass,ber_stellar_mass,color='black',zorder=3,ls='--',lw=1.56,label='1-1')

lims = [
    np.min([ax1.get_xlim(), ax1.get_ylim()]),  # min of both axes
    np.max([ax1.get_xlim(), ax1.get_ylim()]),  # max of both axes
]

ax1.plot(lims,lims,color='saddlebrown',zorder=3,ls='--',lw=3,label='$y=x$')

norm = matplotlib.colors.Normalize(redshifts.min(), redshifts.max())
colors = scat.cmap(norm(redshift_bins))


for counter, redshift in enumerate(redshift_bins):
    bin_centers = np.sqrt(binned_stats[counter][1][1:] * binned_stats[counter][1][:-1])                      
    ax1.plot(bin_centers,binned_stats[counter][0],linewidth=4,color='black', solid_capstyle='round')
    ax1.plot(bin_centers,binned_stats[counter][0],linewidth=3,
             label = '{} $\leq$ z $<$ {}'.format(redshift-2,redshift), solid_capstyle='round',
             color=colors[counter])


# divider = make_axes_locatable(ax1)
# cax = divider.append_axes('right', size='5%', pad=0.05)
# fig.colorbar(scat, cax=cax, orientation='vertical')
# cax.set_ylabel('Redshift z')
cbar = fig.colorbar(scat)
cbar.set_label('z',rotation = 90, fontsize=16)
ax1.set_xlim(10.5**7,(10**12) - 7*10**11)
ax1.set_ylim(10**6,10**13)
ax1.set_ylabel('Modelled M$_*$ / M$_\odot$')
ax1.set_xlabel('Behroozi M$_*$ / M$_\odot$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend(loc='upper left',ncol=2)

plt.tight_layout()
plt.savefig('output graphs/C_Beh.png', dpi= 400)
# plt.show()

