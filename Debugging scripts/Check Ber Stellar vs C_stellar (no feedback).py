# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 16:17:49 2021

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable


hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['subhalo_data']

ber_stellar_mass = dataset['stellar_mass'][:]
c_stellar_mass = dataset['C_stellar_mass'][:]

redshift = dataset['redshift'][:]


stellar_mask = ((ber_stellar_mass > 0) & (c_stellar_mass>0))



redshift = redshift[stellar_mask]
ber_stellar_mass = ber_stellar_mass[stellar_mask]
c_stellar_mass = c_stellar_mass[stellar_mask]



# bins = [10e6,10e7,10e8,10e9,10e10,10e11]

bins = np.logspace(np.log10(10e6), np.log10(10e11), 20)



inds = np.digitize(ber_stellar_mass, bins)

bin_centers = []

plot_array = np.empty((6,len(np.unique(inds))))


for bin_number in np.unique(inds):
    bin_number = bin_number-1
    
    ind_ber_mass = ber_stellar_mass[inds == bin_number]
        
    ind_c_mass = c_stellar_mass[inds == bin_number]
        
    redshift_this_slice = redshift[inds==bin_number]
    
    bin_centers.append((bins[bin_number] + bins[bin_number+1])/2)
        
    for i in range(0,12,2):
        
        
        redshift_mask = ((redshift_this_slice >= i) & (redshift_this_slice < i+1))    
        i = int(i/2)
        median_for_slice = np.median(ind_c_mass[redshift_mask])
        print(median_for_slice)
        plot_array[i,bin_number] = median_for_slice
    




fig, ax1 = plt.subplots(figsize=(8,5))
scat = ax1.scatter(ber_stellar_mass, c_stellar_mass,s=1.5,zorder=-1,alpha=0.5)#, c=redshift,s=1.5)
ax1.plot(ber_stellar_mass,ber_stellar_mass,color='black',zorder=3,ls='--',lw=1.56,label='1-1')
colors = ['red','green','brown','orange','purple']

for i in range(0,5):
    ax1.plot(bin_centers,plot_array[i,:],color=colors[i],label='{} <= z < {}'.format(i*2,(i*2)+2),
             ls='-',marker='x')
            


# divider = make_axes_locatable(ax1)
# cax = divider.append_axes('right', size='5%', pad=0.05)
# fig.colorbar(scat, cax=cax, orientation='vertical')
# cax.set_ylabel('Redshift z')

ax1.set_ylabel('C stellar mass / M$_o$')
ax1.set_xlabel('Behroozi stellar mass / M $_o$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('Modelled Stellar Mass vs Behroozi Formula')
plt.legend()
plt.savefig('output graphs/C_Stellar vs Ber Stellar.png', dpi= 400)
plt.show()

print()
