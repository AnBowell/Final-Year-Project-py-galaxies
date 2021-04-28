import numpy as np
import matplotlib.pyplot as plt

import h5py as h5

import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies')) # This is needed to import C routines.
import Parameters, Monitor, HDF_graph_properties, Halos


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
# plt.style.use('classic')
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)
plt.rcParams['legend.title_fontsize'] = 'x-large'


hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['subhalo_data']

SFR =  dataset['SFR'][:]

zero_mask = SFR > 1

dm_mass = dataset['DM_mass'][:]
redshifts = dataset['redshift'][:]
star_mass = dataset['C_stellar_mass'][:]


zero_mask = star_mass > 0

fig, ax1 = plt.subplots(figsize=(7,4))

order = np.argsort(redshifts[zero_mask])
img = ax1.scatter(dm_mass[zero_mask][order],
            star_mass[zero_mask][order],
            c =redshifts[zero_mask][order],
            s=1)

ax1.set_yscale('log')
ax1.set_xscale('log')


halo_masses = np.arange(10**10, 10**15, 10**11)

# redshifts = [0,1,2,3,4,5,6,7,8]
redshifts = [0,2,4,6,8]


stellar_masses = np.empty((len(halo_masses), len(redshifts)))


for z_counter, z in enumerate(redshifts):
    
    a = 1 / (1 + z)
    
    for mass_counter, halo_mass  in enumerate(halo_masses):
    
        stellar_masses[mass_counter,z_counter] = Halos.HaloProperties.Behroozi_formula(a, z, halo_mass)
        
        

for counter in range(0,len(redshifts)):
    
    ax1.plot(halo_masses, stellar_masses[:,counter],lw=3,label = 'z = {:.1f}'.format(redshifts[counter]))
    
    
cbar = fig.colorbar(img)
cbar.set_label('z',rotation = 90, fontsize=16)
ax1.set_xlabel('M$_\mathrm{DM}$')
ax1.set_ylabel('M$_*$')
ax1.grid(True,linestyle='--')
ax1.set_axisbelow(True)
plt.tight_layout()
plt.savefig('output graphs/Stellar Mass Vs Behroozi.png',dpi=400)
plt.show()