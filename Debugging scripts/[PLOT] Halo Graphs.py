# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:57:04 2021

@author: Andrew
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:04:19 2021

@author: Andrew
"""

import numpy as np
import matplotlib.pyplot as plt

import h5py as h5


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)

hdf_filepath = '../../Input Data/mega_graph_mini_09-11.hdf5'

the_file = h5.File(hdf_filepath,'r')


header_data = the_file['Header']
part_mass = header_data.attrs['part_mass']
list_of_nparts = []



# nprogs = np.empty(0)
# snapshots = np.empty(0)
# x_sep = np.empty(0)

# for i in range(0,976):
    
    
graph = the_file[str(975)]

x_sep = graph['mean_pos'][:][:,0]
nprogs = graph['nparts'][:]*part_mass

snapshots = graph['snapshots'][:]
    
    

fig, ax1 = plt.subplots(figsize=(7,4))
ax1.scatter(x_sep,snapshots)
# ax1.set_xlim(0,3)
# ax1.set_yscale('log')
# ax1.set_ylabel('N')
# ax1.set_xlabel('Mass/$M_\odot$')
# ax1.legend()
# plt.tight_layout()
# plt.savefig('output graphs/Halo Mass Distribution.png',dpi=400)
# plt.savefig('output graphs/Halo_dist.eps')
plt.show()
