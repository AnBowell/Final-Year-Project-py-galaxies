# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 16:17:49 2021

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5



hdf_filepath = '../../Output Data/SMT13Output_class_based.hdf5'

the_file = h5.File(hdf_filepath,'r')

dataset = the_file['subhalo_data']


ber_stellar_mass = dataset['stellar_mass'][:]

cold_gas_mass =  dataset['cold_gas_mass'][:]



stellar_mask = ber_stellar_mass > 0



ber_stellar_mass = ber_stellar_mass[stellar_mask]





cold_gas_mass = cold_gas_mass[stellar_mask]

gas_mask = cold_gas_mass > 0

ber_stellar_mass = ber_stellar_mass[gas_mask]
cold_gas_mass = cold_gas_mass[gas_mask]



fig, ax1 = plt.subplots(figsize=(8,5))
ax1.scatter(ber_stellar_mass, cold_gas_mass)

ax1.set_ylabel('Cold gas mass / M_o')
ax1.set_xlabel('Behroozi stellar mass / M_o')

ax1.set_title('On average, cold gas has {:.2f} more mass than Behroozi stellar'.format((cold_gas_mass/ber_stellar_mass).mean()))

plt.savefig('output graphs/Stellar2 vs cold gas.png', dpi= 400)
plt.show()

print()
