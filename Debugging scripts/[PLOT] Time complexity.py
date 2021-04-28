# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:15:04 2021

@author: Andrew
"""
import sys
import os
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies'))
import numpy as np
import Monitor

import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
# plt.style.use('classic')
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)
plt.rcParams['legend.title_fontsize'] = 'x-large'

filepath_to_class_based = '../../Timing And Memory Tests/TimeTakenPerGraphClassBased.npz'
filepath_to_array_based = '../../Timing And Memory Tests/TimeTakenPerGraphArrayBased.npz'


array_based_data = np.load(filepath_to_array_based)
filepath_to_class_based = np.load(filepath_to_class_based)


fig, ax1 = plt.subplots(figsize = (6,4))

ax1.set_ylabel('Processing time / s')

ax1.set_xlabel('Number of Halos in graph')

ax1.scatter(array_based_data['amount_of_halos'],
            array_based_data['time_taken'],label='Array based')

ax1.scatter(filepath_to_class_based['amount_of_subhalos'],
            filepath_to_class_based['time_taken'],label='Class based')


ax1.grid(True,alpha=0.5, linestyle='--')
ax1.set_xscale('linear')
ax1.set_xlim(0,15000)

ax1.legend()
plt.savefig('../../Timing And Memory Tests/Output Figures/TimeTakenPerGraph.png',
            dpi=450)


plt.plot()