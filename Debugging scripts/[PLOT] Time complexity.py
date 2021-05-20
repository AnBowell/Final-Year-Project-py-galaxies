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
plt.style.use('classic')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)
plt.rcParams['legend.title_fontsize'] = 'x-large'

filepath_to_class_based = '../../Timing And Memory Tests/TimeTakenPerGraphClassBased.npz'
filepath_to_array_based = '../../Timing And Memory Tests/TimeTakenPerGraphArrayBased.npz'


array_based_data = np.load(filepath_to_array_based)
class_based_data = np.load(filepath_to_class_based)



p = np.poly1d(np.polyfit(array_based_data['amount_of_halos'],array_based_data['time_taken'],
               1))
p2 = np.poly1d(np.polyfit(class_based_data['amount_of_halos'],class_based_data['time_taken'],
               1))



to_fit_x = np.linspace(0,13000,500)




fig, ax1 = plt.subplots(figsize = (7,4))

ax1.set_ylabel('Processing Time / s')

ax1.set_xlabel('Number Of Haloes In Graph')

ax1.scatter(array_based_data['amount_of_halos'],
            array_based_data['time_taken'],label='Array-based\n$y=${:.1e}$x+${:.1e}'.format(p[1],p[0]),color='blue')

ax1.plot(to_fit_x, p(to_fit_x))
ax1.plot(to_fit_x, p2(to_fit_x))
ax1.scatter(class_based_data['amount_of_halos'],
            class_based_data['time_taken'],label='Class-based\n$y=${:.1e}$x+${:.1e}'.format(p2[1],p2[0]),color='green')




ax1.grid(True, linestyle='--')
ax1.set_axisbelow(True)
ax1.set_xscale('linear')
ax1.set_xlim(0,13500)
ax1.set_ylim(0,2.2)
ax1.minorticks_on()
# ax.tick_params(axis='both', which='minor',bottom=True)
ax1.legend(loc='upper left')
plt.tight_layout()
plt.savefig('output graphs/TimeTakenPerGraph.png',
            dpi=450)
plt.savefig('output graphs/TimeTakenPerGraph.eps')
plt.plot()