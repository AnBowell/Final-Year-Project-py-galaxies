# -*- coding: utf-8 -*-
"""
Created on Sat May  1 10:45:05 2021

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import time

plt.style.use('classic')

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=14)

dt = np.dtype([('redshift', np.float32), ('pos', np.float64, (3,)),
               ('mass', np.float32), ('vel', np.float64, (3,))])


test_array = np.full(10000,1,dtype=dt)



time_taken_1 = []

time_taken_2 = []


how_many_times = [10,100,1000,10000,100000]



for times in how_many_times:
    
    

    start_time = time.perf_counter()
    for counter in range(0,times):
        y = test_array["redshift"][1]
    end_time = time.perf_counter()
    
    time_taken_1.append(end_time - start_time)

    
    
    
    start_time = time.perf_counter()
    
    for counter in range(0,times):
        x = test_array[1]["redshift"]
    
    end_time = time.perf_counter()
    
    time_taken_2.append(end_time - start_time)
    


fig,ax1 = plt.subplots(figsize = (7,4))


ax1.plot(how_many_times,time_taken_1,label='[Column][Row]',linewidth=2.5)
ax1.plot(how_many_times,time_taken_2,label='[Row][Column]',linewidth=2.5)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Processing Time / s')
ax1.set_xlabel('Number Of Times Indexed ')
ax1.legend(loc='upper left')
ax1.grid(True,linestyle='--')
ax1.set_axisbelow(True)
plt.tight_layout()
plt.savefig('output graphs/Array_based_index.png',dpi =400,bbox_inches =0)
plt.savefig('output graphs/array_index.eps',bbox_inches =0)
plt.show()