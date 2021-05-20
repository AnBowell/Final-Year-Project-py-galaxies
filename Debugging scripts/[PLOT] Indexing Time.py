# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:27:59 2021

@author: Andrew
"""
import numpy as np
import time
import matplotlib.pyplot as plt

plt.style.use('classic')

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large',direction='in')
plt.rc('ytick', labelsize='x-large',direction='in')

# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=14)


class Test:
    def __init__(self):
        
        self.val = 1.0

x = time_taken_array = np.full(100,1, dtype=np.dtype([("Test", np.float64)]))

y = [Test() for i in range(0,100)]


time_taken_array = []

time_taken_class = []


lengths = [10,20,50,100,1000,10e3,10e4,10e5,10e6]

for length in lengths:
    
    start_time = time.perf_counter()
    for i in range(0,int(length)):
        z = x['Test'][1]
    end_time = time.perf_counter()

    time_taken_array.append(end_time- start_time)
    
    
    start_time = time.perf_counter()
    for i in range(0,int(length)):
        z2 = y[1].val
    end_time = time.perf_counter()
    time_taken_class.append(end_time - start_time)

    print(length,' Done')

fig,ax1 = plt.subplots(figsize = (7,4))


ax1.plot(lengths,time_taken_array,label='NumPy Structured Array',linewidth=2.5)
ax1.plot(lengths,time_taken_class,label='List Of Classes',linewidth=2.5)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Processing Time / s')
ax1.set_xlabel('Number Of Times Indexed ')
ax1.legend(loc='upper left')
ax1.grid(True,linestyle='--')
ax1.set_axisbelow(True)
plt.tight_layout()
plt.savefig('output graphs/Indexing Time.png',dpi =400,bbox_inches =0)
plt.savefig('output graphs/Index_time.eps',bbox_inches =0)
plt.show()