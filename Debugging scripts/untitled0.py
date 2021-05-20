# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:27:59 2021

@author: Andrew
"""
import numpy as np
import time
import matplotlib.pyplot as plt

class Test:
    def __init__(self):
        
        self.val = 1.0

x = time_taken_array = np.full(100,1, dtype=np.dtype([("Test", np.float64)]))

y = [Test() for i in range(0,100)]


time_taken_array = []

time_taken_class = []


lengths = [10,20,50,100,1000,10e3,10e4,10e5]

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


ax1.plot(lengths,time_taken_array)
ax1.plot(lengths,time_taken_class)
ax1.set_yscale('log')
ax1.set_xscale('log')

plt.show()