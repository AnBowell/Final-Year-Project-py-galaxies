# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:57:21 2021

@author: Andrew
"""
import numpy as np
import time
import matplotlib 
import matplotlib.pyplot as plt
import sys
from gc import get_referents

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
        self.val1 = 1.0
        self.val2 = 1.0
        self.val3 = 1.0
        self.val4 = 1.0
        self.val5 = np.empty(10000,dtype=float)


# x  = np.empty(1000, dtype=np.dtype([("Val1", float),
#                                       ("Val2", float),
#                                       ("Val3", float),
#                                       ("Val4", float),
#                                       ("Val5", float, (100,))]))



# y = [Test() for i in range(0,1000)]


# print(sys.getsizeof(y),sys.getsizeof(x))


# def getsize(obj):
#     """sum size of object & members."""
#     seen_ids = set()
#     size = 0
#     objects = obj
#     while objects:
#         need_referents = []
#         for obj in objects:
#             if id(obj) not in seen_ids:
#                 seen_ids.add(id(obj))
#                 size += sys.getsizeof(obj)
#                 need_referents.append(obj)
#         objects = get_referents(*need_referents)
#     return size


def getsize(obj):
    """sum size of object & members."""
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


# print(getsize(y[0]),getsize(x[0]))

# FACTER of two offeset is didficult to understand. Python uses memoery poorly. 

time_taken_array = []

time_taken_class = []


lengths = [1,10,20,50,100,1000,10e3,10e4]

for length in lengths:
    
    length = int(length)

    x  = np.empty(length, dtype=np.dtype([("Val1", float),
                                          ("Val2", float),
                                          ("Val3", float),
                                          ("Val4", float),
                                          ("Val5", float, (10000,))]))
    
    time_taken_array.append(getsize(x))
    
    
   
    y = [Test() for i in range(0,length)]

    time_taken_class.append(getsize(y))

    print(length,' Done')




fig,ax1 = plt.subplots(figsize = (7,4))
ax1.plot(lengths,np.array(time_taken_array)/(10**6),label='NumPy Structured Array',linewidth=2.5)
ax1.plot(lengths,np.array(time_taken_class)/(10**6),label='List Of Classes',linewidth=2.5)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Size Of Object / MB')
ax1.set_xlabel('Length Of Object')
ax1.legend(loc='upper left')
ax1.grid(True,linestyle='--')
ax1.set_axisbelow(True)
plt.tight_layout()
plt.savefig('output graphs/Store_mem.eps',bbox_inches =0)
plt.savefig('output graphs/Store_mem.png',dpi =400,bbox_inches =0)

plt.show()