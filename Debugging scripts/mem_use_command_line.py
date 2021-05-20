# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:57:21 2021

@author: Andrew
"""

import sys
from gc import get_referents
import numpy

big_list = [i for i in range(10000)]

class Test:

    def __init__(self):
        self.val1 = 1.0
        self.val2 = 1.0
        self.val3 = 1.0
        self.val4 = 1.0
        self.val5 = big_list
        # self.val5 = np.empty(10000,dtype=float)


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


time_taken_array = []

time_taken_class = []


lengths = [1,10,20,50,100,1000,10e3,10e4]

for length in lengths:
    
    length = int(length)
   
    y = [Test() for i in range(0,length)]

    time_taken_class.append(getsize(y))

    print(length,' Done')



print(time_taken_class)