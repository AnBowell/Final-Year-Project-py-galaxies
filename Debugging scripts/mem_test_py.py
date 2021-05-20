# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:03:45 2021

@author: Andrew
"""
import sys
from gc import get_referents

big_list = [i for i in range(10000000)]

class Test:

    def __init__(self):
        self.val1 = 1.0
        self.val2 = 1.0
        self.val3 = 1.0
        self.val5 = 1.0
        self.val4 = 1.0
        self.list = big_list

def getsize(obj, loop):
    """sum size of object & members."""
    seen_ids = set()
    size = 0
    objects = [obj]
    for loop_counter in range(loop):
        need_referents = []
        for obj in objects:
            if id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


print(getsize(Test(),1))
print(getsize(Test(),2))
print(getsize(Test(),3))
print(getsize(Test(),4))
print(getsize(Test(),5))
print(getsize(Test(),6))