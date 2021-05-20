# -*- coding: utf-8 -*-
"""
Created on Wed May  5 09:49:27 2021

@author: Andrew
"""

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
from gc import get_referents



# dtype = np.dtype([("Val1", float), ("Val2", float),
#                                           ("Val3", float),
#                                           ("Val7",float),
#                                           ("Val4", float),
#                                           ("Val5", float, (1000,))])


small_array = np.empty((3), dtype=float)

large_array = np.empty((1000), dtype=float) 


print('Total memory used (sys.getsizeof())')
print('-----------------------------------')
print('Small array (100 elements): {} bytes'.format(sys.getsizeof(small_array)))
print('Large array (10,000 elements): {} bytes'.format(sys.getsizeof(large_array)))

print('\n\nMemory size of objectes stored (array.nbytes)')
print('-------------------------------------------------')
print('Small array (100 elements): {} bytes'.format(small_array.nbytes))
print('Large array (10,000 elements): {} bytes'.format(large_array.nbytes))


print('\nDifference')
print('-----------')
print('Small array (100 elements): {} bytes'.format(sys.getsizeof(small_array)-small_array.nbytes))
print('Large array (10,000 elements): {} bytes'.format(sys.getsizeof(large_array)-large_array.nbytes))
