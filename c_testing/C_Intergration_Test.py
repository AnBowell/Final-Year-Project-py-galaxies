# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:04:02 2020

@author: Andrew
"""
import c_test

from time import perf_counter

def py_prime_number(n):
    
    for i in range(2, int(n/2)):
        
        if n % i == 0:
            
            return 1
    
    return 0
    
    

number_to_test = 0




start_time = perf_counter()

x = c_test.checkPrimeNumber_c(number_to_test)

end_time = perf_counter()

print('Is number prime? : {}'.format(x))

print('C Code took {}s'.format(end_time - start_time))



start_time = perf_counter()
y = py_prime_number(number_to_test)

end_time = perf_counter()

print('Is number prime? : {}'.format(y))
print('Py Code took {}s'.format(end_time - start_time))