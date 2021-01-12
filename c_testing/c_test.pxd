# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 10:00:45 2020

@author: Andrew
"""
cdef extern from "c_test.h":
    int checkPrimeNumber (int n)