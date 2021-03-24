# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 08:49:41 2021

@author: Andrew
"""
import sys
import os
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies'))
import numpy as np
import Monitor


filepath_to_class_based = '../../Timing And Memory Tests/CumulativeMemoryUsedClassBased.npz'
filepath_to_array_based = '../../Timing And Memory Tests/CumulativeMemoryUsedArrayBased.npz'

analyse = Monitor.AnalyseSavedData(filepath_to_class_based,
                                   filepath_to_array_based)


analyse.plot_cum_mem()