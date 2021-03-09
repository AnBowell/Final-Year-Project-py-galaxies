# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:15:04 2021

@author: Andrew
"""
import sys
import os
sys.path.insert(0,os.path.join(os.getcwd(),'../'))
sys.path.insert(0,os.path.join(os.getcwd(),'../PyGalaxies'))
import numpy as np
import Monitor


filepath_to_class_based = '../../Timing And Memory Tests/TimeTakenPerGraphClassBased.npz'
filepath_to_array_based = '../../Timing And Memory Tests/TimeTakenPerGraph.npz'

analyse = Monitor.AnalyseSavedData(filepath_to_class_based,
                                   filepath_to_array_based)


analyse.plot_time_against_halos()