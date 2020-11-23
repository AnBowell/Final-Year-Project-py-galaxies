# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 12:54:59 2020

@author: Andrew
"""



import numpy as np
from glob import glob
import rasterio 

file_list = sorted(glob('Somefilepath to some files here'))

amount_of_obs = len(file_list)


# Get shape and then find nice divisible strips
shape_of_sat_images = rasterio.open(file_list[0]).read(1).shape

# Find a decent number to divided by. If image 1000 x 1000. you could do 100 
# strips of 10 x 1000 pixels.

window_x_size = 10
Window_y_size = 1000


for strip_x in range(0,100):
    
    storage_array = np.empty(amount_of_obs,10,1000)

    
    for file_counter, file in enumerate(file_list):
        
        # Think I have that windowed read in the right direction.
        
        storage_array[file_counter, :, :] = rasterio.open(file).read(1,
                                            window = Window(strip_x*10,
                                                            0,
                                                            10,
                                                            1000)
                                            
    # Now do smoothing on your storage array
    
    smoothed_array = np.empty(amount_of_obs,10,1000)
    
    for x in range(0, 10):
        for y in range(0,1000):
                
            data_to_smooth = storage_array[:,x,y]
            #Do some smoothing
            
            smoothed_array[:,x,y] = some_smoothed_data
            
            
    # save smoothed_array 
            
    
    
    
    
                                            


