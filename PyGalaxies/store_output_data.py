# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:53:31 2021

@author: Andrew
"""
import numpy
import h5py as h5
import gc




def create_output_file(HDF_output_filepath, halo_properties_dtype,
                        subhalo_properties_dtype):

     # Open file --> this is slow if HDF5 file has remained open on previous
     # run due to an error.
     try:
         output_file = h5.File(HDF_output_filepath, "w")
     except OSError:
         for obj in gc.get_objects():   # Browse through ALL objects
            if isinstance(obj, h5.File):   # Just HDF5 files
                try:
                    obj.close()
                except:
                    pass # Was already closed
         output_file = h5.File(HDF_output_filepath, "w")
                    
        
             
     # Create datasets
     halo_output_dataset = output_file.create_dataset("halo_data", (0,), 
                                                 maxshape=(None,),
                                                 dtype = halo_properties_dtype,
                                                 compression="gzip")
        

     subhalo_output_dataset = output_file.create_dataset("subhalo_data", (0,),
                                             maxshape=(None,),
                                             dtype = subhalo_properties_dtype,
                                             compression="gzip")
     
     
     return output_file
     
     
     
def output_data_this_graph(output_file, array_of_halo_properties,
                        array_of_subhalo_properties):
     

        # Save data to datasets
        
        halo_output_dataset = output_file["halo_data"]
        subhalo_output_dataset = output_file["subhalo_data"]
        
        
        current_length_of_halo_dataset = halo_output_dataset.shape[0]
  
        halo_output_dataset.resize((current_length_of_halo_dataset + 
                                   array_of_halo_properties.shape[0],))
        
        current_length_of_subhalo_dataset = subhalo_output_dataset.shape[0]
        
        subhalo_output_dataset.resize((subhalo_output_dataset.shape[0] + 
                                      array_of_subhalo_properties.shape[0],))
        
        
        halo_output_dataset[current_length_of_halo_dataset:] = array_of_halo_properties
        
        subhalo_output_dataset[current_length_of_subhalo_dataset:] = array_of_subhalo_properties
                  
            

def close_output_file(output_file):
    
    output_file.close()