# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 09:05:48 2021

@author: Andrew
"""

import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'PyGalaxies')) # This is needed to import C routines.
from PyGalaxies import Parameters, Monitor, HDF_graph_properties, Halos
import time
import numpy as np



model_param_filepath='Input_Params/input_params.yml'
debug_flag = False
verbosity = 1 
time_code = True


model_params = Parameters.ModelParams(model_param_filepath,verbosity,debug_flag)
#model_params.precalculate_common_constants()
model_params.load_paramters_to_C()
model_params.read_in_data_tables_from_c()
model_params.output_params() 



# Open the input HDF5 file containing graph groups and create an output file.
HDF_properties = HDF_graph_properties.HDFProperties(model_params)

    
start_time = time.time()


for graph_ID in range(0,HDF_properties.no_of_graphs)[:]:
    
    if HDF_properties.nsubhalos_in_graph[graph_ID] < 1:
        continue
    
    # Read in data from the graph
    
    graph_properties = HDF_graph_properties.GraphProperties(graph_ID,
                                       HDF_properties.graph_input_file,
                                       model_params,
                                       HDF_properties.part_mass)
 

    # Loop over and intialise the halo and sub-halo arrays.  

    list_of_halo_properties = [Halos.HaloProperties(str(graph_ID), halo_ID,
                               graph_properties, HDF_properties.subhalo_dtype) for 
                               halo_ID in range(0, HDF_properties.nhalos_in_graph[graph_ID])]
    
    
    list_of_subhalo_properties = [Halos.SubhaloProperties(graph_ID, host_halo_ID, 
                                  subhalo_ID,graph_properties) for host_halo_ID,subhalo_ID in zip(
                                  graph_properties.host_halos, 
                                  range(0, HDF_properties.nsubhalos_in_graph[graph_ID]))]
    
    
    
    for snap_ID in graph_properties.generation_id:
        
        if snap_ID == HDF_properties.no_data_int:
            continue
        

        # Be careful with this. Will error on very first generation. Needs fixing.
        # Try except would be optimi
        
        try:
            dt = model_params.snap_times[snap_ID] - model_params.snap_times[snap_ID-1]
        except IndexError:
            dt = model_params.snap_times[snap_ID]

        this_generation_halo = graph_properties.graph_halo_ids[
                               graph_properties.generation_start_index[snap_ID]:
                               graph_properties.generation_start_index[snap_ID] + 
                               graph_properties.generation_length[snap_ID]]
            
            

        
            
        for halo_ID in this_generation_halo:
            
            halo = list_of_halo_properties[halo_ID]
            
            halo.snap_ID = snap_ID


            halo.sum_baryon_and_topup(model_params.halo_descend_attrs,
                                      model_params.subhalo_descend_attrs,
                                      list_of_subhalo_properties,
                                      model_params.f_baryon)
            
            
            halo.calculate_SFR_hot_gas_used(dt, list_of_subhalo_properties)



            halo.calculate_hot_gas_temp(model_params.H0, model_params.mu, 
                                        model_params.m_p, model_params.k_B)
            
            halo.calculate_metal_dependent_cooling_rate()
            
            
            halo.calc_mass_of_hot_gas_cooled(model_params.mu, model_params.m_p, 
                                             model_params.G, model_params.beta_prof_ratio,
                                             model_params.beta_prof_ratio_arctan,
                                             dt, list_of_subhalo_properties)


            halo.calc_halo_props_descend(HDF_properties.part_mass, list_of_halo_properties,
                                         model_params.halo_descend_attrs)
            
            
            
            
            
            for subhalo_ID in halo.sub_graph_halo_ids:
                
                subhalo = list_of_subhalo_properties[subhalo_ID]
                
                subhalo.snap_ID = snap_ID
                
                subhalo.calc_subhalo_props_descend(list_of_subhalo_properties,
                                           model_params.subhalo_descend_attrs, 
                                           list_of_halo_properties)
                
                
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





            # halo.halo_descendent_props()

            # halo.calc_halo_DM_descend(HDF_properties.part_mass)

            # halo.calc_halo_attrs_descend(HDF_properties.part_mass, array_of_halo_properties,
            #                              HDF_properties.halo_descend_attrs)


            
            # halo.calc_subhalo_attrs_descend(graph_properties, HDF_properties, array_of_halo_properties)
            
            # halo.add_halo_baryon_mass_then_topup(HDF_properties)
            
            # halo.set_baryon_fraction(array_of_halo_properties, model_params.f_baryon)
            
            # halo.calculate_hot_gas_temp(model_params.H0, model_params.mu, 
            #                             model_params.m_p, model_params.k_B)
            
            # halo.calculate_metal_dependent_cooling_rate()
            
            # halo.add_stellar_mass() # Adds on Berhoozi 
            
            halo.done = True
            
            # array_of_halo_properties[halo_ID] = halo

            HDF_properties.n_halo +=1 
            
        
            
            

    # Save output

    HDF_properties.output_halos(list_of_halo_properties)
    HDF_properties.output_subhalos(list_of_subhalo_properties,
                                   model_params.subhalo_output_list)



# Close input file
HDF_properties.close_graph_io(HDF_properties.graph_input_file)


# Close output files
if HDF_properties.halo_output_iRec > 0: 

    halo_output_iRec = HDF_properties.flush_output()
    
HDF_properties.close_graph_io(HDF_properties.halo_output_file)


end_time = time.time()
print('This took {} seconds'.format(end_time-start_time))