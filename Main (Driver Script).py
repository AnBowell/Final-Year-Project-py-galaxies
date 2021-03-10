# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 17:46:42 2021

@author: Andrew
"""
import os, sys
sys.path.insert(0,os.path.join(os.getcwd(),'PyGalaxies')) # This is needed to import C routines.

from PyGalaxies import parameters, read_graph_data, halo, store_output_data, Monitor
import time
import numpy as np




model_param_filepath='Input_Params/input_params.yml'
debug_flag = False
verbosity = 1 
time_code = True

#The baryonic properties that descend through the graph. These are the 
halo_properties_that_descend = ['mass_hot_gas', 'mass_cold_gas',
                            'mass_ejected_gas', 'mass_intracluster_light']

subhalo_properties_that_descend = ['mass_cold_gas', 'mass_stellar']

model_params = parameters.ModelParams(model_param_filepath,verbosity,debug_flag)
model_params.load_paramters_to_C()
model_params.read_in_data_tables_from_c()
# model_params.output_params() 

start_time = time.perf_counter()

# output must be created before input read in --> Error handling closes any
# open instances of h5py files to ensure no error is thrown when trying to create
# a h5py object that is already open.
output_file = store_output_data.create_output_file(model_params.output_filepath,
                                                   model_params.halo_properties_dtype,
                                                   model_params.subhalo_properties_dtype)

HDF_file, amount_of_graphs_in_file = read_graph_data.open_HDF_file(model_params.input_filepath)



list_of_array_data_structure = read_graph_data.get_datastruct_info(HDF_file,
                                                                   amount_of_graphs_in_file)

if time_code:
    time_monitor = Monitor.TimeMonitor(amount_of_graphs_in_file,
                                       amount_of_timers_required=1)


for graph_ID in range(0, amount_of_graphs_in_file):
    
    if time_code:
        time_monitor.start_timer()
    
    open_graph_data, graph_keys = read_graph_data.load_in_graph_data(HDF_file, graph_ID)
    
    
    # Skip over any graph that does not have sub halos. 
    if 'sub_nparts' not in  graph_keys:
        continue




    data_struct = read_graph_data.create_storage_structure(open_graph_data,
                                                                  graph_ID,
                                                                  list_of_array_data_structure)

    
    generation_data = data_struct['Generation Data']
    halo_data = data_struct['Halo Data']
    halo_contrib_data = data_struct['Halo Contribution Data']
    subhalo_data = data_struct['Subhalo Data']
    subhalo_contrib_data = data_struct['Subhalo Contribution Data']
    
    
    array_of_halo_properties = np.zeros(model_params.nhalos_in_graph[graph_ID],
                                        dtype = model_params.halo_properties_dtype)
    
    array_of_subhalo_properties = np.zeros(model_params.nsubhalos_in_graph[graph_ID],
                                        dtype = model_params.subhalo_properties_dtype)
    
    
    # Modifies arrays in place --> memory reference. This is the same
    # for the rest of the functions.
    halo.initialize_subhalo_properties(array_of_subhalo_properties, 
                                       subhalo_properties_that_descend)
    
    halo.initialize_halo_properties(array_of_halo_properties, 
                                      halo_properties_that_descend)
   
    
    #~~~~~~~~~~ Do vectorized calculations which do not need the loop ~~~~~~~~#
    
    halo.store_essential_halo_data(array_of_halo_properties, halo_data,
                                  graph_ID, model_params.part_mass) 
    
    halo.store_essential_subhalo_data(array_of_subhalo_properties, subhalo_data,
                                      model_params.part_mass, graph_ID) 
    
    # Calculate the temperature of the hot gas.
    
    halo.calculate_hot_gas_temp(array_of_halo_properties,
                                halo_data["rms_radius"],halo_data["redshifts"],
                                model_params.H0, model_params.mu, 
                                model_params.m_p, model_params.k_B)
    

    # Model that relates halo DM mass and redshift to expected stellar mass.
    halo.calculate_Berhoozi_stellar_mass(array_of_halo_properties,
                                         halo_data["redshifts"])
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    

    for snap_ID in generation_data['generation_id']:
        

        
        if snap_ID == model_params.no_data_int:

            continue
        
        try:
            dt = model_params.snap_times[snap_ID] - model_params.snap_times[snap_ID-1]
        except IndexError:
            dt = model_params.snap_times[snap_ID]
        

        
        # The Halo IDs are just 0-n_halos. Array is just created by numpy arange for now
        # as it's a wast of space to store this in the HDF5 file.
        
        this_generation_halo = np.arange(generation_data["generation_start_index"][snap_ID],
                                         generation_data["generation_start_index"][snap_ID] + 
                                         generation_data["generation_length"][snap_ID], 1)
        
        for halo_ID in this_generation_halo:
            
            array_of_halo_properties['snap_ID'][halo_ID] = snap_ID

            n_subhalo = halo_data["nsubhalos"][halo_ID]
            subhalo_start_index = halo_data["subhalo_start_index"][halo_ID]
            
            if n_subhalo > 0:

                # This func finds the halo's central subhalo.
                
                halo.find_central_subhalo(array_of_halo_properties,
                                    array_of_subhalo_properties, halo_data, 
                                    subhalo_data, halo_ID, model_params.no_data_int,
                                    n_subhalo, subhalo_start_index)
                
            else:
                    
                # If there are no sub_halos --> fill with no data int.
                array_of_halo_properties["central_subhalo_ID"][halo_ID] = model_params.no_data_int
            
            
            # Sums up all the baryon mass in the halo and adds any additional
            # infalling hot gas. L-galaxies reionisation function is also
            # included within this function. (Reduces the baryon fraction)
            
            halo.sum_baryon_and_topup(array_of_halo_properties, halo_ID, 
                                      model_params.f_baryon,
                                      halo_properties_that_descend,
                                      subhalo_properties_that_descend,
                                      array_of_subhalo_properties,
                                      halo_data)
 
    
            
            halo.calculate_SFR_hot_gas_used(array_of_halo_properties, 
                               array_of_subhalo_properties, dt, halo_ID, 
                               model_params.no_data_int)
        
    
    
    
            # Calculates metal dependent cooling rate of the hot gas. Inside the 
            # loop as it uses the Halo's metalicity which can change. But the
            # temperature of the hot gas is a vectorized calculation executed
            # before the loop. 
            
            halo.calculate_metal_dependent_cooling_rate(array_of_halo_properties,
                                                        halo_ID)
        
        
            # Use the cooling rate to calculate the mass of hot gas cooled onto
            # the central subhalo. 
        
            halo.calc_mass_of_hot_gas_cooled(model_params.mu, model_params.m_p, 
                                             model_params.G, model_params.beta_prof_ratio,
                                             model_params.beta_prof_ratio_arctan, dt, model_params.no_data_int,
                                             array_of_halo_properties,
                                             array_of_subhalo_properties, halo_ID,
                                             halo_data)
        
        
        
        
            # At end of the halo processing. Push poperties down to descendent 
            # halos.
        
            
            halo.calc_halo_props_descend(array_of_halo_properties, 
                                 halo_properties_that_descend,
                                 halo_data, halo_contrib_data, halo_ID)
            
        
            # Run this function as the last one in the loop. It pushes all properties
            # accumulated by the subhalo down to it's main descendent.
                        
            for subhalo_ID in np.arange(subhalo_start_index,
                            subhalo_start_index + n_subhalo): # Using np.arange as subhalo index 0-->nsubhalo
                # array of halo properties also included. If there is no descendent
                # subhalo to push to, the stellar mass goes into the intracluster medium (halo).
                
                array_of_subhalo_properties['snap_ID'][subhalo_ID] = snap_ID
                
                halo.calc_subhalo_props_descend(array_of_halo_properties,
                                                array_of_subhalo_properties,
                                                subhalo_ID,subhalo_data,
                                                subhalo_contrib_data,
                                                subhalo_properties_that_descend)


    if time_code:
        time_monitor.stop_timer(graph_ID)
        
    
    if graph_ID % 100 == 0:
        print('Graph {} is done'.format(graph_ID))
        
        
    # Save output for each graph. 
    
    store_output_data.output_data_this_graph(output_file,
                                            array_of_halo_properties,
                                            array_of_subhalo_properties)
    
    

store_output_data.close_output_file(output_file)
            
        
end_time = time.perf_counter()


print('This took {}s'.format(end_time - start_time))


if time_code:
   time_monitor.save_all_timing_data_to_disk(HDF_file['nhalos_in_graph'][:],
                                             HDF_file['sub_nhalos_in_graph'][:])