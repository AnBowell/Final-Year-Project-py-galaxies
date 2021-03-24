# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:56:29 2020

@author: Andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import psutil
import time



class TimeMonitor:
    
    def __init__(self, amount_of_graphs, amount_of_timers_required):
        
        self.amount_of_graphs = amount_of_graphs

        self.timer_start_times = []
    

        self.time_taken_array = np.zeros(amount_of_graphs,
                                          dtype = np.dtype([(
                                              'Time Taken / s',
                                                         np.float32)]))

    def start_timer(self):
        
        self.timer_start_times.append(time.perf_counter())
        
    def stop_timer(self, graph_no):
        
        end_time = time.perf_counter()
        
        start_time = self.timer_start_times.pop()
        
        time_taken = end_time - start_time

        
        self.time_taken_array[graph_no] = time_taken
        

    def save_all_timing_data_to_disk(self,amount_of_halos, amount_of_subhalos):
        
        
        no_subhalo_mask = [amount_of_halos > 0]
        
        amount_of_halos = amount_of_halos[no_subhalo_mask]
        amount_of_subhalos = amount_of_subhalos[no_subhalo_mask]
        
        self.time_taken_array = self.time_taken_array[no_subhalo_mask]
        

        np.savez('../Timing And Memory Tests/TimeTakenPerGraphClassBased.npz',
                 amount_of_halos = amount_of_halos, 
                 amount_of_subhalos = amount_of_subhalos,
                 time_taken = self.time_taken_array['Time Taken / s'])

    
        print('Data successfully saved')
        
        
class MemoryMonitor:
    
    def __init__(self, amount_of_graphs):
        
        self.process  = psutil.Process()
        
        self.mem_store = np.zeros(amount_of_graphs,
                                  dtype = np.dtype([(
                                      'memory [MB]',
                                                 np.float32)]))


    def record_mem(self, graph_no):
        
        new_mem = self.process.memory_info().rss
        
        self.mem_store[graph_no] = new_mem / (1024**2)
        

        
    def save_all_mem_data_to_disk(self, amount_of_halos, amount_of_subhalos):
                
        
        no_subhalo_mask = [amount_of_subhalos > 0]
        
        amount_of_halos = amount_of_halos[no_subhalo_mask]
        amount_of_subhalos = amount_of_subhalos[no_subhalo_mask]
        
        self.mem_store = self.mem_store[no_subhalo_mask]

        

        np.savez('../Timing And Memory Tests/CumulativeMemoryUsedClassBased.npz',
                 amount_of_halos = amount_of_halos, 
                 amount_of_subhalos = amount_of_subhalos,
                 mem_used =self.mem_store['memory [MB]'])

    
        print('Data successfully saved')
        

    
        
        
        
        
        
        
        
        
        
        
        


class AnalyseSavedData:
    
    def __init__(self, filepath_to_class_based, filepath_to_array_based):
        
        self.filepath_to_array_based = filepath_to_array_based
        self.filepath_to_class_based = filepath_to_class_based

        self.open_files_extract_data()
    
    def open_files_extract_data(self):
        
        self.array_based_data = np.load(self.filepath_to_array_based)
        self.filepath_to_class_based = np.load(self.filepath_to_class_based)
        
    def plot_time_against_halos(self):
        
        fig, ax1 = plt.subplots(figsize = (6,4))
        
        ax1.set_ylabel('Processing time / s')
        
        ax1.set_xlabel('Number of Halos in graph')
        
        ax1.scatter(self.array_based_data['amount_of_halos'],
                    self.array_based_data['time_taken'])
        
        ax1.scatter(self.filepath_to_class_based['amount_of_subhalos'],
                    self.filepath_to_class_based['time_taken'])
        
        
        ax1.grid(True,alpha=0.5, linestyle='--')
        ax1.set_xscale('linear')
        ax1.set_xlim(0,15000)
    
        ax1.legend()
        plt.savefig('../../Timing And Memory Tests/Output Figures/TimeTakenPerGraph.png',
                    dpi=450)
    

        plt.plot()
        

class Monitor:
    
    def __init__(self, no_of_graphs, amount_of_funcs):
        
        self.amount_of_funcs = amount_of_funcs
        
        self.function_descs ="""
        1. Initialise Halo properties class (does also find most central halo). \n
        2. Calculate mass in common. \n
        3. Gather DM mass from progenitors. \n
        4. Give hot gas etc to descendents. \n
        5. Calc and set baryon fraction. \n
        6. Collect galaxy progenitor info.This includes giving stellar mass to descendents. \n
        7. Simple function to add bary mass to halo from subhalo stars. \n
        8. Calculates top-up from berhoozi et al. \n
        9 Actually adds the top up from berhoozi et al. \n
        """

        self.time_storage_array = np.empty((no_of_graphs, amount_of_funcs, 2))
        
        
        self.memory_storage_array = np.empty((no_of_graphs,
                                               amount_of_funcs,
                                               2))
        
        self.process  = psutil.Process()
                                           
                                         
    def graph_timer_setup(self, n_halos):    

        self.temp_time_store_array = np.empty((n_halos, self.amount_of_funcs))
        
        self.temp_mem_store_array = np.empty((n_halos, self.amount_of_funcs,2))
        
        return None

    def start_timer(self):
        
        self.graph_x_start_mem = self.process.memory_info().rss
        
        self.graph_x_start_time = time.perf_counter_ns()
        
    def store_func_time(self, func_no, halo_no):
        
        
        self.temp_time_store_array[halo_no, func_no-1] = (time.perf_counter_ns() -
                                                          self.graph_x_start_time)
        
        new_mem = self.process.memory_info().rss
        
        
        self.temp_mem_store_array[halo_no, func_no-1,0] = new_mem
        
        self.temp_mem_store_array[halo_no, func_no-1,1] = (new_mem
                                                           - self.graph_x_start_mem)
        
        #self.temp_mem_store_array[halo_no, func_no-1,2] = (self.process.memory_percent())
        
        self.graph_x_start_mem = new_mem
    
        # Down here to avoid mem_processing
        self.graph_x_start_time  = time.perf_counter_ns() 
        return None
    

    def store_average_graph_times(self, graph_no):
        
        self.time_storage_array[graph_no, :, 0] = np.mean(self.temp_time_store_array,
                                                          axis = 0)
        
        self.time_storage_array[graph_no, :, 1] = np.std(self.temp_time_store_array,
                                                          axis = 0)
        
        self.memory_storage_array[graph_no, :, 0] = self.temp_mem_store_array[-1, :, 0]
        
        self.memory_storage_array[graph_no, :, 1] = np.mean(
                                                    self.temp_mem_store_array[:,:,1]
                                                    , axis = 0)


    def save_timing_stats(self, output_file_path, file_name):
        
        file_name = file_name.split('.hdf')[0].split('/')[-1]

        np.save('{}Timing Data {}'.format(output_file_path, file_name),
                self.time_storage_array)
        
        
        # Calculate the mem percentage and save memory data
        
        mem_percentage = (self.memory_storage_array[:,:,0] / 
                          psutil.virtual_memory().total) * 100
        
        shape = list(self.memory_storage_array.shape)
        
        shape[-1] = shape[-1] + 1
        
        z = np.zeros(tuple(shape))
        
        z[:,:,:2] = self.memory_storage_array
        
        z[:,:,-1] = mem_percentage
        
        np.save('{}Memory Data {}'.format(output_file_path, file_name), z)
        
        
        
       
    def plot_timing_barchart(self,output_file_path, file_name, save_file_path):
                
        file_name = file_name.split('.hdf')[0].split('/')[-1]
        
        plotting_data = np.load('{}Timing Data {}.npy'.format(output_file_path, 
                                                           file_name))
        
        
        labels = np.arange(1,len(plotting_data[0,:,0])+1,1)

        
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        ax1.bar(labels,np.nanmean(plotting_data[:,:,0],axis= 0))
        ax1.set_xticklabels(labels)
        ax1.set_xticks(labels)
        ax1.set_ylabel('Mean time taken [ns] over all graphs with sub-halos.')
        ax1.grid(True, alpha = 0.6)
        ax1.set_axisbelow(True)
        ax1.set_xlabel('Function Measured')
        ax2.set_axis_off()
        ax2.text(0.1,0.1,self.function_descs)
        plt.tight_layout()
        plt.savefig('{} Bar graph for {}.jpg'.format(save_file_path, 
                                                           file_name), dpi = 600)
        plt.show()
        
        
    def plot_memory_barchart(self,output_file_path, file_name, save_file_path):
        
        file_name = file_name.split('.hdf')[0].split('/')[-1]
        
        plotting_data = np.load('{}Memory Data {}.npy'.format(output_file_path, 
                                                           file_name))
        
        
        labels = np.arange(1,len(plotting_data[0,:,1])+1,1)

        plotting_data = np.where((plotting_data[:,:,1] < 0),np.nan,plotting_data[:,:,1])


        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        ax1.bar(labels,np.nanmean(plotting_data,axis= 0))
        ax1.set_xticklabels(labels)
        ax1.set_xticks(labels)
        
        ax1.set_ylabel('Mean memory used [B] over all graphs with sub-halos.')
        ax1.grid(True, alpha = 0.6)
        ax1.set_axisbelow(True)
        ax1.set_xlabel('Function Measured')
        ax2.set_axis_off()
        ax2.text(0.1,0.1,self.function_descs)
        plt.tight_layout()
        plt.savefig('{} Memory Bar graph for {}.jpg'.format(save_file_path, 
                                                           file_name), dpi = 600)
        plt.show()