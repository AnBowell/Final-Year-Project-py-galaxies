# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 18:01:51 2021

@author: Andrew
"""
import numpy as np
import h5py as h5
from collections import Counter


def open_HDF_file(HDF_graph_filepath):
    
    graph_input_file = h5.File(HDF_graph_filepath, "r")
    
    amount_of_graphs = len(graph_input_file["/nhalos_in_graph/"])
    
    return graph_input_file, amount_of_graphs


# Needs some optimization and cleaning up - pretty gnarly function.
def get_datastruct_info(HDF_file, amount_of_graphs_in_file):
    
    for graph_ID in range(0, amount_of_graphs_in_file):
    
        dataset_lengths_dtypes = [(dataset,len(HDF_file[str(graph_ID)][dataset]),
                                  HDF_file[str(graph_ID)][dataset].dtype,
                                  int(HDF_file[str(graph_ID)][dataset].size /
                                  len(HDF_file[str(graph_ID)][dataset])))
                                  for dataset in 
                                  HDF_file[str(graph_ID)].keys()]
        
        
        if len(set([x[1] for x in dataset_lengths_dtypes])) != 5:
            continue
        
        dataset_lengths_dtypes.sort(key=lambda x:x[1])
        
        # dtypes = [HDF_file[str(graph_ID)][dataset].dtype for dataset in 
        #                     HDF_file[str(graph_ID)].keys()]
        
        datalength_occurrences = sorted(Counter(list(zip(*dataset_lengths_dtypes))[1]).items())
        
        list_of_array_data_structure = []
        start_index = 0
        
        for length_occurrence in datalength_occurrences:
        
            same_length_datasets = dataset_lengths_dtypes[start_index:start_index + 
                                                    length_occurrence[1]]
            
            
            # Getting around a future numpy depreciation where a numpy dtype
            # can't have a dimension of 1.
            list_to_append = [(x[0],x[2],x[3]) if x[3] > 1 else
                              (x[0],x[2]) for x in same_length_datasets]
            
            list_of_array_data_structure.append(list_to_append)
            
            
            
            start_index += length_occurrence[1]
            
        return list_of_array_data_structure
    
    print("""Error, this is unlikely to happen. The datasets do not have 5 unique
          lengths in any graph. This means the five arrays could not be made
          (Generation IDs, halo properties, sub_halo properties, halo desc and
           prog info, sub_halo desc and prog info.""")
    

def create_storage_structure(open_graph_data, graph_ID, list_of_array_data_structure):
    
    
 
    
    dict_of_numpy_structured_arrays = {}
    
    for dataset_names_and_others in list_of_array_data_structure:
        
        dataset_names = [item[0] for item in dataset_names_and_others]
    
        
 
        if 'generation_id' in dataset_names:
            array_type = 'Generation Data'
        elif 'host_halos' in dataset_names:
            array_type = 'Subhalo Data'
        elif 'sub_direct_desc_contribution' in dataset_names:
            array_type = 'Subhalo Contribution Data'
        elif 'nparts' in dataset_names:
            array_type = 'Halo Data'
        elif 'direct_desc_contribution' in dataset_names:
             array_type = 'Halo Contribution Data'
        else:
            print('Something has gone wrong')
            
        
        
        # Open first dataset and get length. All in this list will be of same length
        
        dataset_length = len(open_graph_data[dataset_names[0]])
        
        data_storage_array = np.empty(dataset_length, 
                                      np.dtype(dataset_names_and_others))
    
        
    
        for dataset_name_and_other in dataset_names_and_others:
            # dataset_name_and_other can have a name, dtype and length. or just
            # name and dtype. This is due to numpy not liking 1 as a dimension
            # for the strcutured array record. name always first so it can always
            # be found with an index of 0
            
            dataset = open_graph_data[dataset_name_and_other[0]][:]
            
            data_storage_array[dataset_name_and_other[0]] = dataset
        
        dict_of_numpy_structured_arrays[array_type] = data_storage_array
    
    return dict_of_numpy_structured_arrays



def load_in_graph_data(HDF_file, graph_ID):
    
    open_graph_data = HDF_file[str(graph_ID)]
    
    graph_keys = open_graph_data.keys()
    
    return open_graph_data, graph_keys
 