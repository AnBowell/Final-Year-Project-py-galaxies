# -*- coding: utf-8 -*-
""" L-galaxies python interface

This code is designed to be clean therefore there are three main classes. 
HDFProperties,GraphProperties, and HaloProperties. All other classes added to 
this file will be part of different physics routines. 

Attributes 
----------
n_halo_max : float
    Development limiter.

Created on Fri Oct  9 15:13:55 2020

@author: Andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import yaml
import time
import gc

n_halo_max = np.inf


class ModelParams:
    """Read in yml parameters and store them.
    
    Simple class to read in and store parameters from the yml file. Simple
    methods included to print out the parameters etc.
    
    Attributes
    ----------
    yml_filepath : str
        Filepath to the yml file containing the model parameters.
    verbosity : int
        The level of detail needed for debugging messages.
    debug_flag : bool
        To print out debugging messages, true or false. 
    param_dict : dictionary 
        Dictionary containing contents of yml file.
    input_filepath : str
        The filepath to the input HDF5 file.
    output_filepath : str
        The filepath to the output HDF5 file.
    f_baryon : float
        Baryon fraction as allocated in yml file.
    io_nRec : int
        IO buffer size.
    sub_halo : bool
        Whether sub_halo data should be included or not.
    omega_m : float
        Density parameter. 
    
    """
    
    def __init__(self,yml_filepath,verbosity,debug_flag):
        """ Key parameters for the model
        
        Parameters
        ----------
        yml_filepath : str
            Filepath to the yml file containing the model parameters.
        verbosity : int
            The level of detail needed for debugging messages.
        debug_flag : bool
            To print out debugging messages, true or false. 
            
        """
        self.yml_filepath = yml_filepath
        self.verbosity = verbosity
        self.debug_flag = debug_flag
        self.param_dict = yaml.load(open(yml_filepath),Loader=yaml.Loader)
        
        self.input_filepath = self.param_dict['input_files']['graph_file']
        self.output_filepath = self.param_dict['output_files']['halo_file']
        self.f_baryon = self.param_dict['cosmology']['f_baryon']['Value']
        self.io_nRec = self.param_dict['performance']['io_nRec']['Value']
        self.sub_halo = self.param_dict['model_switches']['sub_halo']['Value']
        self.omega_m = self.param_dict ['cosmology']['omega_m']['Value']

    def output_params(self):
        """ Simple method to print out parameters.
        
        May be replaced for inbuilt method Peter mentioned.
            
        Returns 
        -------
        None
        
        """
        for item in self.param_dict:
            print("{:20s}: {}".format(item,self.param_dict[item]))
            
        return None

class HDFProperties:
    """This class contains the attributes and methods for the HDF5 files.
    
    This class is used to store attributes that are applied accross all graphs
    such as the dark matter particle mass and the total amount of graphs. It 
    was also seen useful to store the input/output HDF5 files here alongside 
    methods to open/write/close files.
    
    Attributes
    ----------
    HDF_input_filepath : str
        The filepath for the HDF input file. Usually in yml parameters file.
    HDF_output_filepath : str
        The filepath for the HDF output file. Usually in yml parameters file.
    graph_input_file : :obj: 'File'
        The open HDF5 input file. 
    halo_output_file: :obj: 'File'
        The open output HDF5 file to which the halo data is to be stored.
    dtype_halo : dtype
        Numpy dtype describing the types of the different class attributes.
    halo_output_data : ndarry
        Numpy array of type custom defined dtype_halo. 
    halo_output_dataset : :obj: 'Dataset'
        HDF5 dataset to which the halo data should be stored.
    n_halo : int
        A simple counter for the amount of halos processed.
    no_of_graphs : int
        The amount of graphs within the input HDF5 file.
    part_mass : float
        The dark matter particle mass. Universal for all graphs.
    dtype_subhalo_stores : dtype
        Numpy dtype for the structured array in which subhalo attributes are 
        stored. 

    """
    def __init__(self,model_params):
        """Key variables about the HDF files.
        
        Parameters
        ----------
        HDF_input_filepath : str
            The filepath for the HDF input file. Usually in yml parameters file.
        HDF_output_filepath : str
            The filepath for the HDF output file. Usually in yml parameters file.
        
        """
        self.model_params = model_params
        self.HDF_input_filepath = model_params.input_filepath
        self.HDF_output_filepath = model_params.output_filepath
        self.graph_input_file = None
        self.halo_output_file = None
        
        self.dtype_halo = np.dtype([
                ('graph_ID',np.int32),
                ('snap_ID',np.int32),
                ('halo_ID',np.int32),
                ('catalog_ID',np.int64),
                ('mass',np.float32),
                ('mass_baryon',np.float32),
                ('mass_from_progenitors',np.float32)
                ])
        
        self.dtype_subhalo_stores = np.dtype([
                ('sub_graph_ids', np.int32),
                ('prog_DM_mass', np.float64),
                ('stellar_mass', np.float64),
                ('descended', bool),
                ('SFR', np.float64)
                ])

        self.open_graph_input()
        self.halo_output_data = np.empty(self.model_params.io_nRec,
                                         dtype=self.dtype_halo)
        self.halo_output_dataset = None
        self.open_halo_output()
        self.halo_output_iRec = 0
        self.n_halo = 0
        self.no_of_graphs = len(self.graph_input_file['/nhalos_in_graph/'][:])
        self.part_mass = self.graph_input_file['Header'].attrs["part_mass"]

    def open_graph_input(self):
        """Opens input graph file (HDF5) using h5py.
        
        Returns
        -------
        None
        
        """
        self.graph_input_file =  h5.File(self.HDF_input_filepath,'r')

    
    @staticmethod
    def close_graph_io(open_HDF5_file):
        """Closes an open HDF5 file. Static method to allow general use.
        
        Parameters
        ----------
        open_HDF5_file : :obj: 'File'
            An open HDF5 file.
        
        Returns
        -------
        None
        
        """
        open_HDF5_file.close()
        
        return None


    def open_halo_output(self):
        """ Open halo output HDF5 file and create dataset.
        
        Returns
        -------
        None
        
        """

        # for obj in gc.get_objects():   # Browse through ALL objects
        #     if isinstance(obj, h5.File):   # Just HDF5 files
        #         try:
        #             self.halo_output_file = h5.File(self.HDF_output_filepath,
        #                                             'w')
        #         except:
        #             pass # Was already closed
        
        self.halo_output_file = h5.File(self.HDF_output_filepath,'w')
        self.halo_output_dataset = self.halo_output_file.create_dataset(
            'Halos',(0,),maxshape=(None,),dtype=self.dtype_halo,
            compression='gzip')
        
        return None

    def output_halos(self,halos):
        """Stores attributes of HaloProperties class to structered array.
        
        Parameters
        ----------
        halos : ndarry
            Numpy array of type object. Each element contains a HaloProperties
            class.
        
        Returns
        -------
        None
        
        """
        
        for halo in halos:
            self.halo_output_data[self.halo_output_iRec][
                'graph_ID'] = halo.graph_ID
            self.halo_output_data[self.halo_output_iRec][
                'snap_ID'] = halo.snap_ID
            self.halo_output_data[self.halo_output_iRec][
                'halo_ID'] = halo.halo_ID
            self.halo_output_data[self.halo_output_iRec][
                'catalog_ID'] = halo.catalog_ID
            self.halo_output_data[self.halo_output_iRec]['mass'] = halo.mass
            self.halo_output_data[self.halo_output_iRec][
                'mass_baryon']= halo.mass_baryon
            self.halo_output_data[self.halo_output_iRec][
                'mass_from_progenitors']=halo.mass_from_progenitors
            self.halo_output_iRec+=1
            
            if self.halo_output_iRec == self.model_params.io_nRec: 
                self.flush_output()
                
                
        return None

    def flush_output(self):
        """Resizes the HDF5 dataset and saves generated data.

        Returns
        -------
        None
        
        """
        
        self.halo_output_dataset.resize((self.halo_output_dataset.shape[0]+
                                         self.halo_output_iRec,))
        self.halo_output_dataset[-self.halo_output_iRec:] = \
            self.halo_output_data[:self.halo_output_iRec]
            
        self.halo_output_iRec=0
        return None

    # Not yet implementedv parameters['outputFiles']['galaxyFile']
    @staticmethod
    def open_galaxy_output(filepath):
        """ Opens the output HDF5 file in write mode.
        
        Parameters
        ----------
        filepath : str
            Filepath file should be created at.
        
        Returns
        -------
        :obj: 'File'
            HDF5 file open in write mode.
        
        """
        return h5.File(filepath,'w')

class GraphProperties:
    """A container for all data contained within a single graph
    
    This class consists of data gathered from a single graph. This documentat-
    tion will essentially be copied from Will's. 
    
    
    Attributes
    ----------
    desc_start_index : ndarray of type 'int'
        The starting index (pointer) for each halo’sentries in all descendant 
        halo arrays (i.e.direct_desc_ids,direct_desc_contribution,​ etc.). 
        Entries containing 2**30 have nodescendants.
    direct_desc_contribution : ndarray of type 'int'
        The number of dark matter particles contributed ​to​ each direct
        descendent ​from the halo.
    direct_desc_ids : ndarray of type 'int'
        The descendent halo IDs, extracted usingdesc_start_index​ and ​ndesc​
    direct_prog_contribution : ndarray of type 'int'
        The number of dark matter particlescontributed ​by​ each direct 
        progenitor ​to​ thehalo.
    direct_prog_ids : ndarray of type 'int'
        The progenitor halo IDs, extracted usingprog_start_index​ and ​nprog​.
    generation_id : ndarray of type 'int'
        The ID (or number) associated with eachgeneration. Counting starts
        from the earliest snapshot.
    generation_length : ndarray of type 'int'
        The number of halos in each generation
    generation_start_index : ndarray of type 'int'
        The starting index (pointer) for each host halo generation.
    graph_halo_ids : ndarray of type 'int'
        The internal graph halo ID assigned to ahalo. NOTE: These differ from
        the halocatalog. These run from 0 - N​host​ with theindex equal to the
        value.
    halo_catalog_halo_ids : ndarray of type 'int'
        The halo catalog ID assigned to each halo.
    mean_pos : ndarray of type 'float'
        The mean position of the particles in the halo.
    ndesc : ndarray of type 'int'
        The number of descendants for each halo.
    nparts : ndarray of type 'int'
        The number of dark matter particles in each halo.Well,
    nprog : ndarray of type 'int'
        The number of progenitors for each halo.
    prog_start_index : ndarray of type 'int'
        The starting index (pointer) for each halo’sentries in all progenitor
        halo arrays (i.e.direct_prog_ids,direct_prog_contribution,​ etc.).
        Entries containing 2**30 have nodescendants.
    redshifts : ndarray of type 'float'
        The redshift for each halo in the graph.
    snapshots : ndarray of type 'int'
        The index of the snapshot in the snapshottext file dictated in the 
        param file, for eachhalo.
    sub_desc_start_index : ndarray of type 'int'
        The starting index (pointer) for eachsubhalo’s entries in all
        descendant subhaloarrays (i.e. ​sub_direct_desc_ids,
        sub_direct_desc_contribution,​ etc.).Entries containing 2**30 have 
        nodescendants.
    sub_direct_desc_contribution : ndarray of type 'int'
        The number of dark matter particlescontributed ​to​ each direct
        descendent ​fromthe subhalo.
    sub_direct_desc_ids : ndarray of type 'int'
        The descendent subhalo IDs, extractedusing ​sub_desc_start_index​ 
        and sub_ndesc​.
    sub_direct_prog_contribution : ndarray of type 'int'
        The number of dark matter particlescontributed ​by​ each direct
        progenitor ​to​ thesubhalo.
    sub_direct_prog_ids : ndarray of type 'int'
        The progenitor subhalo IDs, extracted using sub_prog_start_index​ and 
        ​sub_nprog​.
    sub_generation_id : ndarray of type 'int'
        The ID (or number) associated with eachgeneration. Counting starts
        from the earliest snapshot.
    sub_generation_length : ndarray of type 'int'
        The number of subhalos in each generation.
    sub_generation_start_index : indarray of type 'int'nt
        The starting index (pointer) for each subhalogeneration.
    sub_graph_halo_ids : ndarray of type 'int'
        The internal graph subhalo ID assigned to asubhalo. NOTE: These differ 
        from the halocatalog. These run from 0 - N​sut​ with theindex equal to 
        the value.
    sub_mean_pos : ndarray of type 'float'
        The mean position of the particles in the subhalo.
    sub_ndesc : ndarray of type 'int'
        The number of descendents for eachsubhalo.
    sub_nparts : ndarray of type 'int'
        The number of dark matter particles in eachhalo.
    sub_nprog : ndarray of type 'int'
        The number of progenitors for each halo.
    sub_prog_start_index : ndarray of type 'int'
        The starting index (pointer) for eachsubhalo’s entries in all 
        progenitor subhaloarrays (i.e. ​sub_direct_prog_ids,
        sub_direct_prog_contribution,​ etc.).Entries containing 2**30 have
        no descendants.
    sub_redshifts : ndarray of type 'float'
        The redshift for each subhalo in the graph.
    sub_snapshots : ndarray of type 'int'
        The index of the snapshot in the snapshottext file dictated in the
        param file, for each subhalo.
    subhalo_catalog_halo_ids : ndarray of type 'int'
        The subhalo catalog ID assigned to each subhalo.
    velocity_dispersion_3D : ndarray of type 'float'
        The velocity dispersion for each halo in the graph.
    half_mass_radius : ndaarray of type 'float'
        The half mass radius for each halo in the graph.
    half_mass_velocity_radius : ndaarray of type 'float'
        The half mass velocity radius for each halo in the graph.
    mean_vel : ndaarray of type 'float'
        An array 3 by N_halos long. It contains the x,y,z velocities for each
        halo in the graph.
    rms_radius : ndarray of type 'float'
        The rms radius of each halo within the graph.
    v_max : ndarray of type 'float'
        The maximum velocity of the halo within the graph.
    sub_halo : bool
        From yml input file. Sub halos or no sub halos basically. 
    sub_velocity_dispersion_3D : ndarray of type 'float'
        The velocity dispersion for each subhalo in the graph.
    sub_half_mass_radius : ndaarray of type 'float'
        The half mass radius for each subhalo in the graph.
    sub_half_mass_velocity_radius : ndaarray of type 'float'
        The half mass velocity radius for each subhalo in the graph.
    sub_mean_vel : ndaarray of type 'float'
        An array 3 by N_halos long. It contains the x,y,z velocities for each
        subhalo in the graph.
    sub_rms_radius : ndarray of type 'float'
        The rms radius of each subhalo within the graph.
    sub_v_max : ndarray of type 'float'
        The maximum velocity of the subhalo within the graph.
    subhalo_start_index : ndarray of type 'int'
        Starting index for the subhalo within the graph. Indexed with the halo
        ID.
    host_halos : ndarray of type 'int'
        The host halo ID for the sub halos. n_subhalos long,indexed with sub-
        halo ID.

    """
    
    def __init__(self,graph_ID,open_HDF_file,model_params,part_mass):
        """ Opening HDF datasets for a graph and saving them to the class.
        
        Parameters
        ----------
        graph_ID : int
            The ID of the graph to be opened.
        open_HDF_file : :obj: 'File'
            Open HDF5 file that is to be read from.
        model_params : obj: 'Class'
            ModelParams class object.  
        part_mass : float
            The mass of the dark matter particles. Attribute in HDF.
        
        """
        self.graph_ID = graph_ID
        
        open_HDF_group = open_HDF_file[str(graph_ID)]
        
        self.desc_start_index = open_HDF_group['desc_start_index'][:]
        self.direct_desc_contribution = \
            open_HDF_group['direct_desc_contribution'][:]*part_mass
        self.direct_desc_ids = open_HDF_group['direct_desc_ids'][:]
        self.direct_prog_contribution = \
            open_HDF_group['direct_prog_contribution'][:]*part_mass
        self.direct_prog_ids = open_HDF_group['direct_prog_ids'][:]
        
        # This is a temporary fix until Will reruns his code. Swap after
        self.generation_length = open_HDF_group['generation_length'][:]
        self.generation_id = np.arange(len(self.generation_length),
                                       dtype=int) #open_HDF_group['generation_id'][:]
    
        self.generation_start_index = \
            open_HDF_group['generation_start_index'][:]
        

        self.halo_catalog_halo_ids = open_HDF_group['halo_catalog_halo_ids'][:]
        self.mean_pos = open_HDF_group['mean_pos'][:]
        self.ndesc = open_HDF_group['ndesc'][:]
        self.nparts = open_HDF_group['nparts'][:]
        self.nprog = open_HDF_group['nprog'][:]
        self.prog_start_index = open_HDF_group['prog_start_index'][:]
        self.redshifts = open_HDF_group['redshifts'][:]
        self.snapshots = open_HDF_group['snapshots'][:]
        
        self.velocity_dispersion_3D = open_HDF_group['3D_velocity_dispersion'][:]
        self.half_mass_radius = open_HDF_group['half_mass_radius'][:]
        self.half_mass_velocity_radius = \
            open_HDF_group['half_mass_velocity_radius'][:]
        self.mean_vel = open_HDF_group['mean_vel'][:]
        self.rms_radius = open_HDF_group['rms_radius'][:]
        self.v_max = open_HDF_group['v_max'][:]
        
        self.sub_halo = model_params.sub_halo
        
        self.graph_halo_ids = np.arange(len(self.nparts))

        self.n_halos_in_graph = open_HDF_group.attrs["nhalos_in_graph"]
        
        # Add to docs
        
        
        
        
        if self.sub_halo:
            # This is temporary until Will uses sub halo stuff. 
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            try:
                self.n_subhalos = open_HDF_group['nsubhalos'][:]
                self.sub_desc_start_index = \
                    open_HDF_group['sub_desc_start_index'][:]
                self.sub_direct_desc_contribution = \
                    open_HDF_group['sub_direct_desc_contribution'][:]*part_mass
                self.sub_direct_desc_ids = open_HDF_group['sub_direct_desc_ids'][:]
    
                self.sub_direct_prog_contribution = \
                    open_HDF_group['sub_direct_prog_contribution'][:]*part_mass
                self.sub_direct_prog_ids = open_HDF_group['sub_direct_prog_ids'][:]
                self.sub_generation_id = open_HDF_group['sub_generation_id'][:]
                self.sub_generation_length = \
                    open_HDF_group['sub_generation_length'][:]
                self.sub_generation_start_index = \
                    open_HDF_group['sub_generation_start_index'][:]
                
                self.sub_mean_pos = open_HDF_group['sub_mean_pos'][:]
                self.sub_ndesc = open_HDF_group['sub_ndesc'][:]
                self.sub_nparts = open_HDF_group['sub_nparts'][:]
                self.sub_nprog = open_HDF_group['sub_nprog'][:]
                self.sub_prog_start_index = \
                    open_HDF_group['sub_prog_start_index'][:]
                self.sub_redshifts = open_HDF_group['sub_redshifts'][:]
                self.sub_snapshots = open_HDF_group['sub_snapshots'][:]
                self.subhalo_catalog_halo_ids = \
                    open_HDF_group['subhalo_catalog_halo_ids'][:]
    
                self.sub_velocity_dispersion_3D = \
                    open_HDF_group['sub_3D_velocity_dispersion'][:]
                self.sub_half_mass_radius = \
                    open_HDF_group['sub_half_mass_radius'][:]
                self.sub_half_mass_velocity_radius = \
                    open_HDF_group['sub_half_mass_velocity_radius'][:]
                self.sub_mean_vel = open_HDF_group['sub_mean_vel'][:]
                self.sub_rms_radius = open_HDF_group['sub_rms_radius'][:]
                self.sub_v_max = open_HDF_group['sub_v_max'][:]
                self.subhalo_start_index = open_HDF_group['subhalo_start_index'][:]    
                self.host_halos = open_HDF_group['host_halos'][:]
                
                self.sub_graph_halo_ids = np.arange(len(self.sub_nparts))
                
            except KeyError:
                self.n_subhalos = 2**30
                self.sub_desc_start_index = 2**30
                self.sub_direct_desc_contribution = 2**30
                self.sub_direct_desc_ids = 2**30
                self.sub_direct_prog_contribution = 2**30
                self.sub_direct_prog_ids = 2**30
                self.sub_generation_id = 2**30
                self.sub_generation_length = 2**30
                self.sub_generation_start_index = 2**30
                self.sub_graph_halo_ids = 2**30
                self.sub_mean_pos = 2**30
                self.sub_ndesc = 2**30
                self.sub_nparts = 2**30
                self.sub_nprog = 2**30
                self.sub_prog_start_index = 2**30
                self.sub_redshifts = 2**30
                self.sub_snapshots = 2**30
                self.subhalo_catalog_halo_ids = 2**30
                self.sub_velocity_dispersion_3D = 2**30
                self.sub_half_mass_radius = 2**30
                self.sub_half_mass_velocity_radius = 2**30
                self.sub_mean_vel = 2**30
                self.sub_rms_radius = 2**30
                self.sub_v_max = 2**30
                self.subhalo_start_index = 2**30
                self.host_halos = 2**30
            
class HaloProperties:
    """A container for the properties needed for each halo.
   
    No sophisticated methods, it just truncates the GraphProperites class to 
    ensure data from the current generation is selected. A simple method is 
    included to  truncate the correct properties for subhalos and store them in
    a structured numpy array.
    
    Attributes
    ----------
    graph_ID : str
        The graph_ID (from HDF5 group).
    snap_ID : int
         The snapshot ID currently being processed.
    halo_ID : int
        The halo ID of currently being processed.
    catalog_ID : int
        The ID of the halo corresponding to the original catalog.
    mass : int
        Mass of halo. Amount of dark matter particles * mass of particle.
    nprog : int
        The amount of progenitors the halo has.
    prog_start : int
        The index at which this halo's progenitors start.
    prog_end : int 
        The index at which this halo's progenitors end.
    prog_ids : narray
        Numpy array of the progenitor IDs for the halo.
    prog_mass : ndarry
        Numpy array of the progenitor mass contributions.
    ndesc : int
        Amount of descendent halos.
    desc_start : int
        The index at which this halo's descendents start.
    desc_end : int
        The index at which this halo's descendents end.
    desc_ids : ndarray of type 'int'
        Numpy array of the halo's descendent IDs of type.
    desc_mass : ndarray of type 'int'
        Numpy array of the descendent mass contributions.
    mass_baryon : float
        Mass of Baryons within the halo.
    mass_from_progenitors : float 
        Total mass of all the progenitor halos.
    mass_baryon_from_progenitors : float
        Total mass of all the Baryons contained within the progenitor halos.
    inclusive_contribution : float
        The amount of mass that goes 'missing' when a halo descends.
    done : bool
        Whether or not the halo has been processed.
    subhalo_start_index : int
        The index at which the subhalos, corresponding to the main halo, start.
    n_subhalo : int
        The amount of subhalos in the host halo.
    sub_graph_halo_ids : ndarray of type 'int'
        The subhalo IDs that are contained within the main halo.
    sub_halo_attrs : ndarray of type 'dtype_subhalo_stores'
        The subhalo properties. A structured numpy arrary to store all the 
        properties of a single subhalo.
    """   
    def __init__(self,graph_ID, snap_ID, halo_ID, graph_properties, part_mass,
                 dtype_subhalo_stores):
       """Clipping graph properties to the correct generation for halo use.
    
       Parameters
       ----------
       graph_ID : str
           The graph_ID (from HDF5 group).
       snap_ID : int
           The snapshot ID currently being processed.
       halo_ID : int
           The halo ID of currently being processed.
       graph_properties : :obj: 'Class'
           An instance of GraphProperties.
       part_mass : int
           The dark matter particle mass.
           
       """
       self.graph_ID = graph_ID
       self.snap_ID = snap_ID
       self.halo_ID = halo_ID
       self.catalog_ID = graph_properties.halo_catalog_halo_ids[halo_ID]
       self.mass = graph_properties.nparts[halo_ID]*part_mass
       
      
       self.nprog = graph_properties.nprog[halo_ID]
       self.prog_start = graph_properties.prog_start_index[halo_ID]
       self.prog_end = self.prog_start + self.nprog
       self.prog_ids = graph_properties.direct_prog_ids[self.prog_start:
                                                        self.prog_end]
       self.prog_mass =  graph_properties.direct_prog_contribution[
                                       self.prog_start:self.prog_end]
           
       self.ndesc = graph_properties.ndesc[halo_ID]
       self.desc_start = graph_properties.desc_start_index[halo_ID]
       self.desc_end = self.desc_start + self.ndesc
       self.desc_ids = graph_properties.direct_desc_ids[self.desc_start:
                                                        self.desc_end]
       self.desc_mass = graph_properties.direct_desc_contribution[
                                       self.desc_start:self.desc_end]
       
       self.mean_pos = graph_properties.mean_pos[halo_ID]
       
       # Add in to documentation 
       
       self.velocity_dispersion_3D = graph_properties.velocity_dispersion_3D[
                                     halo_ID]
       self.mean_vel = graph_properties.mean_vel[halo_ID]
       
       self.rms_radius = graph_properties.rms_radius[halo_ID]
       
       #add to docs 
       self.redshift = graph_properties.redshifts[halo_ID]
       self.stellar_mass = 0.
       
       self.mass_baryon = 0.
       self.mass_from_progenitors = 0.
       self.mass_baryon_from_progenitors = 0.
       self.inclusive_contribution = 0.
       self.done = False
       
       
       
       
       if graph_properties.sub_halo and type(graph_properties.n_subhalos) \
                                                                   != int:

          
           self.subhalo_start_index = graph_properties.subhalo_start_index[
                                                                   halo_ID]
           
           self.n_subhalo = graph_properties.n_subhalos[halo_ID]
           
           
           self.sub_halo_attrs = np.empty(self.n_subhalo, 
                                          dtype = dtype_subhalo_stores)
           
           self.initialize_sub_halo_attrs() # Currently starts all halos with 1.
           
           
           self.sub_graph_halo_ids = graph_properties.sub_graph_halo_ids[
                       self.subhalo_start_index: 
                       self.subhalo_start_index + 
                       self.n_subhalo]  
               
           self.sub_halo_attrs[:]['sub_graph_ids'] = self.sub_graph_halo_ids  
               
           if graph_properties.n_subhalos[self.halo_ID] > 0:
               

               central_sub_halo_ID_pos, central_sub_halo_dist  = \
                                          self.find_central_galaxy(self.mean_vel,
                                          self.velocity_dispersion_3D,
                                          self.mean_pos,
                                          self.rms_radius,
                                          graph_properties.sub_mean_vel[
                                              self.sub_graph_halo_ids],
                                          graph_properties.sub_velocity_dispersion_3D[
                                              self.sub_graph_halo_ids],
                                          graph_properties.sub_mean_pos[
                                              self.sub_graph_halo_ids],
                                          graph_properties.sub_rms_radius[
                                              self.sub_graph_halo_ids]
                                          )
                                          
               # The central sub-halo dist will be used as a variable.
               # Keep in the code.
            
               self.central_galaxy_ID = self.sub_graph_halo_ids[
                                                       central_sub_halo_ID_pos]

    
    @staticmethod
    def find_central_galaxy(vel,rms_vel,pos,rms_radius,sub_vel,sub_rms_vel,
                          sub_pos,sub_rms_radius):
        """ Using phase space finds the closest galaxy to main halo.
        
        This function follows equation 14 in Will's paper:
            
        https://arxiv.org/pdf/2003.01187.pdf.
        

        Parameters
        ----------
        vel : ndarray of float64
            1x3 array of velocities pertaining to the host halo.
        rms_vel : float
            The root-mean-squared velocity of the host halo.
        pos : ndarry of float64
            1x3 array of the mean position of the host halo.
        rms_radius : float
            The root-mean-squared radius of the host halo.
        sub_vel : ndarray of float64
            1x3 array of velocities pertaining to the subhalo.
        sub_rms_vel : float
            The root-mean-squared velocity of the subhalo.
        sub_pos : ndarry of float64
            1x3 array of the mean position of the subhalo.
        sub_rms_radius : float
            The root-mean-squared radius of the subhalo.

        Returns
        -------
        central_subhalo_ID_pos : int
            Array index of the central subhalo's ID
        float
            The phase space distance between the subhalo and host halo center.

        """
        pos_dists = np.sqrt(np.sum((pos - sub_pos) ** 2,axis = 1))/(rms_radius +
                                                                    sub_rms_radius)
        
        vel_dists = np.sqrt(np.sum((vel - sub_vel) ** 2,axis = 1))/(rms_vel +
                                                                   sub_rms_vel)
        
        total = pos_dists + vel_dists
        
        central_subhalo_ID_pos = np.argmin(total)
        
        return central_subhalo_ID_pos, total[central_subhalo_ID_pos]
    
    def initialize_sub_halo_attrs(self):
        """ Assigns stellar mass of 1 to all new sub-halos.
        
        Gives sub-halos initial stellar mass and a bool flag to check whether
        each galaxy has given its mass to a descendent. 
        
        Returns
        -------
        None.
        
        """
        self.sub_halo_attrs[:]['stellar_mass'] = 1
        self.sub_halo_attrs[:]['descended'] = False
        self.sub_halo_attrs[:]['SFR'] = 0.
        
    
    def collect_galaxy_prog_info(self, graph_properties, part_mass,
                                 array_of_halo_properties):
        
        """ Collects infomation from galaxy progenitors
        
        This method finds the galaxy progenitors and calculates how much 
        stellar mass etc is to be ibtained from them.
        

        Parameters
        ----------
        graph_properties : :obj: 'Class'
            Class of the graph properties. Galaxy information to be extracted
            from here.
        part_mass : float
            The dark matter particle mass.
        array_of_halo_properties : ndarry of type 'Class'
            Numpy array containing the halo-classes previously processed. This
            is so that galaxy stellar mass and other properties can be obtained
            from progenitors.

        Returns
        -------
        None.

        """
        
        for i_galaxy, galaxy_ID in enumerate(self.sub_graph_halo_ids):

            nprog = graph_properties.sub_nprog[galaxy_ID]

            if nprog > 0 :
            
                prog_start_index = graph_properties.sub_prog_start_index[
                                                                    galaxy_ID]
                
                prog_end_index = prog_start_index + nprog
                
                prog_sub_ids = graph_properties.sub_direct_prog_ids[prog_start_index:
                                                                prog_end_index]
                    
                prog_masses, prog_total_mass  = gather_prog_contrib_mass(
                                                    graph_properties, 
                                                    prog_start_index,
                                                    prog_end_index,
                                                    part_mass)
                
                self.sub_halo_attrs[i_galaxy]['prog_DM_mass'] = prog_total_mass 
                

                prog_host_halo_IDs = graph_properties.host_halos[prog_sub_ids]

                prog_halo_properties = array_of_halo_properties[
                                                prog_host_halo_IDs]
   
    
                for prog_halo, prog_sub_id in zip(prog_halo_properties,
                                                  prog_sub_ids):
    
    
    
                    row_indicies = np.where((prog_halo.
                                           sub_halo_attrs['sub_graph_ids'] ==
                                           prog_sub_id))
                    
                    
                    if prog_halo.sub_halo_attrs[row_indicies]['descended']:
                        continue
                    
                    prog_sub_mass = prog_halo.sub_halo_attrs[
                                                row_indicies]['stellar_mass']
                    
                    desc_start_index = graph_properties.sub_desc_start_index[
                                                                 prog_sub_id]
                    
                    main_desc_id = graph_properties.sub_direct_desc_ids[desc_start_index]
        
            
                    
        
                    row_index_this_snap = np.where((self.sub_halo_attrs[
                                                    'sub_graph_ids'] 
                                                    == galaxy_ID)) # main_desc_id
                    
           
                    total_mass = self.sub_halo_attrs[row_index_this_snap
                                             ]['stellar_mass'] + prog_sub_mass
                    
         
        
                    self.sub_halo_attrs[int(row_index_this_snap[0])]['stellar_mass'] = \
                            total_mass
                
                   
                    prog_halo.sub_halo_attrs[int(row_indicies[0])]['descended'] = True

                        
            else: 
                
                self.sub_halo_attrs[i_galaxy]['prog_DM_mass'] = 0. 
                
                
    def add_baryon_mass(self):
        
        self.mass_baryon += np.sum(self.sub_halo_attrs[:]['stellar_mass'])    

        return None
    

    
    @staticmethod
    def Behroozi_formula(a, z, halo_mass):
        
        log_Mh = np.log10(halo_mass)
        
        nu = np.exp(-4 * (a ** 2))
        
        #characteristic  stellar/halo mass ratio
        log_epsilon = -1.777 + (-0.006 * (a - 1) + (-0.000) * z) * nu
        
        #characteristic mass
        log_M1  = 11.514 + (-1.793 * (a - 1) + (-0.251) * z) * nu
        
        log_mass_frac = log_Mh - log_M1
    
        
        log_stellar_mass = log_epsilon + log_M1 + \
                           HaloProperties.powerlaw_Behroozi(log_mass_frac, nu, z, a) + \
                           HaloProperties.powerlaw_Behroozi(0, nu, z, a)
                           
        stellar_mass = 10 ** log_stellar_mass                  
        
        return stellar_mass
    
    @staticmethod
    def powerlaw_Behroozi(x, nu, z, a):
        
        alpha = -1.412 + (0.731* (a - 1)) * nu
        
        delta = 3.508 + (2.608 * (a - 1) + -0.043 * z) * nu
        
        gamma = 0.316 + (1.319 * (a - 1) + 0.279  * z) * nu
    
        f_x = -np.log10((10 ** (alpha * gamma)) + 1) + \
              delta*(((np.log10(1 + np.exp(x))) ** gamma)/  (1 + np.exp(10 ** 
                                                                        (-x))))
    
        return f_x
        
    def add_stellar_mass(self):
        
        if self.sub_graph_halo_ids.size > 0:
        
            z = self.redshift
            
            a = 1 / (1 + z)
            
            Behroozi_stellar_mass = self.Behroozi_formula(a, z, self.mass)
            
            index = np.where((self.sub_halo_attrs['sub_graph_ids'] == 
                    self.central_galaxy_ID))
            
            old_star_mass = self.sub_halo_attrs[index]['stellar_mass']
            
            print(Behroozi_stellar_mass, old_star_mass)
            
            stellar_mass = max(Behroozi_stellar_mass, old_star_mass)
            
            star_mass_delta = stellar_mass - old_star_mass
            
            #Maybe change the 1 in future. Time bwteen snapshots
            self.sub_halo_attrs[int(index[0])]['SFR'] =  star_mass_delta / 1 
        
            self.sub_halo_attrs[int(index[0])]['stellar_mass'] = old_star_mass \
                + star_mass_delta

            self.stellar_mass += (old_star_mass + star_mass_delta)
        
        return None
        
    def calculate_infall(self, f_baryon):
        
        # Sums total bary mass. Only stellar for now but more to add.
        
        halo_old_bary_mass = self.mass_baryon
        
        mass_baryon = max(f_baryon * self.mass, self.mass_baryon)
        
        baryon_mass_delta = mass_baryon - halo_old_bary_mass
        
        if baryon_mass_delta > 0:
            
            self.mass_baryon += baryon_mass_delta
        
        return None
            
        
        
        
    
        
                
class PlotHalos: 
    
    def __init__(self,yml_filepath):
        """Load in parameters from the yml graph file.
        

        Parameters
        ----------
        yml_filepath : str
            Filepath to the yml parameters file for plotting.

        Returns
        -------
        None.

        """
        
        self.param_dict = yaml.load(open(yml_filepath),Loader=yaml.Loader)
        
        self.halo_file = h5.File(self.param_dict['input_files']['halo_file'],
                                 'r')
        self.halo_data = self.halo_file['Halos'][:]
        self.baryon_fraction = self.param_dict['cosmo']['baryon_fraction']
        self.plots = self.param_dict['plots']
        
    
        
        self.graph_min = self.param_dict['graphs']['graph_min']
        self.graph_max = self.param_dict['graphs']['graph_max']
        self.snap_min = self.param_dict['snapshots']['snap_min']
        self.snap_max = self.param_dict['snapshots']['snap_max']
        
        self.filter_halos()
    
    def filter_halos(self):
        """Filter selected halosusing graph and snap IDs.
    
        Returns
        -------
        None.

        """
        self.halo_data = self.halo_data[(self.graph_min <= 
                                         self.halo_data['graph_ID']) &
                                        (self.halo_data['graph_ID'] <=
                                         self.graph_max) &
                                        (self.snap_min <= 
                                         self.halo_data['snap_ID']) &
                                        (self.halo_data['snap_ID'] <=
                                         self.snap_max)] 

    def generate_plots(self):
        """Simple method to create plots.
        
        Method to generate plots. Asthetics are set-up here using the 
        parameters read in earlier.
    
        Returns
        -------
        None.

        """
        for plot_kind in self.plots:   
            
            if self.plots[plot_kind]['show']:   
                
                fig,ax1 = plt.subplots(figsize=self.plots[plot_kind]
                                       ['figsize'])
                ax1.set_title(plot_kind.replace('_', ' ').capitalize())
                ax1.set_yscale(self.plots[plot_kind]['yscale'])
                ax1.set_xscale(self.plots[plot_kind]['xscale'])
                ax1.set_ylabel(self.plots[plot_kind]['ylabel'])
                ax1.set_xlabel(self.plots[plot_kind]['xlabel'])
                
                x,y = self.generate_x_y_data(plot_kind)
    
                ax1.plot(x,y,self.plots[plot_kind]['marker_type'])
                ax1.grid('--',alpha=0.6)
                if self.plots[plot_kind]['save']:
                    plt.savefig(self.plots[plot_kind]['save_path'],
                                dpi = self.plots[plot_kind]['dpi'] )
                    
                plt.show()
                
            
    def generate_x_y_data(self,plot_kind):
        """ Generates the data using formulas.
        
        Based on the requested plots from the yml file, this method uses 
        different formulas to create the desired plot.
        
        Parameters
        ----------
        plot_kind : str
            The type of plot requested form the yml file. E.g. 'baryon_fraction'

        Returns
        -------
        x : ndarray of type 'float'
            The x-axis data for the plot. E.g. mass.
        y : ndarray of type 'float'
            The y-axis data for the plot. Generally created with a simple 
            formula.

        """
        if plot_kind == 'baryon_fraction':
        
            x = self.halo_data['mass']
            y = self.halo_data['mass_baryon']/self.halo_data['mass']
    
            return x,y
    
        return None

    # When more plots are added to them here. Keeping it tidy
        
        
            


def set_baryon_fraction(halo,array_of_halo_properties,f_baryon):
    """ Caclulates the mass of baryons in the halo.
    
    Uses the global f_baryon (baryon fraction) variable to calculate the total
    mass provided by baryons.
    
    Parameters
    ----------
    halo : :obj: 'Class'
        HaloProperties class object
    array_of_halo_properties : ndarry
        Numpy array of HaloProperties classes for all halos this snapshot.
    f_baryon : float
        The baryon fractions from the input parameters.
    
    Returns
    -------
    None
        Halo class object is updated
        
    """
    this_inclusive_contribution = 0.
    
    for prog in halo.prog_ids:
        
        this_prog_contribution = array_of_halo_properties[
                                        prog].inclusive_contribution
        
        this_prog_desc_ids = array_of_halo_properties[prog].desc_ids
        

        this_inclusive_contribution += array_of_halo_properties[prog].mass_baryon \
            * this_prog_contribution[this_prog_desc_ids==halo.halo_ID]
        
    
    halo.mass_baryon = max(f_baryon*halo.mass,this_inclusive_contribution)

    return None

def gather_progenitors(halo,part_mass):
    """Sum the mass from the progenitors.

    Parameters
    ----------
    halo : :obj: 'Class'
        HaloProperties class object
    
    Returns
    -------
    None
        Halo class object is updated
    
    """
    halo.mass_from_progenitors = np.sum(halo.prog_mass)*part_mass

    return None


def inclusive_mass_contribution(halo):
    """Calculate the mass of particles that don't descend.
    
    This functions calculates the true mass contributions to the descendents.
    This is needed as not 100% of the mass is passed down.
    
    Parameters
    ----------
    halo : :obj: 'Class'
        HaloProperties class object
        
    Returns
    -------
    None
    
    """
    # Get progenitor information
    mass = halo.mass
    
    desc_conts = halo.desc_mass  # get the number of particles going to a desc
    
    desc_cont_frac = desc_conts / mass  # convert to fraction of this 
    
    # get the particles that don't move forward
    extra_cont = (mass - np.sum(desc_conts)) * desc_cont_frac  
    
    # combine particles that move forward and don't
    desc_truemass_conts = (desc_conts + extra_cont) / mass  
    
    halo.inclusive_contribution = desc_truemass_conts
        
    return None



def gather_prog_contrib_mass(graph_properties, prog_start_index, 
                             prog_end_index, part_mass):
    
    """ Calcualtes the DM mass inherited from galaxy progenitors 
    

    Parameters
    ----------
    graph_properties : :obj: 'Class'
        Graph properties class containing all information for a graph.
    prog_start_index : int
        Progenitor start index for current galaxy.
    prog_end_index : int
        Progenitor end index for current galaxy.
    part_mass : float
        Dark matter particle mass.

    Returns
    -------
    prog_masses : ndarry 
        The DM mass contributions from the progenitors.
    prog_total_mass : float
        The total DM mass inherited from all progenitors.

    """
    
    prog_masses =  graph_properties.sub_direct_prog_contribution[
                   prog_start_index:prog_end_index]
    
    prog_total_mass = np.sum(prog_masses) * part_mass
    
    
        
    return prog_masses, prog_total_mass
        

