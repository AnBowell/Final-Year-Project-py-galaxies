# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:58:24 2020

@author: Andrew
"""

import numpy as np
import h5py as h5
import gc


class HDFProperties:
    """Class containing the attributes and methods for the HDF5 files.

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

    def __init__(self, model_params):
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
       
        self.dtype_halo = np.dtype(
            [
                
                ("graph_ID", np.int32),
                ("snap_ID", np.int32),
                ("halo_ID", np.int32),
                ("catalog_ID", np.int64),
                ("central_subhalo_ID", np.int32),
                ("mean_pos", np.float32, 3),
                ("redshift", np.float32),
                ("metalicity_hot_gas", np.float32),
                ("metalicity_cooling_rate", np.float32),
                ("temperature_hot_gas", np.float32),
                ("velocity_virial", np.float32),
                ("mass_DM", np.float32),
                ("mass_baryon", np.float32),
                ("mass_hot_gas", np.float32),
                ("mass_cold_gas", np.float32),
                ("mass_ejected_gas", np.float32),
                ("mass_intracluster_light", np.float32)
            ]
        )

        self.subhalo_dtype = np.dtype(
            [
                ("graph_ID", np.int32),
                ("snap_ID", np.int32),
                ("host_halo_ID", np.int32),
                ("subhalo_ID", np.int32),
                ("mean_pos", np.float32, 3),
                ("redshift", np.float32),
                ("SFR", np.float32),
                ("DM_mass", np.float32),
                ("stellar_mass", np.float32),
                ("cold_gas_mass", np.float32)
            ]
        )


        self.sub_halo_descend_attrs = ["stellar_mass", "AGN_mass"]
        self.halo_output_dataset = None
        self.open_halo_output() # Make sure this is done first as it closes HDF files.
  
        self.open_graph_input()
        
        self.halo_output_data = np.empty(
            self.model_params.io_nRec, dtype=self.dtype_halo
        )
        
        self.subhalo_output_data = np.empty(
            self.model_params.io_nRec, dtype=self.subhalo_dtype
        )
        
       
        self.halo_output_iRec = 0
        self.subhalo_output_iRec = 0
        self.n_halo = 0
        self.n_subhalo = 0
        self.no_of_graphs = len(self.graph_input_file["/nhalos_in_graph/"][:])
        self.snap_redshifts, self.snap_times = self.read_input_snapshot_times()
        
        self.load_in_HDF_metadata()
        

    def open_graph_input(self):
        """Open input graph file (HDF5) using h5py.

        Returns
        -------
        None

        """
        self.graph_input_file = h5.File(self.HDF_input_filepath, "r")

    @staticmethod
    def close_graph_io(open_HDF5_file):
        """Close an open HDF5 file. Static method to allow general use.

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
    
    
    @staticmethod
    def read_input_snapshot_times():
        
        filepath = 'Input_Params/snapshot_info.txt'
        data = np.loadtxt(filepath).T
        
        snap_redshifts = data[2,:]
        
        snap_times = data[4,:]
        
        return snap_redshifts, snap_times


    def open_halo_output(self):
        """Open halo output HDF5 file and create dataset.

        Returns
        -------
        None

        """
  
        
        try:
            self.halo_output_file = h5.File(self.HDF_output_filepath, "w")

        except OSError:
            for obj in gc.get_objects():   # Browse through ALL objects
               if isinstance(obj, h5.File):   # Just HDF5 files
                   try:
                       obj.close()
                   except:
                       pass # Was already closed    
            self.halo_output_file = h5.File(self.HDF_output_filepath, "w")

        self.halo_output_dataset = self.halo_output_file.create_dataset(
            "halo_data", (0,), maxshape=(None,), dtype=self.dtype_halo, compression="gzip"
        )
        
        self.subhalo_output_dataset = self.halo_output_file.create_dataset(
            "subhalo_data", (0,), maxshape=(None,), dtype=self.subhalo_dtype, compression="gzip"
        )
       
        return None

    def output_halos(self, halos):
        """Store attributes of HaloProperties class to structered array.

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
            self.halo_output_data[self.halo_output_iRec]["graph_ID"] = halo.graph_ID
            self.halo_output_data[self.halo_output_iRec]["snap_ID"] = halo.snap_ID
            self.halo_output_data[self.halo_output_iRec]["halo_ID"] = halo.halo_ID
            self.halo_output_data[self.halo_output_iRec]["catalog_ID"] = halo.catalog_ID
            self.halo_output_data[self.halo_output_iRec]["central_subhalo_ID"] = halo.central_galaxy_ID
            self.halo_output_data[self.halo_output_iRec]["mean_pos"] = halo.mean_pos
            self.halo_output_data[self.halo_output_iRec]["redshift"] = halo.redshift
            self.halo_output_data[self.halo_output_iRec]["metalicity_hot_gas"] = halo.gas_metalicity
            self.halo_output_data[self.halo_output_iRec]["metalicity_cooling_rate"] = halo.metal_dependent_cooling_rate
            self.halo_output_data[self.halo_output_iRec]["temperature_hot_gas"] = halo.hot_gas_temp
            self.halo_output_data[self.halo_output_iRec]["velocity_virial"] = halo.Vvir
            self.halo_output_data[self.halo_output_iRec]["mass_DM"] = halo.mass
            self.halo_output_data[self.halo_output_iRec]["mass_baryon"] = halo.total_halo_baryon_mass
            self.halo_output_data[self.halo_output_iRec]["mass_hot_gas"] = halo.hot_gas_mass
            self.halo_output_data[self.halo_output_iRec]["mass_cold_gas"] = halo.cold_gas
            self.halo_output_data[self.halo_output_iRec]["mass_ejected_gas"] = halo.ejected_gas
            self.halo_output_data[self.halo_output_iRec]["mass_intracluster_light"] = halo.intracluster_stellar_mass
            
            self.halo_output_iRec += 1

            if self.halo_output_iRec == self.model_params.io_nRec:
          
                self.flush_output()

        return None

    def output_subhalos(self, subhalos, subhalo_output_list):
        
        for subhalo in subhalos:
            

           

                
           for output_item in subhalo_output_list:
               
               self.subhalo_output_data[self.subhalo_output_iRec][output_item] = \
                   getattr(subhalo,output_item)
                   
           
           self.subhalo_output_iRec += 1
           if self.subhalo_output_iRec == self.model_params.io_nRec:
               
               self.flush_subhalo_output()

        return None


    def flush_subhalo_output(self):
        self.subhalo_output_dataset.resize(
            (self.subhalo_output_dataset.shape[0] + self.subhalo_output_iRec,))
        
        self.subhalo_output_dataset[-self.subhalo_output_iRec:] = \
            self.subhalo_output_data[:self.subhalo_output_iRec]

        self.subhalo_output_iRec = 0
     
        
        return None

    def flush_output(self):
        """Resize the HDF5 dataset and saves generated data.

        Returns
        -------
        None

        """

        self.halo_output_dataset.resize(
            (self.halo_output_dataset.shape[0] + self.halo_output_iRec,)
        )
        self.halo_output_dataset[-self.halo_output_iRec :] = self.halo_output_data[
            : self.halo_output_iRec
        ]
        
        

        self.halo_output_iRec = 0
        return None

    # Not yet implementedv parameters['outputFiles']['galaxyFile']
    @staticmethod
    def open_galaxy_output(filepath):
        """Open the output HDF5 file in write mode.

        Parameters
        ----------
        filepath : str
            Filepath file should be created at.

        Returns
        -------
        :obj: 'File'
            HDF5 file open in write mode.

        """
        return h5.File(filepath, "w")

    def load_in_HDF_metadata(self):
    
        header_data = self.graph_input_file['Header']
        
        self.part_mass = header_data.attrs['part_mass']
        self.no_data_float = header_data.attrs['NO_DATA_FLOAT']
        self.no_data_int = header_data.attrs['NO_DATA_INT']
        
        self.nhalos_in_graph = self.graph_input_file['nhalos_in_graph'][:]
        self.nsubhalos_in_graph = self.graph_input_file['sub_nhalos_in_graph'][:]
        


class GraphProperties:
    """A container for all data contained within a single graph.

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

    def __init__(self, graph_ID, open_HDF_file, model_params, part_mass):
        """Opening HDF datasets for a graph and saving them to the class.

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

        self.desc_start_index = open_HDF_group["desc_start_index"][:]
        self.direct_desc_contribution = (
            open_HDF_group["direct_desc_contribution"][:] * part_mass
        )
        self.direct_desc_ids = open_HDF_group["direct_desc_ids"][:]
        self.direct_prog_contribution = (
            open_HDF_group["direct_prog_contribution"][:] * part_mass
        )
        self.direct_prog_ids = open_HDF_group["direct_prog_ids"][:]

        # This is a temporary fix until Will reruns his code. Swap after
        self.generation_length = open_HDF_group["generation_length"][:]
        self.generation_id = np.arange(
            len(self.generation_length), dtype=int
        )  # open_HDF_group['generation_id'][:]

        self.generation_start_index = open_HDF_group["generation_start_index"][:]

        self.halo_catalog_halo_ids = open_HDF_group["halo_catalog_halo_ids"][:]
        self.mean_pos = open_HDF_group["mean_pos"][:]
        self.ndesc = open_HDF_group["ndesc"][:]
        self.mass = open_HDF_group["nparts"][:] * part_mass
        self.nprog = open_HDF_group["nprog"][:]
        self.prog_start_index = open_HDF_group["prog_start_index"][:]
        self.redshifts = open_HDF_group["redshifts"][:]
        self.snapshots = open_HDF_group["snapshots"][:]

        self.velocity_dispersion_3D = open_HDF_group["3D_velocity_dispersion"][:]
        self.half_mass_radius = open_HDF_group["half_mass_radius"][:]
        self.half_mass_velocity_radius = open_HDF_group["half_mass_velocity_radius"][:]
        self.mean_vel = open_HDF_group["mean_vel"][:]
        self.rms_radius = open_HDF_group["rms_radius"][:]
        self.v_max = open_HDF_group["v_max"][:]

        #self.sub_halo = model_params.sub_halo

        self.graph_halo_ids = np.arange(len(self.mass))

        self.n_halos_in_graph = open_HDF_group.attrs["nhalos_in_graph"]

        # Add to docs

      
        self.n_subhalos = open_HDF_group["nsubhalos"][:]
        self.sub_desc_start_index = open_HDF_group["sub_desc_start_index"][:]
        self.sub_direct_desc_contribution = (
            open_HDF_group["sub_direct_desc_contribution"][:] * part_mass
        )
        self.sub_direct_desc_ids = open_HDF_group["sub_direct_desc_ids"][:]

        self.sub_direct_prog_contribution = (
            open_HDF_group["sub_direct_prog_contribution"][:] * part_mass
        )
        self.sub_direct_prog_ids = open_HDF_group["sub_direct_prog_ids"][:]
        self.sub_generation_id = open_HDF_group["sub_generation_id"][:]
        self.sub_generation_length = open_HDF_group["sub_generation_length"][:]
        self.sub_generation_start_index = open_HDF_group[
            "sub_generation_start_index"
        ][:]

        self.sub_mean_pos = open_HDF_group["sub_mean_pos"][:]
        self.sub_ndesc = open_HDF_group["sub_ndesc"][:]
        self.sub_nparts = open_HDF_group["sub_nparts"][:]
        self.sub_mass = self.sub_nparts * part_mass
        self.sub_nprog = open_HDF_group["sub_nprog"][:]
        self.sub_prog_start_index = open_HDF_group["sub_prog_start_index"][:]
        self.sub_redshifts = open_HDF_group["sub_redshifts"][:]
        self.sub_snapshots = open_HDF_group["sub_snapshots"][:]
        self.subhalo_catalog_halo_ids = open_HDF_group[
            "subhalo_catalog_halo_ids"
        ][:]

        self.sub_velocity_dispersion_3D = open_HDF_group[
            "sub_3D_velocity_dispersion"
        ][:]
        self.sub_half_mass_radius = open_HDF_group["sub_half_mass_radius"][:]
        self.sub_half_mass_velocity_radius = open_HDF_group[
            "sub_half_mass_velocity_radius"
        ][:]
        self.sub_mean_vel = open_HDF_group["sub_mean_vel"][:]
        self.sub_rms_radius = open_HDF_group["sub_rms_radius"][:]
        self.sub_v_max = open_HDF_group["sub_v_max"][:]
        self.subhalo_start_index = open_HDF_group["subhalo_start_index"][:]
        self.host_halos = open_HDF_group["host_halos"][:]

        self.sub_graph_halo_ids = np.arange(len(self.sub_nparts))

