# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:55:02 2020

@author: Andrew
"""

import yaml
import h5py as h5
import numpy as np
from astropy import constants as const
from LGalaxies.L_Galaxies import C_update_c_model_params, C_read_cooling_functions, C_read_reionization, C_check_if_first_call


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

    def __init__(self, yml_filepath, verbosity, debug_flag):
        """Key parameters for the model.

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
        self.param_dict = yaml.load(open(yml_filepath), Loader=yaml.Loader)

        self.input_filepath = self.param_dict["input_files"]["graph_file"]
        self.output_filepath = self.param_dict["output_files"]["halo_file"]
        self.f_baryon = self.param_dict["cosmology"]["f_baryon"]["Value"]
        self.io_nRec = self.param_dict["performance"]["io_nRec"]["Value"]
        self.sub_halo = self.param_dict["model_switches"]["sub_halo"]["Value"]
        self.omega_m = self.param_dict["cosmology"]["omega_m"]["Value"]
        self.omega_lambda = self.param_dict["cosmology"]["omega_lambda"]["Value"]
        self.omega_gamma = self.param_dict["cosmology"]["omega_gamma"]["Value"]
        self.H0 = self.param_dict["cosmology"]['H0']["Value"]
        self.zr_reionization = self.param_dict["cosmology"]['zr_reionization']["Value"]
        self.z0_reionization = self.param_dict["cosmology"]['z0_reionization']["Value"]
        self.mu = self.param_dict["cosmology"]["mu"]["Value"]
        self.beta_prof_ratio = self.param_dict["cosmology"]["beta_prof_ratio"]["Value"]
        self.beta_prof_ratio_atan = np.arctan(self.beta_prof_ratio)
        
        self.timing = self.param_dict["Monitoring"]["Timing"]["Value"]
        self.timing_graph_save_path = self.param_dict["Monitoring"]["Timing"]["graph_save_path"]
        self.timing_data_save_path = self.param_dict["Monitoring"]["Timing"]["timing_data_save_path"]
        

        
        self.reionize_model = self.param_dict["model_switches"]["reionization_model"]["Value"]
        
        self.omega = self.omega_lambda + self.omega_m + self.omega_gamma
        
        self.G = const.G.value
        self.c = const.c.value
        self.k_B = const.k_B.value
        self.m_p = const.m_p.value
        
        self.halo_properties_list()
        self.subhalo_properties_list()
        self.read_input_snapshot_times()
        self.load_in_HDF_metadata()
        
    def output_params(self):
        """Short method to print out parameters.

        May be replaced for inbuilt method Peter mentioned.

        Returns
        -------
        None

        """
        for item in self.param_dict:
            print("{:20s}: {}".format(item, self.param_dict[item]))

        return None
    
    def load_paramters_to_C(self):
        """ Load in global variables to C (L-Galaxies code).
        
        Annoyingly I don't think there is any better way to do this. Parameters
        need to be inputted into the function C_update_c_model_params in the 
        order described below. This stores them in a global data structure to 
        be accessed by any of the C routines. 
        
        Omega, OmegaLambda, Hubble, G, ReionizationModel, zr_recomb, z0_recomb
        
        This also updates which reionisation model to run
            
        0 for reionization recipie described in Gnedin (2000),
        with the fitting from Okamoto et al. 2008 -> Qi(2010)
        
        1 for reionization recipie described in Gnedin (2000),
        using the fitting formulas given by Kravtsov et al (2004) Appendix B,
        used after Delucia 2004
        

        Returns
        -------
        None.

        """
        
        
        
        C_update_c_model_params(self.omega_m, self.omega_lambda, self.H0, 
                                self.G, self.reionize_model, 
                                self.zr_reionization, self.z0_reionization)
        
        return None
    
    
    def read_input_snapshot_times(self):
        
        filepath = 'Input_Params/snapshot_info.txt'
        data = np.loadtxt(filepath).T
        
        self.snap_redshifts = data[2,:]
        
        self.snap_times = data[4,:]
        
        return None
    
    
    def read_in_data_tables_from_c(self):
        """ Read in data tables for c routines.
        
        Method to read in the data tables for the L-galaxies C routines using C.

        
        
        Returns
        -------
        None.

        """
        
        
        if C_check_if_first_call() == 1:
            C_read_cooling_functions()
            C_read_reionization()
        else: 
            print('Static variables in C already assigned. Data not read in again.')
        
        return None
    
    
    def load_in_HDF_metadata(self):
        
        graph_file = h5.File(self.input_filepath, "r")
        
        header_data = graph_file['Header']
        
        self.part_mass = header_data.attrs['part_mass']
        self.no_data_float = header_data.attrs['NO_DATA_FLOAT']
        self.no_data_int = header_data.attrs['NO_DATA_INT']
        
        self.nhalos_in_graph = graph_file['nhalos_in_graph']
        self.nsubhalos_in_graph = graph_file['sub_nhalos_in_graph']
        
        return None
    
    
    def halo_properties_list(self):
    
        self.halo_properties_dtype = np.dtype(
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
                ("mass_berhoozi_stellar", np.float32),
                ("mass_intracluster_light", np.float32)
            ]
        )
        
        
        
        return None
    
    def subhalo_properties_list(self):
    
        self.subhalo_properties_dtype = np.dtype(
            [
                ("host_halo_ID", np.int32),
                ("subhalo_ID", np.int32),
                ("mean_pos", np.float32, 3),
                ("redshift", np.float32),
                ("phase_space_dist_halo_centre", np.float32),
                ("SFR", np.float32),
                ("mass_DM", np.float32),
                ("mass_cold_gas", np.float32),
                ("mass_stellar", np.float32)
            ]
        )
        
        
        
        return None