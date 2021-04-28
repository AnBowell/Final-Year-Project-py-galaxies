# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:55:02 2020

@author: Andrew
"""

import yaml
import numpy as np
from astropy import constants as const
from astropy import units  as u
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
        self.beta_prof_ratio_arctan = np.arctan(self.beta_prof_ratio)
        self.timing = self.param_dict["Monitoring"]["Timing"]["Value"]
        self.timing_graph_save_path = self.param_dict["Monitoring"]["Timing"]["graph_save_path"]
        self.timing_data_save_path = self.param_dict["Monitoring"]["Timing"]["timing_data_save_path"]
        self.reionize_model = self.param_dict["model_switches"]["reionization_model"]["Value"]
        
        self.omega = self.omega_lambda + self.omega_m + self.omega_gamma
        
        
        astropy_G = const.G
        
        cosmo_units_G = (u.Mpc * (u.km**2))/ (u.solMass * (u.s**2))
        
        
        self.SolMass_Mpc_G  = astropy_G.to(cosmo_units_G).value
        
        
        self.G = astropy_G.value
        
        
        self.c = const.c.value
        self.k_B = const.k_B.value # Be explicit in Units.
        self.m_p = const.m_p.value
        
        self.SFR_efficiency = self.param_dict["star_params"]["SFR_efficiency"]["Value"]

        # This is in 10^10 Msun --> convert
        self.SFR_cold_crit = self.param_dict["star_params"]["SFR_cold_crit"]["Value"]
        self.SFR_cold_crit *= 10 ** 10
        
        
        # Convert to internal units as per L-gal.
        self.EnergySN = float(self.param_dict["SN_params"]["EnergySN"]["Value"])
        
        #1.989e+33 Msun --> g, 3.08568e+24 is Mpc --> cm   3.08568e+24/100000 is Unit time ins  
        # Converting to Msun.Mpc^2(Mpc/Km/s)^-2)
        self.EnergySN =  self.EnergySN / (1.989e+33 * ((3.08568e+24)**2) / ((3.08568e+24/100000)**2))
    
        #Already in correct units.
        self.EtaSN = float(self.param_dict["SN_params"]["EtaSN"]["Value"])

        
        
        self.feedback_reheating_epsilon = self.param_dict["SN_params"]["FeedbackReheatingEpsilon"]["Value"]
        self.reaheat_pre_vel = self.param_dict["SN_params"]["ReheatPreVelocity"]["Value"]
        self.reaheat_slope = self.param_dict["SN_params"]["ReheatSlope"]["Value"]
        
        
        self.halo_descend_attrs = ["hot_gas_mass", "ejected_gas","intracluster_stellar_mass"]
        
        self.subhalo_descend_attrs = ['cold_gas_mass', 'C_stellar_mass']
        
        self.subhalo_output_list = ['graph_ID', 'snap_ID', 'host_halo_ID',
                                    'subhalo_ID','mean_pos','redshift','SFR',
                                    'DM_mass','stellar_mass','cold_gas_mass',
                                    'C_stellar_mass']
        
        self.read_input_snapshot_times()
        
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
                                self.zr_reionization, self.z0_reionization,
                                self.SFR_efficiency, self.SFR_cold_crit,
                                self.EnergySN, self.EtaSN, self.feedback_reheating_epsilon,
                                self.reaheat_pre_vel, self.reaheat_slope)
        
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
    

    def read_input_snapshot_times(self):
        
        try:
            filepath = 'Input_Params/snapshot_info.txt'
            data = np.loadtxt(filepath).T
            
            self.snap_redshifts = data[2,:]
            
            self.snap_times = data[4,:]
        except OSError:
            print('Snapshot info file has not been found. Program will try and continue')
        return None