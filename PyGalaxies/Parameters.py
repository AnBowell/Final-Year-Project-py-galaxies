# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:55:02 2020

@author: Andrew
"""

import yaml

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

        self.timing = self.param_dict["Monitoring"]["Timing"]["Value"]
        self.timing_graph_save_path = self.param_dict["Monitoring"]["Timing"]["graph_save_path"]
        self.timing_data_save_path = self.param_dict["Monitoring"]["Timing"]["timing_data_save_path"]

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