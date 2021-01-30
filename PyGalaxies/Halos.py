# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:57:31 2020

@author: Andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import yaml


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

    def __init__(
        self,
        graph_ID,
        snap_ID,
        halo_ID,
        graph_properties,
        part_mass,
        dtype_subhalo_stores,
    ):
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
        self.mass = graph_properties.nparts[halo_ID] * part_mass

        self.nprog = graph_properties.nprog[halo_ID]
        self.prog_start = graph_properties.prog_start_index[halo_ID]
        self.prog_end = self.prog_start + self.nprog
        self.prog_ids = graph_properties.direct_prog_ids[
            self.prog_start : self.prog_end
        ]
        self.prog_mass = graph_properties.direct_prog_contribution[
            self.prog_start : self.prog_end
        ]

        self.ndesc = graph_properties.ndesc[halo_ID]
        self.desc_start = graph_properties.desc_start_index[halo_ID]
        self.desc_end = self.desc_start + self.ndesc
        self.desc_ids = graph_properties.direct_desc_ids[
            self.desc_start : self.desc_end
        ]
        self.desc_mass = graph_properties.direct_desc_contribution[
            self.desc_start : self.desc_end
        ]
        

        self.mean_pos = graph_properties.mean_pos[halo_ID]

        # Add in to documentation

        self.velocity_dispersion_3D = graph_properties.velocity_dispersion_3D[halo_ID]
        self.mean_vel = graph_properties.mean_vel[halo_ID]

        self.rms_radius = graph_properties.rms_radius[halo_ID]

        # add to docs
        self.redshift = graph_properties.redshifts[halo_ID]
        self.total_halo_stellar_mass = 0.0
        self.hot_gas = 1.0
        self.cold_gas = 2.0
        self.ejected_gas = 3.0 

        self.mass_baryon = 0.0
        self.mass_from_progenitors = 0.0
        self.mass_baryon_from_progenitors = 0.0
        self.inclusive_contribution = 0.0
        self.done = False

        if graph_properties.sub_halo and type(graph_properties.n_subhalos) != int:

            self.subhalo_start_index = graph_properties.subhalo_start_index[halo_ID]

            self.n_subhalo = graph_properties.n_subhalos[halo_ID]

            self.sub_halo_attrs = np.empty(self.n_subhalo, dtype=dtype_subhalo_stores)

            self.initialize_sub_halo_attrs()  # Currently starts all halos with 1.

            self.sub_graph_halo_ids = graph_properties.sub_graph_halo_ids[
                self.subhalo_start_index : self.subhalo_start_index + self.n_subhalo
            ]

            self.sub_halo_attrs[:]["sub_graph_ids"] = self.sub_graph_halo_ids

            if graph_properties.n_subhalos[self.halo_ID] > 0:

                (
                    central_sub_halo_ID_pos,
                    central_sub_halo_dist,
                ) = self.find_central_galaxy(
                    self.mean_vel,
                    self.velocity_dispersion_3D,
                    self.mean_pos,
                    self.rms_radius,
                    graph_properties.sub_mean_vel[self.sub_graph_halo_ids],
                    graph_properties.sub_velocity_dispersion_3D[
                        self.sub_graph_halo_ids
                    ],
                    graph_properties.sub_mean_pos[self.sub_graph_halo_ids],
                    graph_properties.sub_rms_radius[self.sub_graph_halo_ids],
                )

                # The central sub-halo dist will be used as a variable.
                # Keep in the code.

                self.central_galaxy_ID = self.sub_graph_halo_ids[
                    central_sub_halo_ID_pos
                ]

    @staticmethod
    def find_central_galaxy(
        vel, rms_vel, pos, rms_radius, sub_vel, sub_rms_vel, sub_pos, sub_rms_radius
    ):
        """Use phase space finds the closest galaxy to main halo.

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
        pos_dists = np.sqrt(np.sum((pos - sub_pos) ** 2, axis=1)) / (
            rms_radius + sub_rms_radius
        )

        vel_dists = np.sqrt(np.sum((vel - sub_vel) ** 2, axis=1)) / (
            rms_vel + sub_rms_vel
        )

        total = pos_dists + vel_dists

        central_subhalo_ID_pos = np.argmin(total)

        return central_subhalo_ID_pos, total[central_subhalo_ID_pos]

    def initialize_sub_halo_attrs(self):
        """Assign stellar mass of 1 to all new sub-halos.

        Gives sub-halos initial stellar mass and a bool flag to check whether
        each galaxy has given its mass to a descendent.

        Returns
        -------
        None.

        """
        self.sub_halo_attrs[:]["stellar_mass"] = 1
        self.sub_halo_attrs[:]["descended"] = False
        self.sub_halo_attrs[:]["SFR"] = 0.0

    def set_baryon_fraction(self, array_of_halo_properties, f_baryon):
        """Caclulates the mass of baryons in the halo.

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
        this_inclusive_contribution = 0.0

        for prog in self.prog_ids:

            this_prog_contribution = array_of_halo_properties[
                prog
            ].inclusive_contribution

            this_prog_desc_ids = array_of_halo_properties[prog].desc_ids

            this_inclusive_contribution += (
                array_of_halo_properties[prog].mass_baryon
                * this_prog_contribution[this_prog_desc_ids == self.halo_ID]
            )

        self.mass_baryon = max(f_baryon * self.mass, this_inclusive_contribution)

        return None

    def gather_progenitors(self, part_mass):
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
        self.mass_from_progenitors = np.sum(self.prog_mass) * part_mass

        return None

    def halo_descendent_props(self):
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
        mass = self.mass
        
    
        desc_conts = self.desc_mass  # get the number of particles going to a desc

        desc_cont_frac = desc_conts / mass  # convert to fraction of this
        

        # get the particles that don't move forward
        extra_cont = (mass - np.sum(desc_conts)) * (desc_cont_frac /
                                                    np.sum(desc_cont_frac))

        # combine particles that move forward and don't
        desc_truemass_conts = (desc_conts + extra_cont) / mass


        self.inclusive_contribution = desc_truemass_conts

        return None

    def calc_halo_DM_descend(self, part_mass):

        if self.nprog > 0:

            self.gather_progenitors(part_mass)

        else:

            self.mass_from_progenitors = 0.0

    def calc_halo_attrs_descend(
        self, part_mass, array_of_halo_properties, halo_descend_attrs
    ):

        for halo_property in halo_descend_attrs:

            this_inclusive_contribution = 0.0

            for prog in self.prog_ids:

                this_prog_contribution = array_of_halo_properties[
                    prog
                ].inclusive_contribution

                this_prog_desc_ids = array_of_halo_properties[prog].desc_ids

                this_inclusive_contribution += (
                    getattr(array_of_halo_properties[prog], halo_property)
                    * this_prog_contribution[this_prog_desc_ids == self.halo_ID]
                )

            setattr(
                self,
                halo_property,
                this_inclusive_contribution + getattr(self, halo_property),
            )

    def collect_galaxy_prog_info(
        self, graph_properties, HDF_properties, array_of_halo_properties,
    ):
        """Collect infomation from galaxy progenitors.

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
        part_mass = HDF_properties.part_mass
        
        for i_galaxy, galaxy_ID in enumerate(self.sub_graph_halo_ids):

            nprog = graph_properties.sub_nprog[galaxy_ID]

            if nprog > 0:

                prog_start_index = graph_properties.sub_prog_start_index[galaxy_ID]

                prog_end_index = prog_start_index + nprog

                prog_sub_ids = graph_properties.sub_direct_prog_ids[
                    prog_start_index:prog_end_index
                ]

                prog_masses, prog_total_mass = gather_prog_contrib_mass(
                    graph_properties, prog_start_index, prog_end_index, part_mass
                )

                self.sub_halo_attrs[i_galaxy]["prog_DM_mass"] = prog_total_mass

                prog_host_halo_IDs = graph_properties.host_halos[prog_sub_ids]

                prog_halo_properties = array_of_halo_properties[prog_host_halo_IDs]

                for prog_halo, prog_sub_id in zip(prog_halo_properties, prog_sub_ids):

                    row_indicies = np.where(
                        (prog_halo.sub_halo_attrs["sub_graph_ids"] == prog_sub_id)
                    )

                    if prog_halo.sub_halo_attrs[row_indicies]["descended"]:
                        continue

                    prog_sub_mass = prog_halo.sub_halo_attrs[row_indicies][
                        "stellar_mass"
                    ]

                    # desc_start_index = graph_properties.sub_desc_start_index[
                    #     prog_sub_id
                    # ]

                    # main_desc_id = graph_properties.sub_direct_desc_ids[
                    #     desc_start_index
                    # ]

                    row_index_this_snap = np.where(
                        (self.sub_halo_attrs["sub_graph_ids"] == galaxy_ID)
                    )  # main_desc_id

                    total_mass = (
                        self.sub_halo_attrs[row_index_this_snap]["stellar_mass"]
                        + prog_sub_mass
                    )

                    self.sub_halo_attrs[int(row_index_this_snap[0])][
                        "stellar_mass"
                    ] = total_mass

                    prog_halo.sub_halo_attrs[int(row_indicies[0])]["descended"] = True

            else:

                self.sub_halo_attrs[i_galaxy]["prog_DM_mass"] = 0.0
                
            HDF_properties.n_subhalo += 1

    def add_baryon_mass(self):

        self.mass_baryon += np.sum(self.sub_halo_attrs[:]["stellar_mass"])

        return None

    @staticmethod
    def Behroozi_formula(a, z, halo_mass):

        log_Mh = np.log10(halo_mass)

        nu = np.exp(-4 * (a ** 2))

        # characteristic  stellar/halo mass ratio
        log_epsilon = -1.777 + (-0.006 * (a - 1) + (-0.000) * z) * nu

        # characteristic mass
        log_M1 = 11.514 + (-1.793 * (a - 1) + (-0.251) * z) * nu

        log_mass_frac = log_Mh - log_M1

        log_stellar_mass = (
            log_epsilon
            + log_M1
            + HaloProperties.powerlaw_Behroozi(log_mass_frac, nu, z, a)
            + HaloProperties.powerlaw_Behroozi(0, nu, z, a)
        )

        stellar_mass = 10 ** log_stellar_mass

        return stellar_mass

    @staticmethod
    def powerlaw_Behroozi(x, nu, z, a):

        alpha = -1.412 + (0.731 * (a - 1)) * nu

        delta = 3.508 + (2.608 * (a - 1) + -0.043 * z) * nu

        gamma = 0.316 + (1.319 * (a - 1) + 0.279 * z) * nu

        f_x = -np.log10((10 ** (alpha * gamma)) + 1) + delta * (
            ((np.log10(1 + np.exp(x))) ** gamma) / (1 + np.exp(10 ** (-x)))
        )

        return f_x
    
    

    def add_stellar_mass(self):

        if self.sub_graph_halo_ids.size > 0:

            z = self.redshift

            a = 1 / (1 + z)
        

            Behroozi_stellar_mass = self.Behroozi_formula(a, z, self.mass)

            index = np.where(
                (self.sub_halo_attrs["sub_graph_ids"] == self.central_galaxy_ID)
            )

            old_star_mass = self.sub_halo_attrs[index]["stellar_mass"]

            stellar_mass = max(Behroozi_stellar_mass, old_star_mass)

            star_mass_delta = stellar_mass - old_star_mass

            # Maybe change the 1 in future. Time bwteen snapshots
            self.sub_halo_attrs[int(index[0])]["SFR"] = star_mass_delta / 1

            self.sub_halo_attrs[int(index[0])]["stellar_mass"] = (
                old_star_mass + star_mass_delta
            )

            self.total_halo_stellar_mass += old_star_mass + star_mass_delta

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
    def __init__(self, yml_filepath):
        """Load in parameters from the yml graph file.

        Parameters
        ----------
        yml_filepath : str
            Filepath to the yml parameters file for plotting.

        Returns
        -------
        None.

        """
        self.param_dict = yaml.load(open(yml_filepath), Loader=yaml.Loader)

        self.halo_file = h5.File(self.param_dict["input_files"]["halo_file"], "r")
        self.halo_data = self.halo_file["Halos"][:]
        self.baryon_fraction = self.param_dict["cosmo"]["baryon_fraction"]
        self.plots = self.param_dict["plots"]

        self.graph_min = self.param_dict["graphs"]["graph_min"]
        self.graph_max = self.param_dict["graphs"]["graph_max"]
        self.snap_min = self.param_dict["snapshots"]["snap_min"]
        self.snap_max = self.param_dict["snapshots"]["snap_max"]

        self.filter_halos()

    def filter_halos(self):
        """Filter selected halosusing graph and snap IDs.

        Returns
        -------
        None.

        """
        self.halo_data = self.halo_data[
            (self.graph_min <= self.halo_data["graph_ID"])
            & (self.halo_data["graph_ID"] <= self.graph_max)
            & (self.snap_min <= self.halo_data["snap_ID"])
            & (self.halo_data["snap_ID"] <= self.snap_max)
        ]

    def generate_plots(self):
        """Generate plots of chosen variable.

        Method to generate plots. Asthetics are set-up here using the
        parameters read in earlier.

        Returns
        -------
        None.

        """
        for plot_kind in self.plots:

            if self.plots[plot_kind]["show"]:

                fig, ax1 = plt.subplots(figsize=self.plots[plot_kind]["figsize"])
                ax1.set_title(plot_kind.replace("_", " ").capitalize())
                ax1.set_yscale(self.plots[plot_kind]["yscale"])
                ax1.set_xscale(self.plots[plot_kind]["xscale"])
                ax1.set_ylabel(self.plots[plot_kind]["ylabel"])
                ax1.set_xlabel(self.plots[plot_kind]["xlabel"])

                x, y = self.generate_x_y_data(plot_kind)

                ax1.plot(x, y, self.plots[plot_kind]["marker_type"])
                ax1.grid("--", alpha=0.6)
                if self.plots[plot_kind]["save"]:
                    plt.savefig(
                        self.plots[plot_kind]["save_path"],
                        dpi=self.plots[plot_kind]["dpi"],
                    )

                plt.show()

    def generate_x_y_data(self, plot_kind):
        """Generate the data using formulas.

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
        if plot_kind == "baryon_fraction":

            x = self.halo_data["mass"]
            y = self.halo_data["mass_baryon"] / self.halo_data["mass"]

            return x, y

        return None

    # When more plots are added to them here. Keeping it tidy


def gather_prog_contrib_mass(
    graph_properties, prog_start_index, prog_end_index, part_mass
):
    """Calcualtes the DM mass inherited from galaxy progenitors.

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
    prog_masses = graph_properties.sub_direct_prog_contribution[
        prog_start_index:prog_end_index
    ]

    prog_total_mass = np.sum(prog_masses) * part_mass

    return prog_masses, prog_total_mass
