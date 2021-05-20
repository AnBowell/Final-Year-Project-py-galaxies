# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:57:31 2020.

@author: Andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import yaml
from LGalaxies.L_Galaxies import C_get_metaldependent_cooling_rate, C_do_reionization, C_star_formation


class HaloProperties:
    """A container for the properties and methods needed for each halo.
    
    This class extracts the data needed for a single halo from the 
    GraphProperties class. It also contains the methods needed to process
    the haloes such as infall and cooling routines. Each method within has
    its own documentation.

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
    done : bool
        Whether or not the halo has been processed.
    subhalo_start_index : int
        The index at which the subhalos, corresponding to the main halo,
        start.
    n_subhalo : int
        The amount of subhalos in the host halo.
    sub_graph_halo_ids : ndarray of type 'int'
        The subhalo IDs that are contained within the main halo.
    mean_pos : ndarray of type 'float'
        The mean position (x,y,z) of the halo.
    velocity_dispersion_3D : float
        The 3D velocity dispersion of the halo.
    mean_vel : ndarry of type 'float'
        The mean velocity (V_x,V_y,V_z) of the halo.
    rms_radius : float
        The root mean squared radius of the halo
    redshift : float
        The redshift of the halo.
    total_halo_baryon_mass : float
        The total amount of baryons within the halo. E.g. hot gas + 
        stellar and cold gas mass from any subhaloes within.
    hot_gas_mass : float
        The mass of hot gas within the halo.
    hot_gas_temp : float
        The temperature of the hot gas within the halo.
    gas_metalicity : float
        The metalicity of the gas. Currently constant at 1e-4.
    cold_gas : float
        The mass of cold gas 
    metal_dependent_cooling_rate : float
        The cooling rate of the hot gas in ergs cm^3 / s. 
    intracluster_stellar_mass : float
        The stellar mass within the halo that does not belong to any
        subhaloes. This occurs when a subhalo has no descendents.
    R_200 : float
        The radius of the halo's sphere which encloses an overdensity
        of 200 times that of the critical density.
    central_galaxy_ID : int
        The ID of the central subhalo if there is one. If not,
        it is filled with 2**30, the fill value for the model.
    """

    def __init__(
        self,
        graph_ID,
        halo_ID,
        graph_properties,
        dtype_subhalo_stores,
        model_params
    ):
        """Extract halo properties from GraphProperties for halo halo_ID

        Parameters
        ----------
        graph_ID : str
            The graph_ID (from HDF5 group).
        halo_ID : int
            The halo ID of currently being processed.
        graph_properties : :obj: 'Class'
            An instance of GraphProperties.
        dtype_subhalo_stores : numpy dtype
            The dark matter particle mass.
        model_params : Instance of ModelParams Class
            An instance of the model params class which contains all of the 
            global parameters.
        """
        self.graph_ID = graph_ID
        self.snap_ID = None             # Not filled in intialisation 
        self.halo_ID = halo_ID
        
        self.catalog_ID = graph_properties.halo_catalog_halo_ids[halo_ID]
        self.mass = graph_properties.mass[halo_ID] 
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
        self.velocity_dispersion_3D = graph_properties.velocity_dispersion_3D[halo_ID]
        self.mean_vel = graph_properties.mean_vel[halo_ID]

        self.rms_radius = graph_properties.rms_radius[halo_ID]
        self.redshift = graph_properties.redshifts[halo_ID]

        self.total_halo_baryon_mass = 0.0
        self.hot_gas_mass = 0.0
        self.hot_gas_temp = 0.0
        self.gas_metalicity = 1e-4

        self.metal_dependent_cooling_rate = 0.0
        self.intracluster_stellar_mass = 0.0

        self.subhalo_start_index = graph_properties.subhalo_start_index[halo_ID]
        self.n_subhalo = graph_properties.n_subhalos[halo_ID]
        self.sub_graph_halo_ids = graph_properties.sub_graph_halo_ids[
            self.subhalo_start_index : self.subhalo_start_index + self.n_subhalo
        ]
        self.done = False
        
        # Find the radius of the haloes overdensity.
        self.calc_radius_200(self.redshift, model_params)

        # If there are subhalos, find the central one
        if self.n_subhalo > 0:

            central_sub_halo_ID_pos, central_sub_halo_dist = self.find_central_galaxy(
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
            # The central sub-halo dist will be used as a variable
            # in future versions of the code.
            # If close enough declare it a central halo -
            # currently 2**30 as a fill 
            
            if central_sub_halo_dist < 2**30:
                    
                self.central_galaxy_ID = self.sub_graph_halo_ids[
                        central_sub_halo_ID_pos]
            else:
                self.central_galaxy_ID = 2**30
                
        else:
            self.central_galaxy_ID = 2**30
   

    @staticmethod
    def find_central_galaxy(
        vel, rms_vel, pos, rms_radius, sub_vel, sub_rms_vel, sub_pos, sub_rms_radius
    ):
        """Use phase space finds the closest galaxy to main halo.

        This function follows equation 14 in Roper et al. Paper:

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


    def calc_radius_200(self, z, model_params):
        """ Calculate the spherical overdensity radius of the halo
        
        Using equation S1 from the L-Galaxies 2015 model description,
        https://lgalaxiespublicrelease.github.io/Hen15_doc.pdf ,
        calculate the radius of the halo.

        Parameters
        ----------
        z : float
            Redshift of the halo.
        model_params : Instance of class ModelParams
            An instance of the class ModelParams where the global cosmology 
            parameters have been stored.

        Returns
        -------
        None.
        
        """
        Hz = model_params.H0 * (((model_params.omega_m * ((1+z) ** 3)) +
                                 model_params.omega_lambda) ** (1/2))
        
        
        self.R_200 = ((self.mass * model_params.SolMass_Mpc_G) /
                      (100 * (Hz**2))) ** (1 / 3)

        return None
            
    # Halo properties descend in proportion to their mass
    
    def calc_halo_props_descend(self, part_mass, list_of_halo_properties,
                                halo_descend_attrs):
         
        
        if self.ndesc > 0:

            desc_cont_frac = self.desc_mass  / self.mass  # convert to fraction of this
            
    
            # get the particles that don't move forward
            extra_cont = (self.mass - np.sum(self.desc_mass )) * (desc_cont_frac /
                                                        np.sum(desc_cont_frac))
    
            # combine particles that move forward and don't
            proportional_contribution = (self.desc_mass  + extra_cont) / self.mass
            
            for (descendent_id, prop_contrib) in zip(self.desc_ids, 
                                                     proportional_contribution):
    
                for halo_property in halo_descend_attrs:
                
                    to_descend = getattr(self, halo_property) * prop_contrib
                    
                    
                    desc_halo = list_of_halo_properties[descendent_id]
                    
                    setattr(desc_halo, halo_property,
                            to_descend + getattr(desc_halo, halo_property))
            
        

    
    
 

    

    def sum_baryon_and_topup(self, halo_descend_attrs, sub_halo_descend_attrs,
                             list_of_subhalo_properties,f_baryon):
        
#        old_total_baryonic_mass = self.total_halo_baryon_mass

        for halo_baryon_prop in halo_descend_attrs:       
             self.total_halo_baryon_mass += getattr(self, halo_baryon_prop)
        
        for subhalo_baryon_prop in sub_halo_descend_attrs:
            for subhalo_ID in self.sub_graph_halo_ids:
                self.total_halo_baryon_mass += getattr(list_of_subhalo_properties[subhalo_ID], 
                                                       subhalo_baryon_prop)
                

                
                
        old_total_baryonic_mass = self.total_halo_baryon_mass

        true_baryon_frac = f_baryon * self.calculate_reionization_frac(
                                                             self.mass,
                                                             self.redshift)
        
        
        self.total_halo_baryon_mass = max(true_baryon_frac * self.mass, 
                                           self.total_halo_baryon_mass)
        


        # if self.total_halo_baryon_mass > self.mass:
        #     print(self.halo_ID)
        
        change_in_baryonic_mass = self.total_halo_baryon_mass - old_total_baryonic_mass
         
        self.hot_gas_mass += change_in_baryonic_mass

         
         
    def calculate_SFR_hot_gas_used(self, dt, list_of_subhalo_properties):

        if self.central_galaxy_ID != 2**30:
            
            main_subhalo = list_of_subhalo_properties[self.central_galaxy_ID]
         
            z = self.redshift

            a = 1 / (1 + z)

            Behroozi_stellar_mass = self.Behroozi_formula(a, z, self.mass)

            old_star_mass = main_subhalo.stellar_mass

            stellar_mass = max(Behroozi_stellar_mass, old_star_mass)


            main_subhalo.stellar_mass = stellar_mass
            
            #self.hot_gas_mass -= star_mass_delta
        
        return None


    @staticmethod
    def Behroozi_formula(a, z, halo_mass):
        """ Fitted equation from Behroozi et al. 2013 describing star formation.
        
        https://arxiv.org/abs/1207.6105
        This function takes in the halo mass, redshift, and scale factor,
        and then returns the expected stellar mass in the halo. This is an
        equation that has been fitted using real world data so the 'magic
        numbers' found below are arbitrary and can be found in the paper.

        Parameters
        ----------
        a : float
            Scale factor. 1/(1+z)
        z : float
            Redshift.
        halo_mass : float
            Dark matter mass of the halo in solar masses.

        Returns
        -------
        stellar_mass : float
            The expected stellar mass (in solar masses) for a halo of a 
            halo_mass size
        """
        log_Mh = np.log10(halo_mass)
        nu = np.exp(-4 * (a ** 2))
        
        # Characteristic stellar/halo mass ratio
        log_epsilon = -1.777 + ((-0.006 * (a - 1) + (-0.000) * z) * nu) + (-0.119 * (a - 1))
        # Characteristic mass
        log_M1 = 11.514 + (-1.793 * (a - 1) + ((-0.251) * z)) * nu

        log_mass_frac = log_Mh - log_M1

        log_stellar_mass = (
            log_epsilon
            + log_M1
            + HaloProperties.powerlaw_Behroozi(log_mass_frac, nu, z, a)
            - HaloProperties.powerlaw_Behroozi(0, nu, z, a)
        )
        stellar_mass = 10 ** log_stellar_mass

        return stellar_mass

    @staticmethod
    def powerlaw_Behroozi(x, nu, z, a):
        """ Powerlaw formula from behroozi et al. 
        
        This function is used many times in the main berhoozi formula - look
        at Behroozi_formula function documentation for details.
        
        Parameters
        ----------
        x : float
            Log mass fraction.
        nu : float
            Constant from main Behroozi formula.
        z : float
            Redshift of the halo.
        a : float
            Scale factor 1/(1+z).

        Returns
        -------
        f_x : float
            Some constant to be used in the Behroozi formula.

        """

        alpha = -1.412 + (0.731 * (a - 1)) * nu

        delta = 3.508 + (2.608 * (a - 1) + (-0.043 * z)) * nu

        gamma = 0.316 + (1.319 * (a - 1) + (0.279 * z)) * nu

        f_x = -np.log10((10 ** (alpha * x)) + 1) + delta * (
            ((np.log10(1 + np.exp(x))) ** gamma) / (1 + np.exp(10 ** (-x)))
        )

        return f_x
    
    

    
    @staticmethod
    def calculate_reionization_frac(mass, redshift):
        """ Calculates modification fraction for baryon fraction.
        
        L-Galaxies routine. Recipe from Gnedin 2000, with parameters/eqs from
        either Okamoto er al 2008 or Kravtsov et al 2004. This depends upon
        the model chosen (0 or 1 in parameters).      
        
        Parameters
        ----------
        mass : float
            Mass of the halo
        redshift : float
            Redshift of the halo.

        Returns
        -------
        reionization_fraction : float
            Fraction to multiply baryon fraction by.

        """
        
        reionization_fraction = C_do_reionization(mass, redshift)
    
        return reionization_fraction
    
    
    
    
    
    
    @staticmethod
    def calculate_virial_circular_velocity(z, H0, omega_m, omega_lambda, r0):
        """Calculates the cirular vellocity when at equilibrium.
        
        This formula is taken from White and Frenk 1991
        https://ui.adsabs.harvard.edu/abs/1991ApJ...379...52W/abstract
        
        This will then allow for the calculation of the viral temperature.

        Parameters
        ----------
        z : float
            Redshift of the halo.
        H0 : float
            Global constant - Hubble parameter at t0
        r0 : float
            Currently the root mean square radius. Subject to change 
            Also not 100% sure on the units right now. Maybe Mpc - when multiplied
            by H0 this gives us kms^(-1)

        Returns
        -------
        Vc : float
            Viral circular velocity.
        """
        

        
        Hz =  H0 * (((omega_m * ((1+z) ** 3)) + omega_lambda) ** (1/2))

        Vc = 10 * Hz * r0
        
        return Vc
    
    @staticmethod
    def calculate_virial_temperature(Vc, mu, m_p, k_B):
        """ Calculates the virial temperature of the hot gas.
        
        Takes in the virial circular velocity calculated from White and Frenk
        1991 and then returns the virial temperature using an estimation
        from L-Galaxies 2020 - In the model description. From springel 2001 too.

        
        µm_p\3kB --> is 35.9 ?
        
        Vc in this instant is in Mpc
        
        Parameters
        ----------
        Vc : float
            Virial circular velocity.
        mu : float
            Average amount of particles
        m_p : float
            Mass of a proton - const from parameters.
        k_B : float
            Boltzman constant - const from parameters.
        Returns
        -------
        T_virial : float
            Virial temperature.

        """
        
        virial_temperature_const_ms = ((0.5 * mu * m_p) / k_B)

        
        virial_temperature_const_kms = (virial_temperature_const_ms * 1e6)
        
        T_virial = virial_temperature_const_kms * (Vc ** 2)        

        return T_virial
    
    
    def calculate_hot_gas_temp(self, model_params):
        """Set the hot gas temperature equal to that of the virial temp.
    
        Parameters
        ----------
        H0 : float
            Hubble parameter - constant.
        mu : float
            Average amount of particles
        m_p : float
            Mass of a proton - const from parameters.
        k_B : float
            Boltzman constant - const from parameters.
            
        Returns
        -------
        None.

        """
        
        Vc = self.calculate_virial_circular_velocity(self.redshift,
                                                     model_params.H0, 
                                                     model_params.omega_m,
                                                     model_params.omega_lambda,
                                                     self.R_200)

        T_virial = self.calculate_virial_temperature(Vc, model_params.mu,
                                                     model_params.m_p, 
                                                     model_params.k_B)

        
        self.hot_gas_temp = T_virial
        self.V_200 = Vc

        return None
    
    
    def calculate_metal_dependent_cooling_rate(self):
        """Calculate metal dependent cooling rate from c routine.
        
        Returns
        -------
        None.

        """
        log_metalicity = np.log10(self.gas_metalicity)
        log_temp = np.log10(self.hot_gas_temp)
        
        cooling_rate = C_get_metaldependent_cooling_rate(log_temp,
                                                         log_metalicity)
        
        self.metal_dependent_cooling_rate = cooling_rate # This in in ergs s^-1 cm ^3
        #x = PROTONMASS * BOLTZMANN * temp / lambda; // now this has units sec g/cm^3
        
        
        
        
        return None
    
    def calc_mass_of_hot_gas_cooled(self, mu, m_p, G, beta_prof_ratio,
                                         beta_prof_ratio_arctan,dt,
                                         list_of_subhalo_properties):
        """
        Parameters
        ----------
        mu : TYPE
            DESCRIPTION.
        m_p : TYPE
            DESCRIPTION.
        G : TYPE
            DESCRIPTION.
        beta_prof_ratio : TYPE
            DESCRIPTION.
        beta_prof_ratio_arctan : TYPE
            DESCRIPTION.
        dt : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if self.central_galaxy_ID != 2**30:
            #LAMBDA_RATIO is (n_e/n)*(n_i/n) = 0.25  From L-Gal
            
            # Metal dependent cooling in erg cm^3 / s. 
            
            # Convert measurements to match
            
            rms_radius_km = self.R_200 * 3.086e19  #Mpc --> km
            
            rms_radius_cgs = self.R_200 * 3.086e24 #Mpc --> cm
            
            G_cgs = G * 10e2 # SI --> cgs
    

            m_p_cgs = m_p *  (10e2) #--> kg to g
            
            dt_seconds = dt* 3.154e+7 #Years --> seconds

            dynamical_time_at_edge = rms_radius_km/ self.V_200
         
       
            lambda_ratio = 0.25
            
            f = beta_prof_ratio - beta_prof_ratio_arctan
            

            tau_cool_P = ((20. * G_cgs * ((mu * m_p_cgs * rms_radius_cgs) ** 2) /
                         (lambda_ratio * self.metal_dependent_cooling_rate)) * 
                         ((f ** 2)/(beta_prof_ratio ** 3)))
            
     
        
            
            fg0 = self.hot_gas_mass/self.mass
            
           
            
            dt_ratio =  dt_seconds / dynamical_time_at_edge
            
           
            
            tau_ratio  = dynamical_time_at_edge * fg0 / tau_cool_P
            

            if tau_ratio <= 1:
                fg = fg0/ (1 + tau_ratio * dt_ratio)
            else:
                teq_ratio = np.log(tau_ratio)
                
                if dt_ratio <= teq_ratio:
                    
                    fg = fg0 * np.exp(-dt_ratio)
                    
                else: 
                    
                    fg = fg0 / (tau_ratio * (1 + (dt_ratio - teq_ratio)))
                
        
            cooling_gas = (fg0 - fg) * self.mass

            
            list_of_subhalo_properties[self.central_galaxy_ID].cold_gas_mass +=\
                cooling_gas
       
            if cooling_gas > self.hot_gas_mass:
                print('Cooling gas greater than hot gas in halo {}'.format(self.haloID))

            self.hot_gas_mass -= cooling_gas
            
        return None
    

    
    
    
    
    
    
    
    
     
class SubhaloProperties:
    """ A container for any subhalo properties 
    
    A short class containing the attributes of subhalos as well as
    any methods that act on the class. For example, star formation.
    
    Attributes
    ----------
    
    graph_ID : int
        The graph ID the subhalo belongs to.
    snap_ID : int
        The snapshot the subhalo belongs to.
    host_halo_ID : int
        The host halo the subhalo is within.
    subhalo_ID : int
        The indentification number of the subhalo.
    mean_pos : ndarray of type 'float'
        An array containing the (x,y,z) position of the subhalo.
    redshift : float
        The redshift of the subhalo.
    DM_mass : float
        The dark matter mass of the subhalo
    ndesc : int
        The number of descendents the subhalo has.
    desc_start_index : int
        The index at which this subhalo's descendents start.
    desc_ids : ndarry of type 'int'
        The IDs of any descendants of this subhalo.
    C_stellar_mass : float
        The stellar mass within the subhalo as modelled using the
        L-Galaxies C-routines.
    SFR : float
        The start formation rate of the subhalo. Amount of stars
        formed / the time between snapshots.
    stellar_mass : float
        The stellar mass from the Behroozi calculations.
    cold_gas_mass : float
        The mass of cold gas which has been cooled onto the subhalo.
    
    """
    def __init__(self, graph_ID, host_halo_ID, subhalo_ID, 
                 graph_properties):
        """ Extract the relevant subhalo data from graph properties.
        
        Parameters
        ----------
        graph_ID : int
            The current graph ID.
        host_halo_ID : int
            The ID of the halo the current subhalo resides in.
        subhalo_ID : int
            The ID of the current subhalo.
        graph_properties : Instance of class GraphProperties
            An instance of the graph properties class.

        Returns
        -------
        None.

        """
        self.graph_ID = graph_ID
        self.snap_ID = None # Not in initialisation
        self.host_halo_ID = host_halo_ID
        self.subhalo_ID = subhalo_ID
        self.mean_pos = graph_properties.sub_mean_pos[subhalo_ID]
        self.redshift = graph_properties.sub_redshifts[subhalo_ID]
        self.DM_mass = graph_properties.sub_mass[subhalo_ID]
        
        self.ndesc = graph_properties.sub_ndesc[subhalo_ID]
    
        self.desc_start_index = graph_properties.sub_desc_start_index[subhalo_ID]

        self.desc_ids = graph_properties.sub_direct_desc_ids[
             self.desc_start_index : self.desc_start_index +  self.ndesc
        ]
        
        self.C_stellar_mass  = 0.0
        self.SFR = 0.0
        self.stellar_mass = 0.0
        self.cold_gas_mass = 0.0
        
    def calc_subhalo_props_descend(self,list_of_subhalo_properties,
                                   subhalo_descend_attrs, 
                                   list_of_halo_properties):
        """ Pass on all baryons to a single descendent subhalo.
        
        This function iterates over a list of properties that need to
        be passed down to descendents. It utilises getters and setters 
        to extract values from class attributes while only using a 
        list. This means only the list has to be changed when new 
        properties are added.

        Parameters
        ----------
        list_of_subhalo_properties : list of Subhalo Classes
            The main list of classes for subhalo properties.
        subhalo_descend_attrs : list of type 'str'
            List of baryonic properties (attribute names) that are to
            be passed down to descendents.
        list_of_halo_properties : list of Halo Classes
            The main list of classes for halo properties.

        Returns
        -------
        None.

        """
        if self.ndesc > 0 :
            
            largest_desc_ID = self.desc_ids[0]


            descendent_subhalo = list_of_subhalo_properties[largest_desc_ID]
    
            for subhalo_property in subhalo_descend_attrs:
                
                to_descend = getattr(self, subhalo_property)
        
                setattr(descendent_subhalo, subhalo_property,
                        getattr(descendent_subhalo, subhalo_property) 
                                + to_descend)
    
        # If no descendents stars go to host halo as intracluster light.
        elif self.ndesc ==0:
            
            list_of_halo_properties[self.host_halo_ID].intracluster_stellar_mass += \
                self.stellar_mass
                
        
        
        
    def calculate_stars_formed(self, halo, dt):
        """ Calculate the mass of stars formed using C routine.

        Utilises the Cython C routine to calculate the mass of stars 
        formed in the central subhalo.

        Parameters
        ----------
        halo : Instance of class HaloProperties.
            The halo the host halo resides in. This is currently 
            needed as the galactic disk is an approximation based on
            halo size.
        dt : float
            Time between the last snapshot and the current one.

        Returns
        -------
        None.

        """
        # Calculate the mass of cold gas which has been reheated, and
        # the mass that has been used to form stars.
        return_dict = C_star_formation(halo.V_200, halo.R_200/10, self.cold_gas_mass,
                                        0, dt ,1)

        
        stars_formed = return_dict['Returnstars']
        reheated_mass = return_dict['Return_reheated_mass']
    
        
        self.C_stellar_mass += stars_formed
    
        # Take away the stellar mass from the cold gas.
        self.cold_gas_mass -= (stars_formed + reheated_mass)


        if self.cold_gas_mass < 0:
            print('Subhalo {} has used up too much cold gas'.format(self.subhalo_ID))
        
        star_mass_delta = stars_formed 

        self.SFR = star_mass_delta / dt
        
        
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
