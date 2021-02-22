# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 17:47:33 2021

@author: Andrew
"""

import numpy as np
from time import perf_counter
from LGalaxies.L_Galaxies import C_get_metaldependent_cooling_rate, C_do_reionization

def store_halo_data():
    return None


def initialize_halo_properties(array_of_halo_properties,
                          properties_that_descend):
        """Assign stellar mass of 1 to all new sub-halos.

        Gives sub-halos initial stellar mass and a bool flag to check whether
        each galaxy has given its mass to a descendent.

        Returns
        -------
        None.

        """
        for halo_property in properties_that_descend:
    
            array_of_halo_properties[halo_property][:] = 1.
            
        array_of_halo_properties["metalicity_hot_gas"][:] = 1e-4
        #Inplace edit of array
        return 
    
    
def initialize_subhalo_properties(array_of_subhalo_properties,
                          properties_that_descend):
    """Assign stellar mass of 1 to all new sub-halos.

    Gives sub-halos initial stellar mass and a bool flag to check whether
    each galaxy has given its mass to a descendent.

    Returns
    -------
    None.

    """
    for halo_property in properties_that_descend:

        array_of_subhalo_properties[halo_property][:] = 1.
        
    #Inplace edit of array
    return None

    
    
def store_essential_halo_data(array_of_halo_properties, halo_data, graph_ID,
                              part_mass):
    
    array_of_halo_properties['graph_ID'][:] = graph_ID
    array_of_halo_properties['halo_ID'][:] = np.arange(0,len(halo_data))
    array_of_halo_properties['catalog_ID'][:] = halo_data['halo_catalog_halo_ids']
    array_of_halo_properties['mean_pos'][:] = halo_data['mean_pos'] 
    array_of_halo_properties['redshift'][:] = halo_data['redshifts']
    array_of_halo_properties['mass_DM'][:] = halo_data['nparts'] * part_mass
    
    
    return None

def store_essential_subhalo_data(array_of_subhalo_properties, subhalo_data,
                                 part_mass):
    
    array_of_subhalo_properties['subhalo_ID'][:] = np.arange(0,len(subhalo_data))
    array_of_subhalo_properties['host_halo_ID'][:] = subhalo_data["host_halos"]
    array_of_subhalo_properties['mean_pos'][:] = subhalo_data['sub_mean_pos'] 
    array_of_subhalo_properties['redshift'][:] = subhalo_data['sub_redshifts']
    array_of_subhalo_properties['mass_DM'][:] = subhalo_data['sub_nparts'] * part_mass
    
    
    return None


def halo_descendent_props(mass, desc_mass, array_of_halo_properties, halo_ID):
    """Calculate the mass of particles that don't descend.

    This functions calculates the true mass contributions to the descendents.
    This is needed as not 100% of the mass is passed down.

    Parameters
    ----------
    mass : 
        desc
    desc_mass : 
        desc

    Returns
    -------
    None

    """
 

    desc_cont_frac = desc_mass / mass  # convert to fraction of this
    
    
    # get the particles that don't move forward
    extra_cont = (mass - np.sum(desc_mass)) * (desc_cont_frac /
                                                np.sum(desc_cont_frac))

    # combine particles that move forward and don't and save

    array_of_halo_properties['proportional_contribution'][halo_ID] = (desc_mass 
                                                                      + extra_cont) / mass 

    return None



# Optimized compared to the class based code. Runs at end of the loop and pushes
# material directly down to descendents.
def calc_halo_props_descend(array_of_halo_properties, properties_that_descend,
                            halo_data, halo_contrib_data, halo_ID):

    ndesc = halo_data["ndesc"][halo_ID]
         
    if ndesc > 0:
        
        desc_start_index = halo_data["desc_start_index"][halo_ID]
        
        desc_IDs = halo_contrib_data["direct_desc_ids"][desc_start_index:
                                                        desc_start_index +
                                                        ndesc]
        halo_mass = halo_data['nparts'][halo_ID]
        
        desc_mass = halo_contrib_data['direct_desc_contribution'][desc_start_index :
                                                                  desc_start_index +
                                                                  ndesc]
            
        desc_cont_frac = desc_mass / halo_mass  # convert to fraction of this
    
        # get the particles that don't move forward
        extra_cont = (halo_mass - np.sum(desc_mass)) * (desc_cont_frac /
                                                    np.sum(desc_cont_frac))
    
        # combine particles that move forward and don't and save
    
        proportional_contribution = (desc_mass + extra_cont) / halo_mass
    
        for halo_property in properties_that_descend:
            
            array_of_halo_properties[halo_property][desc_IDs] += \
                (proportional_contribution *
                 array_of_halo_properties[halo_property][halo_ID])
        




def sum_baryon_and_topup(array_of_halo_properties, halo_ID, f_baryon,
                         halo_properties_that_descend):
    
    
    total_halo_baryon_mass = sum(array_of_halo_properties[halo_ID][
                                        halo_properties_that_descend])
    
    # WILL NEED TO ADD A LINE HERE TO ADD UP SUB-HALO BAYRONS.
    
    
    mass_halo_DM = array_of_halo_properties["mass_DM"][halo_ID]
    
    
    true_f_baryon = f_baryon * calculate_reionization_frac(mass_halo_DM, 
                                                           array_of_halo_properties[
                                                           "redshift"][halo_ID])
    

    
    # Accrete more baryons to equal the cosmic average.
    new_total_baryon_mass = max(array_of_halo_properties["mass_baryon"][halo_ID],
                                true_f_baryon * mass_halo_DM)
    
    
    # Find the mass of new baryons accumulated - this is the hot gas mass.
    # This cannot be negative.
    change_in_baryonic_mass = new_total_baryon_mass - total_halo_baryon_mass

    
    array_of_halo_properties["mass_baryon"][halo_ID] = new_total_baryon_mass
    
    array_of_halo_properties["mass_hot_gas"][halo_ID] += change_in_baryonic_mass
    
    
    
def calculate_reionization_frac(mass_DM_halo, redshift):
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
        
        reionization_fraction = C_do_reionization(mass_DM_halo, redshift)
    
        return reionization_fraction
    

def calculate_hot_gas_temp(array_of_halo_properties, rms_radii, redshifts, H0,
                           mu, m_p, k_B):
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
    
    
    Vc = calculate_virial_circular_velocity(redshifts, H0, rms_radii)
    
    T_virial = calculate_virial_temperature(Vc, mu, m_p, k_B)
    
    array_of_halo_properties["velocity_virial"][:] = Vc
    array_of_halo_properties["temperature_hot_gas"][:] = T_virial
    
    
    return None


def calculate_virial_circular_velocity(z, H0, r0):
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
    
    
    Vc = 1.67 * ((1 + z) ** (1 / 2)) * H0 * r0
    
    return Vc

    
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


def calculate_metal_dependent_cooling_rate(array_of_halo_properties, halo_ID):
    """Calculate metal dependent cooling rate from c routine.
    
    Returns
    -------
    None.
    
    """
    log_metalicity = np.log(array_of_halo_properties["metalicity_hot_gas"][halo_ID])
    log_temp = np.log(array_of_halo_properties["temperature_hot_gas"][halo_ID])
    
    cooling_rate = C_get_metaldependent_cooling_rate(log_temp,
                                                     log_metalicity)
    
    array_of_halo_properties["metalicity_cooling_rate"] = cooling_rate # This in in ergs s^-1 cm ^3
    
    #x = PROTONMASS * BOLTZMANN * temp / lambda; // now this has units sec g/cm^3
    
    # Density 1/ (r^2) --> Assume core --> integrate, mean density 
    
    return None



def calculate_Berhoozi_stellar_mass(array_of_halo_properties, z):
    """ Fitted equation from Behroozi et al 2013 describing star formation.
    
    https://arxiv.org/abs/1207.6105
    
    This function takes in the halo mass, redhsift and scale factor and
    returns the expected solar mass within the halo. This is a function
    that has been fitted from real world data so the magic numbers found 
    below are arbitrary and can be found in the paper.

    Parameters
    ----------
    a : float
        Scale factor.
    z : float
        Redshift.
    halo_mass : float
        Dark matter mass of the halo.

    Returns
    -------
    stellar_mass : float
        The expected stellar mass for a halo of a halo_mass size.

    """

    a = 1 / (1 + z)
    
    halo_mass_DM = array_of_halo_properties["mass_DM"]
    
    log_Mh = np.log10(halo_mass_DM)

    nu = np.exp(-4 * (a ** 2))

    # characteristic  stellar/halo mass ratio
    log_epsilon = -1.777 + ((-0.006 * (a - 1) + (-0.000) * z) * nu) + (-0.119 * (a - 1))

    # characteristic mass
    log_M1 = 11.514 + (-1.793 * (a - 1) + ((-0.251) * z)) * nu

    log_mass_frac = log_Mh - log_M1

    log_stellar_mass = (
        log_epsilon
        + log_M1
        + powerlaw_Behroozi(log_mass_frac, nu, z, a)
        - powerlaw_Behroozi(0, nu, z, a)
    )

    stellar_mass = 10 ** log_stellar_mass
    
    array_of_halo_properties["mass_berhoozi_stellar"][:] = stellar_mass

    return None


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
        Scale factor of the halo 1/(1+z).

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

    

def calculate_SFR_hot_gas_used(array_of_halo_properties, 
                               array_of_subhalo_properties, dt,
                               halo_ID, no_data_int):
    """ Function to add Behroozi expected stellar mass to existing mass
    
    This function ensures that the stellar mass expected from the Behroozi
    formula is added onto what is already there. There is a check to ensure
    this cannt decrease - i.e if there is already more stellar mass than 
    expected it cannot dissapear. 
    
    The star formation rate is also calculated here and the mass generated
    if any, is taken from the hot gas mass of the halo.
    
    Parameters
    ----------
    dt : float
        Difference between the current and previous snapshot in years.
    
    Returns
    -------
    None.

    """

    central_subhalo_ID = array_of_halo_properties["central_subhalo_ID"][halo_ID]
        
    
    if central_subhalo_ID != no_data_int:
        
        old_mass_stellar = array_of_subhalo_properties["mass_stellar"][central_subhalo_ID]
        
        Berhoozi_mass_stellar = array_of_halo_properties["mass_berhoozi_stellar"][halo_ID]

        mass_stellar = max(Berhoozi_mass_stellar, old_mass_stellar)

        star_mass_delta = mass_stellar - old_mass_stellar

        array_of_subhalo_properties["SFR"][central_subhalo_ID] = star_mass_delta / dt


        array_of_subhalo_properties["mass_stellar"][central_subhalo_ID] = mass_stellar
    
        array_of_halo_properties["mass_hot_gas"][halo_ID] -= star_mass_delta

    return None
    


def find_central_subhalo(array_of_halo_properties, array_of_subhalo_properties,
                        halo_data, subhalo_data, halo_ID, no_data_int, n_subhalo, 
                        subhalo_start_index
   
):
    """Use phase space finds the closest subhalo to main halo.

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
    #vel, rms_vel, pos, rms_radius, sub_vel, sub_rms_vel, sub_pos, sub_rms_radius

    
    vel = halo_data["mean_vel"][halo_ID]
    rms_vel = halo_data["3D_velocity_dispersion"][halo_ID]
    pos = array_of_halo_properties["mean_pos"][halo_ID]
    rms_radius = halo_data["rms_radius"][halo_ID]

    
    # Again, as the IDs are just 0 --> nsubhalos just use np.arange.
    subhalo_ids = np.arange(subhalo_start_index,
                            subhalo_start_index + n_subhalo)
        
    
        
    
    sub_vel = subhalo_data["sub_mean_vel"][subhalo_start_index:
                                           subhalo_start_index + n_subhalo]
        
    sub_rms_vel = subhalo_data["sub_3D_velocity_dispersion"][subhalo_start_index:
                                                             subhalo_start_index +
                                                             n_subhalo]
        
    sub_pos =subhalo_data["sub_mean_pos"][subhalo_start_index:
                                          subhalo_start_index + n_subhalo]
        
    sub_rms_radius = subhalo_data["sub_rms_radius"][subhalo_start_index:
                                                    subhalo_start_index +
                                                    n_subhalo]
          
        
    
    pos_dists = np.sqrt(np.sum((pos - sub_pos) ** 2, axis=1)) / (
        rms_radius + sub_rms_radius
    )

    vel_dists = np.sqrt(np.sum((vel - sub_vel) ** 2, axis=1)) / (
        rms_vel + sub_rms_vel
    )

    total_distance_from_centre = pos_dists + vel_dists
 
    central_subhalo_ID_pos = np.argmin(total_distance_from_centre)
    
    array_of_subhalo_properties["phase_space_dist_halo_centre"][subhalo_ids] = \
        total_distance_from_centre
        
    
    shortest_distance = total_distance_from_centre[central_subhalo_ID_pos]
        

    #If smaller than 2**30 --> fill valuie for now. Will be a cut off for
    # the centre sub halo soon.
    
    if shortest_distance < no_data_int:
                
        array_of_halo_properties["central_subhalo_ID"][halo_ID] = \
            subhalo_ids[central_subhalo_ID_pos]
            
    else:
        
         array_of_halo_properties["central_subhalo_ID"][halo_ID] = \
             no_data_int

    
    return None



def calc_subhalo_props_descend(array_of_halo_properties, 
                               array_of_subhalo_properties, subhalo_ID, subhalo_data,
                               subhalo_contrib_data, subhalo_properties_that_descend):
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

    
    ndesc = subhalo_data["sub_ndesc"][subhalo_ID]
         
    
    if ndesc > 0:
        
        desc_start_index = subhalo_data["sub_desc_start_index"][subhalo_ID]
        
        # To get desc IDs you slice desc_start_index:desc_start_index+ndesc
        # but in this case we just want the first descendent as it's the one
        # with the most particles in common.
        largest_desc_ID = subhalo_contrib_data["sub_direct_desc_ids"][desc_start_index]
        
        for subhalo_property in subhalo_properties_that_descend:
            
            array_of_subhalo_properties[subhalo_property][largest_desc_ID] += \
                 array_of_subhalo_properties[subhalo_property][subhalo_ID]
        

    # New --> Adds any stars that have no galaxy descendent to intracluster mass.
    elif ndesc == 0:
        
        host_halo_ID = array_of_subhalo_properties["host_halo_ID"][subhalo_ID]
        
        array_of_halo_properties["mass_intracluster_light"][host_halo_ID] += \
        array_of_subhalo_properties["mass_stellar"][subhalo_ID]
      
    
   # If ndesc = -1 this is the final generation. There is no need to push the 
   # subhalo properties anywhere
    

        
    