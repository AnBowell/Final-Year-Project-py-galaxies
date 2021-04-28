cimport L_galaxies


def C_update_c_model_params(double Omega, double OmegaLambda, double Hubble, 
                            double G,int ReionizationModel, double zr_recomb,
                            double z0_recomb,float SfrEfficiency, double SfrColdCrit,
                            double EnergySN, double EtaSN, double feedback_reheating_epsilon,
                            double reaheat_pre_vel, double reaheat_slope):
    
    L_galaxies.update_c_model_params(Omega, OmegaLambda, Hubble, G,
                                     ReionizationModel, zr_recomb, z0_recomb,
                                     SfrEfficiency, SfrColdCrit, EnergySN,
                                     EtaSN, feedback_reheating_epsilon,
                                     reaheat_pre_vel, reaheat_slope)
    
    return None

def C_read_cooling_functions():
    
    L_galaxies.read_cooling_functions()
    
    return None

def C_get_metaldependent_cooling_rate(double logTemp, double logZ):

    return L_galaxies.get_metaldependent_cooling_rate(logTemp, logZ)
    

def C_read_reionization():
    
    L_galaxies.read_reionization()
    
    return None

def C_do_reionization(float Mvir, double Zcurr):
    
    return L_galaxies.do_reionization(Mvir, Zcurr)

def C_check_if_first_call():
    
    return L_galaxies.check_if_first_call()


def C_star_formation(double Vc, double rad_cold_gas, double mass_cold_gas, 
                     double time, double dt, int nstep):
    
    return L_galaxies.starformation(Vc, rad_cold_gas, mass_cold_gas, 
                                    time, dt, nstep)