cimport L_galaxies


def C_update_c_model_params(double Omega, double OmegaLambda, double Hubble, 
                            double G,int ReionizationModel, double zr_recomb,
                            double z0_recomb):
    
    L_galaxies.update_c_model_params(Omega, OmegaLambda, Hubble, G,
                                     ReionizationModel, zr_recomb, z0_recomb)
    
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