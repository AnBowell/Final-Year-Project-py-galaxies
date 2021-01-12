# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:20:49 2020

@author: Andrew

"""
    



import L_Galaxies

Omega, OmegaLambda, Hubble, G, ReionizationModel = 1, 0.7, 70, 6.67e-11, 1

zr_recomb, z0_recomb = 15, 6

L_Galaxies.C_update_c_model_params(Omega, OmegaLambda, Hubble, 
                                   G, ReionizationModel,
                                   zr_recomb, z0_recomb)

# Run once
L_Galaxies.C_read_cooling_functions() # Reset kernel everytime this is called += in code.

temp = 100
log_z = 1
print(L_Galaxies.C_get_metaldependent_cooling_rate(temp, log_z))
    
    
    
Mvir,Zcurr = 1e24,1

#Run once
L_Galaxies.C_read_reionization()

reionize_frac = L_Galaxies.C_do_reionization(Mvir, Zcurr)

print(reionize_frac)