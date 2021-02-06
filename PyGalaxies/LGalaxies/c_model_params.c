#include <model_params.h>

void update_c_model_params(double Omega, double OmegaLambda, double Hubble, 
                           double G,int ReionizationModel, double zr_recomb,
                           double z0_recomb){
                            
modelParams.Omega = Omega;
modelParams.OmegaLambda = OmegaLambda;
modelParams.Hubble = Hubble;
modelParams.G = G;
modelParams.ReionizationModel = ReionizationModel;
modelParams.recombination_zr = zr_recomb;
modelParams.recombination_z0 = z0_recomb;     
};  
