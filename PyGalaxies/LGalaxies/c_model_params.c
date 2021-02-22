#include <model_params.h>

first_call = 1;


void update_c_model_params(double OmegaM, double OmegaLambda, double Hubble, 
                           double G,int ReionizationModel, double zr_recomb,
                           double z0_recomb){
                            
modelParams.OmegaM = OmegaM;
modelParams.OmegaLambda = OmegaLambda;
modelParams.Hubble = Hubble;
modelParams.G = G;
modelParams.ReionizationModel = ReionizationModel;
modelParams.recombination_zr = zr_recomb;
modelParams.recombination_z0 = z0_recomb;     


};  




int check_if_first_call() {

  int to_return;

  if (first_call == 1) {

    first_call = 0;
    to_return = 1;
    
  } else {
    to_return = 0;
  }


return to_return;

}