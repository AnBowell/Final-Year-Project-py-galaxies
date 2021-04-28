
void update_c_model_params(double OmegaM, double OmegaLambda, double Hubble, 
                           double G,int ReionizationModel,double zr_recomb,
                           double z0_recomb,float SfrEfficiency, double SfrColdCrit,
                           double EnergySN, double EtaSN, double feedback_reheating_epsilon,
                           double reaheat_pre_vel, double reaheat_slope);
                            
void read_cooling_functions (void);         
                            
double get_metaldependent_cooling_rate (double logTemp, double logZ);
             
void read_reionization (void);

double do_reionization(float Mvir, double Zcurr);

int check_if_first_call (void);

struct ReturnTuple starformation(double Vc, double rad_cold_gas, double mass_cold_gas, 
                     double time, double dt, int nstep);  