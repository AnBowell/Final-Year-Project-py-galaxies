struct allModelParameters{

double OmegaM;
double OmegaLambda;
double Hubble;
double G;
int ReionizationModel;
double recombination_zr;
double recombination_z0;
float SfrEfficiency;
double SfrColdCrit;
double EnergySN;
double EtaSN;
double feedback_reheating_epsilon;
double reaheat_pre_vel;
double reaheat_slope;
};

struct allModelParameters modelParams;

struct ReturnTuple{
    double Returnstars, Return_reheated_mass;
};

static int first_call;
