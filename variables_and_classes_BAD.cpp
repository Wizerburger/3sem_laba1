#define L 128
#define N 1
#define PI 3.14159265
#define EXP 2.71828183
#define K_BOLZ 1.0///1.380649///E-23
///#define T_krit

int amount_of_configurations = 100;
int mcsteps = 1000;
int mcs_averaging = 300;
double T_start = 1.5;
double T_dispersion = 2.0;
double T_step = 0.05;

double J = 1.0;

mt19937 mt_generator(time(nullptr));

struct Particle
{
    double x;
    double y;
    double z;

    int spin;
};

///For L = 8 is enough 1000 mcs
///For L = 16 is pretty enough about pure 5000 monte carlo steps
///For L = 32 needed not less than 10000 for highest temperatures
