#define L 24
#define N 1
#define PI 3.14159265
#define EXP 2.71828183
#define K_BOLZ 1.380649///E-23
///#define T_krit

int amount_of_configurations = 1;
int mcsteps = 100;
double T_start = 0.5;
double T_dispersion = 5.0;

double J = 1.0;

mt19937 mt_generator(time(nullptr));

struct Particle
{
    double x;
    double y;
    double z;

    int spin;
};
