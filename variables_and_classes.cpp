#define L 16
#define N 1
#define PI 3.14159265
#define EXP 2.71828183
///#define T_krit

int amount_of_configurations = 1;
int mcsteps = 1/*10000*/;
double T_start = 0.5;
double T_dispersion = 0.1/*5.0*/;

double J = 1.0;

mt19937 mt_generator(time(nullptr));

struct Particle
{
    double x;
    double y;
    double z;

    int spin;
};
