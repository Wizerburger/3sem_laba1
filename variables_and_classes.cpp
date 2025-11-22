#define L 16
#define N 1
#define PI 3.14159265
#define EXP 2.71828183
///#define T_krit

int amount_of_configurations = 1;
int mcsteps = 10000;
double T_dispersion = 5.0;

//mt19937 gen(time(nullptr));

struct Particle
{
    double x;
    double y;
    double z;

    int spin;
};
