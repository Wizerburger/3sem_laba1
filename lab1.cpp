#include <iostream>
#include <random>
#include <ctime>

using namespace std;

#include "variables_and_classes.cpp"

//mt19937 gen(time(nullptr));

Particle matrix[L][L][N];

void set_spins(Particle matrix[L][L][N]);
bool accept_new_configuration(double z);

double compute_magnetisation (Particle matrix[L][L][N]);
double compute_magnetic_susceptibility (double M, double T);
double compute_heat_capacity (double E, double T);

int main()
{
    for (int configuration = 0; configuration < amount_of_configurations; configuration++)
    {
        set_spins (matrix);
        for (double T_steps = T_start; T_steps < T_start+T_dispersion; T_steps+=0.1)
        {
            for (int mcs = 0; mcs < mcsteps; mcs++)
            {
                double z = 0.0;
                accept_new_configuration(z);
                for (int i = 0; i < L*L*N; i++)
                {
                    cout << i+1 << ": " <<matrix[0][0][i].spin << endl;
                }
            }
        }
    }
    return 0;
}

void set_spins(Particle matrix[L][L][N])
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int k = 0; k < N; k++)
            {
                matrix[i][j][k].spin = 1;
            }
        }
    }
}

bool accept_new_configuration(double z)
{
    uniform_real_distribution<double> real_distribution(0.0, 1.0); ///Устанвливает распределение и диапазон

    double accept_number = (double)real_distribution(mt_generator);

    if (accept_number > z) return false;
    else return true;
}

double compute_magnetisation (Particle matrix[L][L][N])
{
    double M = 0.0;
    for (int i = 0; i < L*L*N; i++)
    {
        M += matrix[0][0][i].spin;
    }

    return M;
}

double compute_magnetic_susceptibility (double M, double T) ///not done
{
    double k;
    double HI = 0.0;

    return HI;
}

double compute_heat_capacity (double E, double T) ///not done
{
    double C = 0.0;

    return C;
}
