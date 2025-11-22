#include <iostream>
#include <random>
#include <ctime>

#include "variables_and_classes.cpp"

using namespace std;

mt19937 gen(time(nullptr));

Particle matrix[L][L][N];

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

double generate_me_double ()

int accept_configuration_of_random(double z)
{
    ///uniform_real_distribution<double> dist_double(0.0, 1.0); ///Устанвливает распределение и диапазон

    double accept_number = (double)gen()/gen.max();

    cout << "Generared value of chetotam: " << accept_number << endl;
    return 0;
}

int main()
{
    for (int configuration = 0; configuration < amount_of_configurations; configuration++)
    {
        set_spins (matrix);
        for (double T_steps = 0.5; T_steps < T_dispersion; T_steps+=0.1)
        {
            for (int mcs = 0; mcs < mcsteps; mcs++)
            {
                accept_configuration_of_random();
                for (int i = 0; i < L*L; i++)
                {

                }
            }
        }
    }
    return 0;
}
