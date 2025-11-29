#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <fstream>
//#include <filesystem>

using namespace std;

#include "variables_and_classes.cpp"

//mt19937 gen(time(nullptr));



Particle matrix[L][L][N];

void set_spins(Particle matrix[L][L][N]);
int random_particle_sequence_number();
bool accept_new_configuration(double z);

double compute_energy (Particle matrix[L][L][N]);
double compute_zet (double H_before, double H_after, double T);
double compute_magnetisation (Particle matrix[L][L][N]);
double compute_magnetic_susceptibility (double M, double T);
double compute_heat_capacity (double E, double T);

int main()
{
    string folder_name = "configurations";

    for (int configuration = 0; configuration < amount_of_configurations; configuration++)
    {
        string output_filename = folder_name + "/configuration_" + to_string(configuration) + ".csv";
        ofstream output_file_stream(output_filename);
        output_file_stream << "T,H,M,HI,C";

        set_spins (matrix);
        for (double T_steps = T_start; T_steps < T_start+T_dispersion; T_steps+=0.1)
        {
            for (int mcs = 0; mcs < mcsteps; mcs++)
            {
                for (int i = 0; i < L*L*N; i++)
                {
                    double H_before = compute_energy (matrix);
                    matrix[0][0][i].spin = -matrix[0][0][i].spin;
                    double H_after = compute_energy (matrix);
                    if (H_after < H_before);
                    else
                    {
                        if (accept_new_configuration(compute_zet(H_before, H_after, T_steps)));
                        else matrix[0][0][i].spin = -matrix[0][0][i].spin;
                    }
                }
            }

            double H = compute_energy(matrix);
            double M = compute_magnetisation(matrix);
            double HI = compute_magnetic_susceptibility(M, T_steps);
            double C = compute_heat_capacity(H, T_steps);
            output_file_stream << showpoint << "\n" << T_steps << "," << H << "," << M << "," << HI << "," << C;


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

bool accept_new_configuration(double z)///pандомайзер с проверкой подходимости
{
    uniform_real_distribution<double> real_distribution(0.0, 1.0); ///USTANAVLIVAET RASPREDELENIE KOTOROE NADO POTOM ISPOLZOVAT'

    double accept_number = (double)real_distribution(mt_generator);

    if (accept_number > z) return false;
    else return true;
}

int random_particle_sequence_number()///pандомайзер на выбор частицы
{
    uniform_int_distribution<int> int_distribution(0, L*L*N-1);

    int random_particle_number = (int)int_distribution(mt_generator);

    return random_particle_number;
}

double compute_energy (Particle matrix[L][L][N])
{
    double H = 0.0;

    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int k = 0; k < N; k++)
            {
                H += J*matrix[i][j][k].spin*matrix[(i+1)%L][j][k].spin; ///(i+1)%L ento granichnie usloviya
                H += J*matrix[i][j][k].spin*matrix[i][(j+1)%L][k].spin;
                #if N != 1
                H += J*matrix[i][j][k].spin*matrix[i][j][(k+1)%N].spin;
                #endif // N
            }
        }
    }

    return H;
}

double compute_zet (double H_before, double H_after, double T)
{
    double z = (double)pow(EXP, (H_after-H_before)/(-K_BOLZ*T));///GLYANUT' CHTO S PEREGRUJENNOST'YU FUNCTSII POW

    return z;
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
    double HI = 0.0;

    return HI;
}

double compute_heat_capacity (double E, double T) ///not done
{
    double C = 0.0;

    return C;
}
