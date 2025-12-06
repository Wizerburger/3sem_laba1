#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <fstream>
//#include <filesystem>

using namespace std;

#include "variables_and_classes.cpp"


Particle matrix[L][L][N];

void set_spins(Particle matrix[L][L][N]);
int random_particle_number();
int random_layer_number();
bool accept_new_configuration(double z);

void print_spins(Particle matrix[L][L][N]);

double compute_energy (Particle matrix[L][L][N]);
double compute_particle_energy(Particle matrix[L][L][N], int i, int j, int k);
double compute_zet (double H_before, double H_after, double T);
double compute_magnetisation (Particle matrix[L][L][N]);
double compute_average_magnetisation (double M);
double compute_magnetic_susceptibility (double M, double T);
double compute_heat_capacity (double E, double T);

double turn_step (Particle matrix[L][L][N]);

int main()
{
    string output_filename = "data_" + to_string (L) + ".dat";
    ofstream output_file_stream(output_filename);
    output_file_stream << "T H M HI C";

    ofstream TEST ("TESTDATA.txt");
    TEST << "MCS T M" << endl;

    cout << "Completing...\n" << endl;

    for (double T_steps = T_start; T_steps < T_start+T_dispersion; T_steps+=T_step)
    {
        cout << "T = " << showpoint << T_steps << " / " << T_start+T_dispersion << endl;

        double H = 0;
        double M = 0;
        double HI = 0;
        double C = 0;
        double M_av = 0;

        for (int configuration = 0; configuration < amount_of_configurations; configuration++)
        {
            set_spins (matrix);

            for (int mcs = 0; mcs < mcsteps; mcs++)
            {
                double M_curr = compute_magnetisation(matrix);
                M += M_curr;
                M_av += compute_average_magnetisation(M_curr);

                for (int a = 0; a < L*L*N; a++)
                {
                    int i = random_particle_number();
                    int j = random_particle_number();
                    int k = random_layer_number();

                    double H_before = compute_particle_energy (matrix, i, j, k);
                    matrix[i][j][0].spin = -matrix[i][j][0].spin;
                    double H_after = compute_particle_energy (matrix, i, j, k);

                    if (H_after <= H_before);
                    else
                    {
                        if (accept_new_configuration(compute_zet(H_before, H_after, T_steps)));
                        else matrix[i][j][0].spin = -matrix[i][j][0].spin;
                    }
                }



                TEST << mcs << " " << T_steps << " " << compute_magnetisation(matrix)/(L*L*N) << "\n";

                H += compute_energy(matrix);
                HI += compute_magnetic_susceptibility(compute_magnetisation(matrix), T_steps);
                C += compute_heat_capacity(compute_energy(matrix), T_steps);
            }
        }
        double H_conf_av = H/(amount_of_configurations*mcsteps*L*L*N*2);
        double M_conf_av = M_av/(amount_of_configurations*mcsteps);
        double HI_conf_av = HI/(amount_of_configurations*mcsteps);
        double C_conf_av = C/(amount_of_configurations*mcsteps);

        output_file_stream << showpoint << "\n" << T_steps << " " << H_conf_av << " " << M_conf_av << " " << HI_conf_av << " " << C_conf_av;
    }

     cout << "\nDone!";

    return 0;
}

double turn_step (Particle matrix[L][L][N])
{
    return 0.0;
}

void print_spins (Particle matrix[L][L][N])
{
    for (int i = 0; i<L; i++)
    {
        for (int j = 0; j<L; j++)
        {
            if (matrix[i][j][0].spin == 1) cout << "+";
            else cout << "-";
        }
        cout << endl;
    }
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

int random_particle_number()///pандомайзер на выбор частицы
{
    uniform_int_distribution<int> int_distribution(0, L-1);

    int random_particle_number = (int)int_distribution(mt_generator);

    return random_particle_number;
}

int random_layer_number()
{
    uniform_int_distribution<int> int_distribution(0, N-1);

    int random_layer_number = (int)int_distribution(mt_generator);

    return random_layer_number;
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
                H += -J*matrix[i][j][k].spin*matrix[(i+1)%L][j][k].spin; ///(i+1)%L ento granichnie usloviya
                H += -J*matrix[i][j][k].spin*matrix[i][(j+1)%L][k].spin;
                #if N != 1
                H += -J*matrix[i][j][k].spin*matrix[i][j][(k+1)%N].spin;
                #endif // N
            }
        }
    }

    return H;
}

double compute_particle_energy (Particle matrix[L][L][N], int i, int j, int k)
{
    double H = 0.0;

    H += -J*matrix[i][j][k].spin*matrix[(i+1+L)%L][j][k].spin;
    H += -J*matrix[i][j][k].spin*matrix[i][(j+1+L)%L][k].spin;
    H += -J*matrix[i][j][k].spin*matrix[(i-1+L)%L][j][k].spin;
    H += -J*matrix[i][j][k].spin*matrix[i][(j-1+L)%L][k].spin;
    #if N != 1
    H += -J*matrix[i][j][k].spin*matrix[i][j][(k+1+L)%L].spin;
    H += -J*matrix[i][j][k].spin*matrix[i][j][(k-1+L)%L].spin;
    #endif // N

    return H;
}

double compute_zet (double H_before, double H_after, double T)
{
    double z = (double)pow(EXP, (H_after-H_before)/(-T));///GLYANUT' CHTO S PEREGRUJENNOST'YU FUNCTSII POW

    return z;
}

double compute_magnetisation (Particle matrix[L][L][N])
{
    double M = 0.0;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int k = 0; k < N; k++)
            {
                M += matrix[i][j][k].spin;
            }
        }
    }

    return M;
}

double compute_average_magnetisation (double M)
{
    double M_av = 0;

    M_av = M/(L*L*N);

    return M_av;
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
