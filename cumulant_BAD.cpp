#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <fstream>
//#include <filesystem>

using namespace std;

#include "variables_and_classes_BAD.cpp"


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
double compute_magnetic_susceptibility (double M_av_pow_2, double M_pow_2_av, double T);
double compute_heat_capacity (double H_av_pow_2, double H_pow_2_av, double T);

double compute_energy_per_particle(double H);
double compute_magnetisation_per_spin (double M);
double compute_magnetic_susceptibility_per_particle(double HI);
double compute_heat_capacity_per_particle(double C);

double compute_cumulant (double M_pow_2_av_pow_2, double M_pow_4_av, double T);

double turn_step (Particle matrix[L][L][N]);

int main()
{
    string output_filename = "cumulant_" + to_string (L) + ".dat";
    ofstream output_file_stream(output_filename);
    output_file_stream << "T cumulant_for_L=" << to_string (L) << endl;

    cout << "Completing...\n" << endl;

    for (double T_current = T_start; T_current < T_start+T_dispersion; T_current+=T_step)
    {
        static int time_start = time(nullptr);

        cout << "T = " << showpoint << T_current << " / " << T_start+T_dispersion << "\t\t" << 100*(T_current-T_start)/T_dispersion << "%\t" << (time(nullptr)-time_start)/60 << " minutes passed..." << endl;

        double M_pow_2 = 0;
        double M_pow_4 = 0;

        for (int configuration = 0; configuration < amount_of_configurations; configuration++)
        {
            set_spins (matrix);

            for (int mcs = 0; mcs < mcsteps+mcs_averaging; mcs++)
            {
                if (mcs >= mcsteps)
                {
                    double M_curr = compute_magnetisation(matrix);
                    M_pow_2 += M_curr*M_curr;
                    M_pow_4 += M_curr*M_curr*M_curr*M_curr;
                }

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
                        if (accept_new_configuration(compute_zet(H_before, H_after, T_current)));
                        else matrix[i][j][0].spin = -matrix[i][j][0].spin;
                    }
                }

                ///TEST << mcs << " " << T_current << " " << compute_magnetisation(matrix)/(L*L*N) << "\n";
            }
        }

        double M_pow_2_av_pow_2 = (M_pow_2/(mcs_averaging*amount_of_configurations))*(M_pow_2/(mcs_averaging*amount_of_configurations));
        double M_pow_4_av = M_pow_4/(mcs_averaging*amount_of_configurations);

        double cumulant = compute_cumulant (M_pow_2_av_pow_2, M_pow_4_av, T_current);

        output_file_stream << showpoint << T_current << " " << cumulant << endl;
    }

    cout << "\nDone!";

    output_file_stream.close();

    return 0;
}

/**double turn_step (Particle matrix[L][L][N]) ///ESLI ZAHOTIM VPIHNUT' VSYO V ODNU FUNTSIYU
{
    return 0.0;
}*/

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



double compute_magnetic_susceptibility (double M_av_pow_2, double M_pow_2_av, double T)
{
    double HI = (M_pow_2_av - M_av_pow_2)/(K_BOLZ*T);

    return HI;
}

double compute_heat_capacity (double H_av_pow_2, double H_pow_2_av, double T)
{
    double C = (H_pow_2_av - H_av_pow_2)/(K_BOLZ*T);

    return C;
}

double compute_energy_per_particle(double H)
{
    double H_per_particle = H/(L*L*N);

    return H_per_particle;
}

double compute_magnetisation_per_spin (double M)
{
    double M_per_spin = M/(L*L*N);

    return M_per_spin;
}

double compute_magnetic_susceptibility_per_particle(double HI)
{
    double HI_per_spin = HI/(L*L*N);

    return HI_per_spin;
}

double compute_heat_capacity_per_particle(double C)
{
    double C_per_spin = C/(L*L*N);

    return C_per_spin;
}

double compute_cumulant (double M_pow_2_av_pow_2, double M_pow_4_av, double T)
{
    double cumulant = 1.0 - (1.0/3.0)*(M_pow_4_av/M_pow_2_av_pow_2);

    return cumulant;
}
