#include <functional>
#include <vector>
#include <string>

using namespace std;

//Parameters
double mu, _mu;
double sigma, _sigma;
double T;
double energy_old, energy_new;
double meps;
double seps;

//Random
Random rnd;

//Simulation
int nstep;
int nblock = pow(10,2);
int L;
  
vector<double> ave(nblock);
vector<double> av2(nblock);
vector<double> sum_prog(nblock);
vector<double> su2_prog(nblock);
vector<double> err_prog(nblock);

//Variables
vector<double> H(nstep);
double integral_psi, norm;

//Metropolis
double range;
double xmax, xmin;
double eps;

//Functions

void Input();
void Minimize();
void Data_Blocking();
void Data_Blocking_sampling();
void my_Metropolis(string);

double Mod_Psi(const double &);
double Psi(const double &);
double Psi_second(const double &);
double V(const double &);

double Integral();
double Integral_point();
double Integral_point_sampling(const double &);
double Integral_psi();

double Boltzmann(double);
double beta();

double Pbc(double);
double error(vector<double> &, vector<double> &, int);
void Print(string);





