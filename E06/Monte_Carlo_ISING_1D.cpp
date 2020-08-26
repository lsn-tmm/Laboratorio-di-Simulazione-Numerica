/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "Monte_Carlo_ISING_1D.h"

void my_system(std::string const &s) { // Utilizzare i comandi system() con string in ingresso
  int a = std::system(s.c_str());
  if (a == 1) std::cerr << "not ok";
  
}

using namespace std;

int main(){
  
  my_system("./clean.sh");

  cout << "metro (1-metropolis, not 1-gibbs): " ;
  cin >> metro;

  
  if (metro == 1){
    my_system("./Metropolis/clean.sh");
    my_system("./Eq_m/clean.sh");
  }
  else{
    my_system("./Gibbs/clean.sh");
    my_system("./Eq_g/clean.sh");
  }
  
  double points = 15.0;
  double delta = (2.0 - 0.5) / points;

  for(double T=0.5;T<2.01;T=T+delta){

    temp = T;
  
    Input(); //Inizialization

    for(int i=0;i<pow(10,4);i++) Move(metro); //extra steps for equilibration
  
    for(int iblk=0; iblk <= nblk; ++iblk) //Simulation
      {
	Reset(iblk);   //Reset block averages
	for(int istep=1; istep <= nstep; ++istep)
	  {
	    Move(metro);
	    if (iblk != 0){ //First block of steps for equilibration
	      Measure();
	      Accumulate(); //Update block averages
	    }
	  }
	Averages(iblk);   //Print results for current block
      }
    ConfFinal(); //Write final configuration
    Print();
    
    ResetAll();

    h = 0.02;

    Restart();

    for(int i=0;i<pow(10,4);i++) Move(metro); //extra steps for equilibration
  
    for(int iblk=0; iblk <= nblk; ++iblk) //Simulation
      {
	Reset(iblk);   //Reset block averages
	for(int istep=1; istep <= nstep; ++istep)
	  {
	    Move(metro);
	    
	    if (iblk != 0){ //First block of steps for equilibration
	      Measure_M();
	      Accumulate(); //Update block averages
	    }
	  }
	Average_M(iblk);   //Print results for current block
      }
    ConfFinal(); //Write final configuration
    Print_M();
  }

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  // ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  //ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Restart(void){
  ifstream input;

  cout << "Read final configuration from file config.final " << endl << endl;
  input.open("config.final");
  for (int i=0; i<nspin; ++i) input >> s[i];
  input.close();

}


void Move(int metro){
  int o;
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    if(metro==1) //Metropolis
      Metropolis(o);
    else //Gibbs sampling
      Gibbs(o);
  }
 
}


void Metropolis(int o){

  double alpha = std::min( 1., exp( beta * Boltzmann(s[o], o) ) );
    
  if (rnd.Rannyu() <= alpha){
    s[o] = -s[o];
    accepted++;
  }

}

void Gibbs(int o) {

  double p = 1. / ( 1. + exp( beta * Boltzmann(s[o], o) / s[o] ) );

  if (rnd.Rannyu() <= p) s[o] = 1;
  else s[o] = -1;
  
  accepted++;

}

double Boltzmann(int sm, int ip) //sm valore dello spin e ip la sua posizione
{
  double ene = 2 * ( -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm );
  return ene;
}

void Measure()
{
  double u = 0.0, m =0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i){
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  
  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[ix] = pow(m,2);

}


void Measure_M()
{
  double m =0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i) m += s[i];
  
  walker[im] = m;

}




void Reset(int iblk) //Reset block averages
{
   
  if(iblk == 1){
    fill(glob_av.begin(), glob_av.end(), 0);
    fill(glob_av2.begin(), glob_av2.end(), 0);
  }

  fill( blk_av.begin(), blk_av.end(), 0);
  
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


void ResetAll(void){

  fill(walker.begin(), walker.end(), 0);
  fill(glob_av.begin(), glob_av.end(), 0);
  fill(glob_av2.begin(), glob_av2.end(), 0);
  fill(blk_av.begin(), blk_av.end(), 0);

  blk_norm = 0;
  attempted = 0;
  accepted = 0;

}


void Accumulate(void) //Update block averages
{

  for(int i=0; i<n_props; ++i)
      blk_av[i] = blk_av[i] + walker[i];
  
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
  ofstream Ene, Heat, Chi;
   const int wd=12;
   std::string dir;

   attempted = iblk * nspin * nstep;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    if (metro == 1) dir = "EQ_m";
    else dir = "EQ_g";
    
    Ene.open(dir + "/output.ene." + to_string(temp),ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open(dir + "/output.heat." + to_string(temp),ios::app);
    stima_c =  pow(beta,2) * ( blk_av[ic]/blk_norm - pow(stima_u * (double)nspin,2) ) / (double)nspin;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Chi.open(dir + "/output.chi." + to_string(temp),ios::app);
    stima_x = beta * blk_av[ix]/blk_norm/(double)nspin;
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}



void Average_M(int iblk) //Print results for current block
{
    
  ofstream Mag;
   const int wd=12;
   std::string dir;

   attempted = iblk * nspin * nstep;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    if (metro == 1) dir = "EQ_m";
    else dir = "EQ_g";
    
    Mag.open(dir + "/output.mag." + to_string(temp),ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Energy
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    cout << "----------------------------" << endl << endl;
}




void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Print(){
  ofstream Ene, Heat, Mag, Chi;
  std::string dir;

  if (metro == 1) dir = "Metropolis";
  else dir = "Gibbs";
  
  Ene.open(dir + "/output.ene.final",ios::app);
  Heat.open(dir + "/output.heat.final",ios::app);
  Chi.open(dir + "/output.chi.final",ios::app);

  Ene << temp << " " << glob_av[iu]/(double)nblk << " " << err_u << endl;
  Heat << temp << " " << glob_av[ic]/(double)nblk << " " << err_c << endl;
  Chi << temp << " " << glob_av[ix]/(double)nblk << " " << err_x << endl;

  Ene.close();
  Heat.close();
  Chi.close();

}

void Print_M(){
  ofstream Mag;
  std::string dir;

  if (metro == 1) dir = "Metropolis";
  else dir = "Gibbs";
 
  Mag.open(dir + "/output.mag.final",ios::app);
   
  Mag << temp << " " << glob_av[im]/(double)nblk << " " << err_m << endl;
    
  Mag.close();

}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
