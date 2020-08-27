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
#include <string>
#include "Monte_Carlo_NVT.h"

using namespace std;

void my_system(std::string const &s) { // Utilizzare i comandi system() con string in ingresso
  int a = std::system(s.c_str());
  if (a == 1) std::cerr << "not ok";
}

int main(){

  fstream read;

  my_system("./clean.sh");

  cout << "Inserire 'solid', 'liquid' or 'gas': " ;
  cin >> file;

  read.open("input." + file);
  if ( read.is_open() ){
    my_system("cp input." + file + " input.dat");
    cout << "Input data initialized. " << endl << endl;
  }
  else {
    cerr << "PROBLEM: File non esistente." << endl;
    return -1;
  }
  read.close();

  Input(); //Inizialization

  Equilibration();
  ConfFinal(); //Write final configuration

  //Correlation();


  my_system("./" + file + "/clean.sh");
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
	//ConfXYZ(nconf);//Write actual configuration in XYZ format
	nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  Print_SI();

  return 0;
}


void Input(void)
{
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

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

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;

  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial

  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

//Evaluate potential energy and virial of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}

void Restart(){
  ifstream ReadConf;

  //Read initial configuration
  cout << "Restart configuration from file config.final " << endl << endl;
  ReadConf.open(file + "/config.final");
  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();


  Measure();

//Print start values for the potential energy and virial
  cout << "Restarting potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;



}


void Equilibration(void){
  ofstream Epot, Pres;
  std::string dir = "eq";
  int steps = 0;

  if (file == "solid") steps = pow(10,3);
  else if (file == "liquid") steps = pow(10,4);
  else if (file == "gas") steps = 5*pow(10,4);

  Epot.open(dir + "/" + dir + "_epot." + file,ios::out);
  Pres.open(dir + "/" + dir + "_press." + file,ios::out);

  for(int i=0; i<steps;i++){ // Equilibration
    Move();

    if( i%(steps/500) == 0){
       Measure();
       //Potential energy per particle
       Epot << i << " " << walker[iv]/(double)npart + vtail << endl;
       //Pressure
       Pres << i << " " << rho * temp + (walker[iw] + ptail * (double)npart) / vol << endl;
    }
  }

  Epot.close();
  Pres.close();

  return;

}


void Correlation(void){
  ofstream Epot, Pres;
  std::string dir = "correlation";

  cout << endl << "Inizio presa dati" << endl;

  Epot.open(dir + "/" + dir + "_epot." + file,ios::out);
  Pres.open(dir + "/" + dir + "_press." + file,ios::out);

  for(int i=0; i<5*pow(10,5);i++){ // Equilibration
    if (i%(int)pow(10,4) == 0) cout << "step:" << i << endl;

    Move();
    Measure();
    //Potential energy per particle
    Epot << i << " " << walker[iv]/(double)npart + vtail << endl;
    //Pressure
    Pres << i << " " << rho * temp + (walker[iw] + ptail * (double)npart) / vol << endl;

  }

  Epot.close();
  Pres.close();

  return;

}


void Move(void)
{
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;


  for(int i=0; i<npart; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);

  //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

  //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

  //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu())
    {
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;

       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

void Measure()
{
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
//INSERT CODE HERE------------------------------------------------------

     walker[ (int)(dr/bin_size) + 2 ] += 2;

//--------------------------------------------------------------------

     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;

  return;
}


void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   double r, gdir;
   ofstream Gofr, Gave, Epot, Pres;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    Epot.open(file + "/output.epot.0",ios::app);
    Pres.open(file + "/output.pres.0",ios::app);
    Gofr.open(file + "/output.gofr.0",ios::app);
    Gave.open(file + "/output.gave.0",ios::app);

    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

    //Potential energy per particle
    Epot  << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
    //Pressure
    Pres  << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

//g(r)
//INSERT CODE HERE--------------------------------------------------------------

    for(int i=2; i<n_props; i++){

      r = bin_size * (i-1);
      double V_r = 4./3. * M_PI * ( pow(r + bin_size,3) - pow(r,3) );

      gdir =  blk_av[i]/blk_norm * pow( rho * npart * V_r ,-1);
      glob_av[i] += gdir;
      glob_av2[i] += gdir*gdir;

      //average value of g(r) in each block
      Gofr << iblk << setw(wd) << r <<  setw(wd)  << glob_av[i]/(double)iblk << endl;

      //final average value of g(𝑟) with statistical uncertainties
      if(iblk == nblk){
	err_gdir=Error(glob_av[i],glob_av2[i],iblk);
	Gave << r << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gdir << endl;
      }

    }

//------------------------------------------------------------------------------

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
}


void Print_SI(void){

  const double K_b = 1.38064852 * pow(10,-23);
  const double eps_A = 120. * K_b;
  const double sigma_A = 0.34;

  vector<string> files = {"epot","pres","gave"};
  fstream input, output_A, output_K;

  for (auto name : files) {

    if (name == "gave") {

      double bin = 0.0, value = 0.0, err = 0.0;

      input.open(file + "/output." + name + ".0", ios::in);
      output_A.open(file + "/argon_" + name + ".out", ios::out);

      for(int i=0; i<nbins; i++){
	input >> bin;
	input >> value;
	input >> err;

	output_A << bin*sigma_A << " " << value*sigma_A << " " << err*sigma_A << endl;
      }
	input.close();
	output_A.close();
    }
    else{
      int blk = 0;
      double value = 0.0, ave_value = 0.0, err = 0.0;

      input.open(file + "/output." + name + ".0", ios::in);
      output_A.open(file + "/argon_" + name + ".out", ios::out);

      for(int i=0; i<nblk; i++){
	input >> blk;
	input >> value;
	input >> ave_value;
	input >> err;

	if (name == "epot"){
	  output_A << blk << " " << value*eps_A << " " << ave_value*eps_A << " " << err*eps_A << endl;
	}
	else if (name == "pres") {
	  output_A << blk << " " << value*eps_A/pow(sigma_A,3) << " " << ave_value*eps_A/pow(sigma_A,3) << " " << err*eps_A/pow(sigma_A,3) << endl;
	}
      }

      input.close();
      output_A.close();

    }

  }

  return;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open(file + "/config.final");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
