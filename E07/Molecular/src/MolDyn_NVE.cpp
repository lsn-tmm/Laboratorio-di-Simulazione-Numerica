/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include <vector>
#include "MolDyn_NVE.h"

void my_system(std::string const &s) { // Utilizzare i comandi system() con string in ingresso
  int a = std::system(s.c_str());
  if (a == 1) std::cerr << "not ok";
}


using namespace std;

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

  my_system("./ " + file + "/clean.sh");

  int eq_step = 0;
  int restart = 0;
  if  (file == "solid") {eq_step = 5*pow(10,3); restart = 3;}
  else if (file == "liquid") {eq_step = 5*pow(10,3); restart = 5;}
  else if (file == "gas") {eq_step = pow(10,4); restart = 5;}

  for (int i=0;i<restart;i++) {
    if (i == 0){
      cout << endl << "FIRST INPUT" << endl;
      Input();    //Inizialization
    }
    else{
      cout << endl << "RESTART N. " << i << endl;
      Restart(); //Equilibration
    }
    int nconf = 1;
    for(int istep=1; istep <= eq_step; ++istep){
      Move();           //Move particles with Verlet algorithm
      if(istep%100 == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
        Measure_T();     //Properties measurement
	//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
      }
    }
    ConfFinal();         //Write final configuration to restart
    ConfOld();           //Write final -1 configuration to restart

  }

  cout << endl << "Temperatura equilibrata: inizio simulazione. " << endl;

  nblk = 20;
  nstep = 5*pow(10,3);

  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
	Measure();
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


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  beta = 1.0/temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

  /*
//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  */

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
  my_system("cp config.fcc config.0");
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  //Evaluate potential energy and virial of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  //cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  //cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  //cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

void Restart(void){
  ifstream ReadOld,ReadConf;
  double Temp;

  //cout << endl << "RESTART" << endl;

  //Read initial configuration with final and old files
  cout << "Read initial configuration from file config.final " << endl;
  cout << "Read initial configuration from file old.final " << endl << endl;
  ReadConf.open(file + "/config.final");
  ReadOld.open(file + "/old.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    ReadOld >> xold[i] >> yold[i] >> zold[i];

    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();
  ReadOld.close();

  //One step Verlet algorithm ( compute v(t) and r(t+dt) )
  Move();

  //Kinetic energy
  double t = 0;
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  Temp = (2.0 / 3.0) * t/(double)npart; //Temperature

  //cout << "Temp : " << Temp << endl;
  if ( T == 0 ){
    if  (file == "solid") T = 0.8;
    else if (file == "liquid") T = 1.1;
    else if (file == "gas") T = 1.2;

    cout << "Equilibrating at the temperature: " << T << endl << endl;
  }


  double fs = sqrt(T/Temp);   // fs = velocity scale factor
  // cout << "fs : " << fs << endl;

  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);

  }

  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}


void Measure_T(){
  double t = 0;
  ofstream Temp;

  Temp.open(file + "/eq_temp.dat",ios::app);

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature

   Temp << stima_temp << endl;
   Temp.close();

   return;

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
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;

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
  fstream input, output_A;

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



void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << endl << "Print final configuration to file config.final " << endl;
  WriteConf.open(file + "/config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOld(void){ //Write final -1 configuration
  ofstream WriteConf;

  cout << "Print old configuration to file old.final " << endl << endl;
  WriteConf.open(file + "/old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
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

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double error(double* AV, double* AV2, int n){

  if (n == 0)
    return 0;
  else
    return sqrt((AV2[n] - pow(AV[n],2))/(n) );

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
