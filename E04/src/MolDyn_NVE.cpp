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

  cout << "Inserire 'solid', 'liquid' or 'gas': ";
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

  my_system("./" + file + "/clean.sh");

  int eq_step = 0;
  int restart = 0;
  if  (file == "solid") {eq_step = 5*pow(10,3); restart = 3;}
  else if (file == "liquid") {eq_step = 5*pow(10,3); restart = 5;}
  else if (file == "gas") {eq_step = pow(10,4); restart = 5;}

  for (int i=0;i< restart;i++) {
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
        Measure_T();     //Temperature measurement
        nconf += 1;
      }
    }
    ConfFinal();         //Write final configuration to restart
    ConfOld();           //Write final -1 configuration to restart

  }


  //Estendo l'ultima fase dell'equilibrazione
  int extra_step = 0;
  if  (file == "solid") extra_step = 5*pow(10,3);
  else if (file == "liquid") extra_step = 5*pow(10,3);
  else if (file == "gas") extra_step = pow(10,4);

  for(int istep=1; istep <= extra_step; ++istep){
    Move();           //Move particles with Verlet algorithm
    if(istep%100 == 0) cout << "Number of time-steps: " << istep << endl;
    if(istep%10 == 0){
      Measure_T();     //Temperature measurement
    }
  }
  ConfFinal();         //Write final configuration to restart
  ConfOld();           //Write final -1 configuration to restart

  cout << endl << "Temperatura equilibrata: inizio simulazione. " << endl;


  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
    Move();           //Move particles with Verlet algorithm
    if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
    if(istep%10 == 0){
      Measure();     //Properties measurement
      //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
      nconf += 1;
    }
  }
  ConfFinal();         //Write final configuration to restart

  Averages(file);        //Compute average values and errors

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

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //Virial
  n_props = 5; //Number of observables

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

void Measure(){ //Properties measurement
  //int bin;
  double v, t, w;
  double vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open(file + "/output_epot.dat",ios::app);
  Ekin.open(file + "/output_ekin.dat",ios::app);
  Temp.open(file + "/output_temp.dat",ios::app);
  Etot.open(file + "/output_etot.dat",ios::app);
  Press.open(file + "/output_press.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 48.0/3.0 * (1.0/pow(dr,12) - 0.5/pow(dr,6));

//Potential energy
       v += vij;
       w += wij;
     }
    }
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_press = rho * temp + w / vol; //Pressure

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();

    return;
}

void Measure_T(){
  double t = 0.;
  ofstream Temp;

  Temp.open(file + "/eq_temp.dat",ios::app);

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature

   Temp << stima_temp << endl;
   Temp.close();

   return;

}

void Averages(std::string file){
  int N,L;
  if (file == "gas"){
    N = nstep/(10 * 50);
  }
  else if (file == "liquid"){
    N = nstep/(10 * 20);
  }
  else {
    N = nstep/(10 * 10);
  }
  L = nstep/ (10 * N);
  double ave[N];
  double av2[N];
  double* sum_prog = new double[N];
  double* su2_prog = new double[N];
  double err_prog[N];

  const double K_b = 1.38064852 * pow(10,-23);
  const double eps = 120. * K_b;
  const double sigma_A = 0.34 * pow(10,-9);

  vector<string> files = {"etot", "ekin", "epot", "temp", "press"};
  fstream input, output, data;

  for(auto name: files){

    //Inizializzazione variabili di appoggio
    for(int i=0;i<N;i++){
      ave[i] = 0;
      av2[i] = 0;
      sum_prog[i] = 0;
      su2_prog[i] = 0;
      err_prog[i] = 0;
    }

    double sum = 0;
    double k = 0;

    input.open(file + "/output_" + name + ".dat", ios::in); //Serie dati della simulazione
    output.open(file + "/ave_" + name + ".out", ios::out); //Data blocking
    data.open(file + "/ave_" + name + ".dat", ios::out); //Conversione SI del data blocking

    for(int i=0;i<N;i++){
      sum = 0;
      for(int j=0;j<L;j++){
	       input >> k;
	       sum += k;
      }
      ave[i] = sum/L;
      av2[i] = pow(ave[i],2);
    }

    for(int i=0;i<N;i++){
      for(int j=0;j<i+1;j++){
	       sum_prog[i] += ave[j];
	       su2_prog[i] += av2[j];
      }
      sum_prog[i] /= (i+1);
      su2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog,su2_prog,i);

      output << sum_prog[i] << " " << err_prog[i] << endl;

      if (name == "temp") data << sum_prog[i]*eps/K_b << " " << err_prog[i]*eps/K_b << endl;
      else if (name == "press") data << sum_prog[i]*eps/pow(sigma_A,3) << " " << err_prog[i]*eps/pow(sigma_A,3) << endl;
      else  data << sum_prog[i]*eps << " " << err_prog[i]*eps << endl;

    }

    input.close();
    output.close();
    data.close();

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
