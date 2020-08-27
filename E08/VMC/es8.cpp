#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <functional>
#include <vector>
#include <string>
#include "random.h"
#include "es8.h"

using namespace std;

int main(){

  Input();

  Minimize();

  Data_Blocking();
  Print("gs_energy.out");

  my_Metropolis("psi.dat");

  return 0;

}


//Funzioni

void Input(){

  ifstream ReadInput;

  //Inizializzazione libreria random
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "RANDOMSEED" ){
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
      }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  ReadInput.open("input.dat");

  ReadInput >> mu;
  ReadInput >> sigma;

  _mu = mu;
  _sigma = sigma;

  ReadInput >> range;
  xmax = range/2.;
  xmin = -range/2.;

  ReadInput >> eps;

  ReadInput >> nblock;
  ReadInput >> nstep;
  L = nstep/nblock;

  ReadInput.close();

}

double Mod_Psi(const double & x){

  return pow(exp( -pow(x-mu,2) / (2.*pow(sigma,2)) ) +  exp( -pow(x+mu,2) / (2.*pow(sigma,2)) ),2);
  //return exp(-x);
}

double Psi(const double & x){

  return exp( -pow(x-mu,2) / (2.*pow(sigma,2)) ) +  exp( -pow(x+mu,2) / (2.*pow(sigma,2)) );

}

double Psi_second(const double & x){

  double exp_m = exp( -pow(x-mu,2) / (2.*pow(sigma,2)) );
  double exp_p = exp( -pow(x+mu,2) / (2.*pow(sigma,2)) );

  return ( (pow(x,2) - 2.*x*mu + pow(mu,2) - pow(sigma,2)) * exp_m + (pow(x,2) + 2.*x*mu + pow(mu,2) - pow(sigma,2)) * exp_p ) / pow(sigma,4);

}

double V(const double & x){

  return pow(x,4) - 2.5 * pow(x,2);

}

double Integral_point(){

  double x = rnd.Rannyu(xmin, xmax);

  return Mod_Psi(x) / integral_psi * (-0.5 * Psi_second(x) + V(x) * Psi(x)) / Psi(x);

}


double Integral_point_sampling(const double & x){

  return 1./ integral_psi * (-0.5 * Psi_second(x) + V(x) * Psi(x)) / Psi(x);

}

double Integral(){

  integral_psi = Integral_psi();
  double sum = 0;
   for(int i=0;i<nstep;i++){
     double y = Integral_point();
     sum += (xmax-xmin)*y;
  }
  return sum/nstep;

}

double Integral_psi(){

 // Computing value of Integral[ |psi_t|^2 ]
  double sum = 0;
  for(int i=0;i<nstep;i++){
    double y = rnd.Rannyu(xmin,xmax);
    sum += (xmax-xmin)*Mod_Psi(y);
  }
  return sum/nstep;

}

void Data_Blocking(){

  fill(ave.begin(), ave.end(), 0);
  fill(av2.begin(), av2.end(), 0);
  fill(sum_prog.begin(), sum_prog.end(), 0);
  fill(su2_prog.begin(), su2_prog.end(), 0);
  fill(err_prog.begin(), err_prog.end(), 0);

  double sum = 0;
  integral_psi = norm;
  mu = _mu;
  sigma = _sigma;

  cout << integral_psi << " " << _mu << " " << _sigma << endl;

  for(int i=0;i<nblock;i++){
    sum = 0;
    for(int j=0;j<L;j++){
      sum += Integral_point();
    }
    ave[i] = (xmax - xmin) * sum/L;
    av2[i] = pow(ave[i],2);
  }

  for(int i=0;i<nblock;i++){
    for(int j=0;j<i+1;j++){
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }
    sum_prog[i] /= (i+1);
    su2_prog[i] /= (i+1);
    err_prog[i] = error(sum_prog,su2_prog,i);

    //cout<< i << " " << sum_prog[i] << " " << err_prog[i] << endl;

  }

}

void Data_Blocking_sampling(){

  fill(ave.begin(), ave.end(), 0);
  fill(av2.begin(), av2.end(), 0);
  fill(sum_prog.begin(), sum_prog.end(), 0);
  fill(su2_prog.begin(), su2_prog.end(), 0);
  fill(err_prog.begin(), err_prog.end(), 0);

  double sum = 0;
  integral_psi = norm;
  mu = _mu;
  sigma = _sigma;

  cout << integral_psi << " " << _mu << " " << _sigma << endl;

  double x = rnd.Rannyu(xmin,xmax);

  for(int i=0;i<nblock;i++){
    sum = 0;
    for(int j=0;j<L;j++){
      double xnew = Pbc( rnd.Rannyu(x-eps,x+eps) ) ;

      double alpha = min( 1., Mod_Psi(xnew)/Mod_Psi(x) );

      if ( rnd.Rannyu() <= alpha){
	x = xnew;
      }
      else x = x;
      sum += Integral_point_sampling(x);
    }
    ave[i] = (xmax - xmin) * sum/L;
    av2[i] = pow(ave[i],2);
  }

  for(int i=0;i<nblock;i++){
    for(int j=0;j<i+1;j++){
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }
    sum_prog[i] /= (i+1);
    su2_prog[i] /= (i+1);
    err_prog[i] = error(sum_prog,su2_prog,i);

    //cout<< i << " " << sum_prog[i] << " " << err_prog[i] << endl;

  }

}

void Minimize(){

  int hit = 0;
  int M = pow(10,2);

  energy_old = Integral();

  T = 0.0025;
  meps = 0.1;
  seps = 0.1;

  for(int k=1;k<=3;k++){

    cout << "Minimalize blk: " << k << endl;

    T -= 1/5 * T;

    // Mettere un while e fattore di convergenza
    for(int i=1;i<M;i++){

      mu = rnd.Rannyu(_mu-meps,_mu+meps);
      sigma = rnd.Rannyu(_sigma-seps,_sigma+seps);

      energy_new = Integral();

      //cout << energy_new << " " << mu << " " << sigma << endl;

      double alpha = min( 1., Boltzmann(energy_new - energy_old) );

      if ( rnd.Rannyu() <= alpha ){
	       _mu = mu;
	       _sigma = sigma;
	       energy_old = energy_new;
	       norm = integral_psi;
	       hit++;
      }

      cout << energy_old << " " << _mu << " " << _sigma << endl;

    }

  }

  //cout << "Minimalize: " << (double)hit/((double)M*3) * 100 << "%" << endl;

}

double Boltzmann(double delta){ return exp( - beta() * delta ); }

double beta(){ return 1./T; }



void my_Metropolis(string file){

  ofstream output;

  int hit = 0;

  output.open(file);

  double x = xmin;

  for(int i=1;i<nstep;i++){

    double xnew = Pbc( rnd.Rannyu(x-eps,x+eps) ) ;

    double alpha = min( 1., Mod_Psi(xnew)/Mod_Psi(x) );

    if ( rnd.Rannyu() <= alpha){
      x = xnew;
      hit++;
    }
    else x = x;

    output << x << endl;

  }

  output.close();

  cout << "Metropolis: " <<(double)hit/(double)nstep * 100 << "%" << endl;

}


double error(vector<double> & AV, vector<double> & AV2, int n){

  if (n == 0)
    return 0;
  else
    return sqrt((AV2[n] - pow(AV[n],2))/(n) );

}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - range * rint(r/range);
}

void Print(string file){

  fstream output;

  output.open(file, ios::out);
  for(int i=0;i<nblock;i++)  output << i+1 << " " << sum_prog[i] << " " << err_prog[i] << endl;
  output.close();

  output.open("parameters.out", ios::out);
  output << _mu << " " << _sigma << " " << norm << endl;
  output.close();

}
