#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <functional>
#include "random.h"

using namespace std;

double error(vector<double> &, vector<double> &, int);
void Set_random(Random &);
void Data_Blocking(const vector<double> &, int, int, string);
double Mod_Psi100(const vector<double> &);
double Mod_Psi210(const vector<double> &);
vector<double> T_uniform(vector<double>, double, Random &);
vector<double> T_gauss(vector<double>, double, Random &);

void my_Metropolis(function<double(const vector<double> &)>, function<vector<double>(vector<double>, double, Random &)>, vector<double>, double, string, Random &);

int main(){

 //Inizializzazione libreria random
  Random rnd;
  Set_random(rnd);


  //my_Metropolis(armonica sferica, T(x|y), posizione iniziale, passo (unità di a_0), nome file, rnd)

  //Origine o vicino
  my_Metropolis(Mod_Psi100, T_uniform, {0.0,0.0,0.0}, 1.35, "psi100_uni.0", rnd);
  my_Metropolis(Mod_Psi100, T_gauss, {0.0,0.0,0.0}, 0.425, "psi100_multi.0", rnd);

  my_Metropolis(Mod_Psi210, T_uniform, {0.1,0.1,0.1}, 3.3, "psi210_uni.0", rnd);
  my_Metropolis(Mod_Psi210, T_gauss, {0.1,0.1,0.1}, 1.05, "psi210_multi.0", rnd);

  //Lontano dall'origine
  my_Metropolis(Mod_Psi100, T_uniform, {10,10,10}, 1.35, "psi100_uni.pos", rnd);
  my_Metropolis(Mod_Psi100, T_gauss, {10,10,10}, 0.425, "psi100_multi.pos", rnd);

  my_Metropolis(Mod_Psi210, T_uniform, {10,10,10}, 3.3, "psi210_uni.pos", rnd);
  my_Metropolis(Mod_Psi210, T_gauss, {10,10,10}, 1.05, "psi210_multi.pos", rnd);

  //Molto lontano dall'origine
  my_Metropolis(Mod_Psi100, T_uniform, {100,100,100}, 1.35, "psi100_uni2.pos", rnd);
  my_Metropolis(Mod_Psi100, T_gauss, {100,100,100}, 0.425, "psi100_multi2.pos", rnd);

  my_Metropolis(Mod_Psi210, T_uniform, {100,100,100}, 3.3, "psi210_uni2.pos", rnd);
  my_Metropolis(Mod_Psi210, T_gauss, {100,100,100}, 1.05, "psi210_multi2.pos", rnd);



  return 0;

}


//Funzioni

double Mod_Psi100(const vector<double> & xyz){

  double r = sqrt( pow(xyz[0],2) + pow(xyz[1],2) + pow(xyz[2],2) );

  return pow( 1/sqrt(M_PI) * exp(-r) ,2);

}

double Mod_Psi210(const vector<double> & xyz){

  double r = sqrt( pow(xyz[0],2) + pow(xyz[1],2) + pow(xyz[2],2) );

  return pow( 1./8. * sqrt(2./M_PI) * r * exp(-r/2.) * xyz[2]/r,2);

}

vector<double> T_uniform(vector<double> v, double eps, Random & rnd){

  v[0] = rnd.Rannyu(v[0]-eps,v[0]+eps);
  v[1] = rnd.Rannyu(v[1]-eps,v[1]+eps);
  v[2] = rnd.Rannyu(v[2]-eps,v[2]+eps);

  return vector<double>(v);
}

vector<double> T_gauss(vector<double> v, double eps, Random & rnd){

  v[0] = rnd.Gauss(v[0],2*eps);
  v[1] = rnd.Gauss(v[1],2*eps);
  v[2] = rnd.Gauss(v[2],2*eps);

  return vector<double>(v);
}


void Data_Blocking(const vector<double> & data, int nstep, int nblock, string file){

  //Inizializzazione variabili incertezza
  int L = nstep/nblock;

  vector<double> ave(nblock);
  vector<double> av2(nblock);
  vector<double> sum_prog(nblock);
  vector<double> su2_prog(nblock);
  vector<double> err_prog(nblock);

  fstream output;

  double sum = 0;
  int k = 0;

  for(int i=0;i<nblock;i++){
    sum = 0;
    for(int j=0;j<L;j++){
      k = j+i*L;
      sum += data[k];
    }
    ave[i] = sum/L;
    av2[i] = pow(ave[i],2);
  }

  output.open("Data/" + file, ios::out);
  for(int i=0;i<nblock;i++){
    for(int j=0;j<i+1;j++){
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }
    sum_prog[i] /= (i+1);
    su2_prog[i] /= (i+1);
    err_prog[i] = error(sum_prog,su2_prog,i);

    output << i << " " << sum_prog[i] << " " << err_prog[i] << endl;

  }
  output.close();

}

double error(vector<double> & AV, vector<double> & AV2, int n){

  if (n == 0)
    return 0;
  else
    return sqrt((AV2[n] - pow(AV[n],2))/(n) );

}


void my_Metropolis(function<double(const vector<double> &)> Psi, function<vector<double>(vector<double>, double, Random &)> T, vector<double> start, double eps, string file, Random & rnd){

  int nstep = pow(10,6);
  int nblock = pow(10,2);
  int eqstep = pow(10,5);

  fstream output;

  vector<double> r(nstep);
  vector<vector<double>> xyz (nstep + eqstep, vector<double>(3));
  int hit = 0;

  xyz[0] = start; //Posizione iniziale

  output.open("Scatter/" + file, ios::out);
  output << xyz[0][0] << " " << xyz[0][1] << " " << xyz[0][2] << endl;

  for(int i=1;i<nstep + eqstep;i++){

    vector<double> v = T(xyz[i-1],eps,rnd);

    double alpha = min( 1., Psi(v)/Psi(xyz[i-1]) );

    if ( rnd.Rannyu() <= alpha){
      xyz[i] = v;
      hit++;
    }
    else xyz[i] = xyz[i-1];

    output << xyz[i][0] << " " << xyz[i][1] << " " << xyz[i][2] << endl;

    if(i >= eqstep){
      r[i-eqstep] = sqrt( pow(xyz[i-eqstep][0],2) + pow(xyz[i-eqstep][1],2) + pow(xyz[i-eqstep][2],2) );
    }

  }

  output.close();
  cout << (double)hit/(double)nstep * 100 << "%" << endl;
  Data_Blocking(r,nstep,nblock,file);


}

void Set_random(Random & rnd){

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

}
