#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <set>
#include "random.h"

using namespace std;

double error(double*, double*, int);
void Set_random(Random &);
double Angle2D(Random &);
double Buffon(int, Random &);

int main (int argc, char *argv[]){

  //Inizializzazione libreria random
  Random rnd;
  Set_random(rnd);

  //Parametri simulazione
  int M = pow(10,7);
  int N = 100;
  int L = M/N;

  //Variabili
  double ave[N];
  double av2[N];
  double sum_prog[N];
  double su2_prog[N];
  double err_prog[N];

  //Inizio simulazione
  fstream output;
  output.open("Data/pi.out", ios::out);
  
  //Inizializzazione variabili
  for(int i=0;i<N;i++){
    ave[i] = 0;
    av2[i] = 0;
    sum_prog[i] = 0;
    su2_prog[i] = 0;
    err_prog[i] = 0;
  }
  
  for(int i=0;i<N;i++){
    ave[i] = Buffon(L, rnd);
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

  }  

  output.close();

  rnd.SaveSeed();
  return 0;
}

double Angle2D(Random & rnd){
  double x,y;

  do{
    x = rnd.Rannyu();
    y = rnd.Rannyu();
  } while( pow(x,2) + pow(y,2) > 1 );

  if (y >= 0) return acos( x / sqrt(pow(x,2) + pow(y,2)) );
  else return 2*M_PI - acos( x / sqrt(pow(x,2) + pow(y,2)) );

}

double Buffon(int N, Random & rnd){
  double L = 0.75;
  double d_mezzi = 0.5;
  int hit = 0;

  for(int i=0;i<N;i++){
    double y = rnd.Rannyu(-d_mezzi,d_mezzi);
    double theta = Angle2D(rnd);

    if ( y+L*sin(theta) < -d_mezzi || y+L*sin(theta) > d_mezzi ) hit++;

  }

  return L*N/((double)hit*d_mezzi);

}

double error(double* AV, double* AV2, int n){
  
  if (n == 0)
    return 0;
  else
    return sqrt( (AV2[n] - pow(AV[n],2))/n );

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
