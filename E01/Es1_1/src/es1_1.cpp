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
double kiquadro(vector<multiset<double>>,int, int);
void Reset(vector<multiset<double>> &);


int main (int argc, char *argv[]){

  //Inizializzazione libreria random
  Random rnd;
  Set_random(rnd);

  //Parametri simulazione
  int M = 1000000;
  int N = 100;
  int L = M/N;

  //Variabili
  vector<double> r(M);
  for(auto &el : r){
    el = rnd.Rannyu();
  }
  double ave[N];
  double av2[N];
  double sum_prog[N];
  double su2_prog[N];
  double err_prog[N];


  //Es1.1
  fstream output;
  output.open("Data/int_1.out", ios::out);
  
  //Reset variabili
  for(int i=0;i<N;i++){
    ave[i] = 0;
    av2[i] = 0;
    sum_prog[i] = 0;
    su2_prog[i] = 0;
    err_prog[i] = 0;
  }
  double sum = 0;
  int k = 0;
  
  for(int i=0;i<N;i++){
    sum = 0;
    for(int j=0;j<L;j++){
      k = j+i*L;
      sum += r[k];
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

  }  

  output.close();

  //Es1.2
  output.open("Data/int_2.out", ios::out);
  
  //Reset variabili
  for(int i=0;i<N;i++){
    ave[i] = 0;
    av2[i] = 0;
    sum_prog[i] = 0;
    su2_prog[i] = 0;
    err_prog[i] = 0;
  }
  sum = 0;
  k = 0;
  
  for(int i=0;i<N;i++){
    sum = 0;
    for(int j=0;j<L;j++){
      k = j+i*L;
      sum += pow(r[k]-0.5,2);
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

  }  

  output.close();

  //1.3
  //Parametri simulazione
  M = 100;
  N = pow(10,4);
  double x = 0;
  vector<multiset<double>> n(M);
  double ki2[M];

  output.open("Data/ki2.out", ios::out);

  for(int i=0;i<100;i++){
    for(int j=0;j<N;j++){
      x = rnd.Rannyu();
      n[(int)(x*100)].insert(x);
    }
    ki2[i] = kiquadro(n,M,N);
    output << i+1 << " " << ki2[i] << endl;
    Reset(n); // Reset dei ogni multiset nel vector
  }

  output.close();

  rnd.SaveSeed();
  return 0;
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

void Reset(vector<multiset<double>> & n){
  for (auto &el:n) el.clear();
}

double kiquadro(vector<multiset<double>> n, int M, int N){
  
  double k = 0;
  
  for(auto el : n){
    k += pow((double)el.size()-N/M,2)*M/N ;
  }

  return k;

}


double error(double* AV, double* AV2, int n){
  
  if (n == 0)
    return 0;
  else
    return sqrt( (AV2[n] - pow(AV[n],2))/n );

}
