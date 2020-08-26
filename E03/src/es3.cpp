#include <iostream>
#include <fstream>
#include <math.h>
#include "random.h"

using namespace std;

double error(double*, double*, int);
void Set_random(Random &);

int main (int argc, char *argv[]){

  //Inizializzazione libreria random
  Random rnd;
  Set_random(rnd);

  //fstream
  fstream output_C, output_P;

  //Parametri simulazione
  int M = pow(10,4);
  int N = pow(10,2);
  int L = M/N;

  double S0 = 100.;
  double K = 100.;
  double T = 1.;
  double r = 0.1;
  double sigma = 0.25;

  //Variabili Call
  double ave_C[N];
  double av2_C[N];
  double sum_prog_C[N];
  double su2_prog_C[N];
  double err_prog_C[N];

  //Variabili Put
  double ave_P[N];
  double av2_P[N];
  double sum_prog_P[N];
  double su2_prog_P[N];
  double err_prog_P[N];
 
  //Diretta

  //Inizializzazione variabili
  for(int i=0;i<N;i++){
    ave_C[i] = 0;
    av2_C[i] = 0;
    sum_prog_C[i] = 0;
    su2_prog_C[i] = 0;
    err_prog_C[i] = 0;

    ave_P[i] = 0;
    av2_P[i] = 0;
    sum_prog_P[i] = 0;
    su2_prog_P[i] = 0;
    err_prog_P[i] = 0;
  }

  double sum_C = 0;
  double sum_P = 0;

  output_C.open("./Data/dirC.dat", ios::out);
  output_P.open("./Data/dirP.dat", ios::out);


  //Data blocking
  for(int i=0;i<N;i++){
    sum_C = 0;
    sum_P = 0;
    for(int j=0;j<L;j++){
      double S_T = S0 * exp( (r - 1./2.*pow(sigma,2)*T) + sigma*rnd.Gauss(0.,T) );
      sum_C += exp(-r*T)*max(0.,S_T-K);
      sum_P += exp(-r*T)*max(0.,K-S_T);
    }
    ave_C[i] = sum_C/L;
    av2_C[i] = pow(ave_C[i],2);
    ave_P[i] = sum_P/L;
    av2_P[i] = pow(ave_P[i],2);
  }

  for(int i=0;i<N;i++){
    for(int j=0;j<i+1;j++){
      sum_prog_C[i] += ave_C[j];
      su2_prog_C[i] += av2_C[j];
      
      sum_prog_P[i] += ave_P[j];
      su2_prog_P[i] += av2_P[j];
    }
    sum_prog_C[i] /= (i+1);
    su2_prog_C[i] /= (i+1);
    err_prog_C[i] = error(sum_prog_C,su2_prog_C,i);

    sum_prog_P[i] /= (i+1);
    su2_prog_P[i] /= (i+1);
    err_prog_P[i] = error(sum_prog_P,su2_prog_P,i);

    output_C << i << " " << sum_prog_C[i] << " " <<  err_prog_C[i]  << endl;
    output_P << i << " " << sum_prog_P[i] << " " <<  err_prog_P[i]  << endl;
    
  }

  output_C.close();
  output_P.close();


  //Discretizzato
  int Step = 100;
  double S[Step];
  double delta = T/Step;
  S[0] = S0;

  //Reset variabili
  for(int i=0;i<N;i++){
    ave_C[i] = 0;
    av2_C[i] = 0;
    sum_prog_C[i] = 0;
    su2_prog_C[i] = 0;
    err_prog_C[i] = 0;

    ave_P[i] = 0;
    av2_P[i] = 0;
    sum_prog_P[i] = 0;
    su2_prog_P[i] = 0;
    err_prog_P[i] = 0;
  }

  output_C.open("./Data/disC.dat", ios::out);
  output_P.open("./Data/disP.dat", ios::out);


  //Data blocking
  for(int i=0;i<N;i++){
    sum_C = 0;
    sum_P = 0;
    for(int j=0;j<L;j++){
      for(int t=1;t<=Step;t++){
	S[t] = S[t-1]*exp( (r - 1./2.*pow(sigma,2)) * (delta) + sigma * rnd.Gauss(0.,1.) * sqrt(delta));
      }
      sum_C += exp(-r*T)*max(0.,S[Step]-K);
      sum_P += exp(-r*T)*max(0.,K-S[Step]);
    }
    ave_C[i] = sum_C/L;
    av2_C[i] = pow(ave_C[i],2);
    ave_P[i] = sum_P/L;
    av2_P[i] = pow(ave_P[i],2);
  }

  for(int i=0;i<N;i++){
    for(int j=0;j<i+1;j++){
      sum_prog_C[i] += ave_C[j];
      su2_prog_C[i] += av2_C[j];
      
      sum_prog_P[i] += ave_P[j];
      su2_prog_P[i] += av2_P[j];
    }
    sum_prog_C[i] /= (i+1);
    su2_prog_C[i] /= (i+1);
    err_prog_C[i] = error(sum_prog_C,su2_prog_C,i);

    sum_prog_P[i] /= (i+1);
    su2_prog_P[i] /= (i+1);
    err_prog_P[i] = error(sum_prog_P,su2_prog_P,i);

    output_C << i << " " << sum_prog_C[i] << " " <<  err_prog_C[i]  << endl;
    output_P << i << " " << sum_prog_P[i] << " " <<  err_prog_P[i]  << endl;
    
  }

  output_C.close();
  output_P.close();
  
  rnd.SaveSeed();
  return 0;
  
}


double error(double* AV, double* AV2, int n){
  
  if (n == 0)
    return 0;
  else
    return sqrt((AV2[n] - pow(AV[n],2))/(n) );

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
