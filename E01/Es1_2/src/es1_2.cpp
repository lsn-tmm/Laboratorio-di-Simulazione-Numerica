#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <set>
#include "random.h"
#include "Exponential.h"
#include "Lorentzian.h"

using namespace std;

void Set_random(Random &);

int main (int argc, char *argv[]){

  //Inizializzazione libreria random
  Random rnd;
  Set_random(rnd);
  
  //Inizio programma
  Exponential exp(1.);
  Lorentzian lor(1.,0.);
  uint N[4] = {1,2,10,100};
  vector<double> sum(3);
  fstream Std,Exp,Lor;

  for (auto el : N){
    Std.open("Data/std" + to_string(el) + ".out", ios::out);
    Exp.open("Data/exp" + to_string(el) + ".out", ios::out);
    Lor.open("Data/lor" + to_string(el) + ".out", ios::out);
    for(uint i=0;i<pow(10,4);i++){
      sum = {0,0,0};
      for(uint j=0;j<el;j++){
	sum[0] += rnd.Rannyu();
	sum[1] += exp.Eval( rnd.Rannyu() );
	double x;
	//Controllo sull'intervallo di generazione per non far fallire l'istogramma
	do x = lor.Eval( rnd.Rannyu() ); while( x < -10 || x > 10); 
	sum[2] += x;
      }
      Std << sum[0]/(double)el << endl;
      Exp << sum[1]/(double)el << endl;
      Lor << sum[2]/(double)el << endl;
    }
    Std.close();
    Exp.close();
    Lor.close();
  }

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
