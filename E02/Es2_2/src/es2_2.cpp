#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <set>
#include "random.h"
#include "Vettore3D.h"

using namespace std;

double error(double, double, int);
void Set_random(Random &);
void Angle3D(double&, double&, Random&);

int main (int argc, char *argv[]){

  //Inizializzazione libreria random
  Random rnd;
  Set_random(rnd);

  //Parametri simulazione
  int N = pow(10,4);
  int Step = pow(10,2);

  //Variabili
  Vettore3D r[N];
  double ave;
  double av2;

  //1. Reticolo
  fstream output;
  output.open("Data/reticolo.out", ios::out);
  output << 0 << " " << 0 << endl; //Punto di partenza del cammino
  
  for(int i=0;i<Step;i++){
    double sum = 0;
    double sum2 = 0;
    for(int j=0;j<N;j++){
      int dice = (int)rnd.Rannyu(0,6);
      double x = 0,y = 0,z = 0;
      switch (dice){
      case 0 : x = 1; break;
      case 1 : y = 1; break;
      case 2 : z = 1; break;
      case 3 : x = -1; break;
      case 4 : y = -1; break;
      case 5 : z = -1; break;
      }

      Vettore3D a(x,y,z);
      r[j] = r[j] + a;
      sum = sum + pow(r[j].Modulo(),2);
      sum2 = sum2 + pow(r[j].Modulo(),4);

    }
    
    ave = sum/N;
    av2 = sum2/N;

    output << sqrt(ave) << " " << error(sqrt(ave), sqrt(av2), N) << endl;
  }

  output.close();

  //Reset variabili
  for (auto &el:r) el = Vettore3D();
  
  //2. Continuo
  output.open("Data/continuo.out", ios::out);
  output << 0 << " " << 0 << endl;// Punto di partenza del cammino
  
  for(int i=0;i<Step;i++){
    double sum = 0;
    double sum2 = 0;
    for(int j=0;j<N;j++){
      double theta = 0;
      double phi = 0;
      Angle3D(theta,phi,rnd);

      Vettore3D a( sin(theta)*cos(phi) , sin(theta)*sin(phi)  , cos(theta) );
      r[j] = r[j] + a;
      sum = sum + pow(r[j].Modulo(),2);
      sum2 = sum2 + pow(r[j].Modulo(),4);

    }
    
    ave = sum/N;
    av2 = sum2/N;

    output << sqrt(ave) << " " << error(sqrt(ave), sqrt(av2), N) << endl;
  }
  
  output.close();
  
  rnd.SaveSeed();
  return 0;
}


double error(double AV, double AV2, int n){

    return sqrt((AV2 - pow(AV,2))/(n-1) );

}

void Angle3D(double& t, double& p, Random& rand){
  t = rand.Rannyu(0.,M_PI);
  p = rand.Rannyu(0,2*M_PI);
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
