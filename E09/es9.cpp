#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <functional>
#include <ctime>
#include <random>
#include <vector>
#include <string>
#include "random.h"
#include "classi.h"

using namespace std;

void Print(Population &, Metropoly*, string);
void Reset(Population &);

int main(){
  
  srand(time(0));
  Population P(10000);

  Square sq;
  Circumference cr;
 
  for(uint i=0; i<gen; i++) P.NewGeneration(&sq);
  Print(P,&sq,"sq");

  Reset(P);

  for(uint i=0; i<gen; i++) P.NewGeneration(&cr);
  Print(P,&cr,"cr");


  return 0;
  
}


void Print(Population & P, Metropoly * metro, string file){

  cout << "Crossover rate: " << hit_c/(P.get_route().size()/2 * gen) << endl;
  cout << "Mutation rate: " << hit_m/(P.get_route().size() * gen) << endl << endl;

  fstream output;

  output.open("Data/" + file + "_bestvalue.dat", ios::out);
  for (auto el:bests) output << el << endl;
  output.close();

  output.open("Data/" + file + "_halfaverage.dat", ios::out);
  for (auto el:averages) output << el << endl;
  output.close();
  
  output.open("Data/" + file + "_bestpath.dat", ios::out);
  for (auto el:P.get_best().get_DNA()) output << metro->get_city()[el].first << " " << metro->get_city()[el].second << endl;
  output << (*metro->get_city().begin()).first << " " << (*metro->get_city().begin()).second << endl;
  output.close();
 
}

void Reset(Population & P){

  P.get_first() = true;
  P.get_best().set_l1(100);

  bests.clear();
  averages.clear();
  hit_m = 0;
  hit_c = 0;

}
