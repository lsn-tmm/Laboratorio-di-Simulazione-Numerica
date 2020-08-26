#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
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
void swap(vector<Individual> &);

int main(){

  mpi::environment env;
  mpi::communicator world;

  int N_m = 10;

  srand(time(0) + world.rank());
  Population P(10000);

  Square sq;

  for(uint i=1; i<=gen; i++){
    if(i%N_m == 0){
      //std::cout << "I am process " << world.rank() << " with best of " << P.get_best().get_l1() << "." << std::endl;
      Individual rec;
      vector<Individual> all_best;
      if (world.rank() == 0){
        gather(world, P.get_best(), all_best, 0);
        for(int i=0; i<100; i++) swap(all_best);
      }
      else{
        gather(world, P.get_best(), 0);
      }
      mpi::scatter(world, all_best, rec, 0);
      P.get_route()[0] = rec;
      P.set_best(rec);
      //std::cout << "I am process " << world.rank() << " with best of " << P.get_best().get_l1() << "." << std::endl;
    }
    P.NewGeneration(&sq);
    Print(P,&sq,"sq" + to_string(world.rank()));
  }

  return 0;

}


void Print(Population & P, Metropoly * metro, string file){

  //cout << "Crossover rate: " << hit_c/(P.get_route().size()/2 * gen) << endl;
  //cout << "Mutation rate: " << hit_m/(P.get_route().size() * gen) << endl << endl;

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

void swap(vector<Individual> & v){
  Individual appo;
  int i = (int)(((double)rand()/RAND_MAX)*(v.size()));
  int j = (int)(((double)rand()/RAND_MAX)*(v.size()));

  appo = v[i];
  v[i] = v[j];
  v[j] = appo;
}
