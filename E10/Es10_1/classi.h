#include <stdlib.h>    
#include <iostream>     
#include <fstream>
#include <ctime>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <vector>
#include <string>
#include "random.h"

using namespace std;

vector<double> bests;
double hit = 0;
uint gen = 2*pow(10,4);
uint annealing = 0;


class Metropoly {

  Random _rnd;
  vector<pair<double,double>> _city;
  
 public:

  Metropoly(){
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
	  _rnd.SetRandom(seed,p1,p2);
	}
      }
      input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  }
  double Rannyu() { return _rnd.Rannyu(); }
  double Rannyu(double min, double max) { return _rnd.Rannyu(min, max); }
  vector<pair<double,double>> & get_city() { return _city; }
  int Pbc(int i) { return i%_city.size();  }
  virtual double distance(pair<double,double> & A, pair<double,double> & B) = 0;

  double L1(vector<int> & dna) {
    double sum = 0.;
    for (uint i=0;i<dna.size();i++) sum +=  distance( _city[dna[i]], _city[dna[Pbc(i+1)]] );
    return sum;
  }

  double L2(vector<int> & dna) {
    double sum = 0.;
    for (uint i=0;i<dna.size();i++) sum +=  pow(distance( _city[dna[i]], _city[dna[Pbc(i+1)]] ),2);
    return sum;
  }
      
}; 

class Square : public Metropoly {
  
public:

  Square(){  for(int i=0;i<32;i++) get_city().push_back( make_pair(Rannyu(), Rannyu()) );  }
  Square(uint N){  for(uint i=0;i<N;i++) get_city().push_back( make_pair(Rannyu(), Rannyu()) );  }
  
  double distance(pair<double,double> & A, pair<double,double> & B) { return ( sqrt( pow( A.first - B.first,2) + pow(A.second - B.second,2) ) ); }
      
};

class Circumference : public Metropoly {
  
 public:
  
  Circumference(){ for(int i=0;i<32;i++) get_city().push_back( make_pair(1., Rannyu(0.,2.*M_PI)) ); }
  Circumference(uint N){ for(uint i=0;i<N;i++) get_city().push_back( make_pair(1., Rannyu(0.,2.*M_PI)) ); }

  double distance(pair<double,double> & A, pair<double,double> & B) { return ( sqrt( pow( A.first*cos(A.second) - B.first*cos(B.second),2) + pow( A.first*sin(A.second) - B.first*sin(B.second),2) ) ); }
      
};



class Individual{

  uint _size;
  double _l1;
  double _l2;
  double _fitness;
  vector<int> _DNA;

 public:

 Individual() : _size(32), _l1(100), _l2(100) { for(uint i=0;i<_size;i++) _DNA.push_back(i); }
  Individual(std::vector<int>::iterator it, std::vector<int>::iterator it2) { _DNA = vector<int>(it,it2); _size = _DNA.size(); _l1=100.; _l2=100.; _fitness=0;}

  vector<int> & get_DNA() { return _DNA; }
  void set_l1(double l1) { _l1 = l1; }
  void set_l2(double l2) { _l2 = l2; }
  void set_fitness(double fitness) { _fitness = fitness; }

  double get_l1() { return _l1; }
  double get_l2() { return _l2; }
  double get_fitness() { return _fitness; }
  uint get_size() { return _size; }

  Individual & new_DNA(){
    for(uint k=0;k<pow(10,4);k++) (*this).swap((int)(1+((double)rand()/RAND_MAX)*(_size-1)),(int)(1+((double)rand()/RAND_MAX)*(_size-1)));
    return *this;
 }

  //Genetic mutations

  Individual & Mutation(){
    int i,j;
    const int m = (int)((double)rand()/RAND_MAX*4.);
    //cout << m << endl;
    switch (m){
    case 0:
      i = (int)(1+((double)rand()/RAND_MAX)*(_size-1));
      j = (int)(1+((double)rand()/RAND_MAX)*(_size-1));
      this->swap(i,j);
      //cout << "Swap: index " << i << " with index " << j << endl;
      break;
    case 1:
      i = (int)(2+((double)rand()/RAND_MAX)*(_size/2-2));
      this->permutation(i);
      //cout << "Permutation: " << i << " contiguous city permuted." << endl;
      break;
    case 2:
      j = (int)(1+((double)rand()/RAND_MAX)*(_size-1));
      i = (int)(1+((double)rand()/RAND_MAX)*(_size-1-j));
      this->shift(i,j);
      //cout << "Shift: " << i << " cities by " << j << " steps." << endl;
      break;
    case 3:
      i = (int)(2+((double)rand()/RAND_MAX)*(_size-2));
      this->reverse(i);
      //cout << "Reverse: " << i << " contiguous city reversed." << endl;
      break;
    }

    return *this;
  }

  
  Individual & swap(int i, int j) {// swap i-th city with j-th
    int appo = this->get_DNA()[i];
    this->get_DNA()[i] = this->get_DNA()[j];
    this->get_DNA()[j] = appo;
    
    return *this;
  }

  Individual & permutation(int m) {// permutation among m cities
    int i = (int)(1+((double)rand()/RAND_MAX)*((_size-1)-(2*m-1)));
    for(int j=i;j<i+m;j++) this->swap(j,j+m);

    return *this;
  }

  Individual & shift(int m, int n) {// shift m cities by n steps
    if(m==1){
      uint i = (int)(1+((double)rand()/RAND_MAX)*((_size-1)-(m+n-1)));
      std::vector<int>::iterator it = this->get_DNA().begin()+i;
      int appo = *it;
      this->get_DNA().erase(it);
      if ( i+n >= _size) this->get_DNA().insert(this->get_DNA().end() , appo);
      else this->get_DNA().insert(it+n, appo);
    }
    else{
      uint i = (int)(1+((double)rand()/RAND_MAX)*((_size-1)-(m+n-1)));
      std::vector<int>::iterator it1 = this->get_DNA().begin()+i;
      std::vector<int>::iterator it2 = this->get_DNA().begin()+(i+m);
      vector<int> appo(it1,it2);
      this->get_DNA().erase(it1,it2);
      if ( i+n >= (_size-appo.size()) ) this->get_DNA().insert( this->get_DNA().end() , appo.begin(), appo.end());
      else this->get_DNA().insert(it1+n, appo.begin(), appo.end());
    }

    return *this;
  }

  Individual & reverse(int m) {// reverse m cities
    int i = (int)(1+((double)rand()/RAND_MAX)*((_size-1)-(m-1)));
    std::vector<int>::iterator it1 = this->get_DNA().begin()+i;
    std::vector<int>::iterator it2 = this->get_DNA().begin()+(i+m);
    vector<int> appo(it1,it2);
    this->get_DNA().erase(it1,it2);
    this->get_DNA().insert(it1, appo.rbegin(), appo.rend()); 

    return *this;
  }
  
};


class Population{

  Random _rnd;
  double _Pc;
  double _Pm;
  double _T;
  vector<Individual> _route;
  vector<Individual> _newgen;
  bool _first;
  Individual _best;

 public:

  Population(){};
 Population(uint N) : _Pc(0.5), _Pm(0.1), _T(1){

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
	  _rnd.SetRandom(seed,p1,p2);
	}
      }
      input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    _route = vector<Individual>(1.);
    
    for(uint i=1;i<N;i++) {
      _route.push_back( _route[0] );
      _route[i].new_DNA();
    }

    _first = true;
 }

 vector<Individual> & get_route() { return _route; }
 vector<Individual> & get_newgen() { return _newgen; }
 Individual & get_best() { return _best; }
 bool & get_first() { return _first; }

 void set_Pc(double p){ _Pc = p; }
 void set_Pm(double p){ _Pm = p; }
 void set_T(double t){ _T = t; }
 
 double Rannyu() { return _rnd.Rannyu(); }
 double Rannyu(double min, double max) { return _rnd.Rannyu(min, max); }
 double Beta(){ return 1./_T; }
 double Boltzmann(double E){ return exp(-Beta()*E); }

 void _sort(){
   std::sort( std::begin(_route), std::end(_route), [](Individual & a,Individual & b) -> bool{  return a.get_fitness() > b.get_fitness(); });
 }

 void Initialize(Metropoly * metro){
   for (auto &el : _route) {
     el.set_l1( metro->L1( el.get_DNA() ) );
     el.set_fitness( Boltzmann(el.get_l1()) );
   }
   _sort();
   _first = false;
 }

 Individual selection(){
   int i = (int)(pow( (double)rand()/RAND_MAX, 2) *_route.size());
   return _route[i];

 }

 void add_to_newgen(Individual n){
   _newgen.push_back(n);
 }

 void crossover(Individual A, Individual B){

   int k = 20;//(int)Rannyu(1.,A.get_size()-2);
   vector<int> Achr(A.get_DNA().begin(),A.get_DNA().begin()+k);
   vector<int> Bchr(B.get_DNA().begin(),B.get_DNA().begin()+k);

   for(uint i=0;i<Achr.size();i++){
     vector<int>::iterator it = find(A.get_DNA().begin(), A.get_DNA().end(), Bchr[i]);
     A.get_DNA().erase(it);
     it = find(B.get_DNA().begin(), B.get_DNA().end(), Achr[i]);
     B.get_DNA().erase(it);
   }

   Achr.insert( Achr.end(), B.get_DNA().begin(), B.get_DNA().end() );
   Bchr.insert( Bchr.end(), A.get_DNA().begin(), A.get_DNA().end() );

   Individual a(Achr.begin(),Achr.end());
   Individual b(Bchr.begin(),Bchr.end());

   add_to_newgen(a);
   add_to_newgen(b);
   
 }
 
 double generation_best_l1(){ return _route[0].get_l1(); }
 double generation_best_l2(){ return _route[0].get_l2(); }

 double average_best_half_l1(){ 
   return std::accumulate(begin(_route), begin(_route)+_route.size()/2, 0., [](double sum, Individual o){ return o.get_l1() + sum; })/((double)_route.size()/2.);
 }

 double average_best_half_l2(){ 
   return std::accumulate(begin(_route), begin(_route)+_route.size()/2, 0., [](double sum, Individual o){ return o.get_l2() + sum; })/((double)_route.size()/2.);
 }


 void NewGeneration(Metropoly * metro){ //Simulated Annealing
  
   if (_first) Initialize(metro);
   annealing++;

   Individual id = _route[0];
   id.Mutation();
   id.set_l1( metro->L1( id.get_DNA() ) );

   if ( annealing%(int)(gen/pow(10,2)) == 0 ) {
     //cout << "Acceptance rate block: " << hit/(get_route().size() * (gen/pow(10,4))) << endl;
     _T -= _T/5.;
     hit = 0;
   }

   double alpha = std::min( 1., Boltzmann(id.get_l1() - _route[0].get_l1()) );

   if ( Rannyu() <= alpha ) {
     add_to_newgen(id);
     hit++;
   }
   else add_to_newgen(_route[0]);

   _route = _newgen;

   for (auto &el : _route){
     el.set_fitness( Boltzmann(el.get_l1()) );
   }

   bests.push_back( generation_best_l1() );

   if ( generation_best_l1() < _best.get_l1() ) _best = _route[0];
   
   _newgen.clear();

 }

 
};

