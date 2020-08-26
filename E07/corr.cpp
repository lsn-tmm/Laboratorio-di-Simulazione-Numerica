#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <functional>
#include <vector>
#include <string>
#include <numeric>

using namespace std;

void Autocorrelation(string, string);

int main(){

  // check //Autocorrelation("MonteCarlo/correlation/correlation_epot.solid", "ac_epot.solid");
  // check //Autocorrelation("MonteCarlo/correlation/correlation_press.solid", "ac_press.solid");

  // check //Autocorrelation("MonteCarlo/correlation/correlation_epot.liquid", "ac_epot.liquid");
  // check //Autocorrelation("MonteCarlo/correlation/correlation_press.liquid", "ac_press.liquid");

  // check //Autocorrelation("MonteCarlo/correlation/correlation_epot.gas", "ac_epot.gas");
  Autocorrelation("MonteCarlo/correlation/correlation_press.gas", "ac_press.gas");
  
  return 0;
  
}

void Autocorrelation(string Input, string Output){

  fstream Read, Write;
  vector<double> data;
  //vector<double> chi;
  double x;

  Read.open(Input, ios::in);

  while (true) {
    Read >> x;
    Read >> x;
    if( Read.eof() ) break;
    data.push_back(x);
  } 

  Read.close();

  Write.open(Output, ios::out);

  double sum = std::accumulate(data.begin(), data.end(), 0.0);
  double mean = sum / data.size();

  double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
  double var = sq_sum / data.size() - mean * mean;

  
  for(int i=0;i<data.size();i++){
    double sum_prog = 0;
    double sum1 = 0;
    double sum2 = 0;
    for(int j=0;j<(data.size()-i);j++){
	sum_prog += data[j] * data[j+i];
	sum1 += data[j];
	sum2 += data[i+j];
    }

    double ave_prog = sum_prog/(data.size()-i);
    double ave1 = sum1/(data.size()-i);
    double ave2 = sum2/(data.size()-i);

    if (i%50000 == 0) cout << i << endl;
  
    Write << (ave_prog - ave1*ave2)/var << endl;
   }

  Write.close();
  
}
