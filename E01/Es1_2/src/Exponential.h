#ifndef _Exponential_h_
#define _Exponential_h_

#include "DistribuzioneBase.h"
#include <cmath>

class Exponential : public DistribuzioneBase {

 private:

  double lambda_;

 public:

  Exponential(double);

  void Set_lambda(double);
  double Get_lambda() const;
  
  double Eval(double) const;

};

#endif
