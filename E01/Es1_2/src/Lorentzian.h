#ifndef _Lorentzian_h_
#define _Lorentzian_h_

#include "DistribuzioneBase.h"
#include <cmath>

class Lorentzian : public DistribuzioneBase {

 private:

  double gamma_;
  double mu_;

 public:

  Lorentzian(double,double);

  void Set_gamma(double);
  void Set_mu(double);
  double Get_gamma() const;
  double Get_mu() const;
  
  
  double Eval(double) const;

};

#endif
