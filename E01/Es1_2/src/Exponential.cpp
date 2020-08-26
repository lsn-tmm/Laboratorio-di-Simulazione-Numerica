#include "Exponential.h"


Exponential::Exponential(double l) : lambda_(l) {}

void Exponential::Set_lambda(double l){ lambda_ = l; }

double Exponential::Get_lambda() const { return lambda_; }

double Exponential::Eval(double x) const{

  return (-1.)/lambda_ * log(1-x);

}

