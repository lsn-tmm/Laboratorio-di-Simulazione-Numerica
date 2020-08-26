#include "Lorentzian.h"


Lorentzian::Lorentzian(double g, double m) : gamma_(g), mu_(m) {}

void Lorentzian::Set_gamma(double g){ gamma_ = g; }

void Lorentzian::Set_mu(double m){ mu_ = m; }

double Lorentzian::Get_gamma() const { return gamma_; }

double Lorentzian::Get_mu() const { return mu_; }

double Lorentzian::Eval(double x) const{

  return gamma_ * tan(M_PI *( x - 0.5) );

}
