#ifndef _Vettore3D_h_
#define _Vettore3D_h_

#include <cmath>

class Vettore3D{

 private:

  double x_;
  double y_;
  double z_;

  double rho_;
  double phi_;
  double theta_;

 public:

  Vettore3D();
  Vettore3D(double,double,double);

  double Get_rho();
  double Get_phi();
  double Get_theta();

  double Get_x() const;
  double Get_y() const;
  double Get_z() const;

  void Set_x(double);
  void Set_y(double);
  void Set_z(double);

  void Clear();
  double Modulo();

  Vettore3D operator+(const Vettore3D&);  



};

#endif
