#include "Vettore3D.h"

Vettore3D::Vettore3D() : x_(0), y_(0), z_(0) {}

Vettore3D::Vettore3D(double x, double y, double z) : x_(x), y_(y), z_(z) {}

double Vettore3D::Get_rho(){return rho_;}
double Vettore3D::Get_phi(){return phi_;}
double Vettore3D::Get_theta(){return theta_;}

double Vettore3D::Get_x() const {return x_;}
double Vettore3D::Get_y() const {return y_;}
double Vettore3D::Get_z() const {return z_;}

void Vettore3D::Set_x(double x){

  x_ = x;  

}

void Vettore3D::Set_y(double y){

  y_ = y;  

}

void Vettore3D::Set_z(double z){

  z_ = z;
  
}

void Vettore3D::Clear(){ x_ = 0; y_ = 0; z_= 0; }

double Vettore3D::Modulo(){ return sqrt(pow(x_,2)+pow(y_,2)+pow(z_,2)); }

Vettore3D Vettore3D::operator+(const Vettore3D& V){

  return Vettore3D(x_+V.Get_x(), y_+V.Get_y(), z_+V.Get_z());

}
