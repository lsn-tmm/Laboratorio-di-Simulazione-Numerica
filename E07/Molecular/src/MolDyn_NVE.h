/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, nbins,bin_size,sd, vtail,ptail;
double T = 0;
double walker[m_props];

// averages
double acc,att;
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double err_gdir, err_pot, err_press, stima_pres;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
std::string file;

// thermodynamical state
int npart;
double energy,beta,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblk;
double delta;

//functions
void Input(void);
void Restart(void);
void Reset(int);
void Accumulate(void);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
void Measure_T(void);
void Averages();
void Averages(int);
void Print_SI(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);
double error(double*, double*, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
