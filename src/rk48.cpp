#include "rk48.hpp"

// +------------------------------+
// | 8th Order Central Difference |
// +------------------------------+
void RK48::grad( const scalar_field& f, scalar_field& dfdx )
{

  static real a =  4 / (5*m_dx);
  static real b = -1 / (5*m_dx);
  static real c =  4 / (105*m_dx);
  static real d = -1 / (280*m_dx);

  int NX = m_nx;

  // left periodic boundary
  int i = 0; 
  dfdx[i] = a*( f[i+1] - f[NX-1] ) + b*( f[i+2] - f[NX-2] ) 
      + c*( f[i+3] - f[NX-3] ) + d*( f[i+4] - f[NX-4] );

  i = 1;
  dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[i+2] - f[NX-1] ) 
      + c*( f[i+3] - f[NX-2] ) + d*( f[i+4] - f[NX-3] );

  i = 2;
  dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[i+2] - f[i-2] ) 
      + c*( f[i+3] - f[NX-1] ) + d*( f[i+4] - f[NX-2] );

  i = 3;
  dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[i+2] - f[i-2] ) 
      + c*( f[i+3] - f[i-3] ) + d*( f[i+4] - f[NX-1] );

  for (i=4; i<NX-4; i++)
  {
    dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[i+2] - f[i-2] ) 
      + c*( f[i+3] - f[i-3] ) + d*( f[i+4] - f[i-4] );
  }

  // right periodic boundary
  i = NX-1;
  dfdx[i] = a*( f[0] - f[i-1] ) + b*( f[1] - f[i-2] ) 
      + c*( f[2] - f[i-3] ) + d*( f[3] - f[i-4] );

  i = NX-2;
  dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[0] - f[i-2] ) 
      + c*( f[1] - f[i-3] ) + d*( f[2] - f[i-4] );

  i = NX-3;
  dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[i+2] - f[i-2] ) 
      + c*( f[0] - f[i-3] ) + d*( f[1] - f[i-4] );

  i = NX-4;
  dfdx[i] = a*( f[i+1] - f[i-1] ) + b*( f[i+2] - f[i-2] ) 
      + c*( f[i+3] - f[i-3] ) + d*( f[0] - f[i-4] );

}

void RK48::solve( real& t, scalar_field& f, real dt )
{

  size_t NX = m_nx;

  static scalar_field k1( NX, real(0) );
  static scalar_field k2( NX, real(0) );
  static scalar_field k3( NX, real(0) );
  static scalar_field k4( NX, real(0) );
  static scalar_field tmp( NX, real(0) );
  static scalar_field flt( NX, real(0) );

  static real b1 = real(1) / 6;
  static real b2 = real(2) / 6;
  static real b3 = real(2) / 6;
  static real b4 = real(1) / 6;

  // stage 1
  eval_rhs( t, f, k1 );

  // stage 2
  double t2 = t + .5*dt;
  for (std::size_t i=0; i<NX; i++) tmp[i] = f[i] + .5*dt*k1[i];
  eval_rhs( t2, tmp, k2 );

  // stage 3
  double t3 = t + .5*dt;
  for (std::size_t i=0; i<NX; i++) tmp[i] = f[i] + .5*dt*k2[i];
  eval_rhs( t3, tmp, k3 );

  // stage 4
  t += dt;
  for (std::size_t i=0; i<NX; i++) tmp[i] = f[i] + dt*k3[i];
  eval_rhs( t, tmp, k4 );

  // put it all together
  for (std::size_t i=0; i<NX; i++)
    f[i] = f[i] + dt*( b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] );

}

/*
void filter( const std::vector<double>& f, std::vector<double>& flt )
{

  double Chi =  -252 / double(1024);
  double a   = 210 / double(1024);
  double b   = - 120 / double(1024);
  double c   =  45 / double(1024);
  double d   =  - 10 / double(1024);
  double e   =   1 / double(1024);

  // left periodic boundary
  int i = 0; 
  flt[i] = Chi*f[i] + a*( f[i+1] + f[NX-1] ) + b*( f[i+2] + f[NX-2] ) 
      + c*( f[i+3] + f[NX-3] ) + d*( f[i+4] + f[NX-4] ) 
      + e*( f[i+5] + f[NX-5] );

  i = 1;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[NX-1] ) 
      + c*( f[i+3] + f[NX-2] ) + d*( f[i+4] + f[NX-3] ) 
      + e*( f[i+5] + f[NX-4] );

  i = 2;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[i+3] + f[NX-1] ) + d*( f[i+4] + f[NX-2] ) 
      + e*( f[i+5] + f[NX-3] );

  i = 3;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[i+3] + f[i-3] ) + d*( f[i+4] + f[NX-1] )
      + e*( f[i+5] + f[NX-2] );

  i = 4;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[i+3] + f[i-3] ) + d*( f[i+4] + f[i-4] )
      + e*( f[i+5] + f[NX-1] );

  for (i=5; i<NX-5; i++)
  {
    flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[i+3] + f[i-3] ) + d*( f[i+4] + f[i-4] ) + e*( f[i+5] + f[i-5] );
  }

  // right periodic boundary
  i = NX-1;
  flt[i] = Chi*f[i] + a*( f[0] + f[i-1] ) + b*( f[1] + f[i-2] ) 
      + c*( f[2] + f[i-3] ) + d*( f[3] + f[i-4] ) 
      + e*( f[4] + f[i-5] );

  i = NX-2;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[0] + f[i-2] ) 
      + c*( f[1] + f[i-3] ) + d*( f[2] + f[i-4] ) 
      + e*( f[3] + f[i-5] );

  i = NX-3;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[0] + f[i-3] ) + d*( f[1] + f[i-4] ) 
      + e*( f[2] + f[i-5] );

  i = NX-4;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[i+3] + f[i-3] ) + d*( f[0] + f[i-4] ) 
      + e*( f[1] + f[i-5] );

  i = NX-5;
  flt[i] = Chi*f[i] + a*( f[i+1] + f[i-1] ) + b*( f[i+2] + f[i-2] ) 
      + c*( f[i+3] + f[i-3] ) + d*( f[i+4] + f[i-4] ) 
      + e*( f[0] + f[i-5] );

}
*/