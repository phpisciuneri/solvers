#include "emac24.hpp"

void EMac24::grad( const scalar_field& f, scalar_field& dfdx )
{
  static int call_count = 3;
  
  call_count++;
  switch( call_count % 4 ) {
  case 0:
    forward_diff( f, dfdx );
    break;
  case 1:
    backward_diff( f, dfdx );
    break;
  case 2:
    backward_diff( f, dfdx );
    break;
  case 3:
    forward_diff( f, dfdx );
    break;
  default:
    std::exit(1); // this is an error, handle better
    break;
  }
}

void EMac24::solve( real& t, scalar_field& U, real dt )
{
  static scalar_field U1( m_nx, real(0) );
  static scalar_field U2( m_nx, real(0) );
  static scalar_field tmp( m_nx, real(0) );

  // predictor
  eval_rhs( t, U, tmp );
  for (std::size_t i=0; i<m_nx; i++) 
    U1[i] = U[i] + dt*tmp[i];

  // corrector
  eval_rhs( t, U1, tmp );
  for (std::size_t i=0; i<m_nx; i++) 
    U2[i] = U1[i] + dt*tmp[i];
  
  // update
  t += dt;
  for (std::size_t i=0; i<m_nx; i++) 
    U[i] = .5*( U[i] + U2[i] );
}

void EMac24::forward_diff( const scalar_field& f, scalar_field& df )
{
  static real b_L     =  real(5) / ( 12*m_dx );
  static real a_L     = -real(8) / ( 12*m_dx );
  static real Upsilon = -real(8) / ( 12*m_dx );
  static real a_R     =  real(8) / ( 12*m_dx );
  static real b_R     =  real(3) / ( 12*m_dx );

  // left periodic boundary
  size_t i = 0;
  df[i] = b_L*f[m_nx-2] + a_L*f[m_nx-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[i+2];

  i = 1;
  df[i] = b_L*f[m_nx-1] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[i+2];

  for (i=2; i<m_nx-2; i++)
    df[i] = b_L*f[i-2] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[i+2];

  // right periodic boundary
  i = m_nx - 2;
  df[i] = b_L*f[i-2] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[0];

  i = m_nx - 1;
  df[i] = b_L*f[i-2] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[0] + b_R*f[1];
}

void EMac24::backward_diff( const scalar_field& f, scalar_field& df )
{
  static real b_L     = -real(3) / ( 12*m_dx );
  static real a_L     = -real(8) / ( 12*m_dx );
  static real Upsilon =  real(8) / ( 12*m_dx );
  static real a_R     =  real(8) / ( 12*m_dx );
  static real b_R     = -real(5) / ( 12*m_dx );

  // left periodic boundary
  size_t i = 0;
  df[i] = b_L*f[m_nx-2] + a_L*f[m_nx-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[i+2];

  i = 1;
  df[i] = b_L*f[m_nx-1] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[i+2];

  for (i=2; i<m_nx-2; i++)
    df[i] = b_L*f[i-2] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[i+2];

  // right periodic boundary
  i = m_nx - 2;
  df[i] = b_L*f[i-2] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[i+1] + b_R*f[0];

  i = m_nx - 1;
  df[i] = b_L*f[i-2] + a_L*f[i-1] + Upsilon*f[i] + a_R*f[0] + b_R*f[1];
}




