#include "mac22.hpp"

void Mac22::grad( const scalar_field& f, scalar_field& dfdx )
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

void Mac22::solve( real& t, scalar_field& U, real dt )
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

void Mac22::forward_diff( const scalar_field& f, scalar_field& df )
{
  // this is the viscous notation - fix
  static real Upsilon = -1 / m_dx;
  static real a_R     =  1 / m_dx;

  // right periodic boundary
  size_t i = m_nx-1; 
  df[i] = Upsilon*f[i] + a_R*f[0];

  for (i=0; i<m_nx-1; i++)
    df[i] = Upsilon*f[i] + a_R*f[i+1];
}

void Mac22::backward_diff( const scalar_field& f, scalar_field& df )
{
  // this is the viscous notation - fix
  static real a_L     = -1 / m_dx;
  static real Upsilon =  1 / m_dx;
  
  // left periodic boundary
  size_t i = 0; 
  df[i] = Upsilon*f[i] + a_L*f[m_nx-1];

  for (i=1; i<m_nx; i++)
    df[i] = Upsilon*f[i] + a_L*f[i-1];
}


