#include <iostream>
#include <fstream>
#include <vector>

#include "field.hpp"
#include "solver.hpp"
#include "rk48.hpp"
#include "mac22.hpp"
#include "emac24.hpp"

const double pi = 3.14159265358979323846264338327;
double tmax = 2*pi;
const int max_iter = 201; //1400;

void init_sin( const Field& field, scalar_field& f )
{
  for (std::size_t i=0; i<field.nx(); ++i) 
    f[i] = std::sin( field.x(i) );
}

void init_hat( const Field& field, scalar_field& f )
{
  for (std::size_t i=0; i<field.nx(); ++i)
  {
    if ( field.x(i) < 1 || field.x(i) > 2 )
      f[i] = 0.;
    else
      f[i] = 1.;
  }
}

void init_smooth_hat( const Field& field, scalar_field& f )
{
  for (std::size_t i=0; i<field.nx(); ++i)
  {
    real x = field.x(i);
    if ( x < 1.5 )
      f[i] = std::tanh( double(15)*(x-1) );
    else
      f[i] = -std::tanh( double(15)*(x-2) );
  }
}

void Solver::eval_rhs( double t, const scalar_field& f, scalar_field& rhs ) 
{

  static scalar_field dfdx( f.size(), real(0) );
  
  grad( f, dfdx );

  for (int i=0; i<f.size(); i++)
    rhs[i] = -dfdx[i];

}

int main()
{

  Field field;
  RK48  rk48( field );
  Mac22 mac22( field );
  EMac24 emac24( field );

  std::size_t NX = field.nx();
 
  // check stability of time step
  real dt = tmax / max_iter;
  /*if ( dt > .163*field.dx() )
  {
    std::cout << "Error: dt is too large." << std::endl;
    std::cout << dt << " > " << .163*field.dx() << std::endl;
    return 0;
  }*/

  // allocate vector
  scalar_field U( NX, real(0) );

  // initial condition
  //init_sin( field, U );
  init_hat( field, U );
  //init_smooth_hat( field, U );

  // output initial condition
  std::ofstream out( "profile.out" );
  for (std::size_t i=0; i<NX; i++) out << field.x(i) << " ";
  out << std::endl;
  for (std::size_t i=0; i<NX; i++) out << U[i] << " ";
  out << std::endl;


  // main time loop
  double t = 0;
  for (int n=0; n<max_iter; n++)
  {

    //rk48.solve( t, U, dt );
    mac22.solve( t, U, dt );
    //emac24.solve( t, U, dt );
    
  } // time loop

  std::cout << t << " | " << 2*pi << std::endl; 

  // output final condition
  for (std::size_t i=0; i<NX; i++) out << U[i] << " ";
  out << std::endl;

  return 0;

}