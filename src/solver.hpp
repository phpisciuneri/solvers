#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include "defs.hpp"

class Solver 
{
public:

  // convective derivative
  virtual void grad(  const scalar_field&, scalar_field& ) = 0;
  
  // viscous derivative
  virtual void gradv( const scalar_field&, scalar_field& ) = 0;
  
  // solve
  virtual void solve( real& t, scalar_field& U, real dt ) = 0;

  void eval_rhs( real t, const scalar_field& U, scalar_field& rhs );

private:


};

#endif // SOLVER_HPP