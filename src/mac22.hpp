#ifndef MAC22_HPP_
#define MAC22_HPP_

#include "defs.hpp"
#include "field.hpp"
#include "solver.hpp"

class Mac22 : public Solver
{
public:

  Mac22( const Field& field ) : m_dx( field.dx() ), m_nx( field.nx() ) {}

  void grad( const scalar_field& f, scalar_field& dfdx );

  void gradv( const scalar_field& f, scalar_field& dfdx ) {}

  void solve( real& t, scalar_field& f, real dt );

private:

  // functions
  void forward_diff( const scalar_field& f, scalar_field& df );
  void backward_diff( const scalar_field& f, scalar_field& df );

  real   m_dx;
  std::size_t m_nx;

};

#endif // MAC22_HPP_