#ifndef EMAC24_HPP_
#define EMAC24_HPP_

#include "defs.hpp"
#include "field.hpp"
#include "solver.hpp"

class EMac24 : public Solver
{
public:

  EMac24( const Field& field ) : m_dx( field.dx() ), m_nx( field.nx() ) {}

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

#endif // EMAC24_HPP_