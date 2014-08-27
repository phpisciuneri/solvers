#ifndef RK48_HPP_
#define RK48_HPP_

#include "defs.hpp"
#include "field.hpp"
#include "solver.hpp"

class RK48 : public Solver
{
public:

  RK48( const Field& field ) : m_dx( field.dx() ), m_nx( field.nx() ) {}

  void grad( const scalar_field& f, scalar_field& dfdx );

  void gradv( const scalar_field& f, scalar_field& dfdx ) {}

  void solve( real& t, scalar_field& f, real dt );

private:

  real   m_dx;
  size_t m_nx;

};

#endif // RK48_HPP_