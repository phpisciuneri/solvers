#ifndef FIELD_HPP_
#define FIELD_HPP_

#include "defs.hpp"

class Field
{
public:

  Field()
  {
    const double pi = 3.14159265358979323846264338327;
    m_xmin = 0;
    m_xmax = 2*pi;
    m_nx = 201;
    m_dx = ( m_xmax - m_xmin ) / m_nx;
  }

  // accessors
  const real   xmin() const { return m_xmin; }
  const real   xmax() const { return m_xmax; }
  const real   dx()   const { return m_dx;   }
  const size_t nx()   const { return m_nx;   }

  const real x( std::size_t i ) const { return m_xmin + i*m_dx; }

private:

  real   m_xmin;
  real   m_xmax;
  real   m_dx;
  size_t m_nx;

};

#endif // FIELD_HPP_