#include "gtest/gtest.h"
/*
#include "../defs.hpp"
#include "../field.hpp"
#include "../solver.hpp"
#include "../mac22.hpp"

void Solver::eval_rhs( double t, const scalar_field& f, scalar_field& rhs ) {}

void init_sin( const Field& field, scalar_field& f )
{
  for (std::size_t i=0; i<field.nx(); ++i) 
    f[i] = std::sin( field.x(i) );
}

namespace {

  // The fixture for testing class Foo.
  class Mac22Test : public ::testing::Test {
  protected:

    Mac22Test() : m_mac22( m_field ) {

      // You can do set-up work for each test here.
      m_f.resize( m_field.nx(), 0. );
      init_sin( m_field, m_f ); 

    }

    // Objects declared here can be used by all tests in the test case for Foo.
    Field        m_field;
    Mac22        m_mac22;
    scalar_field m_f;

  };

  // Tests that Mac22::grad returns the expected derivative
  TEST_F(Mac22Test, Grad) {
    scalar_field dfdx( m_field.nx(), 0. );
    for (std::size_t i=0; i<dfdx.size(); i++)
      ASSERT_NEAR( dfdx[i], std::cos( m_field.x(i) ), .5*m_field.dx() );
  }

}  // namespace
*/

TEST( SampleTest, Simple ) {
  EXPECT_EQ( 1, 1 );
}

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}