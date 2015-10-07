

#include "unittest++/UnitTest++.h"


const double tolerance = 1e-10;

TEST(sanity)
{
  CHECK_EQUAL(1, 2 - 1);
} 
  
  
SUITE(functions)
{
  TEST(H)
  {
    double r = 5.0 + 3;
    r = r+r;
    CHECK_CLOSE( 0, 0, tolerance);
  }
    
}