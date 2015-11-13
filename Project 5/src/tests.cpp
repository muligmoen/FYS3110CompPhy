#include "catch.hpp"

#include "vector.hpp"


TEST_CASE( "SANITIY" )
{
  CHECK( 1 == 1 );
  CHECK( 2 == 3-1 );
}

TEST_CASE( "Vector tests" )
{
  const int N = 4; // Changing this leads to errors (Hard-coded matrix mult)
  Vector<int> vec(N);
  for (int ii=0; ii<N; ii++){
    vec[ii] = ii;
  }
  
  SECTION( "Indexing" )
  {
    for (int ii=0; ii<N; ii++){
      CHECK(vec[ii] == ii);
    }
  }
  
  SECTION( "Multiplying with 2x Identity matrix" )
  {
    auto doubleVec = vec.multiply(2, 0);
    for (int ii=0; ii<N; ii++){
      CHECK(doubleVec[ii] == 2*vec[ii]);
    }
  }
  
  SECTION( "Multiplying with a more advanced matrix" )
  {
    vec[0] = 5;
    vec[1] = 7;
    vec[2] = 8;
    vec[3] = 3;
    
    const int a = 2;
    const int b = 5;
    auto Avec = vec.multiply(2, 5);
  
    CHECK(Avec[0] == 45);
    CHECK(Avec[1] == 79);
    CHECK(Avec[2] == 66);
    CHECK(Avec[3] == 46);
    
  }
  
  
}