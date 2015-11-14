#include "catch.hpp"

#include "vector.hpp"
//#include "integrate.hpp"


TEST_CASE( "SANITIY" )
{
  CHECK( 1 == 1 );
  CHECK( 2 == 3-1 );
}

TEST_CASE( "Simple vector functions" , "[Vector]")
{
  const int N = 4; // Changing this leads to errors (Hard-coded matrix mult)
  Vector<int> vec(N);
  for (int ii=0; ii<N; ii++){
    vec[ii] = ii;
  }
  
  SECTION( "Check size of vector" )
  {
    REQUIRE( N == vec.size() );
  }
  
  SECTION( "Indexing" )
  {
    for (int ii=0; ii<N; ii++){
      CHECK(vec[ii] == ii);
    }
  }
  
  SECTION( "Copy initializing to another vector" )
  {
    auto vec2 = vec; 
    for (int ii=0; ii<N; ii++){
      CHECK( vec2[ii] == vec[ii] );
    } 
  }
  
  SECTION( "Copying from one vector to the other" )
  {
    Vector<int> vec2(N);
    vec2 = vec;
    for (int ii=0; ii<N; ii++){
      CHECK( vec2[ii] == vec[ii] );
    } 
  }
  
  SECTION( "Negative of a vector" )
  {
    auto neg_vec = -vec;
    
    for (int ii=0; ii<N; ii++){
      CHECK( neg_vec[ii] == -vec[ii] );
    }
  }
  
  SECTION( "Add two vectors" )
  {
    auto vec2 = vec + vec;
    for (int ii=0; ii<N; ii++){
      CHECK(vec2[ii] == 2*vec[ii]);
    }
  }
  
}


TEST_CASE( "Initializing vector from function" )
{
  auto f = [](double x){return x*x;};
  
  const int N = 4;
  
  auto vec = init_vector(0, 3, N, f);
  
  for (int ii=0; ii<N; ii++){
    CHECK(ii*ii == Approx(vec[ii]));
  }
}

TEST_CASE( "Advanced vector tests", "[Vector]" )
{
  const int N = 4;
  Vector<int> vec(N);
  vec[0] = 5;
  vec[1] = 7;
  vec[2] = 8;
  vec[3] = 3;
  
  Vector<double> double_vec(N);
  double_vec[0] = vec[0];
  double_vec[1] = vec[1];
  double_vec[2] = vec[2];
  double_vec[3] = vec[3];
  
  const int a = 2;
  const int b = 5;
  
  SECTION( "Multiplying with 2x Identity matrix" )
  {
    auto doubleVec = multiply(vec, 2, 0);
    for (int ii=0; ii<N; ii++){
      CHECK(doubleVec[ii] == 2*vec[ii]);
    }
  }
  
  SECTION( "Multiplying with a more advanced matrix" )
  {
    auto Avec = multiply(vec, a, b);
  
    CHECK(Avec[0] == 45);
    CHECK(Avec[1] == 79);
    CHECK(Avec[2] == 66);
    CHECK(Avec[3] == 46); 
  }
  
  SECTION( "Multiplying inplace" )
  {
    multiply_inplace(vec, a, b);
    
    CHECK(vec[0] == 45);
    CHECK(vec[1] == 79);
    CHECK(vec[2] == 66);
    CHECK(vec[3] == 46); 
    
  }
  
  SECTION( "Inverse matrices")
  {
    const double ad = a;
    const double bd = b;
    Vector<double> sympy(N); // Solution from sympy
    sympy[0] = 0.879765395894428;
    sympy[1] = 0.648093841642229;
    sympy[2] = 0.26099706744868;
    sympy[3] = 0.847507331378299;
    
    SECTION( "Inverse matrix solver" )
    {
      auto Avec = solve(double_vec, ad, bd);
      
      for (int ii=0; ii<N; ii++){
        CHECK( Avec[ii] == Approx(sympy[ii]) );
      }
      
      // Using the multiplication
      auto avev = multiply(Avec, (double)a, (double)b);
      for (int ii=0; ii<N; ii++){
        CHECK( avev[ii] == Approx(double_vec[ii]) );
      }
    }
    
    SECTION( "Inplace matrix solver" )
    {
      auto new_vec = double_vec;
      
      solve_inplace(new_vec, (double)a, (double)b);
      
      for (int ii=0; ii<N; ii++){
        CHECK( new_vec[ii] == Approx(sympy[ii]) );
      }
    }
  }
    
    
}