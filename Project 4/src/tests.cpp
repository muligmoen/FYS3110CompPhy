#include "catch.hpp"

#include "lattice.hpp"

#include <cmath>

/*! \file tests.cpp
 * \brief Contains unit tests of lattice.hpp
 */

//! Testing if catch makes a TEST_CASE
TEST_CASE( "SANITIY" )
{
  CHECK( 1 == 1 );
  CHECK( 2 == 3-1 );
}

//! Testing basic operations on a lattice
TEST_CASE( "Setting lattice, indexing and changing elements" , "[lattice]")
{
  int max_x = 2;
  int max_y = 4;
  Lattice lattice(max_x, max_y, init::up);
  
  SECTION( "Check successful init") {
    for (int ii=0; ii<max_y; ii++){
      for (int jj=0; jj<max_x; jj++){
        REQUIRE( lattice(ii,jj) == 1 );
      }
    }
  }
  
  SECTION( "Change element by index") {
    lattice(0, 1) = -5;
    CHECK( lattice(0,1) == -5 );
  }
  
  SECTION( "Indexing periodic boundary" ) {
    lattice(1, 0) = -5;
    lattice(1, max_y - 1) = 3; // last y-element
    
    SECTION( "Wrapping around x-axis" ){
      CHECK(lattice(1+max_x, 0) == -5);
    }
    SECTION( "Wrapping around y-axis" ){
      CHECK(lattice(1, max_y) == -5);
    }
    SECTION( "Wrapping around negative values" ){
      CHECK(lattice(1, -1) == 3); // wrap around negative
    }
  }
}

//! Testing energy and magnetism of a lattice
TEST_CASE( "Energies and spins in lattice", "[lattice]" )
{
  Lattice lattice(2, 2, init::up);
  
  SECTION( "Check total spins" ) {
    CHECK( lattice.sum_spins() == 4 );
    lattice(1,0) = -5;
    CHECK( lattice.sum_spins() == 3-5 );
  }
  
  SECTION( "Check total energy" ) {
    CHECK( lattice.energy() == -4*4 ); // 4 spins with 4 energy 
  }
  SECTION( "Change in energy by flipping single spin" ) {
    const int deltaE = lattice.dE(1,1);
    const int Ebefore = lattice.energy();
    lattice(1,1) *= -1; //flipping
    const int Eafter = lattice.energy();
    CHECK( deltaE == Eafter - Ebefore );
    CHECK( lattice.dE(1,1) == Ebefore - Eafter ); // Reverse flip
  }
}
    
