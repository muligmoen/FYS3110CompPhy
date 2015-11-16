#ifndef integrate_hpp
#define integrate_hpp

#include "vector.hpp"

/*! \file integrate.hpp
 * 
 * @brief Functions to solve the diffusion problem
 * 
 * These functions solves the diffusion problem with schemes based
 * on tridiagonal matrix problems
 */

//! Contains functions to solve the diffusion equation
/*!
 * The Euler and Crank-Nocolson methods use matrices to compute the next
 * derivative, and there is also a MC method
 */
namespace diffusion{

  //! Does forward euler steps on the initial vector
  /*!
  * returns a vector where 'steps' euler steps has been completed.
  * s is \f$ \frac{\Delta t}{\Delta x^2} \f$.
  * 
  * Boundary conditions are set to 0 on both ends
  * 
  * The internal method is based on a multiplication
  * with a sparse matrix
  */
  Vector<double> forward_euler(const Vector<double> &init_vec,
                              const double s, const int steps);

  //! Does backwards euler steps on the initial vector
  /*!
  * returns a vector where 'steps' euler steps has been completed.
  * s is \f$ \frac{\Delta t}{\Delta x^2} \f$.
  * 
  * Boundary conditions are set to 0 on both ends
  * 
  * The internal method is based on an inverse multiplication
  * with a sparse matrix
  */
  Vector<double> backward_euler(const Vector<double> &init_vec,
                                const double s, const int steps);

  //! Does Crank-Nicholson steps on the initial vector
  /*!
  * 
  * 
  * The internal method is based on an inverse multiplication
  * with a sparse matrix
  */
  Vector<double> Crank_Nicolson(const Vector<double> &init_vec,
                                const double s, const int steps);
  }
#endif