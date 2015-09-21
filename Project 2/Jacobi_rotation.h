#ifndef Jacobi_rotation_H
#define Jacobi_rotation_H

#include <armadillo>


void find_cos_sin(double a_kk, double a_ll, double a_kl, double &c, double &s);

void rotate(arma::Mat<double> &A, double c, double s, int k, int l);

void max_err_offdiag(const arma::Mat<double> &A, int &k, int &l, double &err);

double abs_sum_offdiag(const arma::Mat<double> &A);


#endif