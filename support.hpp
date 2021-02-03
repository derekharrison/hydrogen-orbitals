/*
 * support.hpp
 *
 *  Created on: Jan 30, 2021
 *      Author: d-w-h
 */

#ifndef SUPPORT_HPP_
#define SUPPORT_HPP_

#include "user_types.hpp"

int factorial(int n);
long double P_l(long double x, int l);
long double dm_dxm(long double x, int m, int l);
double legendre_polynomial(int m, int l, double x);
long double laguerre_polynomial(long double x, int a, int k);
void calc_psi(d_data domain_data,
              q_nums quantum_numbers,
              p_params physical_params,
              s_data* solution_data);

#endif /* SUPPORT_HPP_ */
