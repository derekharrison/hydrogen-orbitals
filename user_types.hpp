/*
 * user_types.hpp
 *
 *  Created on: Feb 3, 2021
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include "complex.hpp"

typedef struct domain_data {
    int n_r;
    int n_theta;
    int n_phi;
    double R;
} d_data;

typedef struct quantum_numbers {
    int n;
    int l;
    int m;
} q_nums;

typedef struct physical_parameters {
    double a0;
    double h;
    double mp;
    double kp;
    double ke;
} p_params;

typedef struct solution_data {
    Complex*** psi;
    Complex*** psi_square;
    double* r_p;
    double* theta_p;
    double* phi_p;
    double E;
    double error;

} s_data;

#endif /* USER_TYPES_HPP_ */
