/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include <cmath>
#include "complex.hpp"
#include "mem_ops.hpp"
#include "support.hpp"
#include "user_types.hpp"

int main(int argc, char* argv[]) {

    d_data domain_data;
    q_nums quantum_numbers;
    p_params physical_params;
    s_data solution_data;

    /* Parameters */
    domain_data.n_phi = 20;     //Number of nodes in phi direction
    domain_data.n_theta = 40;   //Number of nodes in theta direction
    domain_data.n_r = 200;      //Number of nodes in r direction
    domain_data.R = 40.0;       //Length of domain

    quantum_numbers.n = 4;      //Principal quantum number
    quantum_numbers.l = 2;      //Azimuthal number
    quantum_numbers.m = 1;      //Magnetic quantum number

    physical_params.a0 = 1.0;    //Bohr radius
    physical_params.h = 1.0;     //Constant
    physical_params.mp = 1.0;    //Mass of particle
    physical_params.kp = 1.0;    //Potential constant

    /* Allocate memory for solution data */
    solution_data.psi = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
    solution_data.psi_square = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
    solution_data.r_p = new double[domain_data.n_r];
    solution_data.theta_p = new double[domain_data.n_theta];
    solution_data.phi_p = new double[domain_data.n_phi];

    /* Compute wave function */
    calc_psi(domain_data,
             quantum_numbers,
             physical_params,
             &solution_data);

    /* Print some results */
    for(int i = 0; i < domain_data.n_r; ++i) {
        printf("psi_top_real[%i]: %f, psi_top_im[%i]: %f, psi_bottom_real[%i]: %f, psi_bottom_im[%i]: %f, real: %f, im: %f\n",
                i, solution_data.psi[i][0][0].a, i, solution_data.psi[i][0][0].b, i, solution_data.psi[i][0][domain_data.n_phi-1].a,
                i, solution_data.psi[i][0][domain_data.n_phi-1].b, solution_data.psi[i][domain_data.n_theta/2][domain_data.n_phi/2].a, solution_data.psi[i][domain_data.n_theta/4][domain_data.n_phi/4].b);
    }

    printf("error in percentage : %f\n", solution_data.error);

    return 0;
}

