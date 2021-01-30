/*
 * mem_ops.cpp
 *
 *  Created on: Jan 3, 2021
 *      Author: d-w-h
 */

#include "complex.hpp"

Complex*** mat3D(int n_r, int n_theta, int n_phi) {

    Complex*** f = new Complex**[n_r];

    for(int i = 0; i < n_r; ++i) {
        f[i] = new Complex*[n_theta];
    }

    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            f[i][j] = new Complex[n_phi];
        }
    }

    return f;
}

void free_mat3D(Complex*** f, int n_r, int n_theta) {

    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            delete [] f[i][j];
        }
    }

    for(int i = 0; i < n_r; ++i) {
        delete [] f[i];
    }

    delete [] f;

}
