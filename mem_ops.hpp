/*
 * mem_ops.hpp
 *
 *  Created on: Jan 3, 2021
 *      Author: d-w-h
 */

#ifndef MEM_OPS_HPP_
#define MEM_OPS_HPP_

#include "complex.hpp"

Complex*** mat3D(int n_r, int n_theta, int n_phi);
void free_mat3D(Complex*** f, int n_r, int n_theta);

#endif /* MEM_OPS_HPP_ */
