/*
 * main.cpp
 *
 *  Created on: Jan 29, 2021
 *      Author: d-w-h
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include "complex.hpp"
#include "mem_ops.hpp"
#include "support.hpp"

int main(int argc, char* argv[]) {

    /* Parameters */
    int n_phi = 20;     //Number of nodes in phi direction
    int n_theta = 40;   //Number of nodes in theta direction
    int n_r = 200;      //Number of nodes in r direction

    double R = 40.0;    //Length of domain

    int n = 4;          //Principal quantum number
    int l = 2;          //Azimuthal number
    int m = 1;          //Magnetic quantum number

    double a0 = 1.0;    //Bohr radius

    /* Allocate memory */
    Complex*** psi = mat3D(n_r, n_theta, n_phi);
    Complex*** psi_square = mat3D(n_r, n_theta, n_phi);
    double *r = new double[n_r];
    double* r_p = new double[n_r];
    double* theta_p = new double[n_theta];
    double* phi_p = new double[n_phi];

    /* Some calculations */
    double dr = R/n_r;
    double dtheta = 2*M_PI/n_theta;
    double dphi = M_PI/(n_phi-1);

    /* Initialize r vector */
    for(int i = 0; i < n_r; ++i) {
        r[i] = i*dr + 0.5*dr;
    }

    /* Initialize r_p vector */
    for(int i = 0; i < n_r; ++i) {
        r_p[i] = i*dr;
    }

    /* Initialize theta_p vector */
    for(int i = 0; i < n_theta; ++i) {
        theta_p[i] = i*dtheta + 0.5*dtheta;
    }

    /* Initialize phi_p vector */
    for(int i = 0; i < n_phi; ++i) {
        phi_p[i] = i*dphi;
    }

    /* Compute psi */
    for(int i = 0; i < 1; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                psi[i][j][k].a = 0.0;
                psi[i][j][k].b = 0.0;
            }
        }
    }

    for(int i = 1; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double rho = 2*r_p[i]/(n*a0);
                double beta = exp(-rho/2) * pow(rho, l) * laguerre_polynomial(rho, 2*l+1, n-l-1) * legendre_polynomial(m, l, cos(phi_p[k]));
                Complex e_im_theta(cos(m*theta_p[j]), sin(m*theta_p[j]));
                psi[i][j][k].a = e_im_theta.a * beta;
                psi[i][j][k].b = e_im_theta.b * beta;
            }
        }
    }

    /* Normalize psi */
    Complex integral(0, 0);

    /* Top pole */
    for(int i = 1; i < n_r; ++i) {
        double dV = 2/3*M_PI*(1 - cos(dphi/2))*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1]);
        Complex psi_conj(psi[i][0][0].a, -psi[i][0][0].b);
        Complex psi_sq = psi[i][0][0] * psi_conj;
        psi_sq.a = psi_sq.a*dV;
        psi_sq.b = psi_sq.b*dV;
        integral = integral + psi_sq;
    }

    /* Central nodes */
    for(int i = 1; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 1; k < n_phi - 1; ++k) {
                double dV = r_p[i]*r_p[i]*sin(phi_p[k])*dphi*dtheta*dr;
                Complex psi_conj(psi[i][j][k].a, -psi[i][j][k].b);
                Complex psi_sq = psi[i][j][k] * psi_conj;
                psi_sq.a = psi_sq.a*dV;
                psi_sq.b = psi_sq.b*dV;
                integral = integral + psi_sq;
            }
        }
    }

    /* Bottom pole */
    for(int i = 1; i < n_r; ++i) {
        double dV = 2/3*M_PI*(1 - cos(dphi/2))*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1]);
        Complex psi_conj(psi[i][0][n_phi-1].a, -psi[i][0][n_phi-1].b);
        Complex psi_sq = psi[i][0][n_phi-1] * psi_conj;
        psi_sq.a = psi_sq.a*dV;
        psi_sq.b = psi_sq.b*dV;
        integral = integral + psi_sq;
    }

    /* Normalize */
    for(int i = 1; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                psi[i][j][k].a = psi[i][j][k].a/sqrt(integral.a);
                psi[i][j][k].b = psi[i][j][k].b/sqrt(integral.a);
            }
        }
    }

    /* Compute probability density */
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                Complex psi_conj(psi[i][j][k].a, -psi[i][j][k].b);
                psi_square[i][j][k] = psi[i][j][k] * psi_conj;
            }
        }
    }

    /* Compute max probability density */
    double max_pd = -1e+8;
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                if(psi_square[i][j][k].a > max_pd) {
                    max_pd = psi_square[i][j][k].a;
                }
            }
        }
    }

    /* Compute max psi_real */
    double max_real = -1e+8;
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                if(psi[i][j][k].a > max_real) {
                    max_real = psi[i][j][k].a;
                }
            }
        }
    }

    /* Export some data */
    std::ofstream myfile;
    std::string file_name = "grid_data.txt";
    myfile.open(file_name);
    myfile << n_r << " " << n_theta << " " <<  n_phi << " " << max_pd << " " << max_real << "\n";
    myfile.close();

    /* Export probability density data */
    std::ofstream myfile_pd;
    file_name = "data.txt";
    myfile_pd.open(file_name);

    /* Export psi central nodes */
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < 1; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double z = r_p[i]*cos(phi_p[k]);
                myfile_pd << x << " " << z << " " << psi_square[i][j][k].a << "\n";
            }
        }
    }

    for(int i = 0; i < n_r; ++i) {
        for(int j = n_theta/2; j < n_theta/2 + 1; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double z = r_p[i]*cos(phi_p[k]);
                myfile_pd << x << " " << z << " " << psi_square[i][j][k].a << "\n";
            }
        }
    }

    myfile_pd.close();

    /* Export psi real data */
    std::ofstream myfile_real;
    file_name = "data_real.txt";
    myfile_real.open(file_name);

    int theta_slice = 0;//n_theta/4; // Should be a value between 0 and n_theta/2 - 2
    /* Export psi central nodes */
    for(int i = 0; i < n_r; ++i) {
        for(int j = theta_slice; j < theta_slice + 1; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double z = r_p[i]*cos(phi_p[k]);
                myfile_real << x << " " << z << " " << psi[i][j][k].a << "\n";
            }
        }
    }

    for(int i = 0; i < n_r; ++i) {
        for(int j = n_theta/2 + theta_slice; j < n_theta/2 + theta_slice + 1; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double z = r_p[i]*cos(phi_p[k]);
                myfile_real << x << " " << z << " " << psi[i][j][k].a << "\n";
            }
        }
    }

    myfile_real.close();

    /* Export psi im data */
    std::ofstream myfile_im;
    file_name = "data_im.txt";
    myfile_im.open(file_name);

    theta_slice = n_theta/8; // Should be a value between 0 and n_theta/2 - 2
    /* Export psi central nodes */
    for(int i = 0; i < n_r; ++i) {
        for(int j = theta_slice; j < theta_slice + 1; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double z = r_p[i]*cos(phi_p[k]);
                myfile_im << x << " " << z << " " << psi[i][j][k].b << "\n";
            }
        }
    }

    for(int i = 0; i < n_r; ++i) {
        for(int j = n_theta/2 + theta_slice; j < n_theta/2 + theta_slice + 1; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double z = r_p[i]*cos(phi_p[k]);
                myfile_im << x << " " << z << " " << psi[i][j][k].b << "\n";
            }
        }
    }

    myfile_im.close();

    /* Export psi real data at n_phi/2 */
    std::ofstream myfile_phi_2;
    file_name = "data_real_phi_2.txt";
    myfile_phi_2.open(file_name);

    /* Export psi central nodes */
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = n_phi/2; k < n_phi/2 + 1; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double y = r_p[i]*sin(phi_p[k])*sin(theta_p[j]);
                myfile_phi_2 << x << " " << y << " " << psi[i][j][k].a << "\n";
            }
        }
    }

    myfile_phi_2.close();

    /* Export psi im data at n_phi/2 */
    std::ofstream myfile_phi_2_im;
    file_name = "data_real_phi_2_im.txt";
    myfile_phi_2_im.open(file_name);

    /* Export psi central nodes */
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = n_phi/2; k < n_phi/2 + 1; ++k) {
                double x = r_p[i]*sin(phi_p[k])*cos(theta_p[j]);
                double y = r_p[i]*sin(phi_p[k])*sin(theta_p[j]);
                myfile_phi_2_im << x << " " << y << " " << psi[i][j][k].b << "\n";
            }
        }
    }

    myfile_phi_2_im.close();

    /* Print results */
    for(int i = 0; i < n_r; ++i) {
        printf("psi_top_real[%i]: %f, psi_top_im[%i]: %f, psi_bottom_real[%i]: %f, psi_bottom_im[%i]: %f, real: %f, im: %f\n",
                i, psi[i][0][0].a, i, psi[i][0][0].b, i, psi[i][0][n_phi-1].a, i, psi[i][0][n_phi-1].b, psi[i][n_theta/2][n_phi/2].a, psi[i][n_theta/4][n_phi/4].b);
    }

    return 0;
}

