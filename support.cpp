/*
 * support.cpp
 *
 *  Created on: Jan 30, 2021
 *      Author: d-w-h
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include "complex.hpp"
#include "mem_ops.hpp"
#include "user_types.hpp"

int factorial(int n) {

	int p = 1;

    for(int i = 1; i < n + 1; ++i) {
        p = p * i;
    }

    return p;
}

long double P_l(long double x, int l) {

	long double P_l_x = 0.0;

	for(int k = 0; k < l/2 + 1; ++k) {
        P_l_x += 1.0/pow(2.0, l) * pow(-1.0, k) * factorial(2*l - 2*k) * pow(x, l - 2*k) /
                (factorial(k) * factorial(l-k) * factorial(l-2*k));
    }

    return P_l_x;
}

long double dm_Pl_dxm(long double x, int m, int l) {

	long double f = 0.0;
    long double dx = 0.0001;

    if(m == 0) {
        return f = P_l(x, l);
    }
    else if(m == 1) {
        return f = (P_l(x+dx, l) - P_l(x, l))/dx;
    }
    else {
        f = (dm_Pl_dxm(x+dx, m-1, l) - dm_Pl_dxm(x, m-1, l))/dx;
    }

    return f;
}

double legendre_polynomial(int m, int l, double x) {

	double P_l_m_x = 0.0;

    if(m > -1) {
        return P_l_m_x = pow(-1.0, m) * pow(1.0 - x*x, 0.5*m) * dm_Pl_dxm(x, m, l);
    }
    else {
        m = -m;
        P_l_m_x = pow(-1.0, m) * factorial(l-m) / factorial(l+m) * legendre_polynomial(m, l, x);
    }

    return P_l_m_x;
}

long double laguerre_polynomial(long double x, int a, int k) {

	long double f = 0.0;

    if(k == 0) {
        return f = 1.0;
    }
    else if(k == 1) {
        return f = 1.0 + a - x;
    }
    else {
        f = ((2.0*(k-1) + 1.0 + a - x) * laguerre_polynomial(x, a, k-1) - (k-1+a) * laguerre_polynomial(x, a, k-2)) / k;
    }

    return f;
}

void init_r_and_phi(double* r, double* phi, d_data domain_data) {

    double dr = domain_data.R/domain_data.n_r;
    double dphi = M_PI/(domain_data.n_phi-1);

    for(int i = 0; i < domain_data.n_r; ++i) {
        r[i] = i*dr + 0.5*dr;
    }

    /* Initialize phi vector */
    for(int i = 0; i < domain_data.n_phi - 1; ++i) {
        phi[i] = i*dphi + 0.5*dphi;
    }
}

void init_domain_coordinates(d_data domain_data, s_data* solution_data) {

	double dr = domain_data.R/domain_data.n_r;
    double dtheta = 2*M_PI/domain_data.n_theta;
    double dphi = M_PI/(domain_data.n_phi-1);

    /* Initialize r_p vector */
    for(int i = 0; i < domain_data.n_r; ++i) {
        solution_data->r_p[i] = i*dr;
    }

    /* Initialize theta_p vector */
    for(int i = 0; i < domain_data.n_theta; ++i) {
        solution_data->theta_p[i] = i*dtheta + 0.5*dtheta;
    }

    /* Initialize phi_p vector */
    for(int i = 0; i < domain_data.n_phi; ++i) {
        solution_data->phi_p[i] = i*dphi;
    }
}

void compute_wave_function(d_data domain_data, p_params physical_params, q_nums quantum_numbers, s_data* solution_data) {

	for(int i = 0; i < 1; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                solution_data->psi[i][j][k].a = 0.0;
                solution_data->psi[i][j][k].b = 0.0;
            }
        }
    }

    for(int i = 1; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double rho = 2*solution_data->r_p[i]/(quantum_numbers.n*physical_params.a0);
                double beta = exp(-rho/2) * pow(rho, quantum_numbers.l) *
                              laguerre_polynomial(rho, 2*quantum_numbers.l+1, quantum_numbers.n-quantum_numbers.l-1) *
                              legendre_polynomial(quantum_numbers.m, quantum_numbers.l, cos(solution_data->phi_p[k]));
                Complex e_im_theta(cos(quantum_numbers.m*solution_data->theta_p[j]), sin(quantum_numbers.m*solution_data->theta_p[j]));
                solution_data->psi[i][j][k].a = e_im_theta.a * beta;
                solution_data->psi[i][j][k].b = e_im_theta.b * beta;
            }
        }
    }
}

void normalize_psi(double* r, d_data domain_data, s_data* solution_data) {

    Complex integral(0, 0);
    double dr = domain_data.R/domain_data.n_r;
    double dtheta = 2*M_PI/domain_data.n_theta;
    double dphi = M_PI/(domain_data.n_phi-1);

    /* Top pole */
    for(int i = 1; i < domain_data.n_r; ++i) {
        double dV = 2/3*M_PI*(1 - cos(dphi/2))*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1]);
        Complex psi_conj(solution_data->psi[i][0][0].a, -solution_data->psi[i][0][0].b);
        Complex psi_sq = solution_data->psi[i][0][0] * psi_conj;
        psi_sq.a = psi_sq.a*dV;
        psi_sq.b = psi_sq.b*dV;
        integral = integral + psi_sq;
    }

    /* Central nodes */
    for(int i = 1; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 1; k < domain_data.n_phi - 1; ++k) {
                double dV = solution_data->r_p[i]*solution_data->r_p[i]*sin(solution_data->phi_p[k])*dphi*dtheta*dr;
                Complex psi_conj(solution_data->psi[i][j][k].a, -solution_data->psi[i][j][k].b);
                Complex psi_sq = solution_data->psi[i][j][k] * psi_conj;
                psi_sq.a = psi_sq.a*dV;
                psi_sq.b = psi_sq.b*dV;
                integral = integral + psi_sq;
            }
        }
    }

    /* Bottom pole */
    for(int i = 1; i < domain_data.n_r; ++i) {
        double dV = 2/3*M_PI*(1 - cos(dphi/2))*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1]);
        Complex psi_conj(solution_data->psi[i][0][domain_data.n_phi-1].a, -solution_data->psi[i][0][domain_data.n_phi-1].b);
        Complex psi_sq = solution_data->psi[i][0][domain_data.n_phi-1] * psi_conj;
        psi_sq.a = psi_sq.a*dV;
        psi_sq.b = psi_sq.b*dV;
        integral = integral + psi_sq;
    }

    /* Normalize */
    for(int i = 1; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                solution_data->psi[i][j][k].a = solution_data->psi[i][j][k].a/sqrt(integral.a);
                solution_data->psi[i][j][k].b = solution_data->psi[i][j][k].b/sqrt(integral.a);
            }
        }
    }

}

void compute_probability_density(d_data domain_data, s_data *solution_data) {

    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                Complex psi_conj(solution_data->psi[i][j][k].a, -solution_data->psi[i][j][k].b);
                solution_data->psi_square[i][j][k] = solution_data->psi[i][j][k] * psi_conj;
            }
        }
    }
}

double compute_max_pd(d_data domain_data, s_data* solution_data) {

    double max_pd = -1e+8;

    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                if(solution_data->psi_square[i][j][k].a > max_pd) {
                    max_pd = solution_data->psi_square[i][j][k].a;
                }
            }
        }
    }

    return max_pd;
}

double compute_max_psi_real(d_data domain_data, s_data* solution_data) {

    double max_real = -1e+8;

    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                if(solution_data->psi[i][j][k].a > max_real) {
                    max_real = solution_data->psi[i][j][k].a;
                }
            }
        }
    }

    return max_real;
}

void export_grid_data(std::string file_name, d_data domain_data, double max_pd, double max_real) {

	std::ofstream myfile;
    myfile.open(file_name);
    myfile << domain_data.n_r << " "
    	   << domain_data.n_theta << " "
    	   <<  domain_data.n_phi << " "
    	   << max_pd << " "
    	   << max_real << "\n";
    myfile.close();
}

void export_probability_density(std::string file_name, d_data domain_data, s_data* solution_data) {

    std::ofstream myfile_pd;
    myfile_pd.open(file_name);

    /* Export psi central nodes */
    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = 0; j < 1; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double z = solution_data->r_p[i]*cos(solution_data->phi_p[k]);
                myfile_pd << x << " " << z << " " << solution_data->psi_square[i][j][k].a << "\n";
            }
        }
    }

    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = domain_data.n_theta/2; j < domain_data.n_theta/2 + 1; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double z = solution_data->r_p[i]*cos(solution_data->phi_p[k]);
                myfile_pd << x << " " << z << " " << solution_data->psi_square[i][j][k].a << "\n";
            }
        }
    }

    myfile_pd.close();
}

void export_psi_real(std::string file_name, d_data domain_data, s_data* solution_data) {

	std::ofstream myfile_real;
    myfile_real.open(file_name);

    int theta_slice = 0; // Should be a value between 0 and n_theta/2 - 2
    /* Export psi central nodes */
    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = theta_slice; j < theta_slice + 1; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double z = solution_data->r_p[i]*cos(solution_data->phi_p[k]);
                myfile_real << x << " " << z << " " << solution_data->psi[i][j][k].a << "\n";
            }
        }
    }

    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = domain_data.n_theta/2 + theta_slice; j < domain_data.n_theta/2 + theta_slice + 1; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double z = solution_data->r_p[i]*cos(solution_data->phi_p[k]);
                myfile_real << x << " " << z << " " << solution_data->psi[i][j][k].a << "\n";
            }
        }
    }

    myfile_real.close();
}

void export_psi_im(std::string file_name, d_data domain_data, s_data* solution_data) {

	std::ofstream myfile_im;
    myfile_im.open(file_name);

    int theta_slice = 0; // Should be a value between 0 and n_theta/2 - 2
    /* Export psi central nodes */
    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = theta_slice; j < theta_slice + 1; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double z = solution_data->r_p[i]*cos(solution_data->phi_p[k]);
                myfile_im << x << " " << z << " " << solution_data->psi[i][j][k].b << "\n";
            }
        }
    }

    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = domain_data.n_theta/2 + theta_slice; j < domain_data.n_theta/2 + theta_slice + 1; ++j) {
            for(int k = 0; k < domain_data.n_phi; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double z = solution_data->r_p[i]*cos(solution_data->phi_p[k]);
                myfile_im << x << " " << z << " " << solution_data->psi[i][j][k].b << "\n";
            }
        }
    }

    myfile_im.close();
}

void export_data_phi_2(std::string file_name, d_data domain_data, s_data* solution_data) {

    std::ofstream myfile_phi_2;
    myfile_phi_2.open(file_name);

    /* Export psi central nodes */
    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = domain_data.n_phi/2; k < domain_data.n_phi/2 + 1; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double y = solution_data->r_p[i]*sin(solution_data->phi_p[k])*sin(solution_data->theta_p[j]);
                myfile_phi_2 << x << " " << y << " " << solution_data->psi[i][j][k].a << "\n";
            }
        }
    }

    myfile_phi_2.close();
}

void export_data_phi_2_im(std::string file_name, d_data domain_data, s_data* solution_data) {

	std::ofstream myfile_phi_2_im;
    myfile_phi_2_im.open(file_name);

    /* Export psi central nodes */
    for(int i = 0; i < domain_data.n_r; ++i) {
        for(int j = 0; j < domain_data.n_theta; ++j) {
            for(int k = domain_data.n_phi/2; k < domain_data.n_phi/2 + 1; ++k) {
                double x = solution_data->r_p[i]*sin(solution_data->phi_p[k])*cos(solution_data->theta_p[j]);
                double y = solution_data->r_p[i]*sin(solution_data->phi_p[k])*sin(solution_data->theta_p[j]);
                myfile_phi_2_im << x << " " << y << " " << solution_data->psi[i][j][k].b << "\n";
            }
        }
    }

    myfile_phi_2_im.close();

}

void verify_solution(d_data domain_data, p_params physical_params, q_nums quantum_numbers, s_data* solution_data, double* r, double* phi) {

    Complex*** ratio = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
    double dr = domain_data.R/domain_data.n_r;
    double dtheta = 2*M_PI/domain_data.n_theta;
    double dphi = M_PI/(domain_data.n_phi-1);
    double alpha = physical_params.h/(2*physical_params.mp);
    double E = physical_params.ke/(quantum_numbers.n*quantum_numbers.n);

    for(int i = 1; i < domain_data.n_r - 1; ++i) {
        for(int j = 1; j < domain_data.n_theta - 1; ++j) {
            for(int k = 1; k < domain_data.n_phi - 1; ++k) {
                double dV1 = solution_data->r_p[i]*solution_data->r_p[i]*dr;
                double dV2 = solution_data->r_p[i]*solution_data->r_p[i]*sin(solution_data->phi_p[k])*sin(solution_data->phi_p[k])*dtheta;
                double dV3 = solution_data->r_p[i]*solution_data->r_p[i]*sin(solution_data->phi_p[k])*dphi;

                double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                double aN = r[i]*r[i]*alpha/(dV1*dr);
                double aW = alpha/(dV2*dtheta);
                double aE = alpha/(dV2*dtheta);
                double aB = sin(phi[k-1])*alpha/(dV3*dphi);
                double aT = sin(phi[k])*alpha/(dV3*dphi);

                Complex psi_s(0,0);
                psi_s.a = -aS*solution_data->psi[i-1][j][k].a;
                psi_s.b = -aS*solution_data->psi[i-1][j][k].b;

                Complex psi_n(0,0);
                psi_n.a = -aN*solution_data->psi[i+1][j][k].a;
                psi_n.b = -aN*solution_data->psi[i+1][j][k].b;

                Complex psi_w(0,0);
                psi_w.a = -aW*solution_data->psi[i][j-1][k].a;
                psi_w.b = -aW*solution_data->psi[i][j-1][k].b;

                Complex psi_e(0,0);
                psi_e.a = -aE*solution_data->psi[i][j+1][k].a;
                psi_e.b = -aE*solution_data->psi[i][j+1][k].b;

                Complex psi_b(0,0);
                psi_b.a = -aB*solution_data->psi[i][j][k-1].a;
                psi_b.b = -aB*solution_data->psi[i][j][k-1].b;

                Complex psi_t(0,0);
                psi_t.a = -aT*solution_data->psi[i][j][k+1].a;
                psi_t.b = -aT*solution_data->psi[i][j][k+1].b;

                Complex psi_p_k_rp(0,0);
                psi_p_k_rp.a = -physical_params.kp/solution_data->r_p[i]*solution_data->psi[i][j][k].a;
                psi_p_k_rp.b = -physical_params.kp/solution_data->r_p[i]*solution_data->psi[i][j][k].b;

                Complex Fp(0,0);
                Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + psi_p_k_rp;

                double ap = aS + aN + aW + aE + aB + aT - E;

                Complex ap_psi_p(0,0);
                ap_psi_p.a = -ap*solution_data->psi[i][j][k].a;
                ap_psi_p.b = -ap*solution_data->psi[i][j][k].b;

                ratio[i][j][k] = Fp / ap_psi_p;
            }
        }
    }

    /* Compute error */
    solution_data->error = 0.0;
    for(int i = 1; i < domain_data.n_r - 1; ++i) {
        for(int j = 1; j < domain_data.n_theta - 1; ++j) {
            for(int k = 1; k < domain_data.n_phi - 1; ++k) {
                if(fabs(ratio[i][j][k].a) > 0.2) {
                    solution_data->error = solution_data->error + fabs(1.0 - ratio[i][j][k].a);
                }
            }
        }
    }

    solution_data->error = 100 * solution_data->error / ((domain_data.n_r - 2)*(domain_data.n_theta - 2)*(domain_data.n_phi - 2));

}

void calc_psi(d_data domain_data,
              q_nums quantum_numbers,
              p_params physical_params,
              s_data* solution_data) {

    /* Allocate memory */
    double *r = new double[domain_data.n_r];
    double* phi = new double[domain_data.n_phi-1];

    /* Initialize r and phi vectors */
    init_r_and_phi(r, phi, domain_data);

    /* Initialize domain coordinate vectors */
    init_domain_coordinates(domain_data, solution_data);

    /* Compute psi */
    compute_wave_function(domain_data, physical_params, quantum_numbers, solution_data);

    /* Normalize psi */
    normalize_psi(r, domain_data, solution_data);

    /* Compute probability density */
    compute_probability_density(domain_data, solution_data);

    /* Compute max probability density */
    double max_pd = compute_max_pd(domain_data, solution_data);

    /* Compute max psi_real */
    double max_real = compute_max_psi_real(domain_data, solution_data);

    /* Export domain data */
    export_grid_data("grid_data.txt", domain_data, max_pd, max_real);

    /* Export probability density data */
    export_probability_density("data.txt", domain_data, solution_data);

    /* Export psi real data */
    export_psi_real("data_real.txt", domain_data, solution_data);

    /* Export psi im data */
    export_psi_im("data_im.txt", domain_data, solution_data);

    /* Export psi real data at n_phi/2 */
    export_data_phi_2("data_real_phi_2.txt", domain_data, solution_data);

    /* Export psi im data at n_phi/2 */
    export_data_phi_2_im("data_im_phi_2.txt", domain_data, solution_data);

    /* Verify solution numerically */
    verify_solution(domain_data, physical_params, quantum_numbers, solution_data, r, phi);

    /* Deallocate memory */
    delete [] r;
    delete [] phi;

}
