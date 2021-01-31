/*
 * support.cpp
 *
 *  Created on: Jan 30, 2021
 *      Author: d-w-h
 */

#include <math.h>

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
