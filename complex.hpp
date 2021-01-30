/*
 * complex.hpp
 *
 *  Created on: Jan 3, 2021
 *      Author: d-w-h
 */

#ifndef COMPLEX_HPP_
#define COMPLEX_HPP_

class Complex {
public:
    double a;
    double b;
    Complex();
    Complex(double a, double b);
    Complex operator+(const Complex& m);
    Complex operator*(const Complex& m);
};

#endif /* COMPLEX_HPP_ */
