/*
 * complex.cpp
 *
 *  Created on: Jan 3, 2021
 *      Author: d-w-h
 */

#include "complex.hpp"

Complex::Complex() {}

Complex::Complex(double a, double b) {
    this->a = a;
    this->b = b;
}

Complex Complex::operator+(const Complex& m) {
    Complex result_addition(0,0);
    result_addition.a = this->a + m.a;
    result_addition.b = this->b + m.b;

    return result_addition;
}

Complex Complex::operator*(const Complex& m) {
    Complex result_multiplication(0,0);
    result_multiplication.a = this->a*m.a - this->b*m.b;
    result_multiplication.b = this->a*m.b + this->b*m.a;

    return result_multiplication;
}

Complex Complex::operator/(const Complex& m) {
    Complex result_division(0,0);
    double c2_d2 = m.a*m.a + m.b*m.b;
    result_division.a = (this->a*m.a + this->b*m.b) / (c2_d2 + 1e-20);
    result_division.b = (this->b*m.a - this->a*m.b) / (c2_d2 + 1e-20);

    return result_division;
}
