#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <iostream>

#include <mcl/bls12_381.hpp>
using namespace mcl::bn;

class Polynomial {
public:
    Polynomial();
    Polynomial(const Fr& constant);
    Polynomial(const std::vector<Fr>& coefficients);
    Polynomial(const Polynomial& other);

    // 0 polynomial -> -1
    // otherwise, the degree of the polynomial
    int getDegree() const;

    Polynomial& operator=(const Polynomial& other);

    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;
    Polynomial operator/(const Polynomial& other) const;
    Polynomial operator%(const Polynomial& other) const;

    Polynomial operator-() const;

    Polynomial& operator+=(const Polynomial& other);
    Polynomial& operator-=(const Polynomial& other);
    Polynomial& operator*=(const Polynomial& other);
    Polynomial& operator/=(const Polynomial& other);
    Polynomial& operator%=(const Polynomial& other);

    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;
    Fr operator()(const Fr& x) const;
    Fr operator[](uint i) const;

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);

private:
    std::vector<Fr> coefficients;
    void removeLeadingZeros();
    void divide(const Polynomial& other, Polynomial& quotient, Polynomial& remainder) const;
};

#endif  // POLYNOMIAL_HPP
