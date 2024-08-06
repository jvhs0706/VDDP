#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <iostream>

#include <mcl/bls12_381.hpp>
using namespace mcl::bn;

class VanishingPolynomial;

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

    Polynomial operator/(const VanishingPolynomial& other) const;
    Polynomial operator%(const VanishingPolynomial& other) const;

    Polynomial operator-() const;

    Polynomial& operator+=(const Polynomial& other);
    Polynomial& operator-=(const Polynomial& other);
    Polynomial& operator*=(const Polynomial& other);
    Polynomial& operator/=(const Polynomial& other);
    Polynomial& operator%=(const Polynomial& other);

    Polynomial& operator/=(const VanishingPolynomial& other);
    Polynomial& operator%=(const VanishingPolynomial& other);

    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;
    Fr operator()(const Fr& x) const;
    Fr operator[](uint i) const;

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);

protected:
    std::vector<Fr> coefficients;
    void removeLeadingZeros();
    void divide(const Polynomial& other, Polynomial& quotient, Polynomial& remainder) const;
    void divide(const VanishingPolynomial& other, Polynomial& quotient, Polynomial& remainder) const;
};

// X^len - 1
class VanishingPolynomial : public Polynomial {
public:
    VanishingPolynomial(uint len);
    VanishingPolynomial(const VanishingPolynomial& other);
    Fr operator()(const Fr& x) const;

    // ban the following operations: +=, -=, *=, /=, %=, 
    // because the result may not be a vanishing polynomial
    Polynomial& operator+=(const Polynomial& other) = delete;
    Polynomial& operator-=(const Polynomial& other) = delete;
    Polynomial& operator*=(const Polynomial& other) = delete;
    Polynomial& operator/=(const Polynomial& other) = delete;
    Polynomial& operator%=(const Polynomial& other) = delete;
    Polynomial& operator/=(const VanishingPolynomial& other) = delete;
    Polynomial& operator%=(const VanishingPolynomial& other) = delete;
    
};

#endif  // POLYNOMIAL_HPP
