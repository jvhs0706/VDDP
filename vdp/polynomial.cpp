#include "polynomial.hpp"
#include <stdexcept>
#include <cassert>
#include <execution>

#include "utils.hpp"

using namespace std;

Polynomial::Polynomial() : coefficients(0) {}

Polynomial::Polynomial(const Fr& constant) : coefficients({constant}) {
    removeLeadingZeros();
}

Polynomial::Polynomial(const std::vector<Fr>& coefficients) : coefficients(coefficients) {
    removeLeadingZeros();
}

Polynomial::Polynomial(const Polynomial& other) : coefficients(other.coefficients) {
    removeLeadingZeros();
}

int Polynomial::getDegree() const {
    return (signed) coefficients.size() - 1;
}

Polynomial& Polynomial::operator=(const Polynomial& other) {
    coefficients = other.coefficients;
    return *this;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    std::vector<Fr> resultCoefficients(std::max(coefficients.size(), other.coefficients.size()));
    for (uint i = 0; i < resultCoefficients.size(); ++i) {
        Fr a = (i < coefficients.size()) ? coefficients[i] : Fr(0);
        Fr b = (i < other.coefficients.size()) ? other.coefficients[i] : Fr(0);
        resultCoefficients[i] = a + b;
    }
    Polynomial out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    std::vector<Fr> resultCoefficients(std::max(coefficients.size(), other.coefficients.size()));
    for (uint i = 0; i < resultCoefficients.size(); ++i) {
        Fr a = (i < coefficients.size()) ? coefficients[i] : Fr(0);
        Fr b = (i < other.coefficients.size()) ? other.coefficients[i] : Fr(0);
        resultCoefficients[i] = a - b;
    }
    Polynomial out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    if (coefficients.empty() || other.coefficients.empty()) {
        return Polynomial();
    }
    auto resultDegree = coefficients.size() + other.coefficients.size() - 1;
    auto log_ntt_size = ceilLog2(resultDegree + 1);
    auto ntt_size = 1 << log_ntt_size;
    auto omega = getRootOfUnity(log_ntt_size);

    auto this_coef_extend = coefficients;
    this_coef_extend.resize(ntt_size, Fr(0));
    auto other_coef_extend = other.coefficients;
    other_coef_extend.resize(ntt_size, Fr(0));

    auto this_ntt = ntt_helper(this_coef_extend, omega);
    auto other_ntt = ntt_helper(other_coef_extend, omega);

    vector<Fr> result_ntt(ntt_size);
    for_each(
        // std::execution::par,
        result_ntt.begin(),
        result_ntt.end(),
        [&](Fr& v)
        {
            auto idx = &v - &result_ntt.front();
            v = this_ntt[idx] * other_ntt[idx];
        }
    );

    auto resultCoefficients = ntt_helper(result_ntt, 1 / omega);
    for_each(
        // std::execution::par,
        resultCoefficients.begin(),
        resultCoefficients.end(),
        [&](Fr& v)
        {
            v /= ntt_size;
        }
    );

    Polynomial out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

Polynomial Polynomial::operator/(const Polynomial& other) const {
    Polynomial quotient, remainder;
    divide(other, quotient, remainder);
    return quotient;
}

Polynomial Polynomial::operator%(const Polynomial& other) const {
    Polynomial quotient, remainder;
    divide(other, quotient, remainder);
    return remainder;
}

Polynomial Polynomial::operator-() const {
    std::vector<Fr> resultCoefficients(coefficients.size());
    for (size_t i = 0; i < coefficients.size(); ++i) {
        resultCoefficients[i] = -coefficients[i];
    }
    return Polynomial(resultCoefficients);
}

Polynomial& Polynomial::operator+=(const Polynomial& other) {
    *this = *this + other;
    return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& other) {
    *this = *this - other;
    return *this;
}

Polynomial& Polynomial::operator*=(const Polynomial& other) {
    *this = *this * other;
    return *this;
}

Polynomial& Polynomial::operator/=(const Polynomial& other) {
    *this = *this / other;
    return *this;
}

Polynomial& Polynomial::operator%=(const Polynomial& other) {
    *this = *this % other;
    return *this;
}

Fr Polynomial::operator()(const Fr& x) const {
    Fr result = Fr(0);
    Fr xPower = Fr(1);
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * xPower;
        xPower *= x;
    }
    return result;
}

Fr Polynomial::operator[](uint i) const {
    assert (i < coefficients.size());
    return coefficients[i];
}

bool Polynomial::operator==(const Polynomial& other) const
{
    for (size_t i = 0; i < coefficients.size(); ++i) {
        if (coefficients[i] != other.coefficients[i]) {
            return false;
        }
    }
    return true;
}

bool Polynomial::operator!=(const Polynomial& other) const
{
    return !(*this == other);
}

void Polynomial::removeLeadingZeros() {
    while (!coefficients.empty() && coefficients.back() == Fr(0)) {
        coefficients.pop_back();
    }
}

void Polynomial::divide(const Polynomial& other, Polynomial& quotient, Polynomial& remainder) const {
    if (other.coefficients.empty()) {
        throw std::invalid_argument("Division by zero");
    }
    else if (other.coefficients.size() == 1){
        quotient.coefficients.resize(coefficients.size());
        for (int i = 0; i < coefficients.size(); i++)
        {
            quotient.coefficients[i] = coefficients[i] / other.coefficients[0];
        }
        remainder = Polynomial();
        return;
    }
    else {
        quotient.coefficients.resize(max(this -> getDegree() - other.getDegree() + 1, 0));
        remainder = *this;

        for (int i = this -> getDegree(); i >= other.getDegree(); i--)
        {
            quotient.coefficients[i - other.getDegree()] = remainder.coefficients[i] / other.coefficients.back();
            const Fr& factor = quotient.coefficients[i - other.getDegree()];
            if (!factor.isZero()) {
                for (int j = 0; j < other.coefficients.size(); j ++)
                {
                    remainder.coefficients[remainder.coefficients.size() - (j+1)] -= factor * other.coefficients[other.coefficients.size() - (j+1)];
                }
            }
            remainder.coefficients.pop_back();
        }

        remainder.removeLeadingZeros();
        return;
    }
}

ostream& operator<< (ostream& os, const Polynomial& p) {
    if (p.coefficients.empty()) {
        os << "0";
        return os;
    }
    else {
        for (int i = p.getDegree(); i >= 0; --i) {
            if (p.coefficients[i].isZero()) continue;
            if (i < p.getDegree()) os << " + ";
            os << p.coefficients[i];
            if (i > 0) os << " * X";
            if (i > 1) os << "^" << i;
        }
        return os;
    }
}

VanishingPolynomial::VanishingPolynomial(uint len)  {
    if (len == 0) {
        throw std::invalid_argument("Vanishing polynomial degree must be positive");
    }
    coefficients.resize(len + 1, Fr(0));
    coefficients[0] = Fr(-1);
    coefficients[len] = Fr(1);
}

VanishingPolynomial::VanishingPolynomial(const VanishingPolynomial& other) : Polynomial(other) {}

Fr VanishingPolynomial::operator()(const Fr& x) const {
    Fr res;
    Fr::pow(res, x, getDegree());
    return res - Fr(1);
}

void Polynomial::divide(const VanishingPolynomial& other, Polynomial& quotient, Polynomial& remainder) const
{
    if (coefficients.empty()) {
        quotient = Polynomial();
        remainder = Polynomial();
        return;
    }
    else {
        quotient.coefficients.resize(max(this -> getDegree() - other.getDegree() + 1, 0));
        remainder = *this;

        for (int i = this -> getDegree(); i >= other.getDegree(); i--)
        {
            const Fr& factor = remainder.coefficients[i];
            if (!factor.isZero()) {
                remainder.coefficients[i - other.getDegree()] += factor;
            }
            remainder.coefficients.pop_back();
            quotient.coefficients[i - other.getDegree()] = factor;
        }

        remainder.removeLeadingZeros();
        return;
    } 
}

Polynomial Polynomial::operator/(const VanishingPolynomial& other) const
{
    Polynomial quotient, remainder;
    divide(other, quotient, remainder);
    return quotient;
}

Polynomial Polynomial::operator%(const VanishingPolynomial& other) const
{
    Polynomial quotient, remainder;
    divide(other, quotient, remainder);
    return remainder;
}

Polynomial& Polynomial::operator/=(const VanishingPolynomial& other)
{
    *this = *this / other;
    return *this;
}

Polynomial& Polynomial::operator%=(const VanishingPolynomial& other)
{
    *this = *this % other;
    return *this;
}

