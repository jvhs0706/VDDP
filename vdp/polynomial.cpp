#include "polynomial.hpp"
#include <stdexcept>
#include <cassert>

using namespace std;

Polynomial::Polynomial() : coefficients() {}

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
    for (size_t i = 0; i < resultCoefficients.size(); ++i) {
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
    for (auto i = 0; i < resultCoefficients.size(); ++i) {
        Fr a = (i < coefficients.size()) ? coefficients[i] : Fr(0);
        Fr b = (i < other.coefficients.size()) ? other.coefficients[i] : Fr(0);
        resultCoefficients[i] = a - b;
    }
    Polynomial out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    std::vector<Fr> resultCoefficients(coefficients.size() + other.coefficients.size() - 1, Fr(0));
    for (size_t i = 0; i < coefficients.size(); ++i) {
        for (size_t j = 0; j < other.coefficients.size(); ++j) {
            resultCoefficients[i + j] += coefficients[i] * other.coefficients[j];
        }
    }
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
        for (int i = 0; i < coefficients.size(); i++)
        {
            quotient.coefficients.push_back(coefficients[i] / other.coefficients[0]);
        }
        remainder = Polynomial();
        return;
    }
    else if (coefficients.size() < other.coefficients.size()) {
        quotient = Polynomial();
        remainder = *this;
        remainder.removeLeadingZeros();
        return;
    }
    else{
        Polynomial temp = *this;

        Fr factor = temp.coefficients.back() / other.coefficients.back();
        if (factor != 0)
        {   
            for (int i = 0; i < other.coefficients.size(); i ++)
            {
                temp.coefficients[temp.coefficients.size() - (i+1)] -= factor * other.coefficients[other.coefficients.size() - (i+1)];
            }
        }
        temp.coefficients.pop_back();

        temp.divide(other, quotient, remainder);
        quotient.coefficients.push_back(factor);
    }
}

ostream& operator<< (ostream& os, const Polynomial& p) {
    if (p.coefficients.empty()) {
        os << "0";
        return os;
    }
    else {
        for (int i = p.getDegree(); i >= 0; --i) {
            if (i == 0) {
                os << p.coefficients[i];
            }
            else {
                os << p.coefficients[i] << " X^" << i << " + ";
            }
        }
        return os;
    }
}
