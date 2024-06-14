#include "polynomial.hpp"
#include <stdexcept>
#include <cassert>

using namespace std;

template <typename T>
Polynomial<T>::Polynomial() : coefficients() {}

template <typename T>
Polynomial<T>::Polynomial(const T& constant) : coefficients({constant}) {
    removeLeadingZeros();
}

template <typename T>
Polynomial<T>::Polynomial(const std::vector<T>& coefficients) : coefficients(coefficients) {
    removeLeadingZeros();
}

template <typename T>
Polynomial<T>::Polynomial(const Polynomial<T>& other) : coefficients(other.coefficients) {
    removeLeadingZeros();
}

template <typename T>
int Polynomial<T>::getDegree() const {
    return (signed) coefficients.size() - 1;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>& other) {
    coefficients = other.coefficients;
    return *this;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator+(const Polynomial<T>& other) const {
    std::vector<T> resultCoefficients(std::max(coefficients.size(), other.coefficients.size()));
    for (size_t i = 0; i < resultCoefficients.size(); ++i) {
        T a = (i < coefficients.size()) ? coefficients[i] : T(0);
        T b = (i < other.coefficients.size()) ? other.coefficients[i] : T(0);
        resultCoefficients[i] = a + b;
    }
    Polynomial<T> out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator-(const Polynomial<T>& other) const {
    std::vector<T> resultCoefficients(std::max(coefficients.size(), other.coefficients.size()));
    for (auto i = 0; i < resultCoefficients.size(); ++i) {
        T a = (i < coefficients.size()) ? coefficients[i] : T(0);
        T b = (i < other.coefficients.size()) ? other.coefficients[i] : T(0);
        resultCoefficients[i] = a - b;
    }
    Polynomial<T> out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator*(const Polynomial<T>& other) const {
    std::vector<T> resultCoefficients(coefficients.size() + other.coefficients.size() - 1);
    for (size_t i = 0; i < coefficients.size(); ++i) {
        for (size_t j = 0; j < other.coefficients.size(); ++j) {
            resultCoefficients[i + j] += coefficients[i] * other.coefficients[j];
        }
    }
    Polynomial<T> out(resultCoefficients);
    out.removeLeadingZeros();
    return out;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator/(const Polynomial<T>& other) const {
    Polynomial<T> quotient, remainder;
    divide(other, quotient, remainder);
    return quotient;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator%(const Polynomial<T>& other) const {
    Polynomial<T> quotient, remainder;
    divide(other, quotient, remainder);
    return remainder;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator-() const {
    std::vector<T> resultCoefficients(coefficients.size());
    for (size_t i = 0; i < coefficients.size(); ++i) {
        resultCoefficients[i] = -coefficients[i];
    }
    return Polynomial<T>(resultCoefficients);
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& other) {
    *this = *this + other;
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& other) {
    *this = *this - other;
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& other) {
    *this = *this * other;
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator/=(const Polynomial<T>& other) {
    *this = *this / other;
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator%=(const Polynomial<T>& other) {
    *this = *this % other;
    return *this;
}

template <typename T>
T Polynomial<T>::operator()(const T& x) const {
    T result = T(0);
    T xPower = T(1);
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * xPower;
        xPower *= x;
    }
    return result;
}

template <typename T>
T Polynomial<T>::operator[](uint i) const {
    assert (i < coefficients.size());
    return coefficients[i];
}

template <typename T>
bool Polynomial<T>::operator==(const Polynomial<T>& other) const
{
    for (size_t i = 0; i < coefficients.size(); ++i) {
        if (coefficients[i] != other.coefficients[i]) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool Polynomial<T>::operator!=(const Polynomial<T>& other) const
{
    return !(*this == other);
}

template <typename T>
void Polynomial<T>::removeLeadingZeros() {
    while (!coefficients.empty() && coefficients.back() == T(0)) {
        coefficients.pop_back();
    }
}

template <typename T>
void Polynomial<T>::divide(const Polynomial<T>& other, Polynomial<T>& quotient, Polynomial<T>& remainder) const {
    if (other.coefficients.empty()) {
        throw std::invalid_argument("Division by zero");
    }
    if (coefficients.size() < other.coefficients.size()) {
        quotient = Polynomial<T>();
        remainder = *this;
        remainder.removeLeadingZeros();
        return;
    }

    Polynomial<T> temp = *this;

    T factor = temp.coefficients.back() / other.coefficients.back();
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

template <typename T>
ostream& operator<< (ostream& os, const Polynomial<T>& p) {
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



#include <mcl/bls12_381.hpp>
template class Polynomial<mcl::bn::Fr>;
template ostream& operator<< (ostream& os, const Polynomial<mcl::bn::Fr>& p);