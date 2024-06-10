#include "polynomial.hpp"
#include <stdexcept>

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
        return;
    }
    std::vector<T> resultCoefficients(coefficients.size() - other.coefficients.size() + 1);
    Polynomial<T> current = *this;
    for (size_t i = resultCoefficients.size(); i-- > 0;) {
        T factor = current.coefficients.back() / other.coefficients.back();
        resultCoefficients[i] = factor;
        std::vector<T> subtractCoefficients(i + other.coefficients.size());
        subtractCoefficients[i] = factor;
        Polynomial<T> subtract(subtractCoefficients);
        current -= subtract * other;
        current.removeLeadingZeros();
    }
    quotient = Polynomial<T>(resultCoefficients);
    quotient.removeLeadingZeros();
    remainder = current;
}

#include <mcl/bls12_381.hpp>
template class Polynomial<mcl::bn::Fr>;