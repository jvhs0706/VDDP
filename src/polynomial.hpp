#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>


template <typename T>
class Polynomial {
    public:
    Polynomial();
    Polynomial(const T& constant);
    Polynomial(const std::vector<T>& coefficients);
    Polynomial(const Polynomial<T>& other);

    // 0 polynomial -> -1
    // otherwise, the degree of the polynomial
    int getDegree() const;

    Polynomial<T>& operator=(const Polynomial<T>& other);

    Polynomial<T> operator+(const Polynomial<T>& other) const;
    Polynomial<T> operator-(const Polynomial<T>& other) const;
    Polynomial<T> operator*(const Polynomial<T>& other) const;
    Polynomial<T> operator/(const Polynomial<T>& other) const;
    Polynomial<T> operator%(const Polynomial<T>& other) const;

    Polynomial<T>& operator+=(const Polynomial<T>& other);
    Polynomial<T>& operator-=(const Polynomial<T>& other);
    Polynomial<T>& operator*=(const Polynomial<T>& other);
    Polynomial<T>& operator/=(const Polynomial<T>& other);
    Polynomial<T>& operator%=(const Polynomial<T>& other);

    T operator()(const T& x) const;

    private:
    std::vector<T> coefficients;
    void removeLeadingZeros();
    void divide(const Polynomial<T>& other, Polynomial<T>& quotient, Polynomial<T>& remainder) const;
};

#endif  // POLYNOMIAL_HPP