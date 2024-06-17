#ifndef UTILS_HPP
#define UTILS_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include "polynomial.hpp"

using namespace mcl::bn;
using namespace std;

Polynomial<Fr> randomPolynomial(uint deg);

Fr getRootOfUnity(uint n);

Polynomial<Fr> ntt(const vector<Fr>& a, Fr& omega, bool inverse = true);

uint ceilLog2(uint n);

#endif // UTILS_HPP