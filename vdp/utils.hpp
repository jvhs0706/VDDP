#ifndef UTILS_HPP
#define UTILS_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include "polynomial.hpp"

using namespace mcl::bn;
using namespace std;

Polynomial randomPolynomial(uint deg);

Fr getRootOfUnity(uint logN);

vector<Fr> ntt_helper(const vector<Fr>& n, const Fr& omega);

Polynomial ntt_vec_to_poly(const vector<Fr>& a, Fr& omega);

Polynomial ntt_vec_to_poly_given_omega(const vector<Fr>& a, const Fr& omega);

uint ceilLog2(uint n);

#endif // UTILS_HPP