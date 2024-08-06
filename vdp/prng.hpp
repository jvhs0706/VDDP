#ifndef PRNG_HPP
#define PRNG_HPP

#include <mcl/bls12_381.hpp>
#include <vector>

using namespace std;
using namespace mcl::bn;

#define LEGENDRE_PRNG_NON_SQ 5

vector<bool> LegendrePRNG(const Fr& key, uint len, vector<Fr>& rt_vec);

#endif // PRNG_HPP