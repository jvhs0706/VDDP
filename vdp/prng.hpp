#ifndef PRNG_HPP
#define PRNG_HPP

#include <mcl/bls12_381.hpp>
#include <vector>

using namespace std;
using namespace mcl::bn;

int LegendreSymbol(const Fr& x);

bool LegendrePRF(const Fr& x, const Fr& key);

vector<bool> LegendrePRNG(const Fr& key, uint len);

#endif // PRNG_HPP