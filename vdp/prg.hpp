#ifndef PRG_HPP
#define PRG_HPP

#include <mcl/bls12_381.hpp>
#include <vector>

using namespace std;
using namespace mcl::bn;

int LegendreSymbol(const Fr& x);

bool LegendrePRF(const Fr& x, const Fr& key);

#endif // PRG_HPP