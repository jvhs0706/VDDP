#include "prng.hpp"
#include <execution>

int LegendreSymbol(const Fr& x)
{
    if (x.isZero()) return 0;
    else 
    {
        Fr temp;
        bool b = Fr::squareRoot(temp, x);
        if (b) return 1;
        else return -1;
    }
}

bool LegendrePRF(const Fr& x, const Fr& key)
{
    return LegendreSymbol(x + key) >= 0;
}

vector<bool> LegendrePRNG(const Fr& key, uint len)
{
    // generate a random bit string of length len using the LegendrePRF function
    vector<uint> idx(len);
    vector<bool> res(len);
    auto* idx_begin = &idx.front();
    for_each(
        std::execution::par,
        idx.begin(),
        idx.end(),
        [&key, &idx_begin, &res](uint& item)
        {
            item = &item - idx_begin;
            res[item] = LegendrePRF(item, key);
        }
    );
    return res;
}