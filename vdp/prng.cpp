#include "prng.hpp"
#include <execution>

vector<bool> LegendrePRNG(const Fr& key, uint len, vector<Fr>& rt_vec)
{
    // if rt_vec.len() != len, then we need to resize rt_vec
    if (rt_vec.size() != len) rt_vec.resize(len);
    auto rt_begin = &rt_vec.front();

    // generate a random bit string of length len using the LegendrePRF function
    vector<bool> res(len);
    for_each(
        std::execution::par,
        rt_vec.begin(),
        rt_vec.end(),
        [&](Fr& rt)
        {   
            auto idx = &rt - rt_begin;
            Fr temp = idx + key;
            res[idx] = Fr::squareRoot(rt, temp);
            if (!res[idx]) {
                if (!Fr::squareRoot(rt, temp * LEGENDRE_PRNG_NON_SQ)) throw std::runtime_error("ERROR: PRNG is not working correctly!");
            }
        }
    );
    return res;
}