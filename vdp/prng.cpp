#include "prng.hpp"
#include <execution>

LegendrePRNGPubParam LegendrePRNGTrustedSetup(uint len)
{
    // len > 0
    assert (len);
    // len is a power of 2
    assert((len & (len - 1)) == 0);
    
    auto pp = trustedSetup(len + 1);
    Fr omega_gen;
    Polynomial F_range = ntt_vec_to_poly(arange(len), omega_gen);
    
    vector<Fr> F_omega_params(len + 1, Fr(0));
    F_omega_params[0] = - 1;
    F_omega_params[len] = 1;
    Polynomial F_omega(F_omega_params);

    G1 com_F_range = commitPoly(F_range, pp.gVec);
    G1 com_F_omega = commitPoly(F_omega, pp.gVec);

    return {len, pp, omega_gen, F_range, F_omega, com_F_range, com_F_omega};
}

vector<Fr> arange(uint len)
{
    vector<Fr> res(len);
    Fr* res_begin = &res.front();
    for_each(
        std::execution::par,
        res.begin(),
        res.end(),
        [&](Fr& rt)
        {
            auto idx = &rt - res_begin;
            rt = idx;
        }
    );
    return res;
}

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

void convertLegendrePRNG(const vector<Fr>& rt_vec, const vector<bool>& res, Polynomial& F_rt, Polynomial& F_res, const LegendrePRNGPubParam& pp)
{
    F_rt = ntt_vec_to_poly_given_omega(rt_vec, pp.omega_gen);
    F_res = ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen);
}

void commitLegendrePRNG(const Polynomial& F_rt, const Polynomial& F_res, const Polynomial& R_rt, const Polynomial& R_res, G1& com_rt, G1& com_res, const LegendrePRNGPubParam& pp)
{
    com_rt = commitPoly(F_rt, R_rt, pp.pp.gVec, pp.pp.hVec);
    com_res = commitPoly(F_res, R_res, pp.pp.gVec, pp.pp.hVec);
}

// (5-4 res ) * input == sqare(rt)
// bool proveLegendrePRNG(const Fr& key, const vector<Fr>& rt_vec, const vector<bool>& res, 
//     const Fr& r_key, const Polynomial& R_rt, const Polynomial& R_res,
//     const G1& com_key, const G1& com_rt, const G1& com_res, 
//     const LegendrePRNGPubParam& pp, Timer& ptimer, Timer& vtimer)
// {   
//     return true;
// }