#ifndef PRNG_HPP
#define PRNG_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include "proof.hpp"
#include "utils.hpp"

using namespace std;
using namespace mcl::bn;

#define LEGENDRE_PRNG_NON_SQ 5

struct LegendrePRNGPubParam
{
    uint len;
    PubParam pp;
    Fr omega_gen;
    Polynomial F_range;
    Polynomial F_omega;
    G1 com_F_range;
    G1 com_F_omega; 
};

LegendrePRNGPubParam LegendrePRNGTrustedSetup(uint len);

vector<Fr> arange(uint len);

vector<bool> LegendrePRNG(const Fr& key, uint len, vector<Fr>& rt_vec);

void convertLegendrePRNG(const vector<Fr>& rt_vec, const vector<bool>& res, Polynomial& F_rt, Polynomial& F_res, const LegendrePRNGPubParam& pp);

void commitLegendrePRNG(const Polynomial& F_rt, const Polynomial& F_res, const Polynomial& R_rt, const Polynomial& R_res, G1& com_rt, G1& com_res, const LegendrePRNGPubParam& pp);

bool proveLegendrePRNG(const Fr& key, const vector<Fr>& rt_vec, const vector<bool>& res, 
    const Fr& r_key, const Polynomial& R_rt, const Polynomial& R_res,
    const G1& com_key, const G1& com_rt, const G1& com_res, 
    const LegendrePRNGPubParam& pp, Timer& ptimer, Timer& vtimer);

#endif // PRNG_HPP