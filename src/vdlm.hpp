#ifndef VDLP_HPP
#define VDLP_HPP

#include <iostream>
#include <vector>
#include "proof.hpp"
#include "prng.hpp"

using namespace std;

vector<bool> probToBits(double p, uint n);

// if a[i] is true, res[i] = b[i], otherwise res[i] = c[i]
// in arithmetic: res[i] = a[i] * b[i] + (1 - a[i]) * c[i] = a[i] * (b[i] - c[i]) + c[i]
// vector<bool> mux(const vector<bool>& a, const vector<bool>& b, const vector<bool>& c);

struct BernoulliStepCache{
    Fr key;
    Fr r_key;
    Polynomial F_prev;
    Polynomial R_prev;
    Polynomial F_cur;
    Polynomial R_cur;
    Polynomial F_prng_rt;
    Polynomial R_prng_rt;
    Polynomial F_prng_res;
    Polynomial R_prng_res;
    G1 com_key;
    G1 com_prev;
    G1 com_cur;
    G1 com_prng_rt;
    G1 com_prng_res;
    bool b;
};

vector<bool> BernoulliStepSample(vector<bool> prev, bool b,
    const Fr& key,
    const LegendrePRNGPubParam& pp, 
    const BernoulliStepCache& cache_prev, BernoulliStepCache& cache_cur);

bool BernoulliStepIP(const BernoulliStepCache& cache, const LegendrePRNGPubParam& pp, Timer& ptimer, Timer& vtimer);

vector<bool> Bernoulli(vector<bool> p, const LegendrePRNGPubParam& pp, Timer& ptimer, Timer& vtimer, G1& com_out);



// vector<bool> sampleBernoulli(const vector<vector<bool>>& r_vec, uint len, double p, uint prec);

// vector<uint> sampleGeometric(const vector<vector<vector<bool>>>& r_vec, uint len, double p, uint log_range, uint prec);

// vector<int> sampleLaplacian(const vector<vector<bool>>& r1, const vector<vector<vector<bool>>>& r2, const vector<bool>& r3,
//     uint len, double t, uint log_range, uint prec);

// struct BernoulliPubParam
// {
//     uint len;
//     vector<bool> p_vec;
//     LegendrePRNGPubParam pp;
//     G1 com_F_all1;
// }



#endif // VDLP_HPP