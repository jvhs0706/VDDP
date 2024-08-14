#include "vdlm.hpp"
#include <execution>
#include <cmath>

vector<bool> probToBits(double p, uint n)
{
    if (p < 0 || p >= 1) throw invalid_argument("p must be in [0, 1]");
    // p = 0.res[0]res[1]...res[n-1]
    vector<bool> res;
    for (uint i = 0; i < n - 1; ++ i)
    {
        p *= 2;
        if (p >= 1)
        {
            res.push_back(true);
            p -= 1;
        }
        else res.push_back(false);
    }
    res.push_back(p >= 0.5);

    return res;
}

vector<bool> sampleBernoulli(const vector<vector<bool>>& r_vec, uint len, double p, uint prec)
{
    auto p_bits = probToBits(p, prec);
    // for (uint i = 0; i < prec; ++ i)
    // {
    //     cout << p_bits[i] << " ";
    // }
    // cout << endl;
    vector<bool> res(len, true);

    // assert (r_vec.size() == prec);
    if (r_vec.size() != prec) throw invalid_argument("r_vec.size() != prec");

    vector<uint> idx(len);
    auto* idx_begin = &idx.front();

    for (int i = prec - 1; i >= 0; -- i)
    {
        auto& r = r_vec[i];
        // assert (r.size() == len);
        if (r.size() != len) throw invalid_argument("r.size() != len");

        for_each(
            std::execution::par,
            idx.begin(),
            idx.end(),
            [&](uint& j)
            {
                j = &j - idx_begin;
                auto t = p_bits[i] != r[j];
                res[j] = t ? !r[j] : res[j];
            }
        );
    }
    return res;
}

vector<uint> sampleGeometric(const vector<vector<vector<bool>>>& r_vec, uint len, double p, uint log_range, uint prec)
{   
    vector<uint> res(len, 0);
    auto* res_begin = &res.front();

    if (r_vec.size() != log_range) throw invalid_argument("r_vec.size() != log_range");

    for (uint i = 0; i < log_range; ++ i)
    {
        double p_pow_inv = pow(p, -static_cast<double>(1 << i));
        double q = 1.0 / (1.0 + p_pow_inv);
        // cout << "p_pow_inv = " << p_pow_inv << endl;
        // cout << "q = " << q << endl;
        auto r = sampleBernoulli(r_vec[i], len, q, prec);

        // for (uint j = 0; j < len; ++ j)
        // {
        //     cout << r[j] << " ";
        // }
        // cout << endl;

        for_each(
            std::execution::par,
            res.begin(),
            res.end(),
            [&](uint& j)
            {
                auto idx = &j - res_begin;
                if (r[idx]) j |= (1 << i);
            }
        );
    }
    return res;
}

vector<int> sampleLaplacian(const vector<vector<bool>>& r1, const vector<vector<vector<bool>>>& r2, const vector<bool>& r3,
    uint len, double t, uint log_range, uint prec)
{   
    double p0 = 0.5; // change this later!!!!!!!!!!!!!!
    double p = exp(-1.0/t);
    auto b = sampleBernoulli(r1, len, p0, prec);
    auto x = sampleGeometric(r2, len, 0.5, log_range, prec);
    auto& s= r3;

    vector<int> res(len, 0);

    auto* res_begin = &res.front();

    for_each(
        std::execution::par,
        res.begin(),
        res.end(),
        [&](int& j)
        {
            auto idx = &j - res_begin;
            int t = b[idx] ? 0 : (x[idx] + 1); 
            j = s[idx] ? t : -t;
        }
    );
    

    return res;    
}

// NOT TESTED YET
vector<bool> BernoulliStepSample(vector<bool> prev, bool b,
    const Fr& key,
    const LegendrePRNGPubParam& pp, 
    const BernoulliStepCache& cache_prev, BernoulliStepCache& cache_cur)
{
    vector<bool> res(pp.len);
    vector<uint> idx(pp.len);
    assert (prev.size() == pp.len);
    vector<Fr> rt_vec(pp.len);
    auto r = LegendrePRNG(key, pp.len, rt_vec);

    auto* idx_begin = &idx.front();
    for_each(
        std::execution::par,
        idx.begin(),
        idx.end(),
        [&](uint& j)
        {
            j = &j - idx_begin;
            auto t = b != r[j];
            res[j] = t ? !r[j] : res[j];
        }
    );

    auto& g = pp.pp.gVec[0];
    auto& h = pp.pp.hVec[0];

    cache_cur.key = key;
    cache_cur.r_key.setByCSPRNG();
    cache_cur.F_prev = cache_prev.F_cur;
    cache_cur.R_prev = cache_prev.R_cur;
    convertLegendrePRNG(rt_vec, r, cache_cur.F_prng_rt, cache_cur.F_prng_res, pp);
    cache_cur.R_prng_rt = randomPolynomial(cache_cur.F_prng_rt.getDegree());
    cache_cur.R_prng_res = randomPolynomial(cache_cur.F_prng_res.getDegree());
    
    cache_cur.F_cur = ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen);
    cache_cur.R_cur = randomPolynomial(cache_cur.F_cur.getDegree());

    commitLegendrePRNG(cache_cur.key, cache_cur.F_prng_rt, cache_cur.F_prng_res,
        cache_cur.r_key, cache_cur.R_prng_rt, cache_cur.R_prng_res,
        cache_cur.com_key, cache_cur.com_prng_rt, cache_cur.com_prng_res, pp);

    cache_cur.com_prev = cache_prev.com_cur;
    cache_cur.com_cur = commitPoly(cache_cur.F_cur, cache_cur.R_cur, pp.pp.gVec, pp.pp.hVec);

    return res;
}