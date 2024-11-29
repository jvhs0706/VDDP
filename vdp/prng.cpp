#include "prng.hpp"
#include <execution>

void determineSeed(Fr& key, Fr& r_key, G1& com_key, const LegendrePRNGPubParam& pp, Timer& ptimer, Timer& vtimer, uint& communication)
{
    Fr pub_coin;
    
    ptimer.start();
    key.setByCSPRNG();
    r_key.setByCSPRNG();
    com_key = pp.pp.gVec[0] * key + pp.pp.hVec[0] * r_key;
    
    vtimer.start();
    pub_coin.setByCSPRNG();
    key += pub_coin;
    com_key += pp.pp.gVec[0] * pub_coin;
    vtimer.stop();
    ptimer.stop();

    communication += (sizeof(Fr) + sizeof(G1));
}

LegendrePRNGPubParam LegendrePRNGTrustedSetup(uint len)
{
    // len > 0
    assert (len);
    // len is a power of 2
    assert((len & (len - 1)) == 0);
    
    auto pp = trustedSetup(len + 1);
    Fr omega_gen;
    Polynomial F_range = ntt_vec_to_poly(arange(len), omega_gen);
    
    VanishingPolynomial F_omega(len);

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
        // std::execution::par,
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

void commitLegendrePRNG(const Polynomial& F_rt, const Polynomial& F_res, 
    const Polynomial& R_rt, const Polynomial& R_res, 
    G1& com_rt, G1& com_res, 
    const LegendrePRNGPubParam& pp)
{
    com_rt = commitPoly(F_rt, R_rt, pp.pp.gVec, pp.pp.hVec);
    com_res = commitPoly(F_res, R_res, pp.pp.gVec, pp.pp.hVec);
}

// (5-4 res ) * input == sq(rt)
bool proveLegendrePRNG(const Fr& key, const Polynomial& F_rt, const Polynomial& F_res,
    const Fr& r_key, const Polynomial& R_rt, const Polynomial& R_res,
    const G1& com_key, const G1& com_rt, const G1& com_res, 
    const LegendrePRNGPubParam& pp, Timer& ptimer, Timer& vtimer, uint& communication)
{
    bool accepted = Binary(F_res, R_res, pp.len, pp.omega_gen, com_res, pp.pp, ptimer, vtimer);


    const auto& F_range = pp.F_range;
    const auto& F_omega = pp.F_omega;


    const G1& g = pp.pp.gVec[0];
    const G1& h = pp.pp.hVec[0];


    ptimer.start();
    Polynomial zero_poly = (F_res * Fr(1 - LEGENDRE_PRNG_NON_SQ) + Fr(LEGENDRE_PRNG_NON_SQ)) * (F_range + key) - F_rt * F_rt;
    Polynomial F_quot, F_rem;
    zero_poly.divide(F_omega, F_quot, F_rem);
    assert (F_rem.getDegree() == -1); // It's the zero polynomial
    auto R_quot = randomPolynomial(pp.len);
    auto com_F_quot = commitPoly(F_quot, R_quot, pp.pp.gVec, pp.pp.hVec);
    ptimer.stop();

    communication += sizeof(G1);

    vtimer.start();
    Fr u;
    u.setByCSPRNG();
    const Fr z_omega = F_omega(u);
    vtimer.stop();

    communication += sizeof(Fr);

    ptimer.start();
    Fr z_range;
    auto z_range_proof = provePolyEval(F_range, u, z_range, pp.pp.gVec);
    ptimer.stop();

    communication += (sizeof(Fr) + sizeof(G1));

    vtimer.start();
    accepted &= verifyPolyEval(z_range, z_range_proof, pp.com_F_range, u, g, pp.pp.g2, pp.pp.g2_tau);
    vtimer.stop();


    ptimer.start();
    const auto z_res = F_res(u), z_rt = F_rt(u);
    Fr r_res, r_rt;
    r_res.setByCSPRNG();
    r_rt.setByCSPRNG();    

    const auto com_res_u = g * z_res + h * r_res;
    const auto com_rt_u = g * z_rt + h * r_rt;

    const auto z_res_ = z_res * Fr(1 - LEGENDRE_PRNG_NON_SQ) + Fr(LEGENDRE_PRNG_NON_SQ);
    const auto key_ = key + z_range;

    const auto r_res_ = r_res * Fr(1 - LEGENDRE_PRNG_NON_SQ);
    const auto& r_key_ = r_key;

    vtimer.start();
    const auto com_res_u_ = com_res_u * Fr(1 - LEGENDRE_PRNG_NON_SQ) + g * Fr(LEGENDRE_PRNG_NON_SQ);
    const auto com_key_ = com_key + g * z_range;
    vtimer.stop();

    const auto z_lhs = z_res_ * key_;
    Fr r_lhs;
    r_lhs.setByCSPRNG();
    const auto com_lhs = g * z_lhs + h * r_lhs;
    ptimer.stop();

    communication += sizeof(G1);

    accepted &= Prod(z_lhs, r_lhs, z_res_, r_res_, key_, r_key_, com_lhs, com_res_u_, com_key_, g, h, ptimer, vtimer, communication);

    ptimer.start();
    auto z_rhs = z_rt * z_rt;
    Fr r_rhs;
    r_rhs.setByCSPRNG();
    const auto com_rhs = g * z_rhs + h * r_rhs;
    ptimer.stop();

    communication += sizeof(G1);

    accepted &= Prod(z_rhs, r_rhs, z_rt, r_rt, z_rt, r_rt, com_rhs, com_rt_u, com_rt_u, g, h, ptimer, vtimer, communication);

    ptimer.start();
    const Fr z_quot = F_quot(u);
    Fr r_quot;
    r_quot.setByCSPRNG();
    const auto com_quot = g * z_quot + h * r_quot;
    ptimer.stop();

    communication += sizeof(G1);

    accepted &= EvalSecret(z_quot, r_quot, u, F_quot, R_quot, com_quot, com_F_quot, pp.pp, ptimer, vtimer, communication);
    
    ptimer.start();
    auto z_quot_ = z_quot * z_omega;
    auto r_quot_ = r_quot * z_omega;
    ptimer.stop();

    accepted &= Equal(z_quot_, r_quot_, z_lhs - z_rhs, r_lhs - r_rhs, 
        com_quot * z_omega,
        com_lhs - com_rhs, g, h, ptimer, vtimer, communication);
    
    return accepted;
    
}

vector<bool> verifiableUniformBits(Polynomial& F, Polynomial& R, G1& com, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& communication)
{
    Fr key, r_key;
    G1 com_key;
    determineSeed(key, r_key, com_key, pp, ptimer, vtimer, communication);

    comp_timer.start();
    vector<Fr> rt_vec;
    auto res = LegendrePRNG(key, pp.len, rt_vec);
    Polynomial F_rt;
    auto R_rt = randomPolynomial(F_rt.getDegree());
    G1 com_rt;

    convertLegendrePRNG(rt_vec, res, F_rt, F, pp);
    R = randomPolynomial(F.getDegree());
    commitLegendrePRNG(F_rt, F, R_rt, R, com_rt, com, pp);
    comp_timer.stop();

    communication += sizeof(G1);

    assert(proveLegendrePRNG(key, F_rt, F, r_key, R_rt, R, com_key, com_rt, com, pp, ptimer, vtimer, communication));

    

    return res; 

}