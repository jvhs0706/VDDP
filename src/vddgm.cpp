#include "vddgm.hpp"
#include <execution>
#include <cmath>
#include <fstream>

vector<vector<bool>> BitDecomp(const vector <uint>& u, const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm,
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    vector<Polynomial>& F_out, vector<Polynomial>& R_out, vector<G1>& com_out)
{
    vector<vector<bool>> res(32, vector<bool>(pp.len, false));
    Polynomial R_31 = R_in / Fr(1u << 31);

    // resize F_out and R_out to hold 32 polynomials
    F_out.resize(32);
    R_out.resize(32);
    com_out.resize(32);
    
    for (uint i = 0; i < 32; ++i)
    {   
        comp_timer.start();
        for (uint j = 0; j < pp.len; ++j)
        {
            res[i][j] = (u[j] >> i) & 1;
        }
        comp_timer.stop();
        
        ptimer.start();
        F_out[i] = ntt_vec_to_poly_given_omega(vector<Fr>(res[i].begin(), res[i].end()), pp.omega_gen);
        R_out[i] = (i < 31) ? randomPolynomial(F_out[i].getDegree()) : R_31;
        if (i < 31){
            R_31 -= R_out[i] / Fr(1u << (31-i));
        }
        com_out[i] = commitPoly(F_out[i], R_out[i], pp.pp.gVec, pp.pp.hVec);
        ptimer.stop();

        assert(Binary(
            F_out[i], R_out[i],
            pp.len, pp.omega_gen,
            com_out[i], pp.pp,
            ptimer, vtimer, comm
        )); // Binary check for each bit decomposition
    }
    return res;
}

vector<vector<bool>> BitDecomp(const vector <int>& u, const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm,
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    vector<Polynomial>& F_out, vector<Polynomial>& R_out, vector<G1>& com_out)
{
    vector<vector<bool>> res(32, vector<bool>(pp.len, false));
    Polynomial R_31 = - R_in / Fr(1u << 31);

    // resize F_out and R_out to hold 32 polynomials
    F_out.resize(32);
    R_out.resize(32);
    com_out.resize(32);
    
    for (uint i = 0; i < 32; ++i)
    {   
        comp_timer.start();
        for (uint j = 0; j < pp.len; ++j)
        {
            res[i][j] = (u[j] >> i) & 1;
        }
        comp_timer.stop();
        
        ptimer.start();
        F_out[i] = ntt_vec_to_poly_given_omega(vector<Fr>(res[i].begin(), res[i].end()), pp.omega_gen);
        R_out[i] = (i < 31) ? randomPolynomial(F_out[i].getDegree()) : R_31;
        if (i < 31){
            R_31 += R_out[i] / Fr(1u << (31-i));
        }
        com_out[i] = commitPoly(F_out[i], R_out[i], pp.pp.gVec, pp.pp.hVec);
        ptimer.stop();

        comm += sizeof(G1);

        assert(Binary(
            F_out[i], R_out[i],
            pp.len, pp.omega_gen,
            com_out[i], pp.pp,
            ptimer, vtimer, comm
        )); // Binary check for each bit decomposition
    }
    return res;
}

vector<uint> Abs(const vector<int>& s, const LegendrePRNGPubParam& pp,
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm,
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    vector<Polynomial> F_bd, R_bd;
    vector<G1> com_bd;
    auto s_bd = BitDecomp(s, pp, comp_timer, ptimer, vtimer, comm, F_in, R_in, com_in, F_bd, R_bd, com_bd);
    comp_timer.start();
    vector<uint> res(pp.len, 0);
    for (uint i = 0; i < pp.len; ++i)
    {
        res[i] = (s[i] < 0) ? -s[i] : s[i];
    }
    F_out = ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen);
    R_out = randomPolynomial(F_out.getDegree());
    com_out = commitPoly(F_out, R_out, pp.pp.gVec, pp.pp.hVec);
    comp_timer.stop();
    comm += sizeof(G1);

    ptimer.start();
    vtimer.start();
    auto F_c = (F_in - F_out) / Fr(2);
    auto R_c = (R_in - R_out) / Fr(2);
    auto com_c = (com_in - com_out) * (Fr(1) / Fr(2));
    vtimer.stop();
    ptimer.stop();

    const auto& F_a = F_bd[31];
    const auto& R_a = R_bd[31];
    const auto& com_a = com_bd[31];

    // Check Hadamard product
    assert(Hadamard(F_c, R_c, F_a, R_a, F_in, R_in, 
            pp.len, pp.omega_gen,
            com_c, com_a, com_in, 
            pp.pp, ptimer, vtimer, comm));

    return res;
}

vector<bool> Bexp_step(vector<bool> u, double r, uint prob_nbit, uint idx,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    vector<bool> res(pp.len, false);
    
    double p = exp(-pow(2.0L, idx) / r);
    if (p >= 0.5L / (1u << prob_nbit)) {
        auto p_bits = probToBits(p, prob_nbit);
        Polynomial F_b_i, R_b_i;
        G1 com_b_i;

        // Generate Bernoulli bits
        auto b_i = Bernoulli(p_bits, pp, comp_timer, ptimer, vtimer, comm, F_b_i, R_b_i, com_b_i);
        
        comp_timer.start();
        for (uint i = 0; i < pp.len; ++ i)
        {
            res[i] = (!u[i]) || b_i[i];
        }
        comp_timer.stop();

        ptimer.start();
        F_out = ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen);
        R_out = randomPolynomial(F_out.getDegree());
        com_out = commitPoly(F_out, R_out, pp.pp.gVec, pp.pp.hVec);
        vtimer.start();
        auto F_c = F_out + F_in - Fr(1);
        auto R_c = R_out + R_in;
        auto com_c = com_out + com_in - pp.pp.gVec[0];
        vtimer.stop();
        ptimer.stop();
        comm += sizeof(G1);

        // Check Hadamard product
        assert(Hadamard(F_c, R_c, F_in, R_in, F_b_i, R_b_i, 
                pp.len, pp.omega_gen,
                com_c, com_in, com_b_i, 
                pp.pp, ptimer, vtimer, comm));
        return res;
    }
    else {
        comp_timer.start();
        for (uint i = 0; i < pp.len; ++ i)
        {
            res[i] = (!u[i]);
        }
        comp_timer.stop();

        ptimer.start();
        vtimer.start();
        F_out = - F_in + Fr(1);
        R_out = - R_in;
        com_out = - com_in + pp.pp.gVec[0];
        ptimer.stop();
        vtimer.stop();
        return res;
    }
}



vector<bool> Bexp(vector<uint>u, double r, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    Polynomial& F_out, Polynomial& R_out, G1& com_out)
{   
    vector<Polynomial> F_b_is;
    vector<Polynomial> R_b_is;
    vector<G1> com_b_is;
    auto u_is = BitDecomp(u, pp, comp_timer, ptimer, vtimer, comm, F_in, R_in, com_in, F_b_is, R_b_is, com_b_is);
    vector<bool> res = Bexp_step(u_is[0], r, prob_nbit, 0, pp, comp_timer, ptimer, vtimer, comm, F_b_is[0], R_b_is[0], com_b_is[0], F_out, R_out, com_out);
    for (uint i = 1; i < 32; ++i)
    {
        Polynomial F_out_step, R_out_step;
        G1 com_out_step;
        auto step_res = Bexp_step(u_is[i], r, prob_nbit, i, pp, comp_timer, ptimer, vtimer, comm, F_b_is[i], R_b_is[i], com_b_is[i], F_out_step, R_out_step, com_out_step);
        comp_timer.start();
        vector<bool> new_res(pp.len, false);
        Polynomial F_new_res, R_new_res;
        G1 com_new_res;
        for (uint j = 0; j < pp.len; ++j)
        {
            new_res[j] = res[j] && step_res[j];
        }
        comp_timer.stop();

        ptimer.start();
        F_new_res = ntt_vec_to_poly_given_omega(vector<Fr>(new_res.begin(), new_res.end()), pp.omega_gen);
        R_new_res = randomPolynomial(F_new_res.getDegree());
        com_new_res = commitPoly(F_new_res, R_new_res, pp.pp.gVec, pp.pp.hVec);
        ptimer.stop();
        comm += sizeof(G1);

        assert(Hadamard(F_new_res, R_new_res,
            F_out, R_out,
            F_out_step, R_out_step,
            pp.len, pp.omega_gen,
            com_new_res, com_out, com_out_step,
            pp.pp, ptimer, vtimer, comm));

        res = new_res;
        F_out = F_new_res;
        R_out = R_new_res;
        com_out = com_new_res;
    }
    return res;
}

vector<int> DLap4DGauss(double t, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    double p_z = (exp(1/t) - 1) / (exp(1/t) + 1);
    auto z_config = probToBits(p_z, prob_nbit);
    double p_g = 1 - exp(-1/t);
    vector<vector<bool>> g_config;
    for (uint i = 0; i < 32; ++i)
    {
        double p_i = 1 / (1 + exp(pow(2.0L, i) / t));
        // cout << "p_i: " << p_i << endl; 
        if (p_i >= 0.5L / (1u << prob_nbit)) {
            g_config.push_back(probToBits(p_i, prob_nbit));
        } else {
            break;
        }
    }
    return DiscreteLaplacianNew(z_config, g_config, pp,
        comp_timer, ptimer, vtimer, comm, F_out, R_out, com_out);
}

pair<vector<int>,vector<bool>> DGaussIter(double sigma, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    Polynomial& F_out, Polynomial& R_out, G1& com_out,
    Polynomial& F_acc, Polynomial& R_acc, G1& com_acc)
{
    if (sigma <= 0.0) {
        throw invalid_argument("Sigma must be positive.");
        return {};
    }
    else if (sigma < 1.0) {
        double t = sigma * sigma * ceil(1 / sigma);
        double r = 2 * t * ceil(1 / sigma);
        auto s = DLap4DGauss(t, prob_nbit, pp, comp_timer, ptimer, vtimer, comm, F_out, R_out, com_out);
        Polynomial F_abs, R_abs;
        G1 com_abs;
        auto abs_s = Abs(s, pp, comp_timer, ptimer, vtimer, comm, F_out, R_out, com_out, F_abs, R_abs, com_abs);
        Polynomial F_v, R_v, F_u, R_u;
        G1 com_v, com_u;

        comp_timer.start();
        vector<int> v(pp.len, 0);
        vector<uint> u(pp.len, 0);
        for (uint i = 0; i < pp.len; ++i)
        {
            v[i] = -1 + (signed) (ceil(1 / sigma) * abs_s[i]);
            u[i] = v[i] * v[i];
        }

        vtimer.start();
        F_v =  F_abs * Fr(int(ceil(1/sigma))) - Fr(1); 
        R_v = R_abs * Fr(int(ceil(1/sigma)));
        com_v = com_abs * Fr(int(ceil(1/sigma))) - pp.pp.gVec[0];
        vtimer.stop();

        F_u = ntt_vec_to_poly_given_omega(vector<Fr>(u.begin(), u.end()), pp.omega_gen);
        R_u = randomPolynomial(F_u.getDegree());
        com_u = commitPoly(F_u, R_u, pp.pp.gVec, pp.pp.hVec);
        comp_timer.stop();

        comm += sizeof(G1) * 3;

        assert(Hadamard(F_u, R_u, F_v, R_v, F_v, R_v, pp.len, pp.omega_gen,
            com_u, com_v, com_v, 
            pp.pp, ptimer, vtimer, comm));

        auto b = Bexp(u, r, prob_nbit, pp, 
            comp_timer, ptimer, vtimer, comm, F_u, R_u, com_u,
            F_acc, R_acc, com_acc);
        
        return {s,b};
    }
    else { // sigma >= 1
        double t = sigma * sigma / round(sigma);
        double r = 2 * sigma * sigma;
        auto s = DLap4DGauss(t, prob_nbit, pp, comp_timer, ptimer, vtimer, comm, F_out, R_out, com_out);
        Polynomial F_abs, R_abs;
        G1 com_abs;
        auto abs_s = Abs(s, pp, comp_timer, ptimer, vtimer, comm, F_out, R_out, com_out, F_abs, R_abs, com_abs);
        Polynomial F_v, F_u, R_u;
        G1 com_v, com_u;

        comp_timer.start();
        vector<int> v(pp.len, 0);
        vector<uint> u(pp.len, 0);
        for (uint i = 0; i < pp.len; ++i)
        {
            v[i] = -round(sigma) + (signed) abs_s[i];
            u[i] = v[i] * v[i];
        }

        vtimer.start();
        F_v =  F_abs - Fr(int(round(sigma))); 
        auto& R_v = R_abs;
        com_v = com_abs - pp.pp.gVec[0] * Fr(int(round(sigma)));
        vtimer.stop();

        F_u = ntt_vec_to_poly_given_omega(vector<Fr>(u.begin(), u.end()), pp.omega_gen);
        R_u = randomPolynomial(F_u.getDegree());
        com_u = commitPoly(F_u, R_u, pp.pp.gVec, pp.pp.hVec);
        comp_timer.stop();

        comm += sizeof(G1) * 3;

        assert(Hadamard(F_u, R_u, F_v, R_v, F_v, R_v, pp.len, pp.omega_gen,
            com_u, com_v, com_v, 
            pp.pp, ptimer, vtimer, comm));

        auto b = Bexp(u, r, prob_nbit, pp, 
            comp_timer, ptimer, vtimer, comm, F_u, R_u, com_u,
            F_acc, R_acc, com_acc);
        
        return {s,b};
    }
}

vector<int> DGauss(double sigma, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    vector<bool> accepted(pp.len, false);
    vector<int> res(pp.len, 0);
    uint num_accepted = 0;

    F_out = Fr(0);
    R_out = Fr(0);
    com_out = pp.pp.gVec[0] * Fr(0); // Initialize with the commitment to zero polynomial

    
    while (num_accepted < pp.len) {
        Polynomial F_cur_out, R_cur_out, F_cur_accept, R_cur_accept; 
        G1 com_cur_out, com_cur_accept;
        auto iter_res = DGaussIter(sigma, prob_nbit, pp, comp_timer, ptimer, vtimer, comm,
            F_cur_out, R_cur_out, com_cur_out, F_cur_accept, R_cur_accept, com_cur_accept);
        auto& cur_out = iter_res.first;
        auto& cur_accept = iter_res.second;
        
        vtimer.start();
        assert(commitPoly(F_cur_accept, R_cur_accept, pp.pp.gVec, pp.pp.hVec) == com_cur_accept);
        vtimer.stop();
        
        vector<bool> new_accepted(pp.len, false);
        vector<int> s_filtered(pp.len, 0);
        for (uint i = 0; i < pp.len; ++i) {
            new_accepted[i] = cur_accept[i] && !accepted[i];
            if (new_accepted[i]) {
                res[i] = cur_out[i];
                accepted[i] = true;
                num_accepted += 1;
                s_filtered[i] = cur_out[i]; // Keep the value of s for accepted indices
            }
        }

        comp_timer.start();
        auto F_new_accepted = ntt_vec_to_poly_given_omega(vector<Fr>(new_accepted.begin(), new_accepted.end()), pp.omega_gen);
        auto R_new_accepted = randomPolynomial(F_new_accepted.getDegree());
        auto com_new_accepted = commitPoly(F_new_accepted, R_new_accepted, pp.pp.gVec, pp.pp.hVec);

        auto F_s_filtered = ntt_vec_to_poly_given_omega(vector<Fr>(s_filtered.begin(), s_filtered.end()), pp.omega_gen);
        auto R_s_filtered = randomPolynomial(F_s_filtered.getDegree());
        auto com_s_filtered = commitPoly(F_s_filtered, R_s_filtered, pp.pp.gVec, pp.pp.hVec);
        comp_timer.stop();

        assert(Hadamard(F_s_filtered, R_s_filtered,
            F_cur_out, R_cur_out,
            F_new_accepted, R_new_accepted,
            pp.len, pp.omega_gen,
            com_s_filtered, com_cur_out, com_new_accepted,
            pp.pp, ptimer, vtimer, comm));

        comm += sizeof(G1) * 2;

        ptimer.start();
        F_out += F_s_filtered;
        R_out += R_s_filtered;
        vtimer.start();
        com_out += com_s_filtered;
        vtimer.stop();
        ptimer.stop();

        // cout << "Accepted " << num_accepted << " out of " << pp.len << " samples." << endl;
    }

    return res;
}