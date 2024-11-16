#include "vddlm.hpp"
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

    // remove the falses at the end
    while (res.size() > 0 && !res.back()) res.pop_back();

    return res;
}

vector<bool> Bernoulli(vector<bool> p, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    // the last bit of p must be 1
    if (p.size() == 0 || !p.back()) {
        throw invalid_argument("the last bit of p must be 1");
        return vector<bool>();
    }
    else if (p.size() == 1) {
        Fr key, r_key;
        G1 com_key;
        vector<Fr> rt_vec(pp.len);
        return verifiableUniformBits(F_out, R_out, com_out, pp, comp_timer, ptimer, vtimer);
    }
    else {
        Polynomial F_out_, R_out_;
        G1 com_out_;
        auto p_ = vector<bool>(p.begin() + 1, p.end());
        auto res_ = Bernoulli(p_, pp, comp_timer, ptimer, vtimer, F_out_, R_out_, com_out_);

        Polynomial F_cur, R_cur;
        G1 com_cur;
        auto cur_bit = verifiableUniformBits(F_cur, R_cur, com_cur, pp, comp_timer, ptimer, vtimer);
        auto res = vector<bool>(pp.len);

        if (p[0])
        {
            for (uint i = 0; i < pp.len; ++ i)
            {
                res[i] = cur_bit[i] || res_[i];
            }
            
            F_out = ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen);
            R_out = randomPolynomial(F_out.getDegree());
            com_out = commitPoly(F_out, R_out, pp.pp.gVec, pp.pp.hVec);
            auto Fc = F_out_ + F_cur - F_out;
            auto Rc = R_out_ + R_cur - R_out;
            assert (Hadamard(Fc, Rc, 
                F_out_, R_out_,
                F_cur, R_cur,
                pp.len, pp.omega_gen,
                com_out_ + com_cur - com_out, com_out_, com_cur,
                pp.pp, ptimer, vtimer));
        }
        else
        {
            for (uint i = 0; i < pp.len; ++ i)
            {
                res[i] = cur_bit[i] && res_[i];
            }
            
            F_out = ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen);
            R_out = randomPolynomial(F_out.getDegree());
            com_out = commitPoly(F_out, R_out, pp.pp.gVec, pp.pp.hVec);
            assert (Hadamard(F_out, R_out, 
                F_out_, R_out_,
                F_cur, R_cur,
                pp.len, pp.omega_gen,
                com_out, com_out_, com_cur,
                pp.pp, ptimer, vtimer));
        }

        return res;
    }

}

vector<uint> Geometric(double p, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    vector<uint> res(pp.len, 0);
    auto res_begin = &res.front();
    for (uint i = 0; i < log_range; ++ i)
    {
        double fail_prob_pow_inv = pow(1.0L - p, -static_cast<double>(1 << i));
        double q = 1.0L / (1.0L + fail_prob_pow_inv);
        auto q_bin = probToBits(q, prec);
        if (!q_bin.size()) break;
        Polynomial step_F_out, step_R_out;
        G1 step_com_out;
        auto r = Bernoulli(q_bin, pp, comp_timer, ptimer, vtimer, step_F_out, step_R_out, step_com_out);

        comp_timer.start();
        for_each(
            // std::execution::par,
            res.begin(),
            res.end(),
            [&](uint& j)
            {
                auto idx = &j - res_begin;
                if (r[idx]) j |= (1 << i);
            }
        );
        comp_timer.stop();

        if (i == 0)
        {
            F_out = step_F_out;
            R_out = step_R_out;
            com_out = step_com_out;
        }
        else
        {
            F_out += step_F_out * Fr(1 << i);
            R_out += step_R_out * Fr(1 << i);
            com_out += step_com_out * Fr(1 << i);
        }
    }

    return res;
}