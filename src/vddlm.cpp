#include "vddlm.hpp"
#include <execution>
#include <cmath>
#include <fstream>

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

void readVDDLMConfig(const string& config_file, vector<bool>& z_config, vector<vector<bool>>& g_config)
{
    // empty z_config and g_config
    z_config.clear();
    g_config.clear();

    ifstream in(config_file, ios::in);
    if (!in.is_open()) throw invalid_argument("Cannot open config file");

    string line;
    uint linecount = 0;
    while (getline(in, line))
    {
        if (linecount == 0){ // z 0101111011
            assert (line[0] == 'z');
            // locate the first ' '
            auto space_pos = line.find(' ');
            if (space_pos == string::npos) throw invalid_argument("Invalid config file");
            for (uint i = space_pos + 1; i < line.size(); ++ i)
            {
                if (line[i] == '0') z_config.push_back(false);
                else if (line[i] == '1') z_config.push_back(true);
                else throw invalid_argument("Invalid config file");
            }
        }
        else { // num 0101111011
            auto space_pos = line.find(' ');
            if (space_pos == string::npos) throw invalid_argument("Invalid config file");
            vector<bool> g;
            for (uint i = space_pos + 1; i < line.size(); ++ i)
            {
                if (line[i] == '0') g.push_back(false);
                else if (line[i] == '1') g.push_back(true);
                else throw invalid_argument("Invalid config file");
            }
            g_config.push_back(g);
        }
        linecount ++;
    }

}

void readVDDLMConfig(const string& config_file, vector<bool>& z_config, vector<vector<bool>>& g_config)
{
    // empty z_config and g_config
    z_config.clear();
    g_config.clear();

    ifstream in(config_file, ios::in);
    if (!in.is_open()) throw invalid_argument("Cannot open config file");

    string line;
    uint linecount = 0;
    while (getline(in, line))
    {
        if (linecount == 0){ // z 0101111011
            assert (line[0] == 'z');
            // locate the first ' '
            auto space_pos = line.find(' ');
            if (space_pos == string::npos) throw invalid_argument("Invalid config file");
            for (uint i = space_pos + 1; i < line.size(); ++ i)
            {
                if (line[i] == '0') z_config.push_back(false);
                else if (line[i] == '1') z_config.push_back(true);
                else throw invalid_argument("Invalid config file");
            }
        }
        else { // num 0101111011
            auto space_pos = line.find(' ');
            if (space_pos == string::npos) throw invalid_argument("Invalid config file");
            vector<bool> g;
            for (uint i = space_pos + 1; i < line.size(); ++ i)
            {
                if (line[i] == '0') g.push_back(false);
                else if (line[i] == '1') g.push_back(true);
                else throw invalid_argument("Invalid config file");
            }
            g_config.push_back(g);
        }
        linecount ++;
    }

}

vector<bool> Bernoulli(vector<bool> p, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, Polynomial& F_out, Polynomial& R_out, G1& com_out)
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
        return verifiableUniformBits(F_out, R_out, com_out, pp, comp_timer, ptimer, vtimer, comm);
    }
    else {
        Polynomial F_out_, R_out_;
        G1 com_out_;
        auto p_ = vector<bool>(p.begin() + 1, p.end());
        auto res_ = Bernoulli(p_, pp, comp_timer, ptimer, vtimer, F_out_, R_out_, com_out_, comm);

        Polynomial F_cur, R_cur;
        G1 com_cur;
        auto cur_bit = verifiableUniformBits(F_cur, R_cur, com_cur, pp, comp_timer, ptimer, vtimer, comm);
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
            comm += sizeof(G1);
            auto Fc = F_out_ + F_cur - F_out;
            auto Rc = R_out_ + R_cur - R_out;
            assert (Hadamard(Fc, Rc, 
                F_out_, R_out_,
                F_cur, R_cur,
                pp.len, pp.omega_gen,
                com_out_ + com_cur - com_out, com_out_, com_cur,
                pp.pp, ptimer, vtimer, comm));
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
            comm += sizeof(G1);
            assert (Hadamard(F_out, R_out, 
                F_out_, R_out_,
                F_cur, R_cur,
                pp.len, pp.omega_gen,
                com_out, com_out_, com_cur,
                pp.pp, ptimer, vtimer, comm));
        }

        return res;
    }

}

vector<uint> Geometric(double p, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, Polynomial& F_out, Polynomial& R_out, G1& com_out)
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
        auto r = Bernoulli(q_bin, pp, comp_timer, ptimer, vtimer, comm, step_F_out, step_R_out, step_com_out);

        comp_timer.start();
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
        comp_timer.stop();

        ptimer.start();
        vtimer.start();
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
        vtimer.stop();
        ptimer.stop();
    }

    return res;
}

vector<int> DiscreteLaplacian(double t, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    Polynomial F1, R1, F2, R2;
    G1 com1, com2;
    double p = 1 - exp(-1.0L / t);
    auto res1 = Geometric(p, log_range, prec, pp, comp_timer, ptimer, vtimer, comm, F1, R1, com1);
    auto res2 = Geometric(p, log_range, prec, pp, comp_timer, ptimer, vtimer, comm, F2, R2, com2);
    
    ptimer.start();
    vtimer.start();
    F_out = F1 - F2;
    R_out = R1 - R2;
    com_out = com1 - com2;
    vtimer.stop();
    ptimer.stop();

    // return res1 - res2, make them signed
    vector<int> out = vector<int>(pp.len, 0);
    for (uint i = 0; i < pp.len; ++ i)
    {
        out[i] = (signed) res1[i] - (signed) res2[i];
    }
    return out;
}

vector<uint> GeometricNew(const vector<vector<bool>>& config, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    vector<uint> res(pp.len, 0);
    auto res_begin = &res.front();
    auto log_range = config.size();
    for (uint i = 0; i < log_range; ++ i)
    {
        const auto& q_bin = config[i];
        assert (q_bin.size());
        assert (q_bin.back());
        Polynomial step_F_out, step_R_out;
        G1 step_com_out;
        auto r = Bernoulli(q_bin, pp, comp_timer, ptimer, vtimer, step_F_out, step_R_out, step_com_out);

        comp_timer.start();
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
        comp_timer.stop();

        ptimer.start();
        vtimer.start();
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
        vtimer.stop();
        ptimer.stop();
    }

    return res;
}

vector<int> DiscreteLaplacianNew(const vector<bool>& z_config, const vector<vector<bool>>& g_config, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out)
{
    Polynomial F_geom_, R_geom_;
    G1 com_geom_;
    auto geom = GeometricNew(g_config, pp, comp_timer, ptimer, vtimer, F_geom_, R_geom_, com_geom_);
    auto F_geom = F_geom_ + Fr(1);
    auto& R_geom = R_geom_;
    auto com_geom = com_geom_ + pp.pp.gVec[0];
    for (auto& g : geom) g += 1;

    Polynomial F_zero, R_zero, F_is_positive, R_is_positive;
    G1 com_zero, com_is_positive;

    auto is_zero = Bernoulli(z_config, pp, comp_timer, ptimer, vtimer, F_zero, R_zero, com_zero); 
    auto is_positive = verifiableUniformBits(F_is_positive, R_is_positive, com_is_positive, pp, comp_timer, ptimer, vtimer);

    vector<int> sign(pp.len, 0);
    for (uint i = 0; i < pp.len; ++ i)
    {
        if (is_zero[i]) sign[i] = 0;
        else if (is_positive[i]) sign[i] = 1;
        else sign[i] = -1;
    }

    Polynomial F_sign = ntt_vec_to_poly_given_omega(vector<Fr>(sign.begin(), sign.end()), pp.omega_gen);
    Polynomial R_sign = randomPolynomial(F_sign.getDegree());
    G1 com_sign = commitPoly(F_sign, R_sign, pp.pp.gVec, pp.pp.hVec);

    auto F_c = F_is_positive * Fr(2) + F_zero - Fr(1) - F_sign;
    auto R_c = R_is_positive * Fr(2) + R_zero - R_sign;
    auto com_c = com_is_positive * Fr(2) + com_zero - com_sign - pp.pp.gVec[0];

    assert (Hadamard(F_c, R_c, 
        F_is_positive * Fr(2), R_is_positive * Fr(2),
        F_zero, R_zero,
        pp.len, pp.omega_gen,
        com_c, com_is_positive * Fr(2), com_zero,
        pp.pp, ptimer, vtimer)); 

    vector<int> out(pp.len, 0);
    for (uint i = 0; i < pp.len; ++ i)
    {
        out[i] = geom[i] * sign[i];
    }
    F_out = ntt_vec_to_poly_given_omega(vector<Fr>(out.begin(), out.end()), pp.omega_gen);
    R_out = randomPolynomial(F_out.getDegree());
    com_out = commitPoly(F_out, R_out, pp.pp.gVec, pp.pp.hVec);

    assert (Hadamard(F_out, R_out, 
        F_geom, R_geom,
        F_sign, R_sign,
        pp.len, pp.omega_gen,
        com_out, com_geom, com_sign,
        pp.pp, ptimer, vtimer));

    return out;

}