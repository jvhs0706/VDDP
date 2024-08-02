#include "vdlm.hpp"

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
    vector<bool> res(len, true);

    // assert (r_vec.size() == prec);
    if (r_vec.size() != prec) throw invalid_argument("r_vec.size() != prec");

    for (int i = prec - 1; i >= 0; -- i)
    {
        auto& r = r_vec[i];
        // assert (r.size() == len);
        if (r.size() != len) throw invalid_argument("r.size() != len");
        for (uint j = 0; j < len; ++ j)
        {
            auto t = p_bits[i] != r[j];
            res[j] = t ? !r[j] : res[j];
        }
        
    }
    return res;
}