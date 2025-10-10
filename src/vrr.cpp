#include "vrr.hpp"

VRRPubParam VRRTrustedSetup(uint num_class, uint A, uint B)
{
    assert (num_class && A && B);
    // num_class is a power of 2
    assert((num_class & (num_class - 1)) == 0);

    uint omega_size = A + (num_class - 1) * B;
    // omega_size is a power of 2
    assert((omega_size & (omega_size - 1)) == 0);

    auto pp = trustedSetup(omega_size + 1);
    const Fr x_gen = getRootOfUnity(ceilLog2(num_class));    
    
    vector<Fr> F_distribution;
    Fr x_power = 1;
    for (uint i = 0; i < num_class; i++)
    {   
        uint num = (i == 0) ? A : B;
        vector<Fr> temp((i == 0) ? A : B, x_power);
        F_distribution.insert(F_distribution.end(), temp.begin(), temp.end());
        x_power *= x_gen;
    }
    assert (F_distribution.size() == omega_size);

    Fr omega_gen;
    Polynomial F = ntt_vec_to_poly(F_distribution, omega_gen);
    // Fr omega_power = 1;
    
    vector<Fr> F_omega_params(omega_size + 1, Fr(0));
    F_omega_params[0] = - 1;
    F_omega_params[omega_size] = 1;
    Polynomial F_omega(F_omega_params);

    G1 com_F = commitPoly(F, pp.gVec);
    G1 com_F_omega = commitPoly(F_omega, pp.gVec);

    return {num_class, omega_size, pp, x_gen, omega_gen, F, F_omega, com_F, com_F_omega};
}

pair<G1, G1> VRRCommit(const Fr& ix, const uint is, const Fr& rx, const Fr& rs, const VRRPubParam& vrrpp)
{
    assert (ix < vrrpp.num_class);
    assert (is < vrrpp.omega_size);
    Fr x, s;
    Fr::pow(x, vrrpp.x_gen, ix);
    Fr::pow(s, vrrpp.omega_gen, is);

    auto& g = vrrpp.pp.gVec[0];
    auto& h = vrrpp.pp.hVec[0];

    return {g * x + h * rx, g * s + h * rs};
}

Fr VRRCompute(const uint ix, const uint is, const uint ir, const VRRPubParam& vrrpp)
{   
    assert (ir < vrrpp.omega_size);
    Fr x, z;
    Fr::pow(x, vrrpp.x_gen, ix);
    Fr::pow(z, vrrpp.omega_gen, (is + ir) % vrrpp.omega_size);
    return x * vrrpp.F(z);
}

bool VRR(const uint ix, const uint is, const Fr& rx, const Fr& rs, 
    const G1& comx, const G1& coms, const Fr& y, const uint ir, 
    const VRRPubParam& vrrpp, Timer &ptimer, Timer &vtimer, uint& communication)
{
    bool accepted = true;
    vtimer.start();
    Fr y_pow;
    Fr::pow(y_pow, y, vrrpp.num_class);
    accepted = accepted && (y_pow == 1);
    vtimer.stop();

    ptimer.start();
    auto it = (is + ir) % vrrpp.omega_size;
    Fr t;
    Fr::pow(t, vrrpp.omega_gen, it);
    vtimer.start();
    Fr r;
    Fr::pow(r, vrrpp.omega_gen, ir);
    auto comt = coms * r;
    vtimer.stop();
    auto rt = rs * r;
    auto z = vrrpp.F(t);
    Fr rz;
    rz.setByCSPRNG();
    auto comz = vrrpp.pp.gVec[0] * z + vrrpp.pp.hVec[0] * rz;
    ptimer.stop();

    vtimer.start();
    Fr alpha;
    alpha.setByCSPRNG();
    vtimer.stop();

    ptimer.start();
    auto F_new = vrrpp.F + vrrpp.F_omega * alpha;
    auto com_F_new = vrrpp.com_F + vrrpp.com_F_omega * alpha;


    accepted &= SecretEval(z, rz, t, rt, F_new, comz, comt, com_F_new, vrrpp.pp, ptimer, vtimer, communication);
    auto& g = vrrpp.pp.gVec[0];
    auto& h = vrrpp.pp.hVec[0];
    ptimer.start();
    Fr x;
    Fr::pow(x, vrrpp.x_gen, ix);
    ptimer.stop();

    accepted &= Prod(y, 0, z, rz, x, rx, 
        g * y, comz, comx, g, h, ptimer, vtimer, communication);

    communication += (sizeof(Fr) + sizeof(G1));

    return accepted;

}