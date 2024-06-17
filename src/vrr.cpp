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
        cout << x_power << endl;
        uint num = (i == 0) ? A : B;
        vector<Fr> temp((i == 0) ? A : B, x_power);
        F_distribution.insert(F_distribution.end(), temp.begin(), temp.end());
        x_power *= x_gen;
    }
    assert (F_distribution.size() == omega_size);

    Fr omega_gen;
    Polynomial<Fr> F = ntt(F_distribution, omega_gen, true);
    cout << "omega_gen = " << omega_gen << endl;
    Fr omega_gen_power;
    Fr::pow(omega_gen_power, omega_gen, omega_size);
    cout << "omega_gen_power = " << omega_gen_power << endl;
    cout << "F[1] = "<< F[1] << endl;
    cout << "F(1) = " << F(1) << endl;
    cout << "F(omega_gen) = " << F(omega_gen) << endl;
    
    vector<Fr> F_omega_params(omega_size + 1, Fr(0));
    F_omega_params[0] = - 1;
    F_omega_params[omega_size] = 1;
    Polynomial<Fr> F_omega(F_omega_params);

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
    const VRRPubParam& vrrpp, Timer &ptimer, Timer &vtimer)
{
    vtimer.start();
    Fr y_pow;
    Fr::pow(y_pow, y, vrrpp.num_class);
    cout << y_pow << endl;
    assert (y_pow == 1);
    vtimer.stop();

    return true;

}