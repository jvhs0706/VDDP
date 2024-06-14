#include "poly-commit.hpp"
using namespace std;
using namespace mcl::bn;

pair<vector<G1>, vector<G1>> trustedSetup(const G1& g, const G1& h, G2& g2, G2& g2_tau, uint maxDeg)
{
    Fr tau;
    tau.setByCSPRNG();
    Fr temp = 1;
    g2_tau = g2 * tau;
    vector<G1> gVec(maxDeg + 1);
    vector<G1> hVec(maxDeg + 1);
    for (uint i = 0; i <= maxDeg; i++)
    {
        gVec[i] = g * temp;
        hVec[i] = h * temp;
        temp *= tau;
    }
    return make_pair(gVec, hVec);
}

G1 commitPoly(const Polynomial<Fr>& F, const std::vector<G1>& gVec)
{
    assert(F.getDegree() >= 0);
    assert(F.getDegree() < gVec.size());

    G1 result = gVec[0] * F[0];
    for (uint i = 1; i <= F.getDegree(); i++)
    {
        result += gVec[i] * F[i];
    }
    return result;

}

G1 commitPoly(const Polynomial<Fr>& F, const Polynomial<Fr>& R, const std::vector<G1>& gVec, const std::vector<G1>& hVec)
{
    assert(F.getDegree() >= 0);
    assert(F.getDegree() < gVec.size());
    assert(R.getDegree() >= 0);
    assert(R.getDegree() < hVec.size());

    G1 result = gVec[0] * F[0] + hVec[0] * R[0];
    for (uint i = 1; i <= F.getDegree(); i++)
    {
        result += gVec[i] * F[i];
    }
    for (uint i = 1; i <= R.getDegree(); i++)
    {
        result += hVec[i] * R[i];
    }
    return result;
}


mcl::bn::G1 provePolyEval(const Polynomial<mcl::bn::Fr>& F, 
    const mcl::bn::Fr& x, mcl::bn::Fr& y, const std::vector<mcl::bn::G1>& gVec)
{
    y = F(x);
    auto diff = Polynomial<Fr>({-x, 1});
    return commitPoly((F - y) / diff, gVec);
}

mcl::bn::G1 provePolyEval(const Polynomial<mcl::bn::Fr>& F, const Polynomial<mcl::bn::Fr>& R, 
    const mcl::bn::Fr& x, 
    mcl::bn::Fr& yF, mcl::bn::Fr& yR,
    const std::vector<mcl::bn::G1>& gVec, const std::vector<mcl::bn::G1>& hVec)
{
    yF = F(x);
    yR = R(x);
    auto diff = Polynomial<Fr>({-x, 1});
    return commitPoly((F - yF) / diff, (R - yR) / diff, gVec, hVec);
}